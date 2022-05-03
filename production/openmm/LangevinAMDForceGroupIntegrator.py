import logging
import re 

import numpy as np
from openmm.unit import *
import openmm.unit as unit
import openmm as mm
from openmm.openmm import LangevinIntegrator

from openmmtools.constants import kB
from openmmtools import respa, utils

logger = logging.getLogger(__name__)

# Energy unit used by OpenMM unit system
_OPENMM_ENERGY_UNIT = kilojoules_per_mole

# ============================================================================================
# BASE CLASSES
# ============================================================================================

class PrettyPrintableIntegrator(object):
    """A PrettyPrintableIntegrator can format the contents of its step program for printing.
    This is a mix-in.
    TODO: We should check that the object (`self`) is a CustomIntegrator or subclass.
    """
    def pretty_format(self, as_list=False, step_types_to_highlight=None):
        """Generate a human-readable version of each integrator step.
        Parameters
        ----------
        as_list : bool, optional, default=False
           If True, a list of human-readable strings will be returned.
           If False, these will be concatenated into a single human-readable string.
        step_types_to_highlight : list of int, optional, default=None
           If specified, these step types will be highlighted.
        Returns
        -------
        readable_lines : list of str
           A list of human-readable versions of each step of the integrator
        """
        step_type_dict = {
            0 : "{target} <- {expr}",
            1: "{target} <- {expr}",
            2: "{target} <- sum({expr})",
            3: "constrain positions",
            4: "constrain velocities",
            5: "allow forces to update the context state",
            6: "if({expr}):",
            7: "while({expr}):",
            8: "end"
        }

        if not hasattr(self, 'getNumComputations'):
            raise Exception('This integrator is not a CustomIntegrator.')

        readable_lines = []
        indent_level = 0
        for step in range(self.getNumComputations()):
            line = ''
            step_type, target, expr = self.getComputationStep(step)
            highlight = True if (step_types_to_highlight is not None) and (step_type in step_types_to_highlight) else False
            if step_type in [8]:
                indent_level -= 1
            if highlight:
                line += '\x1b[6;30;42m'
            line += 'step {:6d} : '.format(step) + '   ' * indent_level + step_type_dict[step_type].format(target=target, expr=expr)
            if highlight:
                line += '\x1b[0m'
            if step_type in [6, 7]:
                indent_level += 1
            readable_lines.append(line)

        if as_list:
            return readable_lines
        else:
            return '\n'.join(readable_lines)

    def pretty_print(self):
        """Pretty-print the computation steps of this integrator."""
        print(self.pretty_format())


class ThermostatedIntegrator(utils.RestorableOpenMMObject, PrettyPrintableIntegrator,
                             mm.CustomIntegrator):
    """Add temperature functions to a CustomIntegrator.
    This class is intended to be inherited by integrators that maintain the
    stationary distribution at a given temperature. The constructor adds a
    global variable named "kT" defining the thermal energy at the given
    temperature. This global variable is updated through the temperature
    setter and getter.
    It also provide a utility function to handle per-DOF constants that
    must be computed only when the temperature changes.
    Notice that the CustomIntegrator internally stored by a Context object
    will loose setter and getter and any extra function you define. The same
    happens when you copy your integrator. You can restore the methods with
    the static method ThermostatedIntegrator.restore_interface().
    Parameters
    ----------
    temperature : unit.Quantity
        The temperature of the integrator heat bath (temperature units).
    timestep : unit.Quantity
        The timestep to pass to the CustomIntegrator constructor (time
        units).
    Examples
    --------
    We can inherit from ThermostatedIntegrator to automatically define
    setters and getters for the temperature and to add a per-DOF constant
    "sigma" that we need to update only when the temperature is changed.
    >>> from openmm import openmm, unit
    >>> class TestIntegrator(ThermostatedIntegrator):
    ...     def __init__(self, temperature=298.0*unit.kelvin, timestep=1.0*unit.femtoseconds):
    ...         super(TestIntegrator, self).__init__(temperature, timestep)
    ...         self.addPerDofVariable("sigma", 0)  # velocity standard deviation
    ...         self.addComputeTemperatureDependentConstants({"sigma": "sqrt(kT/m)"})
    ...
    We instantiate the integrator normally.
    >>> integrator = TestIntegrator(temperature=350*unit.kelvin)
    >>> integrator.getTemperature()
    Quantity(value=350.0, unit=kelvin)
    >>> integrator.setTemperature(380.0*unit.kelvin)
    >>> integrator.getTemperature()
    Quantity(value=380.0, unit=kelvin)
    >>> integrator.getGlobalVariableByName('kT')
    3.1594995390636815
    Notice that a CustomIntegrator loses any extra method after a serialization cycle.
    >>> integrator_serialization = openmm.XmlSerializer.serialize(integrator)
    >>> deserialized_integrator = openmm.XmlSerializer.deserialize(integrator_serialization)
    >>> deserialized_integrator.getTemperature()
    Traceback (most recent call last):
    ...
    AttributeError: type object 'object' has no attribute '__getattr__'
    We can restore the original interface with a class method
    >>> ThermostatedIntegrator.restore_interface(integrator)
    True
    >>> integrator.getTemperature()
    Quantity(value=380.0, unit=kelvin)
    >>> integrator.setTemperature(400.0*unit.kelvin)
    >>> isinstance(integrator, TestIntegrator)
    True
    """
    def __init__(self, temperature, *args, **kwargs):
        super(ThermostatedIntegrator, self).__init__(*args, **kwargs)
        self.addGlobalVariable('kT', kB * temperature)  # thermal energy

    @property
    def global_variable_names(self):
        """The set of global variable names defined for this integrator."""
        return set([ self.getGlobalVariableName(index) for index in range(self.getNumGlobalVariables()) ])

    def getTemperature(self):
        """Return the temperature of the heat bath.
        Returns
        -------
        temperature : unit.Quantity
            The temperature of the heat bath in kelvins.
        """
        # Do most unit conversion first for precision
        conversion = _OPENMM_ENERGY_UNIT / kB
        temperature = self.getGlobalVariableByName('kT') * conversion
        return temperature

    def setTemperature(self, temperature):
        """Set the temperature of the heat bath.
        Parameters
        ----------
        temperature : unit.Quantity
            The new temperature of the heat bath (temperature units).
        """
        kT = kB * temperature
        self.setGlobalVariableByName('kT', kT)

        # Update the changed flag if it exist.
        if 'has_kT_changed' in self.global_variable_names:
            self.setGlobalVariableByName('has_kT_changed', 1)

    def addComputeTemperatureDependentConstants(self, compute_per_dof):
        """Wrap the ComputePerDof into an if-block executed only when kT changes.
        Parameters
        ----------
        compute_per_dof : dict of str: str
            A dictionary of variable_name: expression.
        """
        # First check if flag variable already exist.
        if not 'has_kT_changed' in self.global_variable_names:
            self.addGlobalVariable('has_kT_changed', 1)

        # Create if-block that conditionally update the per-DOF variables.
        self.beginIfBlock('has_kT_changed = 1')
        for variable, expression in compute_per_dof.items():
            self.addComputePerDof(variable, expression)
        self.addComputeGlobal('has_kT_changed', '0')
        self.endBlock()

    @classmethod
    def is_thermostated(cls, integrator):
        """Return true if the integrator is a ThermostatedIntegrator.
        This can be useful when you only have access to the Context
        CustomIntegrator, which loses all extra function during serialization.
        Parameters
        ----------
        integrator : simtk.openmm.Integrator
            The integrator to check.
        Returns
        -------
        True if the original CustomIntegrator class inherited from
        ThermostatedIntegrator, False otherwise.
        """
        global_variable_names = set([ integrator.getGlobalVariableName(index) for index in range(integrator.getNumGlobalVariables()) ])
        if not 'kT' in global_variable_names:
            return False
        return super(ThermostatedIntegrator, cls).is_restorable(integrator)

    @classmethod
    def restore_interface(cls, integrator):
        """Restore the original interface of a CustomIntegrator.
        The function restore the methods of the original class that
        inherited from ThermostatedIntegrator. Return False if the interface
        could not be restored.
        Parameters
        ----------
        integrator : simtk.openmm.CustomIntegrator
            The integrator to which add methods.
        Returns
        -------
        True if the original class interface could be restored, False otherwise.
        """
        restored = super(ThermostatedIntegrator, cls).restore_interface(integrator)
        # Warn the user if he is implementing a CustomIntegrator
        # that may keep the stationary distribution at a certain
        # temperature without exposing getters and setters.
        if not restored:
            if hasattr(integrator, 'getGlobalVariableName'):
                global_variable_names = set([ integrator.getGlobalVariableName(index) for index in range(integrator.getNumGlobalVariables()) ])
                if 'kT' in global_variable_names:
                    if not hasattr(integrator, 'getTemperature'):
                        logger.warning("The integrator {} has a global variable 'kT' variable "
                                       "but does not expose getter and setter for the temperature. "
                                       "Consider inheriting from ThermostatedIntegrator.")
        return restored

    @property
    def kT(self):
        """The thermal energy in simtk.openmm.Quantity"""
        return self.getGlobalVariableByName("kT") * _OPENMM_ENERGY_UNIT


class LangevinAMDForceGroupIntegrator(ThermostatedIntegrator):
    """Integrates Langevin dynamics with a prescribed operator splitting and aMD "Boost potential".
    Example: integrator = LangevinAMDForceGroupIntegrator(2092*kilojoules_per_mole, 14644*kilojoules_per_mole,group=2, temperature=310.0*kelvin, collision_rate=1.0 / picoseconds, splitting="V R O R V")
    One way to divide the Langevin system is into three parts which can each be solved "exactly:"
        - R: Linear "drift" / Constrained "drift"
            Deterministic update of *positions*, using current velocities
            x <- x + v dt
        - V: Linear "kick" / Constrained "kick"
            Deterministic update of *velocities*, using current forces
            v <- v + (f/m) dt
                where f = force, m = mass
        - O: Ornstein-Uhlenbeck
            Stochastic update of velocities, simulating interaction with a heat bath
            v <- av + b sqrt(kT/m) R
                where
                a = e^(-gamma dt)
                b = sqrt(1 - e^(-2gamma dt))
                R is i.i.d. standard normal
    We can then construct integrators by solving each part for a certain timestep in sequence.
    (We can further split up the V step by force group, evaluating cheap but fast-fluctuating
    forces more frequently than expensive but slow-fluctuating forces. Since forces are only
    evaluated in the V step, we represent this by including in our "alphabet" V0, V1, ...)
    When the system contains holonomic constraints, these steps are confined to the constraint
    manifold.
    Examples
    --------
        - VVVR
            splitting="O V R V O"
        - BAOAB:
            splitting="V R O R V"
        - g-BAOAB, with K_r=3:
            splitting="V R R R O R R R V"
        - g-BAOAB with solvent-solute splitting, K_r=K_p=2:
            splitting="V0 V1 R R O R R V1 R R O R R V1 V0"
    Attributes
    ----------
    _kinetic_energy : str
        This is 0.5*m*v*v by default, and is the expression used for the kinetic energy
    shadow_work : unit.Quantity with units of energy
       Shadow work (if integrator was constructed with measure_shadow_work=True)
    heat : unit.Quantity with units of energy
       Heat (if integrator was constructed with measure_heat=True)
    References
    ----------
    [Leimkuhler and Matthews, 2015] Molecular dynamics: with deterministic and stochastic numerical methods, Chapter 7
    """

    _kinetic_energy = "0.5 * m * v * v"

    def __init__(self,
                 alphaGroup,
                 EGroup,
                 group = 2,
                 temperature=298.0 * unit.kelvin,
                 collision_rate=1.0 / unit.picoseconds,
                 timestep=1.0 * unit.femtoseconds,
                 splitting="V R O R V",
                 constraint_tolerance=1e-8,
                 measure_shadow_work=False,
                 measure_heat=False,
                 ):
        """Create a Langevin integrator with the prescribed operator splitting.
        Parameters
        ----------
        splitting : string, default: "V R O R V"
            Sequence of "R", "V", "O" (and optionally "{", "}", "V0", "V1", ...) substeps to be executed each timestep.
            Forces are only used in V-step. Handle multiple force groups by appending the force group index
            to V-steps, e.g. "V0" will only use forces from force group 0. "V" will perform a step using all forces.
            "{" will cause Metropolization, and must be followed later by a "}".
        temperature : np.unit.Quantity compatible with kelvin, default: 298.0*unit.kelvin
           Fictitious "bath" temperature
        collision_rate : np.unit.Quantity compatible with 1/picoseconds, default: 1.0/unit.picoseconds
           Collision rate
        timestep : np.unit.Quantity compatible with femtoseconds, default: 1.0*unit.femtoseconds
           Integration timestep
        constraint_tolerance : float, default: 1.0e-8
            Tolerance for constraint solver
        measure_shadow_work : boolean, default: False
            Accumulate the shadow work performed by the symplectic substeps, in the global `shadow_work`
        measure_heat : boolean, default: False
            Accumulate the heat exchanged with the bath in each step, in the global `heat`
        """

        # Compute constants
        gamma = collision_rate
        self._gamma = gamma

        # assign Boost group
        self._group = group

        # Check if integrator is metropolized by checking for M step:
        if splitting.find("{") > -1:
            self._metropolized_integrator = True
            # We need to measure shadow work if Metropolization is used
            measure_shadow_work = True
        else:
            self._metropolized_integrator = False

        # Record whether we are measuring heat and shadow work
        self._measure_heat = measure_heat
        self._measure_shadow_work = measure_shadow_work

        ORV_counts, mts, force_group_nV = self._parse_splitting_string(splitting)

        # Record splitting.
        self._splitting = splitting
        self._ORV_counts = ORV_counts
        self._mts = mts
        self._force_group_nV = force_group_nV

        # Create a new CustomIntegrator
        super(LangevinAMDForceGroupIntegrator, self).__init__(temperature, timestep)

        # Initialize
        self.addPerDofVariable("sigma", 0)

        # Velocity mixing parameter: current velocity component
        h = timestep / max(1, ORV_counts['O'])
        self.addGlobalVariable("a", np.exp(-gamma * h))

        # Velocity mixing parameter: random velocity component
        self.addGlobalVariable("b", np.sqrt(1 - np.exp(- 2 * gamma * h)))

        # Positions before application of position constraints
        self.addPerDofVariable("x1", 0)

        # Set constraint tolerance
        self.setConstraintTolerance(constraint_tolerance)

        # Set aMD boost parameters
        self.addGlobalVariable("alphaGroup", alphaGroup)
        self.addGlobalVariable("EGroup", EGroup)
        self.addGlobalVariable("groupEnergy", 0)
        self.addGlobalVariable("deltaV", 0)
        self.addPerDofVariable("fg", 0)

        # Add global variables
        self._add_global_variables()

        # Add integrator steps
        self._add_integrator_steps()

    @property
    def _step_dispatch_table(self):
        """dict: The dispatch table step_name -> add_step_function."""
        # TODO use methoddispatch (see yank.utils) when dropping Python 2 support.
        dispatch_table = {
            'O': (self._add_O_step, False),
            'R': (self._add_R_step, False),
            '{': (self._add_metropolize_start, False),
            '}': (self._add_metropolize_finish, False),
            'V': (self._add_V_step, True)
        }
        return dispatch_table

    def _add_global_variables(self):
        """Add global bookkeeping variables."""
        if self._measure_heat:
            self.addGlobalVariable("heat", 0)

        if self._measure_shadow_work or self._measure_heat:
            self.addGlobalVariable("old_ke", 0)
            self.addGlobalVariable("new_ke", 0)

        if self._measure_shadow_work:
            self.addGlobalVariable("old_pe", 0)
            self.addGlobalVariable("new_pe", 0)
            self.addGlobalVariable("shadow_work", 0)

        # If we metropolize, we have to keep track of the before and after (x, v)
        if self._metropolized_integrator:
            self.addGlobalVariable("accept", 0)
            self.addGlobalVariable("ntrials", 0)
            self.addGlobalVariable("nreject", 0)
            self.addGlobalVariable("naccept", 0)
            self.addPerDofVariable("vold", 0)
            self.addPerDofVariable("xold", 0)

    def reset_heat(self):
        """Reset heat."""
        if self._measure_heat:
            self.setGlobalVariableByName('heat', 0.0)

    def reset_shadow_work(self):
        """Reset shadow work."""
        if self._measure_shadow_work:
            self.setGlobalVariableByName('shadow_work', 0.0)

    def reset_ghmc_statistics(self):
        """Reset GHMC acceptance rate statistics."""
        if self._metropolized_integrator:
            self.setGlobalVariableByName('ntrials', 0)
            self.setGlobalVariableByName('naccept', 0)
            self.setGlobalVariableByName('nreject', 0)

    def reset(self):
        """Reset all statistics (heat, shadow work, acceptance rates, step).
        """
        self.reset_heat()
        self.reset_shadow_work()
        self.reset_ghmc_statistics()

    def _get_energy_with_units(self, variable_name, dimensionless=False):
        """Retrive an energy/work quantity and return as unit-bearing or dimensionless quantity.
        Parameters
        ----------
        variable_name : str
           Name of the global context variable to retrieve
        dimensionless : bool, optional, default=False
           If specified, the energy/work is returned in reduced (kT) unit.
        Returns
        -------
        work : unit.Quantity or float
           If dimensionless=True, the work in kT (float).
           Otherwise, the unit-bearing work in units of energy.
        """
        work = self.getGlobalVariableByName(variable_name) * _OPENMM_ENERGY_UNIT
        if dimensionless:
            return work / self.kT
        else:
            return work

    def get_shadow_work(self, dimensionless=False):
        """Get the current accumulated shadow work.
        Parameters
        ----------
        dimensionless : bool, optional, default=False
           If specified, the work is returned in reduced (kT) unit.
        Returns
        -------
        work : unit.Quantity or float
           If dimensionless=True, the protocol work in kT (float).
           Otherwise, the unit-bearing protocol work in units of energy.
        """
        if not self._measure_shadow_work:
            raise Exception("This integrator must be constructed with 'measure_shadow_work=True' to measure shadow work.")
        return self._get_energy_with_units("shadow_work", dimensionless=dimensionless)

    @property
    def shadow_work(self):
        return self.get_shadow_work()

    def get_heat(self, dimensionless=False):
        """Get the current accumulated heat.
        Parameters
        ----------
        dimensionless : bool, optional, default=False
           If specified, the work is returned in reduced (kT) unit.
        Returns
        -------
        work : unit.Quantity or float
           If dimensionless=True, the heat in kT (float).
           Otherwise, the unit-bearing heat in units of energy.
        """
        if not self._measure_heat:
            raise Exception("This integrator must be constructed with 'measure_heat=True' in order to measure heat.")
        return self._get_energy_with_units("heat", dimensionless=dimensionless)

    @property
    def heat(self):
        return self.get_heat()

    def get_acceptance_rate(self):
        """Get acceptance rate for Metropolized integrators.
        Returns
        -------
        acceptance_rate : float
           Acceptance rate.
           An Exception is thrown if the integrator is not Metropolized.
        """
        if not self._metropolized_integrator:
            raise Exception("This integrator must be Metropolized to return an acceptance rate.")
        return self.getGlobalVariableByName("naccept") / self.getGlobalVariableByName("ntrials")

    @property
    def acceptance_rate(self):
        """Get acceptance rate for Metropolized integrators."""
        return self.get_acceptance_rate()

    @property
    def is_metropolized(self):
        """Return True if this integrator is Metropolized, False otherwise."""
        return self._metropolized_integrator

    def _add_integrator_steps(self):
        """Add the steps to the integrator--this can be overridden to place steps around the integration.
        """
        # Integrate
        self.addUpdateContextState()
        self.addComputeTemperatureDependentConstants({"sigma": "sqrt(kT/m)"})

        for i, step in enumerate(self._splitting.split()):
            self._substep_function(step)

    def _sanity_check(self, splitting):
        """Perform a basic sanity check on the splitting string to ensure that it makes sense.
        Parameters
        ----------
        splitting : str
            The string specifying the integrator splitting
        mts : bool
            Whether the integrator is a multiple timestep integrator
        allowed_characters : str, optional
            The characters allowed to be present in the splitting string.
            Default RVO and the digits 0-9.
        """

        # Space is just a delimiter--remove it
        splitting_no_space = splitting.replace(" ", "")

        allowed_characters = "0123456789"
        for key in self._step_dispatch_table:
            allowed_characters += key

        # sanity check to make sure only allowed combinations are present in string:
        for step in splitting.split():
            if step[0]=="V":
                if len(step) > 1:
                    try:
                        force_group_number = int(step[1:])
                        if force_group_number > 31:
                            raise ValueError("OpenMM only allows up to 32 force groups")
                    except ValueError:
                        raise ValueError("You must use an integer force group")
            elif step == "{":
                    if "}" not in splitting:
                        raise ValueError("Use of { must be followed by }")
                    if not self._verify_metropolization(splitting):
                        raise ValueError("Shadow work generating steps found outside the Metropolization block")
            elif step in allowed_characters:
                continue
            else:
                raise ValueError("Invalid step name '%s' used; valid step names are %s" % (step, allowed_characters))

        # Make sure we contain at least one of R, V, O steps
        assert ("R" in splitting_no_space)
        assert ("V" in splitting_no_space)
        assert ("O" in splitting_no_space)

    def _verify_metropolization(self, splitting):
        """Verify that the shadow-work generating steps are all inside the metropolis block
        Returns False if they are not.
        Parameters
        ----------
        splitting : str
            The langevin splitting string
        Returns
        -------
        valid_metropolis : bool
            Whether all shadow-work generating steps are in the {} block
        """
        # check that there is exactly one metropolized region
        #this pattern matches the { literally, then any number of any character other than }, followed by another {
        #If there's a match, then we have an attempt at a nested metropolization, which is unsupported
        regex_nested_metropolis = "{[^}]*{"
        pattern = re.compile(regex_nested_metropolis)
        if pattern.match(splitting.replace(" ", "")):
            raise ValueError("There can only be one Metropolized region.")

        # find the metropolization steps:
        M_start_index = splitting.find("{")
        M_end_index = splitting.find("}")

        # accept/reject happens before the beginning of metropolis step
        if M_start_index > M_end_index:
            return False

        #pattern to find whether any shadow work producing steps lie outside the metropolization region
        RV_outside_metropolis = "[RV](?![^{]*})"
        outside_metropolis_check = re.compile(RV_outside_metropolis)
        if outside_metropolis_check.match(splitting.replace(" ","")):
            return False
        else:
            return True

    def _add_R_step(self):
        """Add an R step (position update) given the velocities.
        """
        if self._measure_shadow_work:
            self.addComputeGlobal("old_pe", "energy")
            self.addComputeSum("old_ke", self._kinetic_energy)

        n_R = self._ORV_counts['R']

        # update positions (and velocities, if there are constraints)
        self.addComputePerDof("x", "x + ((dt / {}) * v)".format(n_R))
        self.addComputePerDof("x1", "x")  # save pre-constraint positions in x1
        self.addConstrainPositions()  # x is now constrained
        self.addComputePerDof("v", "v + ((x - x1) / (dt / {}))".format(n_R))
        self.addConstrainVelocities()

        if self._measure_shadow_work:
            self.addComputeGlobal("new_pe", "energy")
            self.addComputeSum("new_ke", self._kinetic_energy)
            self.addComputeGlobal("shadow_work", "shadow_work + (new_ke + new_pe) - (old_ke + old_pe)")

    def _add_V_step(self, force_group="0"):
        """Deterministic velocity update, using only forces from force-group fg.
        Parameters
        ----------
        force_group : str, optional, default="0"
           Force group to use for this step
        """
        if self._measure_shadow_work:
            self.addComputeSum("old_ke", self._kinetic_energy)

        # update velocities, with _mts, the boost is not yet implemented.
        if self._mts:
            self.addComputePerDof("v", "v + ((dt / {}) * f{} / m)".format(self._force_group_nV[force_group], force_group))
        else:
            self.addComputeGlobal("groupEnergy", "energy"+str(self._group))
            self.addComputeGlobal("deltaV", "modify*(EGroup-groupEnergy)^2/(alphaGroup+EGroup-groupEnergy); modify=step(EGroup-groupEnergy)")
            self.addComputePerDof("fg", "f"+str(self._group))
            self.addComputePerDof("v", "v + (dt / {}) * fprime / m; fprime=fother+fg*((1-modify) + modify*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2); fother=f-fg; modify=step(EGroup-groupEnergy)".format(self._force_group_nV["0"]))

        self.addConstrainVelocities()

        if self._measure_shadow_work:
            self.addComputeSum("new_ke", self._kinetic_energy)
            self.addComputeGlobal("shadow_work", "shadow_work + (new_ke - old_ke)")

    def _add_O_step(self):
        """Add an O step (stochastic velocity update)
        """
        if self._measure_heat:
            self.addComputeSum("old_ke", self._kinetic_energy)

        # update velocities
        self.addComputePerDof("v", "(a * v) + (b * sigma * gaussian)")
        self.addConstrainVelocities()

        if self._measure_heat:
            self.addComputeSum("new_ke", self._kinetic_energy)
            self.addComputeGlobal("heat", "heat + (new_ke - old_ke)")

    def _substep_function(self, step_string):
        """Take step string, and add the appropriate R, V, O step with appropriate parameters.
        The step string input here is a single character (or character + number, for MTS)
        """
        function, can_accept_force_groups = self._step_dispatch_table[step_string[0]]
        if can_accept_force_groups:
            force_group = step_string[1:]
            function(force_group)
        else:
            function()

    def _parse_splitting_string(self, splitting_string):
        """Parse the splitting string to check for simple errors and extract necessary information
        Parameters
        ----------
        splitting_string : str
            The string that specifies how to do the integrator splitting
        Returns
        -------
        ORV_counts : dict
            Number of O, R, and V steps
        mts : bool
            Whether the splitting specifies an MTS integrator
        force_group_n_V : dict
            Specifies the number of V steps per force group. {"0": nV} if not MTS
        """
        # convert the string to all caps
        splitting_string = splitting_string.upper()

        # sanity check the splitting string
        self._sanity_check(splitting_string)

        ORV_counts = dict()

        # count number of R, V, O steps:
        for step_symbol in self._step_dispatch_table:
            ORV_counts[step_symbol] = splitting_string.count(step_symbol)

        # split by delimiter (space)
        step_list = splitting_string.split(" ")

        # populate a list with all the force groups in the system
        force_group_list = []
        for step in step_list:
            # if the length of the step is greater than one, it has a digit after it
            if step[0] == "V" and len(step) > 1:
                force_group_list.append(step[1:])

        # Make a set to count distinct force groups
        force_group_set = set(force_group_list)

        # check if force group list cast to set is longer than one
        # If it is, then multiple force groups are specified
        if len(force_group_set) > 1:
            mts = True
        else:
            mts = False

        # If the integrator is MTS, count how many times the V steps appear for each
        if mts:
            force_group_n_V = {force_group: 0 for force_group in force_group_set}
            for step in step_list:
                if step[0] == "V":
                    # ensure that there are no V-all steps if it's MTS
                    assert len(step) > 1
                    # extract the index of the force group from the step
                    force_group_idx = step[1:]
                    # increment the number of V calls for that force group
                    force_group_n_V[force_group_idx] += 1
        else:
            force_group_n_V = {"0": ORV_counts["V"]}

        return ORV_counts, mts, force_group_n_V

    def _add_metropolize_start(self):
        """Save the current x and v for a metropolization step later"""
        self.addComputePerDof("xold", "x")
        self.addComputePerDof("vold", "v")

    def _add_metropolize_finish(self):
        """Add a Metropolization (based on shadow work) step to the integrator.
        When Metropolization occurs, shadow work is reset.
        """
        self.addComputeGlobal("accept", "step(exp(-(shadow_work)/kT) - uniform)")
        self.addComputeGlobal("ntrials", "ntrials + 1")
        self.beginIfBlock("accept != 1")
        self.addComputePerDof("x", "xold")
        self.addComputePerDof("v", "-vold")
        self.addComputeGlobal("nreject", "nreject + 1")
        self.endBlock()
        self.addComputeGlobal("naccept", "ntrials - nreject")
        self.addComputeGlobal("shadow_work", "0")

    def getAlphaGroup(self):
        """Get the value of alpha for the boosted force group."""
        return self.getGlobalVariable(4)*kilojoules_per_mole

    def setAlphaGroup(self, alpha):
        """Set the value of alpha for the boosted force group."""
        self.setGlobalVariable(4, alpha)

    def getEGroup(self):
        """Get the energy threshold E for the boosted force group."""
        return self.getGlobalVariable(5)*kilojoules_per_mole

    def setEGroup(self, E):
        """Set the energy threshold E for the boosted force group."""
        self.setGlobalVariable(5, E)

    def getDeltaV(self):
        """Return the value of the boost potential deltaV for the system at one instant"""
        return self.getGlobalVariable(7)

    def getEffectiveEnergy(self):
        """Given the actual group energy of the system, return the value of the effective potential.

        Parameters
        ----------
        groupEnergy : energy
            the actual potential energy of the boosted force group

        Returns
        -------
        energy
            the value of the effective potential
        """
        groupEnergy = self.getGlobalVariable(6)
        dV = self.getGlobalVariable(7)
        return groupEnergy+dV


