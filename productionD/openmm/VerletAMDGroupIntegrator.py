from openmm import *

class VerletAMDForceGroupIntegrator(CustomIntegrator):
    """AMDForceGroupIntegrator implements a single boost aMD integration algorithm.

    This is similar to AMDIntegrator, but is applied based on the energy of a single force group
    (typically representing torsions).

    For details, see Hamelberg et al., J. Chem. Phys. 127, 155102 (2007).
    """

    def __init__(self, dt, group, alphaGroup, EGroup):
        """Create a AMDForceGroupIntegrator.

        Parameters
        ----------
        dt : time
            The integration time step to use
        group : int
            The force group to apply the boost to
        alphaGroup : energy
            The alpha parameter to use for the boosted force group
        EGroup : energy
            The energy cutoff to use for the boosted force group
        """
        CustomIntegrator.__init__(self, dt)
        self.addGlobalVariable("alphaGroup", alphaGroup)
        self.addGlobalVariable("EGroup", EGroup)
        self.addGlobalVariable("groupEnergy", 0)
        self.addGlobalVariable("deltaV", 0)
        self.addPerDofVariable("x1", 0)
        self.addPerDofVariable("fprime", 0)
        self.addPerDofVariable("fg", 0)
        self.addUpdateContextState();
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputeGlobal("deltaV", "modify*(EGroup - groupEnergy)^2/(alphaGroup+EGroup-groupEnergy); modify=step(EGroup-groupEnergy)")
        self.addComputePerDof("fg", "f"+str(group))
        self.addComputePerDof("fprime", "fother + fg*((1-modify) + modify*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2); fother=f-fg; modify=step(EGroup-groupEnergy)")
        self.addComputePerDof("v", "v+0.5*dt*fprime/m")
        self.addComputePerDof("x", "x+dt*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputePerDof("fprime", "fother + fg*((1-modify) + modify*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2); fother=f-fg; modify=step(EGroup-groupEnergy)")
        self.addComputePerDof("v", "v+0.5*dt*fprime/m+(x-x1)/dt")
        self.addConstrainVelocities()

    def getAlphaGroup(self):
        """Get the value of alpha for the boosted force group."""
        return self.getGlobalVariable(0)*kilojoules_per_mole

    def setAlphaGroup(self, alpha):
        """Set the value of alpha for the boosted force group."""
        self.setGlobalVariable(0, alpha)

    def getEGroup(self):
        """Get the energy threshold E for the boosted force group."""
        return self.getGlobalVariable(1)*kilojoules_per_mole

    def setEGroup(self, E):
        """Set the energy threshold E for the boosted force group."""
        self.setGlobalVariable(1, E)

    def getDeltaV(self, groupEnergy):
        """Return the value of the boost potential deltaV for the system at one instant"""
        return self.getGlobalVariable(3)*kilojoules_per_mole

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
        groupEnergy = self.getGlobalVariable(2)*kilojoules_per_mole
        deltaV = self.getGlobalVariable(3)*kilojoules_per_mole
        return groupEnergy+deltaV

