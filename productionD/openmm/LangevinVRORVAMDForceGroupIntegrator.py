from openmm.openmm import CustomIntegrator
from openmm.unit import *
import numpy as np

kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
class LangevinVRORVAMDForceGroupIntegrator(CustomIntegrator):
    def __init__(self, alphaGroup, EGroup, group=2, temperature = 310*kelvin,
                  collision_rate=1.0/picoseconds, timestep=2.0*femtoseconds):
        gamma = collision_rate*picoseconds    # stupid "SWIG" thing with numpy
        super(LangevinVRORVAMDForceGroupIntegrator, self).__init__(timestep)
        self.addPerDofVariable("sigma", 0)
        self.addGlobalVariable('kT', kB * temperature)
        self.addGlobalVariable("a", np.exp(-gamma))
        self.addGlobalVariable("b", np.sqrt(1 - np.exp(-2 * gamma)))
        self.addPerDofVariable("x1", 0) # before position restraints
        
        # Set aMD boost parameters
        self.addGlobalVariable("alphaGroup", alphaGroup)
        self.addGlobalVariable("EGroup", EGroup)
        self.addGlobalVariable("groupEnergy", 0)
        self.addGlobalVariable("deltaV", 0)
        self.addPerDofVariable("fg", 0)
        
        self.addUpdateContextState()
        # Compute thermal properties, kT doesn't change here
        self.addComputePerDof("sigma", "sqrt(kT/m)")
        # V step
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputeGlobal("deltaV", "modify*(EGroup-groupEnergy)^2/(alphaGroup+EGroup-groupEnergy); modify=step(EGroup-groupEnergy)")
        self.addComputePerDof("fg", "f"+str(group))
        self.addComputePerDof("v", "v + (dt / 2) * fprime / m; fprime=fother+fg*((1-modify) + modify*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2); fother=f-fg; modify=step(EGroup-groupEnergy)")
        #self.addConstrainVelocities()
        # R Step
        self.addComputePerDof("x", "x + (dt / 2)*v")
        self.addComputePerDof("x1", "x")  # save pre-constraint positions in x1
        self.addConstrainPositions()      # x is now constrained
        self.addComputePerDof("v", "v + (x - x1)/(dt / 2)")
        #self.addConstrainVelocities()        
        # O Step
        self.addComputePerDof("v", "(a * v) + (b * sigma * gaussian)")
        #self.addConstrainVelocities()
        # R Step
        self.addComputePerDof("x", "x + (dt / 2)*v")
        self.addComputePerDof("x1", "x")  # save pre-constraint positions in x1
        self.addConstrainPositions()      # x is now constrained
        self.addComputePerDof("v", "v + (x - x1) / (dt / 2)")
        #self.addConstrainVelocities()
        # V step
        self.addComputeGlobal("groupEnergy", "energy"+str(group))
        self.addComputeGlobal("deltaV", "modify*(EGroup-groupEnergy)^2/(alphaGroup+EGroup-groupEnergy); modify=step(EGroup-groupEnergy)")
        self.addComputePerDof("fg", "f"+str(group))
        self.addComputePerDof("v", "v + (dt / 2) * fprime / m; fprime=fother+fg*((1-modify) + modify*(alphaGroup/(alphaGroup+EGroup-groupEnergy))^2); fother=f-fg; modify=step(EGroup-groupEnergy)")
        self.addConstrainVelocities()

