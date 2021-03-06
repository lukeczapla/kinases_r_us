from openmm.unit import *
from openmm import *
from openmm.app import *

def CreateCVForce(system):

    cvforce = CustomCVForce("0")
    system.addForce(cvforce)
    D1 = CustomBondForce("r")
    D1.addBond(929, 2444)
    print(D1.getNumGlobalParameters())
    #D1.addGlobalParameter("dist", "r")
    D2 = CustomBondForce("r")
    D2.addBond(620, 2444)
    #D2.addGlobalParameter("dist", "r")
    X_Phi = CustomTorsionForce("theta")
    X_Phi.addTorsion(2407,2409,2411,2417)
    X_Psi = CustomTorsionForce("theta")
    X_Psi.addTorsion(2409,2411,2417,2419)
    D_Phi = CustomTorsionForce("theta")
    D_Phi.addTorsion(2417,2419,2421,2430)
    D_Psi = CustomTorsionForce("theta")
    D_Psi.addTorsion(2419,2421,2430,2432)
    F_Phi = CustomTorsionForce("theta")
    F_Phi.addTorsion(2430,2432,2434,2450)
    F_Psi = CustomTorsionForce("theta")
    F_Psi.addTorsion(2432,2434,2450,2452)
    D_Chi1 = CustomTorsionForce("theta")
    D_Chi1.addTorsion(2419,2421,2423,2426)
    D_Chi2 = CustomTorsionForce("theta")
    D_Chi2.addTorsion(2421,2423,2426,2428)
    F_Chi1 = CustomTorsionForce("theta")
    F_Chi1.addTorsion(2432,2434,2436,2439)
    F_Chi2 = CustomTorsionForce("theta")
    F_Chi2.addTorsion(2434,2436,2439,2440)
    cvforce.addCollectiveVariable("D1", D1)
    cvforce.addCollectiveVariable("D2", D2)
    cvforce.addCollectiveVariable("X_Phi", X_Phi)
    cvforce.addCollectiveVariable("X_Psi", X_Psi)
    cvforce.addCollectiveVariable("D_Phi", D_Phi)
    cvforce.addCollectiveVariable("D_Psi", D_Psi)
    cvforce.addCollectiveVariable("F_Phi", F_Phi)
    cvforce.addCollectiveVariable("F_Psi", F_Psi)
    cvforce.addCollectiveVariable("D_Chi1", D_Chi1)
    cvforce.addCollectiveVariable("D_Chi2", D_Chi2)
    cvforce.addCollectiveVariable("F_Chi1", F_Chi1)
    cvforce.addCollectiveVariable("F_Chi2", F_Chi2)
    return cvforce
