#!/usr/bin/env amspython
from scm.plams import *
from ase.optimize import BFGS

def get_atoms():
    # water in a box
    mol = from_smiles('O')  #PLAMS Molecule
    mol.lattice = [[3., 0., 0.,], [0., 3., 0.], [0., 0., 3.]]
    return toASE(mol) # convert PLAMS Molecule to ASE Atoms

def get_settings(properties=True):
    # PLAMS Settings configuring the calculation
    s = Settings() 
    s.input.ams.Task = 'SinglePoint'

    if properties:
        s.input.ams.Properties.Gradients = "Yes"
        s.input.ams.Properties.StressTensor = "Yes"

    #Engine definition
    s.input.ForceField.Type = 'UFF' 

    # run in serial
    s.runscript.nproc = 1 
    return s

def singlepoint():
    settings = get_settings()
    atoms = get_atoms()
    atoms.calc = AMSCalculator(settings=settings, name='SinglePoint')
    print("Singlepint")
    print("Energy (eV):")
    print(atoms.get_potential_energy())
    print("Forces (eV/ang):")
    print(atoms.get_forces())
    print("Stress (eV/ang^3):")
    print(atoms.get_stress())

def geoopt():
    print("Geo opt in normal mode: One results dir is saved for every step")
    settings = get_settings()
    atoms = get_atoms()
    atoms.calc = AMSCalculator(settings=settings, name='GeoOpt')
    dyn = BFGS(atoms)
    dyn.run(fmax=0.05)
    print("Optimized energy (eV):")
    print(atoms.get_potential_energy())

def worker_geoopt():
    print("Geo opt in AMSWorker mode: no output files are saved, minimal overhead")
    settings = get_settings(properties=False) #in AMSWorker mode, do not specify the Properties in the settings.
    atoms = get_atoms()
    atoms.calc = AMSCalculator(settings=settings, name='WorkerGeoOpt', amsworker=True)
    dyn = BFGS(atoms)
    dyn.run(fmax=0.05)
    print("Optimized energy (eV):")
    print(atoms.get_potential_energy())

    
def main():
    init()
    singlepoint()
    geoopt()
    worker_geoopt()
    finish()

if __name__ == '__main__':
    main()

