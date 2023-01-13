#!/usr/bin/env amspython
from scm.plams import *
from ase.calculators.emt import EMT
import os

def get_calculator():
    """ 
    Return an ASE calculator for use with Engine ASE. 
    
    The function must be called "get_calculator" 

    """
    return EMT()

def main():
    init()

    mol = from_smiles('O')  
    mol.lattice = [[3., 0., 0.,], [0., 3., 0.], [0., 0., 3.]]

    s = Settings()
    s.runscript.nproc = 1
    s.input.ams.task = 'GeometryOptimization'
    s.input.ams.GeometryOptimization.Convergence.Gradients = 0.01 # hartree/ang
    s.input.ASE.File = os.path.abspath(__file__)   # get it from the current file's get_calculator() method
    #s.input.ASE.Arguments._1 = "my_argument = my_value" to call get_calculator(my_argument = my_value)

    job = AMSJob(settings=s, molecule=mol, name='ams_with_custom_ase_calculator')
    job.run()

    energy = job.results.get_energy(unit='eV')
    print(f"AMS with custom ASE calculator (Engine ASE), EMT potential: final energy {energy:.3f} eV")

    finish()

if __name__ == '__main__':
    main()

