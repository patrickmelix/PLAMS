#!/usr/bin/env amspython
from scm.plams import *
from ase.calculators.emt import EMT
import os

def get_calculator(params_path, **kwargs):
    """ 
    Return an ASE calculator for use with the MLPotential%Backend = ASE. 
    
    The function must be called "get_calculator" 

    ``params_path`` is a string corresponding to the MLPotential%ParameterFile input (ignored in this case)
    
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
    s.input.MLPotential.backend = 'ASE'
    s.input.MLPotential.ASECalculatorFile = os.path.abspath(__file__)   # get it from the current file's get_calculator() method
    #s.input.MLPotential.ParameterFile = os.path.abspath('some-file.txt') # pass a path to a parameter file as the params_path argument in get_calculator(), not done in this example

    job = AMSJob(settings=s, molecule=mol, name='ams_with_custom_ase_calculator')
    job.run()

    energy = job.results.get_energy(unit='eV')
    print(f"AMS with custom ASE calculator (MLPotential Backend%ASE), EMT potential: final energy {energy:.3f} eV")

    finish()

if __name__ == '__main__':
    main()

