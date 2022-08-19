#!/usr/bin/env plams

"""

    This example calculates the solubility of aspirin in water and ethanol with
    COSMO-RS.

    This works for arbitrary molecules, i.e. the AMS COSMO-RS library is not
    used but all DFT calculations are performed from scratch.

"""

def main():
    write_initial_xyz_molecules() # comment this line if restarting or running postanalysis
    run() # comment this line if running postanalysis
    #postanalysis('plams_workdir') # comment this line if calling run()

def write_initial_xyz_molecules():
    """ Writes water.xyz, ethanol.xyz, and aspirin.xyz """
    water = from_smiles("O")
    water.write("water.xyz")
    ethanol = from_smiles("CCO", nconfs=10, forcefield='uff')[0] # lowest-energy conformer
    ethanol.write("ethanol.xyz")
    aspirin = from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O", nconfs=100, forcefield='uff')[0] # lowest-energy conformer
    aspirin.write("aspirin.xyz")

def run(temperature=298.15):
    """
    Runs all jobs and prints results.

    Note: the files water.xyz, ethanol.xyz, and aspirin.xyz must exist in the current working directory
    """

    settings_ams = Settings()
    settings_ams.input.adf.symmetry = 'nosym'
    
    settings_crs = Settings()
    settings_crs.input.temperature = temperature
    settings_crs.input.property._h = 'solubility'

    water = Molecule('water.xyz')
    ethanol = Molecule('ethanol.xyz')
    aspirin = Molecule('aspirin.xyz')
    
    solvents = [water, ethanol] 
    solutes = aspirin
    
    crs_dict = run_crs_ams(settings_ams, settings_crs, solvents, solutes)
    for key in crs_dict:
        print("\nSystem:", key) 
        res = crs_dict[key].get_results(settings_crs.input.property._h.upper())
        print('Solubility (g/L solvent): {:.6f}'.format(res['solubility g_per_L_solvent'][1][0]))
        #print(res) #for all results

def postanalysis(workdir='plams_workdir'):
    """ Read and print the results from a previously finished plams_workdir directory """
    print("Loading jobs from {}".format(workdir))
    jobs = load_all(workdir)
    for job in jobs.values():
        if job.name.startswith('CRSJob'):
            res = job.results.get_results('SOLUBILITY')
            solubility_g_L = res['solubility g_per_L_solvent'][1][0]
            print('{} Solubility (g/L solvent): {:.2f}'.format(job.name, solubility_g_L))

if __name__ == '__main__':
    main()

