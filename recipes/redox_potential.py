import os, sys
from scm.plams import *
import redox_defaults as defaults



### ==== SETTINGS ==== ###

def get_settings(task                :str        = 'GeometryOptimization',
                 state               :str        = 'neutral',
                 phase               :str        = 'vacuum',
                 frequencies         :bool       = True,
                 use_dftb            :bool       = False,
                 use_COSMORS         :bool       = False,
                    ) -> Settings:
    '''
    Method that generates settings for jobs based on settings provided
    '''

    assert state in ['neutral', 'oxidation', 'reduction'], 'argument "state" must be "neutral", "oxidation" or "reduction"'
    assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'
    assert task  in ['GeometryOptimization', 'SinglePoint'], 'argument "task" must be "GeometryOptimization", "SinglePoint" or "COSMO"'

    if use_COSMORS:
        sett = defaults.DFT()

        solvation_block = {
            'surf': 'Delley',
            'solv': 'name=CRS cav0=0.0 cav1=0.0',
            'charged': 'method=Conj corr',
            'c-mat': 'Exact',
            'scf': 'Var All',
            'radii': {
                'H': 1.30,
                'C': 2.00,
                'N': 1.83,
                'O': 1.72,
                'F': 1.72,
                'Si': 2.48,
                'P': 2.13,
                'S': 2.16,
                'Cl': 2.05,
                'Br': 2.16,
                'I': 2.32
            }}

        sett.input.adf.solvation = solvation_block
    else:
        if use_dftb:
            sett = defaults.DFTB()
        else:
            sett = defaults.DFT()

        #load cosmo solvent
        if phase == 'solvent' and not use_COSMORS:
            sett.soft_update(defaults.COSMO())

    sett.input.ams.task = task

    if frequencies:
        sett.soft_update(defaults.frequencies())

    #handle state, if neutral the settings are already correct
    if state == 'oxidation':
        if use_dftb:
            sett.soft_update(defaults.DFTB_oxidized())
        else:
            sett.soft_update(defaults.DFT_oxidized())
    if state == 'reduction':
        if use_dftb:
            sett.soft_update(defaults.DFTB_reduced())
        else:
            sett.soft_update(defaults.DFT_reduced())

    return sett



### ==== UTILITY FUNCTIONS ==== ###

def check_termination_succes(result :Results) -> bool:
    term = result.readrkf('General', 'termination status', 'ams')
    if term == 'NORMAL TERMINATION':
        return True
    elif 'NORMAL TERMINATION' in term:
        return 'WARNING'
    return 


### ==== CALCULATIONS ==== ###

def COSMORS_property(solvent_path        :str, 
                     solute_path         :str,
                     name                :str,
                     temperature         :float      = 298.15
                ) -> float:
    """This method runs a COSMORS property job to obtain the activity coefficient
    which will also calculate G solute which we need to calculate the redox 
    potential.
    """

    sett = Settings()
    sett.input.property._h = 'ACTIVITYCOEF'
    compounds = [Settings(), Settings()]
    compounds[0]._h    = solvent_path;    compounds[1]._h    = solute_path
    compounds[0].frac1 = 1;               compounds[1].frac1 = 0
    
    sett.input.temperature = str(temperature)
    sett.input.compound = compounds

    res = CRSJob(settings=sett, name=name).run().get_results()
    if res:
        #convert the Gibbs energy to hartree (COSMORS gives in kcal/mol)
        return float(Units.convert(res["G solute"][1], 'kcal/mol', 'hartree'))
    else: return False
            

def calculation_step(molecule            :Molecule,
                     task                :str        = 'GeometryOptimization',
                     state               :str        = 'neutral',
                     phase               :str        = 'vacuum',
                     frequencies         :bool       = True,
                     use_dftb            :bool       = False,
                     use_COSMORS         :bool       = False,
                     name                :str        = None,
                     solvent_path        :str        = None,
                    ) -> dict:
    """ Method used to optimize the geometry of molecule using DFT (by default B3LYP)
    Other settings may be supplied using settings which will be soft-updated using DFT_defaults
    State specifies whether the molecule is neutral, oxidised or reduced
    Phase specifies whether the system is in vacuum or solvated
    if use_dftb, the system will be optimised using DFTB (by default GFN1-xTB) instead of DFT
    """

    assert state in ['neutral', 'oxidation', 'reduction'], 'argument "state" must be "neutral", "oxidation", "reduction"'
    assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'
    assert task  in ['GeometryOptimization', 'SinglePoint'], 'argument "task" must be "GeometryOptimization" or "SinglePoint"'

    if use_COSMORS: phase = 'solvent'

    settings = get_settings(task=task, 
                            use_dftb=use_dftb, 
                            use_COSMORS=use_COSMORS,
                            state=state, 
                            phase=phase, 
                            frequencies=frequencies)

    #summarize job in one string
    task_abbrev = {"GeometryOptimization":"GO", "SinglePoint":"SP"}[task]
    job_desc = f'{task_abbrev}_{state}_{phase}'
    if use_COSMORS:
        job_desc += '_COSMO-RS'
    if use_dftb:
        job_desc += '_DFTB'

    print(f'\nStarting calculation {name + "_" + job_desc}')
    print(f'\ttask                 = {task}')
    print(f'\tuse_dftb             = {use_dftb}')
    print(f'\tuse_COSMORS          = {use_COSMORS}')
    print(f'\tfrequencies          = {frequencies}')
    print(f'\tstate                = {state}')
    print(f'\tphase                = {phase}')

    #run the job
    job = AMSJob(molecule   = molecule,
                 settings   = settings,
                 name       = name + '_' + job_desc)
    res = job.run()

    result_dict = {}
    #pull out results
    if check_termination_succes(res):
        print(f'\tSuccessfull          = {check_termination_succes(res)}') #True or WARNING
        #set some default values
        bond_energy = None 
        gibbs_energy = None

        #If we are doing COSMO calculations then we need to run an additional job to obtain the activity coefficient
        #when calculating the activity coefficient, the G solute is also calculated.
        if use_COSMORS:
            resfile = KFFile(res['adf.rkf'])
            cosmo_data = resfile.read_section('COSMO')
            coskf = KFFile(os.path.join(job.path, job.name + '.coskf'))
            for k, v in cosmo_data.items():
                coskf.write('COSMO', k, v)
            res.collect()
            bond_energy = res.readrkf('AMSResults', 'Energy', 'adf')
            gibbs_energy = COSMORS_property(solvent_path, os.path.join(job.path, job.name + '.coskf'), job.name + '_ACTIVITYCOEF')
        #if we dont use COSMO-RS we can just extract the Gibbs and bonding energies from the regular job
        else:
            if use_dftb:
                bond_energy = res.readrkf('AMSResults', 'Energy', 'dftb')
                if frequencies: 
                    gibbs_energy = res.readrkf('Thermodynamics', 'Gibbs free Energy', 'dftb')
            else:
                bond_energy = res.readrkf('Energy', 'Bond Energy', 'adf')
                if frequencies: 
                    gibbs_energy = res.readrkf('Thermodynamics', 'Gibbs free Energy', 'adf')
        
        print(f'\tResults:')
        if not bond_energy is None: 
            result_dict['bond_energy'] = Units.convert(bond_energy, 'hartree', 'eV')
            print(f'\t\tBond Energy  = {result_dict["bond_energy"]:.4f} eV')
        if not gibbs_energy is None: 
            result_dict['gibbs_energy'] = Units.convert(gibbs_energy, 'hartree', 'eV')
            print(f'\t\tGibbs Energy = {result_dict["gibbs_energy"]:.4f} eV')

        #extract also optimised molecule
        if task == 'GeometryOptimization':
            result_dict['geometry'] = res.get_main_molecule()

        #and if the phase is solvent we also need the solvation gibbs energy change
        if phase == 'solvent':
            dG_solvation = res.readrkf('Energy','Solvation Energy (el)','adf') + res.readrkf('Energy','Solvation Energy (cd)','adf')
            result_dict['dG_solvation'] = Units.convert(dG_solvation, 'hartree', 'eV')
            print(f'\t\tdG_solvation = {result_dict["dG_solvation"]:.4f} eV')

    else:
        print(f'\tSuccessfull          = False')

    return result_dict


### ==== MAIN FUNCTION ==== ###

def redox_potential(molecule            :Molecule, 
                    mode                :str,
                    method              :str        = 'screening',
                    name                :str        = None,
                    COSMORS_solvent_path:str        = None,
                        ) -> float:

    assert mode in ['oxidation', 'reduction']
    assert method in ['DC', 'TC-COSMO', 'TC-COSMO-RS', 'screening'], 'Argument "method" must be "DC", "TC-COSMO", "TC-COSMO-RS" or "screening"'

    #set name
    if name is None:
        name = molecule.properties.name

    print(f'========================================================================')
    print(f'Starting redox potential calculation for molecule {name}:\n')

    print('\nInitial coordinates:')
    print(molecule)
    print('Settings:')
    print(f'\tName:              {name}')
    print(f'\tMode:              {mode}')
    print(f'\tMethod:            {method}')

    if mode == 'oxidation':
        Gelectron = -0.0375
    elif mode == 'reduction':
        Gelectron = 0.0375
    #get on with the actual calculations
    if method == 'DC':
        GO_os    = calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='solvent')
        GO_ns    = calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral',  phase='solvent')

        redoxpot = GO_os['gibbs_energy'] - GO_ns['gibbs_energy'] + Gelectron
        
    elif method == 'TC-COSMO':
        GO_nv    = calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum')
        SP_nv_ns = calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', phase='solvent', frequencies=False)
        GO_ns    = calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='solvent', frequencies=False)
        SP_nv_nv = calculation_step(GO_ns['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum',  frequencies=False)

        GO_ov    = calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='vacuum')
        SP_ov_os = calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state=mode, phase='solvent', frequencies=False)
        GO_os    = calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='solvent', frequencies=False)
        SP_ov_ov = calculation_step(GO_os['geometry'], task='SinglePoint', name=name, state=mode, phase='vacuum',  frequencies=False)

        redox_part = GO_ov['gibbs_energy'] + SP_ov_os['dG_solvation'] + (SP_ov_ov['bond_energy'] - GO_ov['bond_energy'])
        neutral_part  = GO_nv['gibbs_energy'] + SP_nv_ns['dG_solvation'] + (SP_nv_nv['bond_energy'] - GO_nv['bond_energy'])
        redoxpot = redox_part - neutral_part + Gelectron

    elif method == 'TC-COSMO-RS':
        assert os.path.exists(COSMORS_solvent_path), f'Solvent database {COSMORS_solvent_path} does not exist'
        GO_nv    = calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum',   frequencies=True)
        SP_nv_ns = calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', use_COSMORS=True, frequencies=False, solvent_path=COSMORS_solvent_path)
        GO_ns    = calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', use_COSMORS=True, frequencies=False, solvent_path=COSMORS_solvent_path)
        SP_ns_nv = calculation_step(GO_ns['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum',   frequencies=False)

        GO_ov    = calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='vacuum', frequencies=True)
        SP_ov_os = calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state=mode, use_COSMORS=True, frequencies=False, solvent_path=COSMORS_solvent_path)
        GO_os    = calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, use_COSMORS=True, frequencies=False, solvent_path=COSMORS_solvent_path)
        SP_os_ov = calculation_step(GO_os['geometry'], task='SinglePoint', name=name, state=mode, phase='vacuum',   frequencies=False)

        redox_part = GO_ov['gibbs_energy'] + SP_ov_os['dG_solvation'] + (SP_os_ov['bond_energy'] - GO_ov['bond_energy'])
        neutral_part  = GO_nv['gibbs_energy'] + SP_nv_ns['dG_solvation'] + (SP_ns_nv['bond_energy'] - GO_nv['bond_energy'])
        redoxpot = redox_part - neutral_part + Gelectron

    elif method == 'screening':
        assert os.path.exists(COSMORS_solvent_path), f'Solvent database {COSMORS_solvent_path} does not exist'
        GO_nv    = calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum', use_dftb=True, frequencies=False)
        COSMO_nv = calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', use_COSMORS=True, frequencies=False, solvent_path=COSMORS_solvent_path)

        GO_ov    = calculation_step(molecule, task='GeometryOptimization', name=name, state=mode, phase='vacuum', use_dftb=True, frequencies=False)
        COSMO_ov = calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state=mode, use_COSMORS=True, frequencies=False, solvent_path=COSMORS_solvent_path)

        redoxpot = COSMO_ov['gibbs_energy'] - COSMO_nv['gibbs_energy'] + Gelectron

    print(f"\nOxidation potential: {redoxpot:.4f} eV")
    return redoxpot






if __name__ == '__main__':
    mode = 'oxidation'
    job_dir = './Test'
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    COSMORS_solvent_path = os.path.abspath('Dichloromethane.coskf')

    results = {}
    for mol_file in ['./molecules/methane.xyz']:
        mol_name = os.path.basename(mol_file).split('.')[0]
        results[mol_name] = {}
        for method in ['screening', 'TC-COSMO', 'TC-COSMO-RS', 'DC']:
            job_name = None #if set to None, a name will be generated
            
            if job_name is None:
                job_name = mol_name + '_' + method

            #calculation part
            init(path=job_dir, folder=job_name)

            workdir = config.default_jobmanager.workdir
            logfile = open(f'{workdir}/{job_name}_python.log', 'w')

            mol = Molecule(mol_file)
            redoxpot = redox_potential(mol, mode, method=method, COSMORS_solvent_path=COSMORS_solvent_path)
            results[mol_name][method] = redoxpot
            finish()

    print('RedOx Potentials:')
    name_len = max(len('System'), max(len(n) for n in results))
    methods = list(results[list(results.keys())[0]])
    method_lens = [max(9, len(m)) for m in methods]
    print(f'{"System".ljust(name_len)} | {" | ".join([m.ljust(l) for m, l in zip(methods, method_lens)])}')
    for name, res in results.items():
        s = f'{name.ljust(name_len)} | {" | ".join([(str(round(res[m],3)) + " eV").rjust(l) for m, l in zip(methods, method_lens)])}'
        print(s)

    print('\nCalculations coomplete!\a')