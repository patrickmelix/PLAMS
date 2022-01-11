import os, sys, csv
from scm.plams import *



class OxidationPotentialCalculator:
	def __init__(self, 
				 logfile			 	= sys.stdout):

		self.logfile = logfile #file to send prints to
			
		#### DEFAULT SETTINGS these can be changed after creating a new OxidationPotentialCalculator object
		#default settings for pre-optimization using DFTB
		self.pre_optimize_defaults = Settings()
		self.pre_optimize_defaults.input.ams.task = 'GeometryOptimization' 
		self.pre_optimize_defaults.input.ams.Properties.NormalModes           = 'Yes'
		self.pre_optimize_defaults.input.ams.Properties.PESPointCharacter     = 'Yes'
		self.pre_optimize_defaults.input.ams.NormalModes.ReScanFreqRange      = '-1000 0'
		self.pre_optimize_defaults.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = -20
		self.pre_optimize_defaults.input.DFTB
		self.pre_optimize_defaults.input.DFTB.Model           = "DFTB3" 
		self.pre_optimize_defaults.input.DFTB.ResourcesDir    = 'DFTB.org/3ob-3-1'

		#DFTB GO settings
		self.DFTB_defaults = Settings()
		self.DFTB_defaults.input.ams.task       = 'GeometryOptimization'
		self.DFTB_defaults.input.ams.Properties.NormalModes       = 'Yes'
		self.DFTB_defaults.input.ams.Properties.PESPointCharacter = 'Yes'
		self.DFTB_defaults.input.ams.NormalModes.ReScanFreqRange  = '-1000 0'
		self.DFTB_defaults.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = -20
		self.DFTB_defaults.input.DFTB
		self.DFTB_defaults.input.DFTB.Model = "GFN1-xTB" 

		#default settings for optimization and singlepoint using DFT
		self.DFT_defaults = Settings()	
		self.DFT_defaults.input.ams.task             = 'GeometryOptimization'
		# self.DFT_defaults.input.adf.basis.type       = 'TZ2P'
		self.DFT_defaults.input.adf.basis.type       = 'DZP'
		self.DFT_defaults.input.adf.basis.core       = 'None'
		self.DFT_defaults.input.adf.xc.hybrid        = 'B3LYP'
		self.DFT_defaults.input.adf.xc.Dispersion    = 'GRIMME3 BJDAMP'
		self.DFT_defaults.input.adf.Relativity.Level = 'None'
		self.DFT_defaults.input.adf.NumericalQuality         = 'Good'
		self.DFT_defaults.input.adf.NumericalQuality         = 'Basic'
		self.DFT_defaults.input.adf.RIHartreeFock.UseMe      = 'Yes'
		self.DFT_defaults.input.adf.RIHartreeFock.Quality    = 'Normal'  
		self.DFT_defaults.input.adf.Symmetry 				 = 'NOSYM'
		self.DFT_defaults.input.ams.UseSymmetry              = 'No'
		self.DFT_defaults.input.ams.properties.NormalModes   = 'Yes'
		self.DFT_defaults.input.ams.Properties.PESPointCharacter     = 'No'
		self.DFT_defaults.input.ams.NormalModes.ReScanFreqRange      = '-1000 0'
		self.DFT_defaults.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = -20

		#default solvent settings for optimization
		self.COSMO_defaults = Settings()
		self.COSMO_defaults.input.adf.Solvation.Solv = "Name=Dichloromethane"

		#default settings for oxidized molecules
		self.oxidized_defaults = Settings()
		self.oxidized_defaults.input.adf.Unrestricted     = 'Yes'
		self.oxidized_defaults.input.adf.SpinPolarization = '1.0'
		self.oxidized_defaults.input.ams.System.Charge    = '1.0'


	def __call__(self, *args, **kwargs):
		return self.oxidation_potential(*args, **kwargs)


	def set_paths(self, 
				  job_dir:str = os.getcwd()
				  ) -> None:

		#default paths
		self.job_dir = job_dir
		if not os.path.exists(job_dir):					os.makedirs(job_dir)
		self.log(f'Job directory:                             {job_dir}')

		self.geo_preopt_dir 	= os.path.join(job_dir,'Geometries', 'Preoptimization', 'Preopt_Geo')
		if not os.path.exists(self.geo_preopt_dir): 	os.makedirs(self.geo_preopt_dir)
		self.log(f'Preoptimization geometry directory:        {self.geo_preopt_dir}')

		self.not_geo_preopt_dir 	= os.path.join(job_dir,'Geometries', 'Preoptimization', 'Not_Preopt_Geo')
		if not os.path.exists(self.not_geo_preopt_dir): os.makedirs(self.not_geo_preopt_dir)
		self.log(f'Failed preoptimization geometry directory: {self.not_geo_preopt_dir}')

		self.geo_opt_dir 	= os.path.join(job_dir,'Geometries', "Optimized_Geo")
		if not os.path.exists(self.geo_opt_dir): 		os.makedirs(self.geo_opt_dir)
		self.log(f'DFT optimized geometry directory:          {self.geo_opt_dir}')

		self.not_geo_opt_dir = os.path.join(job_dir,'Geometries', "Not_Optimized_Geo")
		if not os.path.exists(self.not_geo_opt_dir): 	os.makedirs(self.not_geo_opt_dir)
		self.log(f'Failed DFT optimized geometry directory:   {self.not_geo_opt_dir}')


	def check_termination_succes(self,
								 result :Results) -> bool:
		term = result.readrkf('General', 'termination status', 'ams')
		if term == 'NORMAL TERMINATION':
			return True
		elif 'NORMAL TERMINATION' in term:
			return 'WARNING'
		return 

	
	def oxidation_potential(self, 
							molecule			:Molecule, 
							method 				:str		= 'DC',
							name				:str 		= None,
							settings			:Settings 	= None, 
							pre_optimize 		:bool		= True,
							job_dir				:str 		= os.getcwd(),
								) -> dict:
		

		assert method in ['DC', 'TC'], 'Argument "method" must be "DC" or "TC"'

		

		#set name
		if name is None:
			name = molecule.properties.name

		self.log(f'========================================================================')
		self.log(f'Starting oxidation potential calculation for molecule {name}:\n')

		#set paths for the jobs:
		self.log('Setting paths ...')
		self.set_paths(job_dir)

		self.log('Initial coordinates:')
		self.log(molecule)
		self.log('Settings:')
		self.log(f'\tMethod:            {method}')

		if settings is not None:
			self.log(f'Extra settings:')
			self.log(settings)

		if pre_optimize:
			self.log('Pre-optimizing the molecule ...')
			molecule = self._pre_optimize(molecule, name=name)
			self.log('New coordinates:')
			self.log(molecule)

		#get on with the actual calculations
		if method == 'DC':
			GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='solvent', settings=settings)
			GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral',  phase='solvent', settings=settings)

			oxpot = GO_os['gibbs_energy'] - GO_ns['gibbs_energy']
			
		if method == 'TC':
			GO_nv    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='vacuum', settings=settings)
			SP_ns_nv = self._calculation_step(GO_nv['geometry'], task='SinglePoint', name=name, state='neutral', phase='solvent', settings=settings)
			GO_ns    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='neutral', phase='solvent', settings=settings)
			SP_nv_ns = self._calculation_step(GO_ns['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum', settings=settings)

			GO_ov    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='vacuum', settings=settings)
			SP_os_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state='oxidized', phase='solvent', settings=settings)
			SP_nv_ov = self._calculation_step(GO_ov['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum', settings=settings)
			GO_os    = self._calculation_step(molecule, task='GeometryOptimization', name=name, state='oxidized', phase='solvent', settings=settings)
			SP_nv_os = self._calculation_step(GO_os['geometry'], task='SinglePoint', name=name, state='neutral', phase='vacuum', settings=settings)

			oxidized_part = GO_nv['gibbs_energy'] + SP_ns_nv['gibbs_energy'] + (GO_nv['bond_energy'] - SP_nv_os['bond_energy'])
			neutral_part  = GO_ov['gibbs_energy'] + SP_os_ov['gibbs_energy'] + (SP_nv_ov['bond_energy'] - SP_nv_os['bond_energy'])
			oxpot = oxidized_part - neutral_part

		self.log(f"Oxidation potential: {oxpot:.4f} eV")
		return oxpot


	def _pre_optimize(self,
				  molecule 	:Molecule,
				  name 		:str 		= None,
					 ) -> Molecule:
		"""Method used to pre-optimize the molecule using DFTB3
		"""
		job = AMSJob(molecule 	= molecule,
					 settings 	= self.pre_optimize_defaults,
					 name 		= name + '_preoptimization')
		res = job.run()

		#if the optimization is succesful we return the new molecule
		if self.check_termination_succes(res):
			molecule = res.get_main_molecule()
			molecule.write(os.path.join(self.geo_preopt_dir, name + '.xyz'))
		else:
			molecule.write(os.path.join(self.not_geo_preopt_dir, name + '.xyz'))
		
		return molecule


	def _calculation_step(self, 
				  molecule			:Molecule,
				  task				:str 		= 'GeometryOptimization',
				  use_dftb			:bool		= False,
				  state				:str		= 'neutral',
				  phase				:str		= 'vacuum',
				  solvation_method	:str		= 'COSMO',
				  name				:str		= None,
				  settings			:Settings	= None,
					) -> dict:
		""" Method used to optimize the geometry of molecule using DFT (by default B3LYP)
		Other settings may be supplied using settings which will be soft-updated using self.DFT_defaults
		State specifies whether the molecule is neutral or oxidised
		Phase specifies whether the system is in vacuum or solvated
		if use_dftb, the system will be optimised using DFTB (by default GFN1-xTB) instead of DFT
		"""

		assert state in ['neutral', 'oxidized'], 'argument "state" must be "neutral" or "oxidized"'
		assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'
		assert task  in ['GeometryOptimization', 'SinglePoint'], 'argument "task" must be "GeometryOptimization" or "SinglePoint"'

		#update settings for the calculation either DFT or DFTB
		if use_dftb:
			defaults = self.DFTB_defaults.copy()
		else:
			defaults = self.DFT_defaults.copy()

		#set appropriate task
		defaults.input.ams.task = task

		#update the provided settings if they were given
		if settings is not None:
			settings.soft_update(defaults)
		else:
			settings = defaults

		#handle solvent models
		if phase == 'solvent':
			if solvation_method == 'COSMO':
				settings.soft_update(self.COSMO_defaults)
			elif solvation_method == 'COSMO-RS':
				settings.soft_update(self.COSMORS_defaults)

		#handle state, if neutral the settings are already correct
		if state == 'oxidized':
			settings.soft_update(self.oxidized_defaults)


		#summarize job in one string
		job_desc = f'GO_{state}_{phase}'
		if phase == 'solvent':
			job_desc += f'_{solvation_method}'
		if use_dftb:
			job_desc += '_DFTB'


		self.log(f'Starting calculation {name + "_" + job_desc}')
		self.log(f'\tname                 = {name}')
		self.log(f'\ttask                 = {task}')
		self.log(f'\tuse_dftb             = {use_dftb}')
		self.log(f'\tstate                = {state}')
		self.log(f'\tphase                = {phase}')
		if phase == 'solvent':
			self.log(f'\tsolvation_method     = {solvation_method}')

		# self.log(settings)

		#run the job
		job = AMSJob(molecule 	= molecule,
					 settings 	= settings,
					 name 		= name + '_' + job_desc)

		res = job.run()

		result_dict = {}

		#pull out results
		if self.check_termination_succes(res):
			self.log(f'\tSuccessfull          = {self.check_termination_succes(res)}') #True or WARNING
			#set some values
			bond_energy = res.readrkf('Energy', 'Bond Energy', 'adf')
			result_dict['bond_energy'] = Units.convert(bond_energy, 'hartree', 'eV')

			gibbs_energy = res.readrkf('Thermodynamics', 'Gibbs free Energy', 'adf')
			result_dict['gibbs_energy'] = Units.convert(gibbs_energy, 'hartree', 'eV')
			
			self.log(f'\tResults:')
			self.log(f'\t\tBond Energy    = {result_dict["bond_energy"]:.4f} eV')
			self.log(f'\t\tGibbs Energy   = {result_dict["gibbs_energy"]:.4f} eV')

			if task == 'GeometryOptimization':
				opt_mol = res.get_main_molecule()
				opt_mol.write(os.path.join(self.geo_opt_dir, name + '_' + job_desc + '.xyz'))
				result_dict['geometry'] = res.get_main_molecule()
				self.log(f'\t\tFinal geometry:')
				self.log(result_dict['geometry'])

		else:
			self.log(f'\tSuccessfull          = False')
			if task == 'GeometryOptimization':
				molecule.write(os.path.join(self.not_geo_opt_dir, name + '_' + job_desc + '.xyz'))

		return result_dict


	def log(self, line):
		if not self.logfile is None: 
			self.logfile.write(str(line) + '\n')
			self.logfile.flush()
		print(line)




if __name__ == '__main__':
	#settings
	job_dir = './OxidationPotential_test'
	for mol_file in ["./molecules/NDI.xyz", "./molecules/NDI44.xyz", "./molecules/NDI55.xyz", "./molecules/NDI54.xyz", "./molecules/PDI.xyz"]:
		job_name = None #if set to None, a name will be generated
		method = 'TC'

		if job_name is None:
			job_name = os.path.basename(mol_file).split('.')[0]

		if not os.path.exists(job_dir):
			os.makedirs(job_dir)


		#calculation part
		init(path=job_dir, folder=job_name)

		workdir = config.default_jobmanager.workdir
		logfile = open(f'{workdir}/python.log', 'w')

		ox_potential_calc = OxidationPotentialCalculator(logfile=logfile)
		mol = Molecule(mol_file)
		ox_potential_calc.oxidation_potential(mol, job_dir=workdir, method=method)

		finish()