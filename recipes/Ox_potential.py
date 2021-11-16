import os, sys, csv
from scm.plams import *



class OxidationPotentialCalculator:
	def __init__(self, 
				 job_dir		:str	= os.getcwd(),
				 logfile			 	= sys.stdout):

		self.logfile = logfile #file to send prints to

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

			
		#default settings for pre-optimization using DFTB
		self.DFTB_defaults = Settings()
		self.DFTB_defaults.input.ams.task = 'GeometryOptimization' 
		self.DFTB_defaults.input.ams.Properties.NormalModes           = 'Yes'
		self.DFTB_defaults.input.ams.Properties.PESPointCharacter     = 'Yes'
		self.DFTB_defaults.input.ams.NormalModes.ReScanFreqRange      = '-1000 0'
		self.DFTB_defaults.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = -20
		self.DFTB_defaults.input.DFTB
		self.DFTB_defaults.input.DFTB.Model           = "DFTB3" 
		self.DFTB_defaults.input.DFTB.ResourcesDir    = 'DFTB.org/3ob-3-1'

		#default settings for optimization and singlepoint using DFT
		self.DFT_defaults = Settings()	
		self.DFT_defaults.input.ams.task             = 'GeometryOptimization'
		self.DFT_defaults.input.adf.basis.type       = 'TZ2P'
		self.DFT_defaults.input.adf.basis.core       = 'None'
		self.DFT_defaults.input.adf.xc.hybrid        = 'B3LYP'
		self.DFT_defaults.input.adf.xc.Dispersion    = 'GRIMME3 BJDAMP'
		self.DFT_defaults.input.adf.Relativity.Level = 'None'
		self.DFT_defaults.input.adf.NumericalQuality         = 'Good'
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


	def check_termination_succes(self,
								 result :Results) -> bool:
		return result.readrkf('General', 'termination status', 'ams') == 'NORMAL TERMINATION'

	
	def oxidation_potential(self, 
							molecule			:Molecule, 
							method 				:str		= 'DC',
							solvation_method	:str		= 'COSMO',
							name				:str 		= None,
							settings			:Settings 	= None, 
								) -> dict:
		
		assert method in ['DC', 'TC'], 'Argument "method" must be "DC" or "TC"'
		assert solvation_method in ['COSMO', 'COSMO-RS'], 'Argument "solvation_method" must be "COSMO" or "COSMO-RS"'


		if method == 'DC':
			GO_os = self._geometry_optimization(molecule, state='oxidized', phase='solvent')
			GO_ns = self._geometry_optimization(molecule, state='neutral', phase='solvent')

			self.log(f"Oxidation potential: {GO_os['gibbs_energy'] - GO_ns['gibbs_energy']:.4f}")


	def _geometry_optimization(self, 
				  molecule			:Molecule,
				  state				:str		= 'neutral',
				  phase				:str		= 'vacuum',
				  solvation_method	:str		= 'COSMO',
				  name				:str		= None,
				  settings			:Settings	= None,
				  pre_optimize		:bool		= True
					) -> dict:
		""" Method used to optimize the geometry of molecule using DFT (by default B3LYP)
		Other settings may be supplied using settings which will be soft-updated using self.DFT_defaults
		State specifies whether the molecule is neutral or oxidised
		Phase specifies whether the system is in vacuum or solvated
		"""

		assert state in ['neutral', 'oxidized'], 'argument "state" must be "neutral" or "oxidized"'
		assert phase in ['vacuum', 'solvent'], 'argument "phase" must be "vacuum" or "solvent"'

		#update settings for the calculation
		if settings is None:
			settings = self.DFT_defaults.copy()
		else:
			settings.soft_update(self.DFT_defaults)

		if phase == 'solvent':
			if solvation_method == 'COSMO':
				settings.soft_update(self.COSMO_defaults)
			elif solvation_method == 'COSMO-RS':
				settings.soft_update(self.COSMORS_defaults)

		if state == 'oxidized':
			settings.soft_update(self.oxidized_defaults)

		#set name
		if name is None:
			name = molecule.properties.name

		#summarize job in one string
		job_desc = f'GO_{state}_{phase}'
		if phase == 'solvent':
			job_desc += f'_{solvation_method}'
				
		self.log(f'Starting geometry optimization {name + "_" + job_desc}')
		self.log(f'\tname                 = {name}')
		self.log(f'\tstate                = {state}')
		self.log(f'\tphase                = {phase}')
		if phase == 'solvent':
			self.log(f'\tsolvation_method = {solvation_method}')


		self.log(settings)
		#run the job
		job = AMSJob(molecule 	= molecule,
					 settings 	= settings,
					 name 		= name + '_' + job_desc)

		res = job.run()

		result_dict = {}

		#pull out results
		if self.check_termination_succes(res):
			#set some values
			opt_mol = res.get_main_molecule()
			opt_mol.write(os.path.join(self.geo_opt_dir, name + '_' + job_desc + '.xyz'))

			bond_energy = res.readrkf('Energy', 'Bond Energy', 'adf')
			result_dict['bond_energy'] = Units.convert(bond_energy, 'hartree', 'eV')

			gibbs_energy = res.readrkf('Thermodynamics', 'Gibbs free Energy', 'adf')
			result_dict['gibbs_energy'] = Units.convert(gibbs_energy, 'hartree', 'eV')

			result_dict['geometry'] = res.get_main_molecule()

		else:
			molecule.write(os.path.join(self.not_geo_opt_dir, name + '_' + job_desc + '.xyz'))

		return result_dict


	def log(self, line):
		self.logfile.write(str(line) + '\n')
		print(line)




if __name__ == '__main__':
	#settings
	job_dir = './OxidationPotential'
	mol_file = "./molecules/PDI.xyz"
	job_name = None #if set to None, a name will be generated

	#get commandline options
	import getopt
	opts, left_args = getopt.getopt(sys.argv[1:], 'hm:s:M:N:', ['method=', 'solvation_method=', 'Molecule=', 'Name='])

	method = 'DC'
	solvation_method = 'COSMO'
	for option, argument in opts:
		if option == '-m':
			method = argument.strip()
		if option == '-s':
			solvation_method = argument.strip()
		if option == '-M':
			mol_file = argument.strip()
		if option == '-N':
			job_name = argument.strip()


	if job_name is None:
		job_name = os.path.basename(mol_file).split('.')[0]

	if not os.path.exists(job_dir):
		os.makedirs(job_dir)


	#calculation part
	init(path=job_dir, folder=job_name)

	workdir = config.default_jobmanager.workdir
	logfile = open(f'{workdir}/python.log', 'w')

	ox_potential = OxidationPotentialCalculator(job_dir=workdir, logfile=logfile)
	mol = Molecule(mol_file)
	ox_potential.oxidation_potential(mol, method=method, solvation_method=solvation_method)

	finish()