from collections import OrderedDict
from ..core.functions import add_to_instance
from ..core.basejob import MultiJob
from ..core.results import Results
from ..core.settings import Settings
from ..mol.molecule import Molecule
from ..mol.atom import Atom
from ..interfaces.adfsuite.ams import AMSJob
from ..tools.units import Units
from ..interfaces.molecule.packmol import packmol_mixture
import numpy as np
from scipy.optimize import curve_fit

__all__ = ['GreenKuboViscosityJob', 'GreenKuboViscosityResults']

class GreenKuboViscosityResults(Results):
    """Results class for GreenKuboViscosityJob
    """
    def get_viscosity(self, max_dt_fs=50000, reuse=False):
        """
        Returns the viscosity in mPa*s, and writes 
        
        * viscosity.txt: viscosity in mPa*s
        * green_kubo_viscosity.txt: viscosity vs. time (curve that should converge to the viscosity as t goes to inifinity)
        * fit.txt : A fit of the form A*(1-exp(-x/tau)) to green_kubo_viscosity.txt; A is the final viscosity.
        * moving_avg_visc.txt : viscosity vs. time (moving average of green_kubo_viscosity.txt)
        * info.txt : Contains the upper time-limit of the fit and the viscosity

        The upper limit for the fit is determined by when the moving average starts to decrease.

        max_dt_fs : float
            Maximum correlation time in femtoseconds

        reuse : bool
            Will return the value of viscosity.txt if it exists (and not recalculate the autocorrelation function/integral)
        """

        def f(x, A, tau):
            return A*(1-np.exp(-x/tau)) # A is the final converged viscosity

        if reuse and os.path.exists(self.job.path+'/viscosity.txt'):
            return float(np.loadtxt(self.job.path+'/viscosity.txt'))

        job = self.job.children['production']
        assert job.ok(), "Production job did not finish correctly."

        pressuretensor = np.load(job.path+'/pressuretensor.npy')
        t, visc = job.results.get_green_kubo_viscosity(max_dt_fs=max_dt_fs, pressuretensor=pressuretensor)

        A = np.stack((t, visc), axis=1)
        np.savetxt(self.job.path+"/green_kubo_viscosity.txt", A, header="Time(fs) Viscosity(mPa*s)")

        window = min(1000, int(t[-1]/3))
        window = max(window, 1) #window must be greater than 0
        moving_avg_t = np.convolve(t, np.ones(window)/window, mode='valid')
        moving_avg_visc = np.convolve(visc, np.ones(window)/window, mode='valid')

        flat_region_from = min(1000, int(t[-1]/10)) # femtoseconds
        flat_region_from_indices = moving_avg_t > flat_region_from
        allowed_t = moving_avg_t[flat_region_from_indices]
        allowed_visc = moving_avg_visc[flat_region_from_indices]
        d_moving_avg_visc = np.diff(allowed_visc)
        zero_indices = np.argwhere(np.diff(np.sign(d_moving_avg_visc))).flatten()
        if len(zero_indices) == 0:
            max_index = len(allowed_visc)-1
        else:
            max_index = zero_indices[0]

        fit_until_t = allowed_t[max_index]

        fit_x_indices = t < fit_until_t
        fit_x = t[fit_x_indices]
        fit_y = visc[fit_x_indices]

        popt, _ = curve_fit(f, fit_x, fit_y, p0=(1.0, 1000))
        prediction = f(t, popt[0], popt[1])
        with open(self.job.path+'/viscosity.txt', 'w') as f:
            f.write("{}\n".format(popt[0]))

        with open(self.job.path+'/info.txt', 'w') as f:
            f.write("FITTED UNTIL t = {} fs \n".format(fit_until_t))
            f.write("VISCOSITY = {} mPa*s\n".format(popt[0]))

        A = np.stack((t, prediction), axis=1)
        np.savetxt(self.job.path+"/fit.txt", A, header="Time(fs) Viscosity(mPa*s)")

        A = np.stack((moving_avg_t, moving_avg_visc), axis=1)
        np.savetxt(self.job.path+"/moving_avg_visc.txt", A, header="Time(fs) Viscosity(mPa*s)")

        return popt[0]

class GreenKuboViscosityJob(MultiJob):
    """A class for calculating the Green-Kubo viscosity of a liquid
    """

    _result_type = GreenKuboViscosityResults

    def _default_engine_settings(self):
        s = Settings()
        s.input.ForceField.Type = 'GAFF'
        s.input.ForceField.AnteChamberIntegration = 'Yes'
        return s

    def _default_runscript_settings(self):
        s = Settings()
        s.runscript.nproc = 1
        return s

    def _create_scan_density_job(self, initial_molecule):
        name = 'scan_density'

        from_density = initial_molecule.get_density() * 1e-3
        orig_length = initial_molecule.lattice[0][0]
        new_length = orig_length * (from_density / self.scan_to_density) ** 0.3333

        s = Settings()
        s.input.ams.task = 'MolecularDynamics'
        s.input.ams.MolecularDynamics.Nsteps = self.nsteps[name]
        s.input.ams.MolecularDynamics.TimeStep = self.timestep
        s.input.ams.MolecularDynamics.Deformation.TargetLength = f'{new_length} {new_length} {new_length}'
        s.input.ams.MolecularDynamics.Deformation.StartStep = min(1000, self.nsteps[name] // 2)
        s.input.ams.MolecularDynamics.Thermostat.Temperature = self.temperature
        s.input.ams.MolecularDynamics.Thermostat.Tau = 10
        s.input.ams.MolecularDynamics.Thermostat.Type = 'Berendsen'
        s.input.ams.MolecularDynamics.InitialVelocities.Temperature = self.temperature
        s.input.ams.MolecularDynamics.Checkpoint.Frequency = 10000000
        s.input.ams.MolecularDynamics.Trajectory.WriteVelocities = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteBonds = 'Yes'
        s.input.ams.MolecularDynamics.Trajectory.WriteMolecules = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteCharges = 'No'

        s += self.engine_settings + self.runscript_settings

        job = AMSJob(settings=s, molecule=initial_molecule, name=name)
        @add_to_instance(job)
        def get_lowest_energy_molecule(self):
            energies = self.results.get_history_property('Energy')
            minindex = np.argmin(energies) + 1
            return self.results.get_history_molecule(minindex)

        self.children[name] = job

        return self.children[name]

    def _create_nvt_pre_eq_job(self, scan_density_job):
        name = 'nvt_pre_eq'
        s = Settings()
        s.input.ams.task = 'MolecularDynamics'
        s.input.ams.MolecularDynamics.Nsteps = self.nsteps[name]
        s.input.ams.MolecularDynamics.TimeStep = self.timestep
        s.input.ams.MolecularDynamics.Thermostat.Tau = 10
        s.input.ams.MolecularDynamics.Thermostat.Type = 'Berendsen'
        s.input.ams.MolecularDynamics.Thermostat.Temperature = self.temperature
        s.input.ams.MolecularDynamics.InitialVelocities.Temperature = self.temperature
        s.input.ams.MolecularDynamics.Checkpoint.Frequency = 10000000
        s.input.ams.MolecularDynamics.Trajectory.WriteVelocities = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteBonds = 'Yes'
        s.input.ams.MolecularDynamics.Trajectory.WriteMolecules = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteCharges = 'No'
        s += self.engine_settings + self.runscript_settings

        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            self.molecule = scan_density_job.get_lowest_energy_molecule()

        self.children[name] = job
        return self.children[name]

    def _create_npt_job(self, nvt_pre_eq_job):
        name = 'npt'
        s = Settings()
        s.input.ams.Task = 'MolecularDynamics'
        s.input.ams.MolecularDynamics.Nsteps = self.nsteps[name]
        s.input.ams.MolecularDynamics.TimeStep = self.timestep
        s.input.ams.MolecularDynamics.Thermostat.Tau = 100
        s.input.ams.MolecularDynamics.Thermostat.Type = 'NHC'
        s.input.ams.MolecularDynamics.Thermostat.Temperature = self.temperature
        s.input.ams.MolecularDynamics.Barostat.Type = 'MTK'
        s.input.ams.MolecularDynamics.Barostat.Pressure = '1.0 [bar]'
        s.input.ams.MolecularDynamics.Barostat.Tau = 1000
        s.input.ams.MolecularDynamics.Barostat.Equal = 'XYZ'
        s.input.ams.MolecularDynamics.Trajectory.WriteVelocities = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteBonds = 'Yes'
        s.input.ams.MolecularDynamics.Trajectory.WriteMolecules = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteCharges = 'No'
        s.input.ams.MolecularDynamics.Checkpoint.Frequency = 10000000
        s += self.engine_settings + self.runscript_settings

        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            self.settings.input.ams.MolecularDynamics.InitialVelocities.Type = 'FromFile'
            self.settings.input.ams.MolecularDynamics.InitialVelocities.File = nvt_pre_eq_job.results.rkfpath()
            self.molecule = nvt_pre_eq_job.results.get_main_molecule()

        @add_to_instance(job)
        def get_equilibrated_molecule(self):
            densities = self.results.get_history_property('Density', 'MDHistory')
            analyze_from = (len(densities)*2) // 3

            # take structure closest to the target density
            avg_density = np.mean(densities[analyze_from:]) # amu/bohr^3
            delta = np.array(densities[analyze_from:]) - avg_density
            delta = np.abs(delta)
            min_index = np.argmin(delta)
            min_index += analyze_from
            mol = self.results.get_history_molecule(min_index+1)
            return mol

        self.children[name] = job
        return self.children[name]

        
    def _create_nvt_eq_job(self, npt_job):
        name = 'nvt_eq'
        s = Settings()
        s.input.ams.task = 'MolecularDynamics'
        s.input.ams.MolecularDynamics.Nsteps = self.nsteps[name]
        s.input.ams.MolecularDynamics.TimeStep = self.timestep
        s.input.ams.MolecularDynamics.Thermostat.Tau = 100
        s.input.ams.MolecularDynamics.Thermostat.Type = 'NHC'
        s.input.ams.MolecularDynamics.Thermostat.Temperature = self.temperature
        s.input.ams.MolecularDynamics.InitialVelocities.Temperature = self.temperature
        s.input.ams.MolecularDynamics.Trajectory.WriteVelocities = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteBonds = 'Yes'
        s.input.ams.MolecularDynamics.Trajectory.WriteMolecules = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteCharges = 'No'
        s.input.ams.MolecularDynamics.Checkpoint.Frequency = 10000000
        s += self.engine_settings + self.runscript_settings

        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            self.molecule = npt_job.get_equilibrated_molecule()

        self.children[name] = job
        return self.children[name]

    def _create_production_job(self, nvt_eq_job):
        name = 'production'
        s = Settings()
        s.input.ams.task = 'MolecularDynamics'
        s.input.ams.MolecularDynamics.Nsteps = self.nsteps[name]
        s.input.ams.MolecularDynamics.TimeStep = self.timestep
        s.input.ams.MolecularDynamics.Thermostat.Tau = 100
        s.input.ams.MolecularDynamics.Thermostat.Type = 'NHC'
        s.input.ams.MolecularDynamics.Thermostat.Temperature = self.temperature
        s.input.ams.MolecularDynamics.Trajectory.WriteVelocities = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteBonds = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteMolecules = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteCharges = 'No'
        s.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 5
        s.input.ams.MolecularDynamics.CalcPressure = 'Yes'
        s.input.ams.MolecularDynamics.Checkpoint.Frequency = 10000000
        s += self.engine_settings + self.runscript_settings
        s.runscript.postamble_lines = ['$AMSBIN/cpkf ams.rkf trim.rkf Molecule InputMolecule General EngineResults MDHistory && mv trim.rkf ams.rkf']

        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            self.settings.input.ams.MolecularDynamics.InitialVelocities.Type = 'FromFile'
            self.settings.input.ams.MolecularDynamics.InitialVelocities.File = nvt_eq_job.results.rkfpath()
            self.molecule = nvt_eq_job.results.get_main_molecule()

        @add_to_instance(job)
        def postrun(self):
            """ extract and save pressure tensor for easier loading later """
            pt = self.results.get_history_property('PressureTensor', 'MDHistory')
            pt = np.array(pt)
            np.save(self.path+"/pressuretensor.npy", pt)

        self.children[name] = job
        return self.children[name]


    def __init__(self, molecule, engine_settings=None, runscript_settings=None, nsteps=None, temperature=300, **kwargs):
        """
        molecule: Molecule
            3D molecule (liquid) with a low density.

        engine_settings: Settings
            Settings for the engine (``s.input.forcefield.type = 'GAFF'``)

        runscript_settings: Settings
            Settings for the runscript (``s.runscript.nproc = 4``)

        nsteps: dict
            Dictionary where the default key-values pairs are. Any keys present in the dictionary will override the default values.

            .. code-block:: python

                nsteps = {
                    'scan_density': 5000,
                    'nvt_pre_eq': 1000,
                    'npt': 100000,
                    'nvt_eq': 10000,
                    'production': 400000
                }

        temperature: float
            Temperature in K.

        kwargs: other options to be passed to the MultiJob constructor (for example the name)
        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        self.scan_to_density = 1.5
        self.timestep = 1.0
        self.temperature = temperature
        self.nsteps = {
            'scan_density': 5000,
            'nvt_pre_eq': 1000,
            'npt': 100000,
            'nvt_eq': 10000,
            'production': 400000
        }
        if nsteps:
            self.nsteps.update(nsteps)

        self.engine_settings = engine_settings or self._default_engine_settings()
        self.runscript_settings = runscript_settings or self._default_runscript_settings()

        scan_density_job = self._create_scan_density_job(molecule)
        nvt_pre_eq_job = self._create_nvt_pre_eq_job(scan_density_job)
        npt_job = self._create_npt_job(nvt_pre_eq_job)
        nvt_eq_job = self._create_nvt_eq_job(npt_job)
        production_job = self._create_production_job(nvt_eq_job)

