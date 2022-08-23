from collections import OrderedDict
from ...core.functions import add_to_instance
from ...core.basejob import MultiJob
from ...core.results import Results
from ...core.settings import Settings
from ...mol.molecule import Molecule
from ...mol.atom import Atom
from ...interfaces.adfsuite.ams import AMSJob
from ...tools.units import Units
from ...interfaces.molecule.packmol import packmol_mixture
import numpy as np
from scipy.optimize import curve_fit

__all__ = ['EquilibrateDensityJob', 'EquilibrateDensityResults']

class EquilibrateDensityResults(Results):
    """Results class for EquilibrateDensityJob
    """
        
    def get_density(self):
        """
        Returns the density in g/cm^3
        """

        return self.get_molecule('equilibrated').get_density() * 1e-3

    def get_molecule(self, index='equilibrated', return_index=False):
        """
        Returns a molecule from the NPT simulation

        index : str
            'equilibrated': returns a molecule as close as possible to the average density from the last third of the simulation. 'final': returns the final Molecule.

        return_index : bool
            If True ,returns a 2-tuple (Molecule, 1-based-index) containing the index in History for which the molecule is returned.
        """

        assert index in ['equilibrated', 'final']
        job = self.job.children['npt']
        if index == 'equilibrated':
            densities = job.results.get_history_property('Density', 'MDHistory')
            analyze_from = (len(densities)*2) // 3
            # take structure closest to the target density
            avg_density = np.mean(densities[analyze_from:]) # amu/bohr^3
            delta = np.array(densities[analyze_from:]) - avg_density
            delta = np.abs(delta)
            min_index = np.argmin(delta)
            min_index += analyze_from
            mol = job.results.get_history_molecule(min_index+1)
        elif index == 'final':
            nEntries = job.results.readrkf('History', 'nEntries')
            mol = job.results.get_history_molecule(nEntries)
            min_index = nEntries

        if return_index:
            return mol, min_index
        else:
            return mol

    def get_velocities(self, index='equilibrated'):
        if index == 'final':
            return self.job.children['npt'].results.readrkf('MDResults', 'EndVelocities')

        mol, min_index = self.get_main_molecule(index, return_index=True)
        vels = job.results.readrkf('MDHistory', f'Velocities({min_index})')
        return vels

    def rkfpath(self):
        """ Returns the path to ams.rkf from the npt job """
        return self.job.children['npt'].rkfpath()

class EquilibrateDensityJob(MultiJob):
    """A class for equilibrating the density at a certain temperature and pressure
    """

    _result_type = EquilibrateDensityResults

    def _default_settings(self):
        s = Settings()
        s.input.ForceField.Type = 'GAFF'
        s.input.ForceField.AnteChamberIntegration = 'Yes'
        return s

    def _create_scan_density_job(self, initial_molecule):
        name = 'scan_density'

        from_density = initial_molecule.get_density() * 1e-3
        orig_length = initial_molecule.cell_lengths()
        density_ratio = from_density / self.scan_density_upper
        new_length = [x *  density_ratio**0.333333 for x in orig_length]

        s = Settings()
        s.input.ams.task = 'MolecularDynamics'
        s.input.ams.MolecularDynamics.Nsteps = self.nsteps[name]
        s.input.ams.MolecularDynamics.TimeStep = self.timestep
        s.input.ams.MolecularDynamics.Deformation.TargetLength = ' '.join([str(x) for x in new_length]) 
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

        s += self.settings

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
        s += self.settings

        job = AMSJob(settings=s, name=name)
        if scan_density_job is not None:
            @add_to_instance(job)
            def prerun(self):
                self.molecule = scan_density_job.get_lowest_energy_molecule()
        else:
            job.molecule = self.initial_molecule

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
        s.input.ams.MolecularDynamics.Barostat.Pressure = f'{self.pressure} [bar]'
        s.input.ams.MolecularDynamics.Barostat.Tau = 1000
        s.input.ams.MolecularDynamics.Barostat.Equal = 'XYZ'
        s.input.ams.MolecularDynamics.Trajectory.WriteVelocities = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteBonds = 'Yes'
        s.input.ams.MolecularDynamics.Trajectory.WriteMolecules = 'No'
        s.input.ams.MolecularDynamics.Trajectory.WriteCharges = 'No'
        s.input.ams.MolecularDynamics.Checkpoint.Frequency = 10000000
        s += self.settings

        job = AMSJob(settings=s, name=name)
        @add_to_instance(job)
        def prerun(self):
            self.settings.input.ams.MolecularDynamics.InitialVelocities.Type = 'FromFile'
            self.settings.input.ams.MolecularDynamics.InitialVelocities.File = nvt_pre_eq_job.results.rkfpath()
            self.molecule = nvt_pre_eq_job.results.get_main_molecule()


        self.children[name] = job
        return self.children[name]
        
    def __init__(self, 
                 molecule, 
                 settings=None, 
                 nsteps=None, 
                 temperature=300, 
                 pressure=1.0, 
                 scan_density=True, 
                 scan_density_upper=1.5, 
                 **kwargs):
        """
        molecule: Molecule
            3D molecule (liquid/gas with multiple molecules). 

        settings: Settings
            All non-AMS-Driver settings, for example (``s.input.forcefield.type = 'GAFF'``, ``s.runscript.nproc = 1``)

        nsteps: dict
            Dictionary where the default key-values pairs are. Any keys present in the dictionary will override the default values.

            .. code-block:: python

                nsteps = {
                    'scan_density': 5000,
                    'nvt_pre_eq': 1000,
                    'npt': 100000,
                }

        temperature: float
            Temperature in K.

        pressure: float
            The pressure in bar.

        kwargs: other options to be passed to the MultiJob constructor (for example the name)
        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        self.scan_density_upper = scan_density_upper
        self.timestep = 1.0
        self.temperature = temperature
        self.pressure = pressure
        self.nsteps = {
            'scan_density': 5000,
            'nvt_pre_eq': 1000,
            'npt': 100000,
            'nvt_eq': 10000,
            'production': 400000
        }
        if nsteps:
            self.nsteps.update(nsteps)

        self.settings = settings.copy() or self._default_settings()

        if scan_density:
            scan_density_job = self._create_scan_density_job(molecule)
        else:
            scan_density_job = None
        nvt_pre_eq_job = self._create_nvt_pre_eq_job(scan_density_job)
        npt_job = self._create_npt_job(nvt_pre_eq_job)

