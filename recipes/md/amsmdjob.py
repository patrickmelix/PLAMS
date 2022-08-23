from ...interfaces.adfsuite.ams import AMSJob, AMSResults
from ...core.settings import Settings
from ...tools.kftools import KFFile
from ...tools.units import Units
from typing import Union
import numpy as np

__all__ = ['NVEJob', 'NVTJob', 'NPTJob']

class AMSMDJob(AMSJob):
    @staticmethod
    def _velocities2settings(velocities):
        s = Settings()
        if isinstance(velocities, int) or isinstance(velocities, float) or velocities is None:
            s.input.ams.MolecularDynamics.InitialVelocities.Type = 'Random'
            s.input.ams.MolecularDynamics.InitialVelocities.Temperature = velocities or 300
        elif isinstance(velocities, tuple):
            # file and frame number
            f = velocities[0]
            frame = velocities[1]
            vels = KFFile(f).read('MDHistory', f'Velocities({frame})')
            vels = np.array(vels).reshape(-1,3) * Units.convert(1.0, 'bohr', 'angstrom') # angstrom/fs
            s.input.ams.MolecularDynamics.InitialVelocities.Type = 'Input'
            values = ""
            for x in vels:
                values += 6*" " + " ".join(str(y) for y in x) + "\n"
            s.input.ams.MolecularDynamics.InitialVelocities.Values._h = '   # From {} frame {}'.format(f, frame)
            s.input.ams.MolecularDynamics.InitialVelocities.Values._1 = values
        else:
            s.input.ams.MolecularDynamics.InitialVelocities.Type = 'FromFile'
            s.input.ams.MolecularDynamics.InitialVelocities.File = velocities
        return s

    def __init__(
        self, 
        velocities=300,
        timestep=0.25,
        samplingfreq=100,
        nsteps=1000,
        checkpointfrequency=1000,
        writevelocities=True,
        writebonds=True,
        writemolecules=True,
        writecharges=True,
        writeenginegradients=False,
        calcpressure=False,
        molecule=None,
        **kwargs
    ):
        if isinstance(molecule, AMSJob):
            molecule = molecule.results.get_main_molecule()
        if isinstance(molecule, AMSResults):
            molecule = molecule.get_main_molecule()
        AMSJob.__init__(self, molecule=molecule, **kwargs)
        self.settings.input.ams.Task = 'MolecularDynamics'
        self.settings.input.ams.MolecularDynamics.TimeStep = timestep
        self.settings.input.ams.MolecularDynamics.Trajectory.SamplingFreq = samplingfreq
        self.settings.input.ams.MolecularDynamics.NSteps = nsteps
        self.settings.input.ams.MolecularDynamics.Trajectory.WriteVelocities = str(writevelocities)
        self.settings.input.ams.MolecularDynamics.Trajectory.WriteBonds = str(writebonds)
        self.settings.input.ams.MolecularDynamics.Trajectory.WriteMolecules = str(writemolecules)
        self.settings.input.ams.MolecularDynamics.Trajectory.WriteEngineGradients = str(writeenginegradients)
        self.settings.input.ams.MolecularDynamics.CalcPressure = str(calcpressure)
        self.settings += self._velocities2settings(velocities)

class NVEJob(AMSMDJob):
    def __init__( self, **kwargs):
        AMSMDJob.__init__(self, **kwargs)

        self.remove_blocks(['thermostat', 'barostat', 'deformation'])

    def remove_blocks(self, blocks=None):
        if blocks:
            for block in blocks:
                if block in self.settings.input.ams.MolecularDynamics:
                    del self.settings.input.ams.MolecularDynamics[block]

    @classmethod
    def restart_from(cls, other_job, frame=None, settings=None, **kwargs):
        other_job, velocities, molecule, extra_settings = cls._get_restart_job_velocities_molecule(other_job, frame, settings)
        return cls(settings=extra_settings, velocities=velocities, molecule=molecule, **kwargs)

    @staticmethod
    def _get_restart_job_velocities_molecule(other_job, frame=None, settings=None):
        """
            other_job: str or some AMSMdJob
        """
        if isinstance(other_job, str):
            other_job = AMSJob.load_external(other_job)
        if frame:
            velocities = (other_job.results.rkfpath(), frame)
            molecule = other_job.results.get_history_molecule(frame)
        else:
            velocities = other_job
            molecule = other_job.results.get_main_molecule()

        if settings:
            extra_settings = settings.copy()
        else:
            extra_settings = other_job.settings.copy()

        if 'InitialVelocities' in extra_settings.input.ams.MolecularDynamics:
            del extra_settings.input.ams.MolecularDynamics.InitialVelocities

        if 'System' in extra_settings.input.ams:
            del extra_settings.input.ams.System

        return other_job, velocities, molecule, extra_settings

class NVTJob(NVEJob):
    def __init__(self, 
        temperature=300,
        velocities=None,
        thermostat='NHC',
        tau=None,
        **kwargs):
        NVEJob.__init__(self, velocities = velocities or temperature, **kwargs)

        self.settings.input.ams.MolecularDynamics.Thermostat.Type = thermostat or 'NHC'
        self.settings.input.ams.MolecularDynamics.Thermostat.Temperature = temperature or 300
        self.settings.input.ams.MolecularDynamics.Thermostat.Tau = tau or float(self.settings.input.ams.MolecularDynamics.TimeStep) * 400

        self.remove_blocks(['barostat', 'deformation'])

    @classmethod
    def restart_from(cls, other_job, 
        settings=None,
        temperature=None,
        thermostat=None,
        tau=None,
        frame=None,
        **kwargs):

        other_job, velocities, molecule, extra_settings = cls._get_restart_job_velocities_molecule(other_job, frame, settings)

        thermostat = thermostat or other_job.settings.input.ams.MolecularDynamics.Thermostat.Type 
        temperature = temperature or other_job.settings.input.ams.MolecularDynamics.Thermostat.Temperature
        tau = tau or other_job.settings.input.ams.MolecularDynamics.Thermostat.Tau

        return cls(molecule=molecule, settings=extra_settings, velocities=velocities, thermostat=thermostat, temperature=temperature, tau=tau, **kwargs)


class NPTJob(NVTJob):
    def __init__(self,
        pressure=1.0,
        barostat='MTK',
        barostat_tau=None,
        scale='XYZ',
        equal='None',
        constantvolume=False,
        velocities=None,
        temperature=300,
        **kwargs
    ):
        NVTJob.__init__(self, velocities=velocities or temperature, **kwargs)
        self.settings.input.ams.MolecularDynamics.Barostat.Type = barostat
        self.settings.input.ams.MolecularDynamics.Barostat.Pressure = str(pressure) + ' [bar]'
        self.settings.input.ams.MolecularDynamics.Barostat.Tau = barostat_tau or float(self.settings.input.ams.MolecularDynamics.TimeStep) * 4000
        self.settings.input.ams.MolecularDynamics.Barostat.Scale = scale or 'XYZ'
        self.settings.input.ams.MolecularDynamics.Barostat.Equal = equal or 'None'
        self.settings.input.ams.MolecularDynamics.Barostat.ConstantVolume = str(constantvolume) if isinstance(constantvolume, (bool, str)) else "False"
        self.settings.input.ams.MolecularDynamics.CalcPressure = 'True'

        self.remove_blocks(['deformation'])

    @classmethod
    def restart_from(cls,
        other_job, 
        frame=None, 
        settings=None,
        pressure=None,
        barostat=None,
        barostat_tau=None,
        scale=None,
        equal=None,
        constantvolume=None,
        temperature=None,
        thermostat=None,
        tau=None,
        **kwargs
    ):
        
        other_job, velocities, molecule, extra_settings = cls._get_restart_job_velocities_molecule(other_job, frame, settings)

        thermostat = thermostat or other_job.settings.input.ams.MolecularDynamics.Thermostat.Type 
        temperature = temperature or other_job.settings.input.ams.MolecularDynamics.Thermostat.Temperature
        tau = tau or other_job.settings.input.ams.MolecularDynamics.Thermostat.Tau

        barostat = barostat or other_job.settings.input.ams.MolecularDynamics.Barostat.Type 
        pressure = pressure or other_job.settings.input.ams.MolecularDynamics.Barostat.Pressure 
        barostat_tau = barostat_tau or other_job.settings.input.ams.MolecularDynamics.Barostat.Tau 
        scale = scale or other_job.settings.input.ams.MolecularDynamics.Barostat.Scale 
        equal = equal or other_job.settings.input.ams.MolecularDynamics.Barostat.Equal 
        constantvolume = constantvolume if constantvolume is not None else other_job.settings.input.ams.MolecularDynamics.Barostat.ConstantVolume 

        return cls(
            molecule=molecule,
            velocities=velocities,
            settings=extra_settings,
            thermostat=thermostat,
            temperature=temperature,
            tau=tau,
            pressure=pressure,
            barostat=barostat,
            scale=scale,
            equal=equal,
            constantvolume=constantvolume,
            **kwargs
        )


