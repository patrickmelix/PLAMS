from collections import OrderedDict

from scm.plams.core.basejob import MultiJob
from scm.plams.core.results import Results
from scm.plams.core.settings import Settings
import numpy as np
from scm.plams.recipes.md.amsmdjob import AMSNVEJob

__all__ = ['AMSNVESpawnerJob', 'AMSNVESpawnerResults']

class AMSNVESpawnerResults(Results):
    """Results class for AMSNVESpawnerJob
    """


class AMSNVESpawnerJob(MultiJob):
    """A class for running multiple NVE simulations with initial structures/velocities taken from an NVT trajectory. The NVT trajectory must contain the velocities!
    """

    _result_type = AMSNVESpawnerResults

    def _default_settings(self):
        s = Settings()
        s.input.ForceField.Type = 'GAFF'
        s.input.ForceField.AnteChamberIntegration = 'Yes'
        return s

    def __init__(self,
                 previous_job,  # needs to be finished
                 n_nve=1,
                 name='nvespawnerjob',
                 **kwargs):
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Must contain velocities (WriteVelocities Yes). Note that the trajectory should have been equilibrated before it starts.

        n_nve : int
            The number of NVE simulations to spawn

        All other settings can be set as for an AMSNVEJob (e.g. ``nsteps``).

        """
        MultiJob.__init__(self, children=OrderedDict(), name=name)

        self.previous_job = previous_job
        self.n_nve = n_nve
        self.nve_constructor_settings = kwargs
        self.nve_jobs = []

    def prerun(self):  # noqa F811
        """
        Constructs the children jobs
        """

        # use previously run previous_job
        assert self.previous_job.status != 'created', "You can only pass in a finished AMSJob"
        try:
            self.previous_job.results.readrkf('MDHistory', 'Velocities(1)')
        except KeyError:
            raise KeyError("Couldn't read velocities from {}".format(self.previous_job.results.rkfpath()))

        nframes_in_history = self.previous_job.results.readrkf('History', 'nEntries')

        if self.n_nve > 0:
            interval = nframes_in_history // self.n_nve
            frames = np.linspace(interval, nframes_in_history, self.n_nve, dtype=int)
            for i, frame in enumerate(frames):
                # the nve jobs only get different inside prerun, so need to do something to make them different to prevent the PLAMS rerun prevention
                name = f'nve{i+1}'
                self.settings.input.ams['#'] = name # add a comment line

                self.children[name] = AMSNVEJob.restart_from(
                    self.previous_job,
                    frame=frame,
                    name=name,
                    use_prerun=True,
                    **self.nve_constructor_settings
                )
                self.children[name].from_frame = frame

                self.nve_jobs.append(self.children[name])


