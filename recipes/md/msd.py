from collections import OrderedDict, defaultdict
from ...core.functions import add_to_instance
from ...core.basejob import MultiJob
from ...core.results import Results
from ...core.settings import Settings
from ...mol.molecule import Molecule
from ...mol.atom import Atom
from ...interfaces.adfsuite.ams import AMSJob
from ...interfaces.adfsuite.amsanalysis import AMSAnalysisJob, AMSAnalysisResults
from ...tools.units import Units
from .amsmdjob import AMSNVEJob
import numpy as np

__all__ = ['AMSMSDJob', 'AMSMSDResults', 'AMSMSDPerRegionJob', 'AMSMSDPerRegionResults']

class AMSMSDResults(AMSAnalysisResults):
    """Results class for AMSMSDJob
    """
    def get_msd(self):
        """
        returns time [fs], msd [ang^2]
        """
        msd_xy = self.get_xy()
        time = np.array(msd_xy.x[0]) # fs
        y = np.array(msd_xy.y) * Units.convert(1.0, 'bohr', 'angstrom')**2

        return time, y

    def get_linear_fit(self, start_time_fit_fs=None, end_time_fit_fs=None):
        """
            Fits the MSD between start_time_fit_fs and end_time_fit_fs

            Returns a 3-tuple LinRegress.result, fit_x (fs), fit_y (ang^2)

            result.slope is given in ang^2/fs

        """
        from scipy.stats import linregress
        time, y = self.get_msd()
        end_time_fit_fs = end_time_fit_fs or max(time)
        start_time_fit_fs = start_time_fit_fs or self.job.start_time_fit_fs

        if start_time_fit_fs >= end_time_fit_fs:
            start_time_fit_fs = end_time_fit_fs / 2


        y = y[time >= start_time_fit_fs]
        time = time[time >= start_time_fit_fs] 
        y = y[time <= end_time_fit_fs] 
        time = time[time <= end_time_fit_fs] 

        result = linregress(time, y)
        fit_x = time
        fit_y = result.slope * fit_x + result.intercept

        return result, fit_x, fit_y

    def get_D(self, start_time_fit_fs=None, end_time_fit_fs=None):
        """
        Returns D in m^2/s
        """
        result, _, _ = self.get_linear_fit(start_time_fit_fs=start_time_fit_fs, end_time_fit_fs=end_time_fit_fs)
        D = result.slope * 1e-20 / (6 * 1e-15) # convert from ang^2/fs to m^2/s, divide by 6 because 3-dimensional (2d)
        return D

        

class AMSMSDJob(AMSAnalysisJob):
    """A class for equilibrating the density at a certain temperature and pressure
    """

    _result_type = AMSMSDResults

    def __init__(self, 
                 previous_job,  # needs to be finished
                 max_correlation_time_fs=10000,
                 start_time_fit_fs=2000,
                 atom_indices=None,

                 **kwargs):
        """
        previous_job: AMSJob
            An AMSJob with an MD trajectory. Note that the trajectory should have been equilibrated before it starts.

        All other settings can be set as for AMS

        """
        AMSAnalysisJob.__init__(self, **kwargs)

        self.previous_job = previous_job
        self.max_correlation_time_fs = max_correlation_time_fs
        self.start_time_fit_fs = start_time_fit_fs
        self.atom_indices = atom_indices

    def prerun(self):

        # use previously run previous_job
        assert self.previous_job.status != 'created', "You can only pass in a finished AMSJob"

        historylength = self.previous_job.results.readrkf('History', 'nEntries')
        max_dt_frames = int(self.max_correlation_time_fs / self.previous_job.results.get_time_step())
        max_dt_frames = min(max_dt_frames, historylength // 2)

        self.settings.input.Task = 'MeanSquareDisplacement'
        self.settings.input.MeanSquareDisplacement.Property = 'DiffusionCoefficient'
        self.settings.input.TrajectoryInfo.Trajectory.KFFileName = self.previous_job.results.rkfpath()
        self.settings.input.MeanSquareDisplacement.StartTimeSlope = self.start_time_fit_fs
        self.settings.input.MeanSquareDisplacement.MaxStep = max_dt_frames
        if self.atom_indices:
            self.settings.input.MeanSquareDisplacement.Atoms.Atom = self.atom_indices

class AMSMSDPerRegionResults(Results):
    def get_D(self, start_time_fit_fs=2000, end_time_fit_fs=10000):
        """
        Returns a dictionary of region_name: D (m^2)
        """
        ret = {}
        for name, job in self.job.children.items():
            ret[name] = job.results.get_D()
        return ret

class AMSMSDPerRegionJob(MultiJob):
    _result_type = AMSMSDPerRegionResults

    def __init__(self, 
                    previous_job:AMSJob,
                    max_correlation_time_fs:float=10000,
                    start_time_fit_fs:float=2000,
                    regions=None,
                    per_element=False,
                    **kwargs):

        """
            Creates an AMSAnalysisJob for every Region in previous_job

            regions : dict
                ``{ 'region-name': [1, 4, 8] }`` to create "new" regions with atom indices. If not given, it will use the regions from the previous_job.

        """

        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        self.previous_job = previous_job
        self.max_correlation_time_fs = max_correlation_time_fs
        self.start_time_fit_fs = start_time_fit_fs
        self.regions_dict = regions
        self.per_element = per_element


    @staticmethod
    def get_regions_dict(molecule, per_element:bool=False):
        regions_dict = defaultdict(lambda: [])
        for i, at in enumerate(molecule, 1):
            regions = set([at.properties.region]) if isinstance(at.properties.region, str) else at.properties.region
            if len(regions) == 0:
                region_name = 'NoRegion' if not per_element else f'NoRegion_{at.symbol}'
                regions_dict[region_name].append(i)
            for region in regions:
                region_name = region if not per_element else f'{region}_{at.symbol}'
                regions_dict[region_name].append(i)
            regions_dict['All'].append(i)
            if per_element:
                regions_dict[f'All_{at.symbol}'].append(i)

        return regions_dict

    def prerun(self):
        timestep = self.previous_job.results.get_time_step()
        mol = self.previous_job.results.get_main_molecule()
        regions_dict = self.regions_dict or self.get_regions_dict(mol, per_element=self.per_element)

        for region in regions_dict:
            self.children[region] = AMSMSDJob(
                previous_job=self.previous_job,
                name=region,
                max_correlation_time_fs=self.max_correlation_time_fs,
                start_time_fit_fs=self.start_time_fit_fs,
                atom_indices = regions_dict[region]
            )

