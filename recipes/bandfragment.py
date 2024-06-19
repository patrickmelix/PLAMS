from scm.plams.core.functions import log
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.core.settings import Settings
from scm.plams.recipes.adffragment import ADFFragmentJob, ADFFragmentResults
from scm.plams.tools.units import Units

from os.path import join as opj
from typing import Dict

# Check if ASE is present as it is not a default dependency of PLAMS.
__all__ = ['BANDFragmentJob', 'BandFragmentResults']


class BandFragmentResults(ADFFragmentResults):

    def get_energy_decomposition(self, unit='kJ/mol') -> Dict[str, float]:
        res=self.job.full.results
        res1=self.job.f2.results
        res2=self.job.f1.results

        res = {}

        pos = 2
        res['E_int'] = Units.convert(float(res.grep_output('E_int')[-1].split()[pos]), 'au', unit)
        res['E_disp'] = Units.convert(float(res.grep_output('E_disp')[-1].split()[pos]), 'au', unit)
        res['E_Pauli'] = Units.convert(float(res.grep_output('E_Pauli')[-1].split()[pos]), 'au', unit)
        res['E_elstat'] = Units.convert(float(res.grep_output('E_elstat')[-1].split()[pos]), 'au', unit)
        res['E_orb'] = Units.convert(float(res.grep_output('E_orb')[-1].split()[pos]), 'au', unit)

        res['E_1'] = res1.get_energy(unit=unit)
        res['E_2'] = res2.get_energy(unit=unit)
        
        return res


class BANDFragmentJob(ADFFragmentJob):
    _result_type = BandFragmentResults

    def create_mapping_setting(self):
        if "fragment" in self.full_settings.input.band:
            log("Fragment already present in full_settings. Assuming that the user has already set up the mapping. Skipping the mapping setup.", level=1)
            return
        set1 = Settings()
        set1.filename = opj(self.f1.path,'band.rkf')
        # if no mapping is saved, generate a basic mapping,
        # first fragment 1 then fragment 2
        set2 = Settings()
        set2.filename = opj(self.f2.path,'band.rkf')
        if self.fragment_mapping:
            set1.atommapping = { i+1: i+1 for i in range(len(self.fragment1)) }
            set2.atommapping = { i+1: i+1+len(self.fragment1) for i in range(len(self.fragment2)) }
        else:
            set1.atommapping = self.fragment_mapping[0].copy()
            set1.atommapping = self.fragment_mapping[1].copy()
        self.full_settings.input.band.fragment = [set1, set2]


    def prerun(self):
        self.f1 = AMSJob(name='frag1', molecule=self.fragment1, settings=self.settings)
        self.f2 = AMSJob(name='frag2', molecule=self.fragment2, settings=self.settings)
        self.children=[self.f1,self.f2]

        # create the correct mapping settings for the full job
        self.create_mapping_setting()

        self.full = AMSJob(name = 'full',
            molecule = self.fragment1 + self.fragment2,
            settings = self.settings + self.full_settings)
        self.full.depend += [self.f1,self.f2]

        self.children.append(self.full)


class NOCVBandFragmentJob(BANDFragmentJob):
    _result_type = BandFragmentResults

    def __init__(self, nocv_settings=None, **kwargs):
        """Create a BANDFragmentJob with NOCV calculation.
        
        Args:
            nocv_settings (Settings, optional): Settings for the NOCV calculation. Defaults to None.
        
        """
        super().__init__(**kwargs)
        self.nocv_settings = nocv_settings or Settings()

    def prerun(self):
        # prerun normal BANDFragmentJob
        super().prerun()
        # add NOCV run
        self.nocv=AMSJob(name="NOCV", molecule= self.fragment1 + self.fragment2,
        settings = self.settings + self.full_settings+self.nocv_settings)
        # make sure we at least have the restart section to put the file name
        if not hasattr(self.nocv.settings, 'input'):
            self.nocv.settings.input = Settings()
        if not hasattr(self.nocv.settings.input, 'band'):
            self.nocv.settings.input.band = Settings()
        if not hasattr(self.nocv.settings.input.band, 'restart'):
            self.nocv.settings.input.band.restart = Settings()
        if not hasattr(self.nocv.settings.input.band.restart, 'file'):
            self.nocv.settings.input.band.Restart.File = opj(self.full.path,'band.rkf')
        self.children.append(self.nocv)