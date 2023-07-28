from scm.plams import (
    AMSJob,
    CRSJob,
    KFFile,
    Molecule,
    MultiJob,
    Results,
    Settings,
)

from scm.conformers import ConformersJob
from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob
from collections import OrderedDict
import os

class ADFCOSMORSConfFilter:
    def __init__(self,max_num_confs=None,max_energy_range=None):
        self.max_num_confs = max_num_confs
        self.max_energy_range = max_energy_range


class ADFCOSMORSConfResults(Results):
    pass


class ADFCOSMORSConfJob(MultiJob):
    _result_type = ADFCOSMORSConfResults

    def __init__(self,
        molecule,
        conf_gen=None,
        first_filter=None,
        additional=None,
        thing=True,
        final_filter = None,
        initial_conformers=500,
        coskf_dir = None,
        coskf_name = None):

        super().__init__(children = {})

        self.job_count = 0

        self.mol = molecule

        self.adf_results   = False
        self.cosmo_results = False

        self.coskf_dir     = coskf_dir
        self.coskf_name    = coskf_name
        
        self.initial_conformers = initial_conformers

        self.filters        = [first_filter]
        self.job_settings   = [None]
        if additional is not None:
            for js, f in additional:
                self.job_settings.append(js)
                self.filters.append(f)

        if not self.has_valid_filter_settings():
            pass

        if conf_gen is None:
            self.conf_gen = self.default_confgen()
        else:
            if not isinstance(conf_gen,ConformersJob):
                print("Wrong type for argument conf_gen.  Expected ConformersJob instance.  Using the default conformer generator.")
                self.conf_gen = self.default_confgen()
            else:
                self.conf_gen = conf_gen

        self.children['job_0'] = self.conf_gen
        self._add_filter(self.children['job_0'].settings)

    def default_confgen(self):

        sett = Settings()
        sett.input.AMS.Generator.RDKit
        sett.input.AMS.Generator.RDKit.InitialNConformers = self.initial_conformers
        return ConformersJob(name="conformers_uff", molecule=self.mol, settings=sett)
        
    def make_intermediate_job(self,settings):

        settings.input.AMS.InputConformersSet = self.children[f'job_{self.job_count-1}'].results
        self._add_filter(settings)
        return ConformersJob(name = f"additional_{self.job_count}",settings=settings)

    def make_adf_job(self):

        sett = ADFCOSMORSCompoundJob.adf_settings(False)
        sett.input.AMS.Task = "Optimize"
        sett.input.AMS.GeometryOptimization.UseAMSWorker = "False"

        sett.input.AMS.InputConformersSet = self.children[f'job_{self.job_count-1}'].results
        
        return ConformersJob(name="adf_conformers", settings=sett)
    
    def make_cosmo_job(self):

        sett = ADFCOSMORSCompoundJob.adf_settings(True,elements=list(set(at.symbol for at in self.mol)))
        sett.input.AMS.Task = "Replay"
        sett.input.AMS.Replay.File = self.children['adf_job'].results["conformers.rkf"]
        sett.input.AMS.Replay.StoreAllResultFiles = "True"
        return AMSJob(name="replay", settings=sett)

    def new_children(self):

        if self.job_count < len(self.job_settings)-1:
            self.job_count += 1
            settings = self.job_settings[self.job_count]
            new_job = self.make_intermediate_job(settings)
            return {f'job_{self.job_count}':new_job}
        
        if 'adf_job' not in self.children:
            self.job_count += 1
            return {'adf_job':self.make_adf_job()}

        if 'cosmo_job' not in self.children:
            self.job_count += 1
            return {'cosmo_job':self.make_cosmo_job()}

        return None

    def postrun(self):
        self._make_coskfs()

    def _make_coskfs(self):
        '''
        Copy the COSMO sections from all the files back to the folder with the conformers
        '''
        base_name = self.coskf_name if self.coskf_name is not None else "conformer"
        if self.coskf_dir is None:
            self.coskf_dir = self.children['adf_job'].path
        for i, E in enumerate(self.children['adf_job'].results.get_energies("Ha")):
            if f"Frame{i+1}.rkf" in self.children['cosmo_job'].results:
                cosmo_section = self.children['cosmo_job'].results.read_rkf_section("COSMO", f"Frame{i+1}")
                cosmo_section["Gas Phase Bond Energy"] = E
                name = f"{base_name}_{i}.coskf"
                coskf = KFFile(os.path.join(self.coskf_dir, name), autosave=False)
                for key, val in cosmo_section.items():
                    coskf.write("COSMO", key, val)
                coskf.save()

    def _add_filter(self,sett):

        filt = self.filters[self.job_count]
        if filt is not None:
            if filt.max_num_confs is not None:
                sett.input.AMS.InputMaxConfs = filt.max_num_confs
            if filt.max_energy_range is not None:
                sett.input.AMS.InputMaxEnergy = filt.max_energy_range


    def has_valid_filter_settings(self):
        for js, f in zip(self.job_settings,self.filters):
            if js is not None and not isinstance(js,Settings):
                return False
            if f is not None and not isinstance(f,ADFCOSMORSConfFilter):
                return False
        return True
