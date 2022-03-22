from ...core.basejob import Job, SingleJob
from ...mol.molecule import Molecule
from .ams import AMSJob

__all__ = ['AMSConformersJob']

class AMSConformersJob(SingleJob):

    def __init__(self, name='amsconformers', molecule=None, **kwargs):
        Job.__init__(self, name=name, **kwargs)
        if not isinstance(molecule, Molecule): raise NotImplementedError("TODO: Implement multiple molecules input.")
        self.molecule = molecule.copy()

    def check(self):
        return True

    def get_input(self):
        return AMSJob(settings=self.settings, molecule=self.molecule).get_input()

    def get_runscript(self):
        ret = ''
        if 'preamble_lines' in self.settings.runscript:
            for line in self.settings.runscript.preamble_lines:
                ret += f'{line}\n'
        ret += 'AMS_JOBNAME="{}" AMS_RESULTSDIR=. $AMSBIN/amsconformers'.format(self.name)
        if 'nproc' in self.settings.runscript:
            ret += ' -n {}'.format(self.settings.runscript.nproc)
        ret += ' <"{}"'.format(self._filename('inp'))
        if self.settings.runscript.stdout_redirect:
            ret += ' >"{}"'.format(self._filename('out'))
        ret += '\n\n'
        return ret
