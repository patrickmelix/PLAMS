from scm.plams import *

def test_Properties():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = 'SinglePoint'
    s.runscript.nproc = 1
    job = AMSCalculator(s, name = 'Properties')
    assert 'forces' not in job.implemented_properties
    job.ensure_property('forces')
    assert 'forces' in job.implemented_properties


if __name__ == '__main__':
    test_Properties()
