from scm.plams import *
from ase import Atoms


#CO test system
d = 1.1
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])



def test_AMSJob():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = 'SinglePoint'
    s.runscript.nproc = 1
    job = AMSCalculator(s, name = 'AMSJob')
    job.calculate(co)
    print('energy =', job.results['energy'])
    job.calculate(properties = 'forces')
    print('force list')
    for force_vec in job.results['forces']:
        print(' '.join(map(str,force_vec)))


def test_AMSPipe():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = 'SinglePoint'
    s.runscript.nproc = 1
    job = AMSCalculator(s, name = 'AMSPipe', amsworker = True)
    job.calculate(co)
    print('energy =', job.results['energy'])

    job.calculate(properties = 'forces')
    print('force list')
    for force_vec in job.results['forces']:
        print(' '.join(map(str,force_vec)))


def test_System():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = 'SinglePoint'
    s.runscript.nproc = 1
    s.input.ams.system.atoms._1 = f'C     0.0    0.0     0.0'
    s.input.ams.system.atoms._2 = f'O     0.0    0.0     {d}'

    job = AMSCalculator(s, name = 'System', amsworker = False, restart = True)
    job.calculate()
    print('energy =', job.results['energy'])

    job.calculate(properties = 'forces')
    print('force list')
    for force_vec in job.results['forces']:
        print(' '.join(map(str,force_vec)))



def test_Properties():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = 'SinglePoint'
    s.runscript.nproc = 1
    job = AMSCalculator(s, name = 'Properties')
    assert 'forces' not in job.implemented_properties
    job.set_properties('forces')
    assert 'forces' in job.implemented_properties


def test_PipeWorker():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = 'SinglePoint'
    s.runscript.nproc = 1
    s.amsworker.quiet = True
    job = AMSCalculator(s, name = 'PipeWorker', amsworker = True, restart = True)
    job.calculate(co)

def test_ase_deepcopy_worker():
    print("Deepcopy of AMSCalculator can be safely used for multiple atoms objects")
    settings = get_settings()
    atoms1 = get_atoms()
    atoms2 = get_atoms()
    #set OH bonds for atoms1 and atoms2 to 0.9 and 1.0 A respectively
    atoms1.set_distance(0,1, 0.9, fix = 0)
    atoms2.set_distance(0,1, 1.0, fix = 0)
    atoms1.set_distance(0,2, 0.9, fix = 0)
    atoms2.set_distance(0,2, 1.0, fix = 0)

    with AMSCalculator(settings=settings, name='ASE_deepcopy', amsworker=True) as calc:
        from copy import deepcopy
        atoms1.calc = deepcopy(calc)
        atoms2.calc = deepcopy(calc)
        e1 = atoms1.get_potential_energy()
        e2 = atoms2.get_potential_energy()
        #if no deepcopy is made, this would result in a new calculation (calc
        e1 = atoms1.get_potential_energy()
        assert calc._counter[calc.name] == 2

if __name__ == '__main__':
    test_AMSJob()
    test_AMSPipe()
    test_System()
    test_Properties()
    test_PipeWorker()
    test_ase_deepcopy_worker()
    finish()
