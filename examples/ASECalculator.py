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
    s.input.ams.Properties.Gradients = True
    job = AMSCalculator(s, name = 'Properties')
    job.calculate(co)
    print('found forces:', 'forces' in job.results )


def test_PipeWorker():
    s = Settings()
    s.input.ForceField
    s.input.ams.Task = 'SinglePoint'
    s.runscript.nproc = 1
    s.amsworker.quiet = True
    job = AMSCalculator(s, name = 'PipeWorker', amsworker = True, restart = True)
    job.calculate(co)


if __name__ == '__main__':
    test_AMSJob()
    test_AMSPipe()
    test_System()
    test_Properties()
    test_PipeWorker()
    finish()
