from scm.plams import *
from ase import Atoms


#CO test system
d = 1.1
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])



def test_AMSJob():
    s = Settings()
    s.input.ADF
    s.Calculator.AMSJob
    s.input.ams.Task = 'SinglePoint'

    job = AMSCalculator(s)
    job.calculate(co)
    print('energy =', job.results['energy'])
    job.calculate(properties = 'forces')
    print('force list')
    for force_vec in job.results['forces']:
        print(' '.join(map(str,force_vec)))


def test_AMSPipe():
    s = Settings()
    s.input.ADF
    s.Calculator.Pipe
    s.input.ams.Task = 'SinglePoint'

    job = AMSCalculator(s)
    job.calculate(co)
    print('energy =', job.results['energy'])

    job.calculate(properties = 'forces')
    print('force list')
    for force_vec in job.results['forces']:
        print(' '.join(map(str,force_vec)))


def test_System():
    s = Settings()
    s.input.ADF
    s.Calculator.Pipe
    s.input.ams.Task = 'SinglePoint'

    s.input.ams.system.atoms._1 = f'C     0.0    0.0     0.0'
    s.input.ams.system.atoms._2 = f'O     1.0    0.0     {d}'
    s.input.ams.system.charge = 1.0

    job = AMSCalculator(s)
    job.calculate()
    print('energy =', job.results['energy'])

    job.calculate(properties = 'forces')
    print('force list')
    for force_vec in job.results['forces']:
        print(' '.join(map(str,force_vec)))



def test_Properties():
    s = Settings()
    s.input.ADF
    s.Calculator.AMSJob
    s.input.ams.Task = 'SinglePoint'
    s.input.ams.Properties.Gradients = True
    job = AMSCalculator(s)
    job.calculate(co)
    print('found forces:', 'forces' in job.results )


def test_PipeWorker():
    s = Settings()
    s.input.ADF
    s.Calculator.Pipe
    s.input.ams.Task = 'SinglePoint'
    s.amsworker.quiet = True
    s.Calculator.Pipe.Worker.use_restart_cache=True
    job = AMSCalculator(s)
    job.calculate(co)


if __name__ == '__main__':
    init(folder = 'workdir_AMSJob')
    test_AMSJob()
    finish()
    init(folder = 'workdir_AMSPipe')
    test_AMSPipe()
    finish()
    init(folder = 'workdir_System')
    test_System()
    finish()
    init(folder = 'workdir_Properties')
    test_Properties()
    finish()
    init(folder = 'workdir_PipeWorker')
    test_PipeWorker()
    finish()

