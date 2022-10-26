#!/usr/bin/env amspython
from ase.io import read
from ase.calculators.socketio import SocketClient
from scm.plams import Settings, AMSCalculator

"""
Example illustrating how to use AMS as a client together with i-PI.

The i-PI configuration in input.xml sets up a short ring-polymer molecular
dynamics (RPMD) simulation. 

The connection between server and client is of type 'unixsocket' (a file in
/tmp). This example will only run on Linux.

The settings below set up AMS to use the UFF force field to calculate the
energy, forces, and stress tensor of whatever system the i-PI server requests
and return them to i-PI.

AMS is run in "AMSWorker" mode (interactive mode), which means that AMS does
**not** shut down and start up again between calculations.

To run this example, 

* modify run-server.sh to provide the correct path to the i-pi executable,
* sh run-server.sh
* In a different terminal: $AMSBIN/amspython run-ase.py
"""

def main():
    use_stress = True
    atoms = read('firstframe.xyz')

    sett = Settings()
    sett.input.ams.Task = 'SinglePoint'
    sett.input.ams.Properties.Gradients = 'True'
    sett.input.ams.Properties.StressTensor = str(use_stress)
    sett.input.ForceField.Type = 'UFF'
    sett.runscript.nproc = 1

    calc = AMSCalculator(settings=sett, amsworker=True)
    atoms.set_calculator(calc)

    client = SocketClient(unixsocket='driver-irpmd-16')

    try:
        client.run(atoms, use_stress=use_stress)
    finally:
        calc.worker.stop()

if __name__ == '__main__':
    main()
