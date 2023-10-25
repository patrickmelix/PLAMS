#!/usr/bin/env amspython
from scm.plams import *
import numpy as np
import os
import matplotlib.pyplot as plt
from typing import Sequence, Tuple, List
from scm.plams.trajectories.analysis import autocorrelation, power_spectrum
from pathlib import Path

""" 

Example illustrating how to calculate an IR spectrum as an average over multiple NVE simulations.

This only works for gasphase molecules for which the AMS engine returns atomic charges.

For more information see the "MD intro with Python" tutorial in AMS

"""


def main():
    init()
    mol = from_smiles("CO", forcefield="uff")

    s = Settings()
    s.input.ReaxFF.ForceField = "CHO.ff"
    s.runscript.preamble_lines = ["export OMP_NUM_THREADS=1"]
    s.runscript.nproc = 1  # run in serial

    nvt_eq_job = AMSNVTJob(
        molecule=mol,
        settings=s,
        name="NVT-equilibration-1",
        temperature=300,
        thermostat="Berendsen",
        tau=200,
        timestep=1.0,
        nsteps=1000,
        samplingfreq=500,
    )
    nvt_eq_job.run()

    nvt_prod_job = AMSNVTJob.restart_from(
        nvt_eq_job,
        name="NVT-production-1",
        nsteps=5000,
        thermostat="NHC",
        samplingfreq=100,
        writecharges=True,
    )

    nvt_prod_job.run()

    nvespawner_job = AMSNVESpawnerJob(
        nvt_prod_job,
        name="nvespawner-" + nvt_prod_job.name,
        n_nve=2,  # the number of NVE simulations to run
        samplingfreq=1,  # write to disk every frame for velocity autocorrelation
        writevelocities=True,
        writecharges=True,
        writebonds=False,
        timestep=0.5,  # ideally use smaller timestep
        nsteps=10000,
    )
    nvespawner_job.run()

    freq, intens = get_average_spectrum(nvespawner_job.children.values())

    xy = np.stack((freq, intens), axis=1)

    np.savetxt("spectrum.txt", xy)
    log("Saving spectrum in spectrum.txt")

    finish()


def get_average_spectrum(jobs: Sequence[AMSJob]) -> Tuple[List[float], List[float]]:
    all_spectra = []
    for job in jobs:
        freq, intens = manual_single_spectrum(job)
        all_spectra.append(intens)

    intens = np.mean(all_spectra, axis=0)

    freq, intens = moving_average(freq, intens, window=20)

    return freq, intens


def moving_average(x, y, window: int) -> Tuple[List[float], List[float]]:
    if not window:
        return x, y
    window = min(len(x) - 1, window)
    if window <= 1:
        return x, y
    ret_x = np.convolve(x, np.ones(window) / window, mode="valid")
    ret_y = np.convolve(y, np.ones(window) / window, mode="valid")
    return ret_x, ret_y


def manual_single_spectrum(job) -> Tuple[List[float], List[float]]:
    coords = job.results.get_history_property("Coords")
    charges = job.results.get_history_property("Charges", "MDHistory")

    time_step = job.results.get_time_step()

    nEntries = len(coords)

    coords = np.array(coords).reshape(nEntries, -1, 3) * Units.convert(
        1.0, "bohr", "angstrom"
    )
    charges = np.array(charges).reshape(nEntries, -1, 1)

    qd = coords * charges

    dipole_moment = np.sum(qd, axis=1)
    dipole_dt = np.diff(dipole_moment, axis=0)

    data = dipole_dt

    acf = autocorrelation(data, max_dt=5000, normalize=False)
    times = np.arange(len(data)) * time_step
    freq, intens = power_spectrum(times, acf, max_freq=4000)

    intens *= 1e5  # approximate scaling / unit fix

    return freq, intens


if __name__ == "__main__":
    main()
