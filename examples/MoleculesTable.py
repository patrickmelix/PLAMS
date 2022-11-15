#!/usr/bin/env amspython
from scm.plams import *
from collections import Counter
import sys
import os

"""
Example printing a table of the AMS Molecule analysis from reactive MD simulations.

For each species and frame, print the number of species in that frame.

Usage:
$AMSBIN/amspython MoleculesTable.py /path/to/ams.rkf
"""

def analyze(ams_rkf_path):
    if not os.path.exists(ams_rkf_path):
        print(f"Couldn't find the file {ams_rkf_path}")
        return 1

    job = AMSJob.load_external(ams_rkf_path)

    try:
        n_molecules = job.results.readrkf('Molecules', 'Num molecules')
    except KeyError:
        print("Couldn't find Molecules section on ams.rkf. You need to enable MolecularDynamics%Trajectory%WriteMolecules (before running the MD simulation).")
        return 1

    # get the names of the molecules (molecular formula)
    molecule_type_range = range(1, n_molecules+1)  # 1, 2, ..., n_molecules
    names = [job.results.readrkf('Molecules', f'Molecule name {i}') for i in molecule_type_range]

    # read the Mols.Type from each History element
    mols_type_list = job.results.get_history_property('Mols.Type') # list of length nFrames
    frames_list = []
    counts_list = []

    # loop over the frames
    # store frame numbers in frames_list
    # store the counts-per-molecule-type in counts_list
    for frame, mols_types in enumerate(mols_type_list, 1):
        frames_list.append(frame)

        counts = Counter(mols_types)
        counts_list.append([counts[x] for x in molecule_type_range])

    # now print the data
    header = "frame " + " ".join(names)
    print(header)
    for frame, counts in zip(frames_list, counts_list):
        line = f"{frame} " + " ".join(str(x) for x in counts)
        print(line)

    return 0

def main():
    if len(sys.argv) != 2:
        print("Usage: $AMSBIN/amspython MoleculesTable.py /path/to/ams.rkf")
        exit(1)

    ams_rkf_path = sys.argv[1]

    exit_status = analyze(ams_rkf_path)
    exit(exit_status)


if __name__ == '__main__':
    main()

