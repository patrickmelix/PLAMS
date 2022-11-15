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

def print_results(names, counts_list):
    """
        This function prints the results in a space-separated table.

        names: list of str
            Molecule names

        counts_list: list of list of int
            List with dimensions nFrames x nMoleculeTypes
    """
    header = "frame " + " ".join(names)
    print(header)

    for frame, counts in enumerate(counts_list, 1):
        line = f"{frame} " + " ".join(str(x) for x in counts)
        print(line)

def analyze(ams_rkf_path):
    """
        ams_rkf_path: str
            Path to an ams.rkf file from a reactive MD simulation

        Returns: 2-tuple (names, counts)
            ``names``: list of length nMoleculesTypes. ``counts``: list of
            length nFrames, each item a list of length nMoleculesTypes
            containing an integer with the number of molecules of that type at
            the particular frame.

    """
    
    if not os.path.exists(ams_rkf_path):
        raise FileNotFoundError(f"Couldn't find the file {ams_rkf_path}")

    job = AMSJob.load_external(ams_rkf_path)

    try:
        n_molecules = job.results.readrkf('Molecules', 'Num molecules')
    except KeyError:
        raise KeyError("Couldn't find Molecules section on ams.rkf. You need to enable MolecularDynamics%Trajectory%WriteMolecules (before running the MD simulation).")

    # get the names of the molecules (molecular formula)
    molecule_type_range = range(1, n_molecules+1)  # 1, 2, ..., n_molecules
    names = [job.results.readrkf('Molecules', f'Molecule name {i}') for i in molecule_type_range]

    # read the Mols.Type from each History element
    mols_type_list = job.results.get_history_property('Mols.Type') # list of length nFrames
    counts_list = []

    # loop over the frames
    # store the counts-per-molecule-type in counts_list
    for mols_types in mols_type_list:
        counts = Counter(mols_types)
        counts_list.append([counts[x] for x in molecule_type_range])

    return names, counts_list

def main():
    if len(sys.argv) != 2:
        print("Usage: $AMSBIN/amspython MoleculesTable.py /path/to/ams.rkf")
        exit(1)

    ams_rkf_path = sys.argv[1]

    try:
        names, counts_list = analyze(ams_rkf_path)
        print_results(names, counts_list)
    
    except Exception as e:
        print(e)
        exit(1)

if __name__ == '__main__':
    main()

