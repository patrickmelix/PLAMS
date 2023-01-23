#!/usr/bin/env plams
import numpy as np

# Perform a geometry optimization of a water molecule and compute
# the vibrational normal modes using GFN1-xTB. If you do not have
# a DFTB license, remove the line with DFTB settings and instead set
# settings.input.ForceField.Type = 'UFF'

# You could also load the geometry from an xyz file: 
# molecule = Molecule('path/my_molecule.xyz')
# or generate a molecule from SMILES:
# molecule = from_smiles('O')
molecule = Molecule()
molecule.add_atom(Atom(symbol='O', coords=(0,0,0)))
molecule.add_atom(Atom(symbol='H', coords=(1,0,0)))
molecule.add_atom(Atom(symbol='H', coords=(0,1,0)))

settings = Settings()
settings.input.ams.Task = 'GeometryOptimization'
settings.input.ams.Properties.NormalModes = 'Yes'
settings.input.dftb.Model = 'GFN1-xTB'

job = AMSJob(molecule=molecule, settings=settings, name='water_optimization')
result = job.run()

# Calculation timing

timings = result.get_timings()

# Results

energy = result.get_energy(unit='kcal/mol')

try:
    homo = result.get_homo_energies(unit='eV')[0]
    lumo = result.get_lumo_energies(unit='eV')[0]

    homo_lumo_gap = result.get_smallest_HOMO_LUMO_gap(unit='eV')

    dipole_moment = np.linalg.norm(np.array(job.results.get_dipolemoment()))
    dipole_moment *= Units.convert(1.0, 'au', 'debye')
except (KeyError, AttributeError):
    homo = lumo = homo_lumo_gap = dipole_moment = None

frequencies = result.get_frequencies(unit='cm^-1')

# Optimized structure

optimized_molecule = result.get_main_molecule()

# Unlike python lists, where the index of the first element is 0, 
# the index of the first atom in the molecule object is 1
bond_angle = optimized_molecule[1].angle(optimized_molecule[2], optimized_molecule[3])

print('== Results ==')
print('Optimized geometry:')
print(optimized_molecule)
print('')
print('Bond angle  : {:.1f} degrees'.format(Units.convert(bond_angle, 'rad', 'degree')))
print('')
print('Calculation time: {:.3f} s'.format(timings['elapsed']))
print('')
print('Energy      : {:.3f} kcal/mol'.format(energy))

if homo is not None:
    print('HOMO        : {:.3f} eV'.format(homo))
    print('LUMO        : {:.3f} eV'.format(lumo))
    print('HOMO-LUMO gap : {:.3f} eV'.format(homo_lumo_gap))

if dipole_moment is not None:
    print('Dipole moment: {:.3f} debye'.format(dipole_moment))
