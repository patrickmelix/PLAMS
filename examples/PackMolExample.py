#!/usr/bin/env amspython
from scm.plams import from_smiles, packmol_liquid, packmol_mixture, packmol_solid_liquid, packmol_solid_liquid_mixture, fromASE
from ase.build import fcc111

def printsummary(mol):
    print(f'{len(mol)} atoms, density = {mol.get_density()*1e-3:.3f} g/cm^3, box = {mol.lattice[0][0]:.3f}, {mol.lattice[1][1]:.3f}, {mol.lattice[2][2]:.3f}, formula = {mol.get_formula()}')

def main():
    water = from_smiles('O')
    acetonitrile = from_smiles('CC#N')

    print('PURE LIQUID\n')

    print('pure liquid from approximate number of atoms and exact density (in g/cm^3), cubic box with auto-determined size')
    out = packmol_liquid(water, n_atoms=193, density=1.0)
    printsummary(out)
    out.write('water-1.xyz')

    print('pure liquid from approximate density (in g/cm^3) and an orthorhombic box')
    out = packmol_liquid(water, density=1.0, box_bounds=[0., 0., 0., 8., 12., 14.])
    printsummary(out)
    out.write('water-2.xyz')

    print('pure liquid with explicit number of molecules and exact density')
    out = packmol_liquid(water, n_molecules=64, density=1.0)
    printsummary(out)
    out.write('water-3.xyz')

    print('pure liquid with explicit number of molecules and box')
    out = packmol_liquid(water, n_molecules=64, box_bounds=[0., 0., 0., 12., 13., 14.])
    printsummary(out)
    out.write('water-4.xyz')

    # MIXTURES
    x_water = 0.666                # mole fraction
    x_acetonitrile = 1-x_water     # mole fraction
    density = (x_water*1.0 + x_acetonitrile*0.76) / (x_water + x_acetonitrile)  # weighted average of pure component densities

    print(f'\nMIXTURES. x_water = {x_water:.3f}, x_acetonitrile = {x_acetonitrile:.3f}, target density = {density:.3f} g/cm^3\n')

    print('2-1 water-acetonitrile from approximate number of atoms and exact density (in g/cm^3), cubic box with auto-determined size')
    out = packmol_mixture(molecules=[water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], n_atoms=200, density=density)
    printsummary(out)
    out.write('water-acetonitrile-1.xyz')

    print('2-1 water-acetonitrile from approximate density (in g/cm^3) and box bounds')
    out = packmol_mixture(molecules=[water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], box_bounds=[0, 0, 0, 13.2, 13.2, 13.2], density=density)
    printsummary(out)
    out.write('water-acetonitrile-2.xyz')

    print('2-1 water-acetonitrile from explicit number of molecules and density, cubic box with auto-determined size')
    out = packmol_mixture(molecules=[water, acetonitrile], n_molecules=[32, 16], density=density)
    printsummary(out)
    out.write('water-acetonitrile-3.xyz')

    print('2-1 water-acetonitrile from explicit number of molecules and box')
    out = packmol_mixture(molecules=[water, acetonitrile], n_molecules=[32, 16], box_bounds=[0, 0, 0, 13.2, 13.2, 13.2])
    printsummary(out)
    out.write('water-acetonitrile-4.xyz')

    print('\nSOLID-LIQUID INTERFACES\n')
    slab = fromASE(fcc111('Al', size=(4,6,3), vacuum=15.0, orthogonal=True, periodic=True))

    print('water surrounding an Al slab, from an approximate density')
    out = packmol_solid_liquid(slab, water, density=1.0)
    printsummary(out)
    out.write('al-water-pure.xyz')

    print('2-1 water-acetonitrile mixture surrounding an Al slab, from mole fractions and an approximate density')
    out = packmol_solid_liquid_mixture(slab, [water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], density=density)
    printsummary(out)
    out.write('al-water-acetonitrile.xyz')

if __name__ == '__main__':
    main()

