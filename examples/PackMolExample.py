#!/usr/bin/env amspython
from scm.plams import from_smiles, packmol, packmol_on_slab, fromASE, refine_lattice
from ase.build import fcc111

def printsummary(mol, details=None):
    s = f'{len(mol)} atoms, density = {mol.get_density()*1e-3:.3f} g/cm^3, box = {mol.lattice[0][0]:.3f}, {mol.lattice[1][1]:.3f}, {mol.lattice[2][2]:.3f}, formula = {mol.get_formula()}'
    if details:
        s+= f' #added molecules per species: {details["n_molecules"]}, mole fractions: {details["mole_fractions"]}'
    print(s)

def main():
    water = from_smiles('O')
    acetonitrile = from_smiles('CC#N')

    print('PURE LIQUID\n')

    print('pure liquid from approximate number of atoms and exact density (in g/cm^3), cubic box with auto-determined size')
    out = packmol(water, n_atoms=194, density=1.0)
    printsummary(out)
    out.write('water-1.xyz')

    print('pure liquid from approximate density (in g/cm^3) and an orthorhombic box')
    out = packmol(water, density=1.0, box_bounds=[0., 0., 0., 8., 12., 14.])
    printsummary(out)
    out.write('water-2.xyz')

    print('pure liquid with explicit number of molecules and exact density')
    out = packmol(water, n_molecules=64, density=1.0)
    printsummary(out)
    out.write('water-3.xyz')

    print('pure liquid with explicit number of molecules and box')
    out = packmol(water, n_molecules=64, box_bounds=[0., 0., 0., 12., 13., 14.])
    printsummary(out)
    out.write('water-4.xyz')

    print('water-5.xyz: pure liquid in non-orthorhombic box (requires AMS2022 or later)')
    # first place the molecules in a cube surrounding the desired lattice
    # then gradually change into the desired lattice using refine_lattice()
    # note that the molecules may become distorted by this procedure
    lattice = [[10., 2., -1.], [-5., 8., 0.], [0., -2., 11.]]
    temp_out = packmol(water, n_molecules=32, box_bounds=[
        0, 0, 0,
        max(lattice[i][0] for i in range(3))-min(lattice[i][0] for i in range(3)),
        max(lattice[i][1] for i in range(3))-min(lattice[i][1] for i in range(3)),
        max(lattice[i][2] for i in range(3))-min(lattice[i][2] for i in range(3))
    ])
    out = refine_lattice(temp_out, lattice=lattice)
    if out is not None:
        out.write('water-5.xyz')

    # MIXTURES
    x_water = 0.666                # mole fraction
    x_acetonitrile = 1-x_water     # mole fraction
    density = (x_water*1.0 + x_acetonitrile*0.76) / (x_water + x_acetonitrile)  # weighted average of pure component densities

    print(f'\nMIXTURES. x_water = {x_water:.3f}, x_acetonitrile = {x_acetonitrile:.3f}, target density = {density:.3f} g/cm^3\n')

    print('2-1 water-acetonitrile from approximate number of atoms and exact density (in g/cm^3), cubic box with auto-determined size')
    out, details = packmol(molecules=[water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], n_atoms=200, density=density, return_details=True)
    printsummary(out, details)
    out.write('water-acetonitrile-1.xyz')

    print('2-1 water-acetonitrile from approximate density (in g/cm^3) and box bounds')
    out, details = packmol(molecules=[water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], box_bounds=[0, 0, 0, 13.2, 13.2, 13.2], density=density, return_details=True)
    printsummary(out, details)
    out.write('water-acetonitrile-2.xyz')

    print('2-1 water-acetonitrile from explicit number of molecules and density, cubic box with auto-determined size')
    out, details = packmol(molecules=[water, acetonitrile], n_molecules=[32, 16], density=density, return_details=True)
    printsummary(out, details)
    out.write('water-acetonitrile-3.xyz')

    print('2-1 water-acetonitrile from explicit number of molecules and box')
    out = packmol(molecules=[water, acetonitrile], n_molecules=[32, 16], box_bounds=[0, 0, 0, 13.2, 13.2, 13.2])
    printsummary(out)
    out.write('water-acetonitrile-4.xyz')

    print('\nSOLID-LIQUID INTERFACES\n')
    slab = fromASE(fcc111('Al', size=(4,6,3), vacuum=15.0, orthogonal=True, periodic=True))

    print('water surrounding an Al slab, from an approximate density')
    out = packmol_on_slab(slab, water, density=1.0)
    printsummary(out)
    out.write('al-water-pure.xyz')

    print('2-1 water-acetonitrile mixture surrounding an Al slab, from mole fractions and an approximate density')
    out = packmol_on_slab(slab, [water, acetonitrile], mole_fractions=[x_water, x_acetonitrile], density=density)
    printsummary(out)
    out.write('al-water-acetonitrile.xyz')

if __name__ == '__main__':
    main()

