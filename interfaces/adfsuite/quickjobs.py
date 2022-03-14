import os
from ...mol.molecule import Molecule
from .amsworker import AMSWorker
from ...core.settings import Settings

__all__ = ['preoptimize', 'refine_density']

def preoptimize(molecule: Molecule, model:str='UFF', settings:Settings=None, nproc:int=1, maxiterations:int=100):
    """
        Returns an optimized Molecule (or list of optimized molecules)

        molecule: Molecule or list of Molecules
            Molecule to optimize
        model: str
            Shorthand for some model, e.g. 'UFF'
        settings: Settings
            Custom engine settings (overrides ``model``)
        nproc: int
            Number of processes
        maxiterations: int
            Maximum number of iterations for the geometry optimization.
    """
    my_settings = _get_quick_settings(model, settings, nproc)
    single_molecule = isinstance(molecule, Molecule)

    input_molecules = [molecule] if single_molecule else molecule
    output_molecules = []
    with AMSWorker(my_settings) as worker:
        for i,mol in enumerate(input_molecules):
            results = worker.GeometryOptimization(name=f'preoptimization_{i}', molecule=mol, maxiterations=maxiterations, pretendconverged=True)
            output_molecules.append(results.get_main_molecule())

    if single_molecule:
        return output_molecules[0]
    else:
        return output_molecules

def refine_density(molecule: Molecule, density:float, step_size=0.05, model:str='UFF', settings:Settings=None, nproc:int=1, maxiterations:int=100):
    """

        Performs a series of geometry optimizations with densities approaching
        ``density``. This can be useful if you want to compress a system to a
        given density, but cannot just use apply_strain() (because
        apply_strain() also scales the bond lengths).

        This function can be useful if for example packmol does not succeed to
        pack molecules with the desired density. Packmol can then generate a
        structure with a lower density, and this function can be used to
        increase the density to the desired value.

        Returns: a Molecule with the requested density.

        molecule: Molecule
            The molecule must have a 3D lattice
        density: float
            Target density in g/cm^3
        step_size: float
            Step size for the density (in g/cm^3). Set step_size to a large number to only use 1 step.
        model: str
            e.g. 'UFF'
        settings: Settings
            Engine settings (overrides ``model``)
        maxiterations: int
            maximum number of iterations for the geometry optimization.


    """

    assert(len(molecule.lattice) == 3)
    tolerance = 1e-3
    my_settings = _get_quick_settings(model, settings, nproc)
    current_density = molecule.get_density() * 1e-3
    output_molecule = molecule.copy()
    counter = 0
    with AMSWorker(my_settings) as worker:
        while (current_density < density-tolerance or current_density > density+tolerance):
            counter += 1
            if current_density < density-step_size:
                new_density = current_density + step_size
            elif current_density > density+step_size:
                new_density = current_density - step_size
            else:
                new_density = density

            strain = (current_density/new_density)**(1/3.0)
            strain = strain - 1.0
            output_molecule.apply_strain([strain, strain, strain, 0, 0, 0], voigt_form=True)

            results = worker.GeometryOptimization(name=f'preoptimization_{counter}', molecule=output_molecule, maxiterations=maxiterations, pretendconverged=True)

            output_molecule = results.get_main_molecule()
            current_density = output_molecule.get_density() * 1e-3

    return output_molecule


def model_to_settings(model:str):
    """
    Returns Settings
    """
    settings = Settings()
    if model == 'UFF':
        settings.input.ForceField.Type = 'UFF'
    elif model == 'GFNFF':
        settings.input.GFNFF
    elif model == 'ANI-2x':
        settings.input.MLPotential.Model = 'ANI-2x'
    else:
        raise ValueError("Unknown model: {}".format(model))

    return settings

def _get_quick_settings(model, settings, nproc):
    if settings is None:
        my_settings = model_to_settings(model)
    else:
        my_settings = settings.copy()

    if nproc:
        my_settings.runscript.nproc = nproc

    return my_settings

