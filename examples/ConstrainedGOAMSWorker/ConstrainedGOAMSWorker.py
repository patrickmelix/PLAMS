#!/usr/bin/env amspython
# coding: utf-8

# ## Constrained geometry optimizations with the AMSWorker
# 
# The `AMSWorker` class allows geometry optimization of multiple molecules without the overhead of AMS start-up for each molecule.
# Here we want to optimize the geometries of three different molecules, with benzene as the common substructure,
# but wit different substituents.

import os
from scm.plams import Molecule
from scm.plams import plot_grid_molecules

path = os.path.join(os.environ["AMSRESOURCES"], "Molecules", "TestMols", "PLAMS")
filenames = ["tbut_benzene.in", "o_di_tbut_benzene.in", "tri_tbut_benzene.in"]
filenames = [os.path.join(path, fn) for fn in filenames]
molecules = [Molecule(fn) for fn in filenames]

plot_grid_molecules(molecules, molsPerRow=3)


# The structures are unoptimized, as the crowded geometry of the third structure demonstrates.

from scm.plams import plot_molecule
plot_molecule(molecules[2]);


# The geometry of the three structures can be optimized with the AMSWorker as follows.

from scm.plams import AMSWorker, Settings

# Create the general settings object and start upt the amsworker
settings = Settings()
settings.input.ForceField.Type = "UFF"
worker = AMSWorker(settings)

stackmol = Molecule()
for i, mol in enumerate(molecules):
    results = worker.GeometryOptimization("go%i"%(i), mol)
    stackmol += results.get_main_molecule()

plot_molecule(stackmol);


# We may prefer to perform the optimization while constraining the positions of the benzene carbon atoms, so that the benzene rings can be stacked directly on top of one another. The constraints need to be set by a call to the `SetConstraints()` method of the `AMSWorker`, and will apply to only a single geometry optimization.

stackmol = Molecule()
for i, mol in enumerate(molecules):
    # Set the constraints each time
    s = Settings()
    s.input.ams.Constraints.Atom = [1, 2, 3, 4, 5, 6]
    worker.SetConstraints(s, mol)

    results = worker.GeometryOptimization("constrained%i"%(i), mol)
    stackmol += results.get_main_molecule()

plot_molecule(stackmol);


# If we set contraints on one molecule, and then run a geometry optimization on another molecule, this will result in an error.

from scm.plams import JobError
from scm.plams import from_smiles
from scm.amspipe import AMSPipeError

s = Settings()
s.input.ams.Constraints.Atom = [1, 2, 3, 4, 5, 6]
worker.SetConstraints(s, mol)

try:
    results = worker.GeometryOptimization("water", from_smiles("O"))
except (JobError, AMSPipeError) as exc:
    print(str(exc))


# If we then run the geometry optimization again, the constraints will no longer apply.

results = worker.GeometryOptimization("water", from_smiles("O"))
plot_molecule(results.get_main_molecule());


worker.stop();

