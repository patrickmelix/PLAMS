import sys
import os
from scm.plams import AMSJob, Settings

# Read from a previous result file
path = sys.argv[1]
job = AMSJob.load_external(path)
charges, types, patch = job.results.get_forcefield_params()
print(charges)
print(types)
print(patch)

# Add the atom types and charges to the molecule
mol = job.molecule
for i, at in enumerate(mol.atoms):
    at.properties.forcefield.charge = charges[i]
    at.properties.forcefield.type = types[i]

# Write a patch file
outfile = open("patch.dat", "w")
outfile.write(str(patch))
outfile.close()

# Set up new job settings
settings = Settings()
settings.input.ForceField.Type = "GAFF"
settings.input.ForceField.ForceFieldPatchFile = "patch.dat"
settings.input.ams.Task = "SinglePoint"

# Create the new job, and write the input file
job = AMSJob(molecule=mol, settings=settings)
text = job.get_input()
outfile = open("ams.in", "w")
outfile.write(text)
outfile.close()
