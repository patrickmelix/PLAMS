from scm.conformers import ConformersJob, ConformersResults

# This example shows how to use the AMS's Conformers tool via PLAMS

def print_conformers(conformers: ConformersResults):
   for i, conf in enumerate(conformers):
      print(f"Conformer {i}: Energy = {conf.properties.energy_kcalmol} [kcal/mol]")
      print(conf)

ethanol = from_smiles('OCC')

# Simple example: generate conformers for ethanol using default settings
# - default Task: Generate
# - default Generator Method: RDKit 
# - default Engine: ForceField (with UFF)

job = ConformersJob(name='conformers_generation', molecule=ethanol)
job.run()

print("Conformers generated using RDKit and UFF:")
print_conformers(job.results.conformers)

# Re-optimize the conformers generated in the previous steps using the GFN1-xTB engine:

sett = Settings()

# In the 'ams' part of the settings input you can specify 
# the input options for the Conformers tool, which are 
# described in the Conformers tool user manual.

sett.input.ams.Task = "Optimize"
sett.input.ams.InputConformersSet = job.results

# You can specify the engine to be used (and the engine options) like you would 
# for an AMSJob. See the AMSJob documentation for more details.

sett.input.DFTB.Model = "GFN1-xTB"

# Note: here we do not specify the input molecule because we are passing the results 
# of a previous ConformersJob (ConformersResults) via the "InputConformersSet" input.

reoptimize_job = ConformersJob(name='optimize_conformers', settings=sett)
reoptimize_job.run()

print("Conformers re-optimized using the more accurate GFN1-xTB method:")
print_conformers(reoptimize_job.results.conformers)
