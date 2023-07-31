from scm.plams.recipes.adfcosmorsconformers import ADFCOSMORSConfJob, ADFCOSMORSConfFilter
from scm.plams import Molecule, from_smiles, init, finish, Settings
from scm.conformers import ConformersJob

init()

mol = from_smiles("CC(=O)O")

conf_sett = Settings()
conf_sett.input.AMS.Generator.RDKit
conf_sett.input.AMS.Generator.RDKit.InitialNConformers = 7
conf_job = ConformersJob(name="conformers_uff", molecule=mol, settings=conf_sett)

dftb_sett = Settings()
dftb_sett.input.AMS.Task = "Optimize"
dftb_sett.input.DFTB

fil1 = ADFCOSMORSConfFilter(5,10)
fil2 = ADFCOSMORSConfFilter(2,10)

a = ADFCOSMORSConfJob(
    mol,
    conf_gen     = conf_job,
    first_filter = fil1,
    additional   = [(dftb_sett,fil2)],
    coskf_name   = "acetic_acid",
    coskf_dir    = "test_coskfs"
    )
a.run()

finish()
