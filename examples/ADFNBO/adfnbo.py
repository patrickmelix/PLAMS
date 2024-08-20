from scm.plams import Settings, Molecule
from scm.plams.recipes.adfnbo import ADFNBOJob

mol = Molecule("methane.xyz")

s = Settings()
s.input.AMS.Task = "SinglePoint"
s.input.ADF.basis.type = "DZP"
s.input.ADF.xc.lda = "SCF VWN"
s.input.ADF.relativity.level = "scalar"
s.adfnbo = ["write", "spherical", "fock"]

j = ADFNBOJob(molecule=mol, settings=s)
r = j.run()

lines = r.get_output_chunk(begin="NATURAL BOND ORBITALS (Summary):", end="Charge unit", inc_begin=True, inc_end=True)
for line in lines:
    print(line)
