import pytest

from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.unit_tests.test_basejob import DummySingleJob
from scm.plams.unit_tests.test_helpers import temp_file_path
from scm.plams.tools.job_analysis import JobAnalysis
from scm.plams.core.settings import Settings


class TestJobAnalysis:

    @pytest.fixture(scope="class")
    def dummy_single_jobs(self):
        # Generate dummy jobs for a selection of molecules and input settings
        smiles = ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", "Sys"]
        jobs = []
        for i, s in enumerate(smiles):
            sett = Settings()
            sett.input.task = "GeometryOptimization" if i % 2 else "SinglePoint"
            sett.input.ams.Properties.NormalModes = "Yes" if i % 3 else "No"
            if i < 5:
                sett.input.ADF.Basis.Type = "TZP"
                if i % 2:
                    sett.input.ADF.xc.gga = "pbe"
            else:
                sett.input.DFTB
            if s == "Sys":
                sett.input.System = [
                    "Ar 0.0000000000       0.0000000000       0.0000000000",
                    "Ar 1.6050000000       0.9266471820       2.6050000000",
                ]
                mol = None
            else:
                mol = from_smiles(s)
            jobs.append(DummySingleJob(wait=i / 100, molecule=mol, settings=sett))

        for j in jobs:
            j.run()
            j.ok()
        return jobs

    def test_init_with_jobs(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        assert len(ja.jobs) == 10

    def test_init_with_paths(self, dummy_single_jobs):
        ja = JobAnalysis(paths=[j.path for j in dummy_single_jobs])
        assert len(ja.jobs) == 10

    def test_default_fields(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.remove_field("Path")
        ja.add_molecule_fields()
        ja.add_timing_fields()

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | CPUTime | SysTime | ElapsedTime |
|--------------|------|-------|----------|---------|--------|---------|---------|-------------|
| plamsjob     | True | True  | None     | C2H6    | CC     | None    | None    | None        |
| plamsjob.002 | True | True  | None     | CH4     | C      | None    | None    | None        |
| plamsjob.003 | True | True  | None     | H2O     | O      | None    | None    | None        |
| plamsjob.004 | True | True  | None     | CH4O    | CO     | None    | None    | None        |
| plamsjob.005 | True | True  | None     | C3H8    | CCC    | None    | None    | None        |
| plamsjob.006 | True | True  | None     | C4H10   | CCCC   | None    | None    | None        |
| plamsjob.007 | True | True  | None     | C3H8O   | CCCO   | None    | None    | None        |
| plamsjob.008 | True | True  | None     | C6H14   | CCCCCC | None    | None    | None        |
| plamsjob.009 | True | True  | None     | C4H10O  | CCCOC  | None    | None    | None        |
| plamsjob.010 | True | True  | None     | None    | None   | None    | None    | None        |"""
        )

        ja.remove_timing_fields()
        ja.add_job_parent_fields()

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | ParentPath | ParentName |
|--------------|------|-------|----------|---------|--------|------------|------------|
| plamsjob     | True | True  | None     | C2H6    | CC     | None       | None       |
| plamsjob.002 | True | True  | None     | CH4     | C      | None       | None       |
| plamsjob.003 | True | True  | None     | H2O     | O      | None       | None       |
| plamsjob.004 | True | True  | None     | CH4O    | CO     | None       | None       |
| plamsjob.005 | True | True  | None     | C3H8    | CCC    | None       | None       |
| plamsjob.006 | True | True  | None     | C4H10   | CCCC   | None       | None       |
| plamsjob.007 | True | True  | None     | C3H8O   | CCCO   | None       | None       |
| plamsjob.008 | True | True  | None     | C6H14   | CCCCCC | None       | None       |
| plamsjob.009 | True | True  | None     | C4H10O  | CCCOC  | None       | None       |
| plamsjob.010 | True | True  | None     | None    | None   | None       | None       |"""
        )

    def test_add_remove_filter_fields(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.add_timing_fields()
        ja.add_field("Wait", lambda j: j.wait)
        ja.add_field("Output", lambda j: j.results.read_file("$JN.out")[:5])
        ja.remove_field("Path")

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | CPUTime | SysTime | ElapsedTime | Wait | Output |
|--------------|------|-------|----------|---------|--------|---------|---------|-------------|------|--------|
| plamsjob     | True | True  | None     | C2H6    | CC     | None    | None    | None        | 0.0  | Dummy  |
| plamsjob.002 | True | True  | None     | CH4     | C      | None    | None    | None        | 0.01 | Dummy  |
| plamsjob.003 | True | True  | None     | H2O     | O      | None    | None    | None        | 0.02 | Dummy  |
| plamsjob.004 | True | True  | None     | CH4O    | CO     | None    | None    | None        | 0.03 | Dummy  |
| plamsjob.005 | True | True  | None     | C3H8    | CCC    | None    | None    | None        | 0.04 | Dummy  |
| plamsjob.006 | True | True  | None     | C4H10   | CCCC   | None    | None    | None        | 0.05 | Dummy  |
| plamsjob.007 | True | True  | None     | C3H8O   | CCCO   | None    | None    | None        | 0.06 | Dummy  |
| plamsjob.008 | True | True  | None     | C6H14   | CCCCCC | None    | None    | None        | 0.07 | Dummy  |
| plamsjob.009 | True | True  | None     | C4H10O  | CCCOC  | None    | None    | None        | 0.08 | Dummy  |
| plamsjob.010 | True | True  | None     | None    | None   | None    | None    | None        | 0.09 | Dummy  |"""
        )

        ja.remove_empty_fields()

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | Formula | Smiles | Wait | Output |
|--------------|------|-------|---------|--------|------|--------|
| plamsjob     | True | True  | C2H6    | CC     | 0.0  | Dummy  |
| plamsjob.002 | True | True  | CH4     | C      | 0.01 | Dummy  |
| plamsjob.003 | True | True  | H2O     | O      | 0.02 | Dummy  |
| plamsjob.004 | True | True  | CH4O    | CO     | 0.03 | Dummy  |
| plamsjob.005 | True | True  | C3H8    | CCC    | 0.04 | Dummy  |
| plamsjob.006 | True | True  | C4H10   | CCCC   | 0.05 | Dummy  |
| plamsjob.007 | True | True  | C3H8O   | CCCO   | 0.06 | Dummy  |
| plamsjob.008 | True | True  | C6H14   | CCCCCC | 0.07 | Dummy  |
| plamsjob.009 | True | True  | C4H10O  | CCCOC  | 0.08 | Dummy  |
| plamsjob.010 | True | True  | None    | None   | 0.09 | Dummy  |"""
        )

        ja.remove_uniform_fields()

        assert (
            ja.to_table()
            == """\
| Name         | Formula | Smiles | Wait |
|--------------|---------|--------|------|
| plamsjob     | C2H6    | CC     | 0.0  |
| plamsjob.002 | CH4     | C      | 0.01 |
| plamsjob.003 | H2O     | O      | 0.02 |
| plamsjob.004 | CH4O    | CO     | 0.03 |
| plamsjob.005 | C3H8    | CCC    | 0.04 |
| plamsjob.006 | C4H10   | CCCC   | 0.05 |
| plamsjob.007 | C3H8O   | CCCO   | 0.06 |
| plamsjob.008 | C6H14   | CCCCCC | 0.07 |
| plamsjob.009 | C4H10O  | CCCOC  | 0.08 |
| plamsjob.010 | None    | None   | 0.09 |"""
        )

        ja.remove_uniform_fields(tol=0.1, ignore_empty=True)
        ja.filter_fields(lambda vals: all([not v or "H" not in v for v in vals]))

        assert (
            ja.to_table()
            == """\
| Formula |
|---------|
| C2H6    |
| CH4     |
| H2O     |
| CH4O    |
| C3H8    |
| C4H10   |
| C3H8O   |
| C6H14   |
| C4H10O  |
| None    |"""
        )

    def test_field_groups(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.add_timing_fields()
        ja.add_field("Wait", lambda j: j.wait, group="w")

        assert ja.field_groups == {
            "job_info": ["Path", "Name", "OK", "Check", "ErrorMsg"],
            "mol": ["Formula", "Smiles"],
            "timing": ["CPUTime", "SysTime", "ElapsedTime"],
            "w": ["Wait"],
        }

        ja.remove_field_group("mol")
        ja.remove_field_group("w")

        assert ja.field_groups == {
            "job_info": ["Path", "Name", "OK", "Check", "ErrorMsg"],
            "timing": ["CPUTime", "SysTime", "ElapsedTime"],
        }

    def test_reorder_fields(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.add_timing_fields()
        ja.add_field("Wait", lambda j: j.wait, group="w")

        assert ja.field_names == [
            "Path",
            "Name",
            "OK",
            "Check",
            "ErrorMsg",
            "Formula",
            "Smiles",
            "CPUTime",
            "SysTime",
            "ElapsedTime",
            "Wait",
        ]

        ja.reorder_fields(["Name", "Wait"])

        assert ja.field_names == [
            "Name",
            "Wait",
            "Path",
            "OK",
            "Check",
            "ErrorMsg",
            "Formula",
            "Smiles",
            "CPUTime",
            "SysTime",
            "ElapsedTime",
        ]

        ja.sort_fields(lambda k: str(k))

        assert ja.field_names == [
            "CPUTime",
            "Check",
            "ElapsedTime",
            "ErrorMsg",
            "Formula",
            "Name",
            "OK",
            "Path",
            "Smiles",
            "SysTime",
            "Wait",
        ]

    def test_settings_fields(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.add_settings_input_fields()
        ja.add_settings_field(("runscript", "shebang"))
        ja.remove_field("Path")

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Input.Task           | Input.Adf.Basis.Type | Input.Ams.Properties.Normalmodes | Input.Adf.Xc.Gga | Runscript.Shebang |
|--------------|------|-------|----------|---------|--------|----------------------|----------------------|----------------------------------|------------------|-------------------|
| plamsjob     | True | True  | None     | C2H6    | CC     | SinglePoint          | TZP                  | No                               | None             | #!/bin/sh         |
| plamsjob.002 | True | True  | None     | CH4     | C      | GeometryOptimization | TZP                  | Yes                              | pbe              | #!/bin/sh         |
| plamsjob.003 | True | True  | None     | H2O     | O      | SinglePoint          | TZP                  | Yes                              | None             | #!/bin/sh         |
| plamsjob.004 | True | True  | None     | CH4O    | CO     | GeometryOptimization | TZP                  | No                               | pbe              | #!/bin/sh         |
| plamsjob.005 | True | True  | None     | C3H8    | CCC    | SinglePoint          | TZP                  | Yes                              | None             | #!/bin/sh         |
| plamsjob.006 | True | True  | None     | C4H10   | CCCC   | GeometryOptimization | None                 | Yes                              | None             | #!/bin/sh         |
| plamsjob.007 | True | True  | None     | C3H8O   | CCCO   | SinglePoint          | None                 | No                               | None             | #!/bin/sh         |
| plamsjob.008 | True | True  | None     | C6H14   | CCCCCC | GeometryOptimization | None                 | Yes                              | None             | #!/bin/sh         |
| plamsjob.009 | True | True  | None     | C4H10O  | CCCOC  | SinglePoint          | None                 | Yes                              | None             | #!/bin/sh         |
| plamsjob.010 | True | True  | None     | None    | None   | GeometryOptimization | None                 | No                               | None             | #!/bin/sh         |"""
        )

        ja.remove_settings_fields()
        ja.add_settings_input_fields(include_system_block=True)
        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | Input.Task           | Input.Adf.Basis.Type | Input.Ams.Properties.Normalmodes | Input.Adf.Xc.Gga | Input.System.0                                        | Input.System.1                                        |
|--------------|------|-------|----------|---------|--------|----------------------|----------------------|----------------------------------|------------------|-------------------------------------------------------|-------------------------------------------------------|
| plamsjob     | True | True  | None     | C2H6    | CC     | SinglePoint          | TZP                  | No                               | None             | None                                                  | None                                                  |
| plamsjob.002 | True | True  | None     | CH4     | C      | GeometryOptimization | TZP                  | Yes                              | pbe              | None                                                  | None                                                  |
| plamsjob.003 | True | True  | None     | H2O     | O      | SinglePoint          | TZP                  | Yes                              | None             | None                                                  | None                                                  |
| plamsjob.004 | True | True  | None     | CH4O    | CO     | GeometryOptimization | TZP                  | No                               | pbe              | None                                                  | None                                                  |
| plamsjob.005 | True | True  | None     | C3H8    | CCC    | SinglePoint          | TZP                  | Yes                              | None             | None                                                  | None                                                  |
| plamsjob.006 | True | True  | None     | C4H10   | CCCC   | GeometryOptimization | None                 | Yes                              | None             | None                                                  | None                                                  |
| plamsjob.007 | True | True  | None     | C3H8O   | CCCO   | SinglePoint          | None                 | No                               | None             | None                                                  | None                                                  |
| plamsjob.008 | True | True  | None     | C6H14   | CCCCCC | GeometryOptimization | None                 | Yes                              | None             | None                                                  | None                                                  |
| plamsjob.009 | True | True  | None     | C4H10O  | CCCOC  | SinglePoint          | None                 | Yes                              | None             | None                                                  | None                                                  |
| plamsjob.010 | True | True  | None     | None    | None   | GeometryOptimization | None                 | No                               | None             | Ar 0.0000000000       0.0000000000       0.0000000000 | Ar 1.6050000000       0.9266471820       2.6050000000 |"""
        )

    def test_get_set_del_item(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        del ja["Path"]
        del ja["Check"]
        ja["OK"] = lambda j: "Yes" if j.ok() else "No"
        ja["Id"] = lambda j: j.name.split(".")[-1]

        assert ja["Name"] == [
            "plamsjob",
            "plamsjob.002",
            "plamsjob.003",
            "plamsjob.004",
            "plamsjob.005",
            "plamsjob.006",
            "plamsjob.007",
            "plamsjob.008",
            "plamsjob.009",
            "plamsjob.010",
        ]
        assert ja["Smiles"] == ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", None]
        assert ja["OK"] == ["Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"]
        assert ja["Id"] == ["plamsjob", "002", "003", "004", "005", "006", "007", "008", "009", "010"]

    def test_to_table(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.remove_field("Path")
        ja.add_settings_input_fields()
        ja.remove_empty_fields()
        ja.remove_uniform_fields()

        assert (
            ja.to_table(max_col_width=10, max_rows=6)
            == """\
| Name          | Formula | Smiles | Input.Task    | Input.Adf.... | Input.Ams.... | Input.Adf.... |
|---------------|---------|--------|---------------|---------------|---------------|---------------|
| plamsjob      | C2H6    | CC     | SinglePoin... | TZP           | No            | None          |
| plamsjob.0... | CH4     | C      | GeometryOp... | TZP           | Yes           | pbe           |
| plamsjob.0... | H2O     | O      | SinglePoin... | TZP           | Yes           | None          |
| ...           | ...     | ...    | ...           | ...           | ...           | ...           |
| plamsjob.0... | C6H14   | CCCCCC | GeometryOp... | None          | Yes           | None          |
| plamsjob.0... | C4H10O  | CCCOC  | SinglePoin... | None          | Yes           | None          |
| plamsjob.0... | None    | None   | GeometryOp... | None          | No            | None          |"""
        )

    def test_to_csv(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.remove_field("Path")
        ja.add_settings_input_fields()
        ja.remove_empty_fields()
        ja.remove_uniform_fields()

        with temp_file_path(".csv") as tfp:
            ja.to_csv_file(tfp)
            with open(tfp) as tf:
                csv = tf.read()

        assert (
            csv
            == """\
Name,Formula,Smiles,Input.Task,Input.Adf.Basis.Type,Input.Ams.Properties.Normalmodes,Input.Adf.Xc.Gga
plamsjob,C2H6,CC,SinglePoint,TZP,No,
plamsjob.002,CH4,C,GeometryOptimization,TZP,Yes,pbe
plamsjob.003,H2O,O,SinglePoint,TZP,Yes,
plamsjob.004,CH4O,CO,GeometryOptimization,TZP,No,pbe
plamsjob.005,C3H8,CCC,SinglePoint,TZP,Yes,
plamsjob.006,C4H10,CCCC,GeometryOptimization,,Yes,
plamsjob.007,C3H8O,CCCO,SinglePoint,,No,
plamsjob.008,C6H14,CCCCCC,GeometryOptimization,,Yes,
plamsjob.009,C4H10O,CCCOC,SinglePoint,,Yes,
plamsjob.010,,,GeometryOptimization,,No,
"""
        )

    def test_to_dataframe(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.remove_field("Path")
        ja.add_settings_input_fields()
        ja.remove_empty_fields()
        ja.remove_uniform_fields()

        df = ja.to_dataframe()

        assert df.shape == (10, 7)
        assert df.columns.to_list() == [
            "Name",
            "Formula",
            "Smiles",
            "Input.Task",
            "Input.Adf.Basis.Type",
            "Input.Ams.Properties.Normalmodes",
            "Input.Adf.Xc.Gga",
        ]
        assert df.Formula.to_list() == [
            "C2H6",
            "CH4",
            "H2O",
            "CH4O",
            "C3H8",
            "C4H10",
            "C3H8O",
            "C6H14",
            "C4H10O",
            None,
        ]
