import pytest
import shutil

from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.core.jobmanager import JobManager
from scm.plams.unit_tests.test_basejob import DummySingleJob
from scm.plams.unit_tests.test_helpers import temp_file_path, skip_if_no_scm_pisa, skip_if_no_scm_libbase
from scm.plams.tools.job_analysis import JobAnalysis
from scm.plams.core.settings import Settings, JobManagerSettings


class TestJobAnalysis:

    @pytest.fixture(scope="class")
    def dummy_single_jobs(self):
        # Generate dummy jobs for a selection of molecules and input settings
        smiles = ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", "Sys"]
        jobs = []
        for i, s in enumerate(smiles):
            sett = Settings()
            sett.input.ams.task = "GeometryOptimization" if i % 2 else "SinglePoint"
            sett.input.ams.Properties.NormalModes = "True" if i % 3 else "False"
            if i < 5:
                sett.input.ADF.Basis.Type = "TZP"
                if i % 2:
                    sett.input.ADF.xc.gga = "pbe"
            else:
                sett.input.DFTB
            if s == "Sys":
                sett.input.ams.System.Atoms = [
                    "Ar 0.0000000000       0.0000000000       0.0000000000",
                    "Ar 1.6050000000       0.9266471820       2.6050000000",
                ]
                mol = None
            else:
                mol = from_smiles(s)
            jobs.append(DummySingleJob(wait=i / 100, molecule=mol, settings=sett, name="dummyjob"))

        for j in jobs:
            j.run()
            j.ok()

        yield jobs

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
| dummyjob     | True | True  | None     | C2H6    | CC     | None    | None    | None        |
| dummyjob.002 | True | True  | None     | CH4     | C      | None    | None    | None        |
| dummyjob.003 | True | True  | None     | H2O     | O      | None    | None    | None        |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | None    | None    | None        |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | None    | None    | None        |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None    | None    | None        |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None    | None    | None        |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None    | None    | None        |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None    | None    | None        |
| dummyjob.010 | True | True  | None     | None    | None   | None    | None    | None        |"""
        )

        ja.remove_timing_fields()
        ja.add_job_parent_fields()

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | ParentPath | ParentName |
|--------------|------|-------|----------|---------|--------|------------|------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | None       | None       |
| dummyjob.002 | True | True  | None     | CH4     | C      | None       | None       |
| dummyjob.003 | True | True  | None     | H2O     | O      | None       | None       |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | None       | None       |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | None       | None       |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None       | None       |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None       | None       |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None       | None       |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None       | None       |
| dummyjob.010 | True | True  | None     | None    | None   | None       | None       |"""
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
| dummyjob     | True | True  | None     | C2H6    | CC     | None    | None    | None        | 0.0  | Dummy  |
| dummyjob.002 | True | True  | None     | CH4     | C      | None    | None    | None        | 0.01 | Dummy  |
| dummyjob.003 | True | True  | None     | H2O     | O      | None    | None    | None        | 0.02 | Dummy  |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | None    | None    | None        | 0.03 | Dummy  |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | None    | None    | None        | 0.04 | Dummy  |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None    | None    | None        | 0.05 | Dummy  |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None    | None    | None        | 0.06 | Dummy  |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None    | None    | None        | 0.07 | Dummy  |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None    | None    | None        | 0.08 | Dummy  |
| dummyjob.010 | True | True  | None     | None    | None   | None    | None    | None        | 0.09 | Dummy  |"""
        )

        ja.remove_empty_fields()

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | Formula | Smiles | Wait | Output |
|--------------|------|-------|---------|--------|------|--------|
| dummyjob     | True | True  | C2H6    | CC     | 0.0  | Dummy  |
| dummyjob.002 | True | True  | CH4     | C      | 0.01 | Dummy  |
| dummyjob.003 | True | True  | H2O     | O      | 0.02 | Dummy  |
| dummyjob.004 | True | True  | CH4O    | CO     | 0.03 | Dummy  |
| dummyjob.005 | True | True  | C3H8    | CCC    | 0.04 | Dummy  |
| dummyjob.006 | True | True  | C4H10   | CCCC   | 0.05 | Dummy  |
| dummyjob.007 | True | True  | C3H8O   | CCCO   | 0.06 | Dummy  |
| dummyjob.008 | True | True  | C6H14   | CCCCCC | 0.07 | Dummy  |
| dummyjob.009 | True | True  | C4H10O  | CCCOC  | 0.08 | Dummy  |
| dummyjob.010 | True | True  | None    | None   | 0.09 | Dummy  |"""
        )

        ja.remove_uniform_fields()

        assert (
            ja.to_table()
            == """\
| Name         | Formula | Smiles | Wait |
|--------------|---------|--------|------|
| dummyjob     | C2H6    | CC     | 0.0  |
| dummyjob.002 | CH4     | C      | 0.01 |
| dummyjob.003 | H2O     | O      | 0.02 |
| dummyjob.004 | CH4O    | CO     | 0.03 |
| dummyjob.005 | C3H8    | CCC    | 0.04 |
| dummyjob.006 | C4H10   | CCCC   | 0.05 |
| dummyjob.007 | C3H8O   | CCCO   | 0.06 |
| dummyjob.008 | C6H14   | CCCCCC | 0.07 |
| dummyjob.009 | C4H10O  | CCCOC  | 0.08 |
| dummyjob.010 | None    | None   | 0.09 |"""
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

    def test_filter_jobs(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.add_timing_fields()
        ja.add_field("Wait", lambda j: j.wait)
        ja.add_field("Output", lambda j: j.results.read_file("$JN.out")[:5])
        ja.remove_field("Path")

        ja.filter_jobs(lambda d: d["Formula"] is None or "O" in d["Smiles"])

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | CPUTime | SysTime | ElapsedTime | Wait | Output |
|--------------|------|-------|----------|---------|--------|---------|---------|-------------|------|--------|
| dummyjob     | True | True  | None     | C2H6    | CC     | None    | None    | None        | 0.0  | Dummy  |
| dummyjob.002 | True | True  | None     | CH4     | C      | None    | None    | None        | 0.01 | Dummy  |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | None    | None    | None        | 0.04 | Dummy  |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None    | None    | None        | 0.05 | Dummy  |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None    | None    | None        | 0.07 | Dummy  |"""
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
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | RunscriptShebang |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | #!/bin/sh        |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | #!/bin/sh        |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | #!/bin/sh        |"""
        )

        ja.remove_settings_fields()
        ja.add_settings_input_fields(include_system_block=True)

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | InputAmsSystemAtoms0                                  | InputAmsSystemAtoms1                                  |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|-------------------------------------------------------|-------------------------------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | None                                                  | None                                                  |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | None                                                  | None                                                  |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | Ar 0.0000000000       0.0000000000       0.0000000000 | Ar 1.6050000000       0.9266471820       2.6050000000 |"""
        )

    def test_get_set_del_item(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        del ja["Path"]
        del ja["Check"]
        ja["OK"] = lambda j: "Yes" if j.ok() else "No"
        ja["Id"] = lambda j: j.name.split(".")[-1]

        assert ja["Name"] == [
            "dummyjob",
            "dummyjob.002",
            "dummyjob.003",
            "dummyjob.004",
            "dummyjob.005",
            "dummyjob.006",
            "dummyjob.007",
            "dummyjob.008",
            "dummyjob.009",
            "dummyjob.010",
        ]
        assert ja["Smiles"] == ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", None]
        assert ja["OK"] == ["Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"]
        assert ja["Id"] == ["dummyjob", "002", "003", "004", "005", "006", "007", "008", "009", "010"]
        with pytest.raises(KeyError):
            ja["Foo"]
        with pytest.raises(KeyError):
            del ja["Bar"]

    def test_get_set_del_attributes(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        del ja.Path
        del ja.Check

        ja.OK = lambda j: "Yes" if j.ok() else "No"
        ja.Id = lambda j: j.name.split(".")[-1]

        assert ja.Name == [
            "dummyjob",
            "dummyjob.002",
            "dummyjob.003",
            "dummyjob.004",
            "dummyjob.005",
            "dummyjob.006",
            "dummyjob.007",
            "dummyjob.008",
            "dummyjob.009",
            "dummyjob.010",
        ]
        assert ja.Smiles == ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", None]
        assert ja.OK == ["Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"]
        assert ja.Id == ["dummyjob", "002", "003", "004", "005", "006", "007", "008", "009", "010"]
        with pytest.raises(AttributeError):
            ja.Foo
        with pytest.raises(AttributeError):
            del ja.Bar

    def test_to_table(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.remove_field("Path")
        ja.add_settings_field(("Input", "AMS", "Properties", "NormalModes"))
        ja.remove_empty_fields()
        ja.remove_uniform_fields()

        assert (
            ja.to_table(max_col_width=10, max_rows=6)
            == """\
| Name          | Formula | Smiles | InputAmsPr... |
|---------------|---------|--------|---------------|
| dummyjob      | C2H6    | CC     | False         |
| dummyjob.0... | CH4     | C      | True          |
| dummyjob.0... | H2O     | O      | True          |
| ...           | ...     | ...    | ...           |
| dummyjob.0... | C6H14   | CCCCCC | True          |
| dummyjob.0... | C4H10O  | CCCOC  | True          |
| dummyjob.0... | None    | None   | False         |"""
        )

    def test_to_csv(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.remove_field("Path")
        ja.add_settings_field(("Input", "AMS", "Properties", "NormalModes"))
        ja.remove_empty_fields()
        ja.remove_uniform_fields()

        with temp_file_path(".csv") as tfp:
            ja.to_csv_file(tfp)
            with open(tfp) as tf:
                csv = tf.read()

        assert (
            csv
            == """\
Name,Formula,Smiles,InputAmsPropertiesNormalmodes
dummyjob,C2H6,CC,False
dummyjob.002,CH4,C,True
dummyjob.003,H2O,O,True
dummyjob.004,CH4O,CO,False
dummyjob.005,C3H8,CCC,True
dummyjob.006,C4H10,CCCC,True
dummyjob.007,C3H8O,CCCO,False
dummyjob.008,C6H14,CCCCCC,True
dummyjob.009,C4H10O,CCCOC,True
dummyjob.010,,,False
"""
        )

    def test_to_dataframe(self, dummy_single_jobs):
        try:
            import pandas  # noqa F401
        except ImportError:
            pytest.skip("Skipping test as cannot find pandas package.")

        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.remove_field("Path")
        ja.add_settings_field(("Input", "AMS", "Properties", "NormalModes"))
        ja.remove_empty_fields()
        ja.remove_uniform_fields()

        df = ja.to_dataframe()

        assert df.shape == (10, 5)
        assert df.columns.to_list() == [
            "Name",
            "Formula",
            "Smiles",
            "InputTask",
            "InputAmsPropertiesNormalmodes",
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


class TestJobAnalysisWithPisa(TestJobAnalysis):

    @pytest.fixture(scope="class")
    def dummy_single_jobs(self):
        skip_if_no_scm_pisa()

        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB, ADF

        # Generate dummy jobs for a selection of molecules and input settings
        smiles = ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", "Sys"]
        jobs = []
        for i, s in enumerate(smiles):
            sett = AMS()
            sett.Task = "GeometryOptimization" if i % 2 else "SinglePoint"
            sett.Properties.NormalModes = "Yes" if i % 3 else "No"
            if i < 5:
                sett.Engine = ADF()
                sett.Engine.Basis.Type = "TZP"
                if i % 2:
                    sett.Engine.XC.GGA = "pbe"
            else:
                sett.Engine = DFTB()
            if s == "Sys":
                sett.System.Atoms = [
                    "Ar 0.0000000000       0.0000000000       0.0000000000",
                    "Ar 1.6050000000       0.9266471820       2.6050000000",
                ]
                mol = None
            else:
                mol = from_smiles(s)

            if i > 5:
                sett = Settings({"input": sett})

            jobs.append(DummySingleJob(wait=i / 100, molecule=mol, settings=sett, name="dummyjob"))

        jm = JobManager(JobManagerSettings())
        for j in jobs:
            j.run(jobmanager=jm)
            j.ok()
        yield jobs

        jm._clean()
        shutil.rmtree(jm.workdir)

    def test_settings_fields(self, dummy_single_jobs):
        ja = JobAnalysis(jobs=dummy_single_jobs)
        ja.add_molecule_fields()
        ja.add_settings_input_fields()
        ja.add_settings_field(("runscript", "shebang"))
        ja.remove_field("Path")

        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | RunscriptShebang |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | #!/bin/sh        |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | #!/bin/sh        |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | #!/bin/sh        |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | #!/bin/sh        |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | #!/bin/sh        |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | #!/bin/sh        |"""
        )

        ja.remove_settings_fields()
        ja.add_settings_input_fields(include_system_block=True)

        # N.B. small discrepancy in system block atoms column name
        assert (
            ja.to_table()
            == """\
| Name         | OK   | Check | ErrorMsg | Formula | Smiles | InputAdfBasisType | InputAmsTask         | InputAmsPropertiesNormalmodes | InputAdfXcGga | InputAmsSystem0Atoms_10                               | InputAmsSystem0Atoms_11                               |
|--------------|------|-------|----------|---------|--------|-------------------|----------------------|-------------------------------|---------------|-------------------------------------------------------|-------------------------------------------------------|
| dummyjob     | True | True  | None     | C2H6    | CC     | TZP               | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.002 | True | True  | None     | CH4     | C      | TZP               | GeometryOptimization | True                          | pbe           | None                                                  | None                                                  |
| dummyjob.003 | True | True  | None     | H2O     | O      | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.004 | True | True  | None     | CH4O    | CO     | TZP               | GeometryOptimization | False                         | pbe           | None                                                  | None                                                  |
| dummyjob.005 | True | True  | None     | C3H8    | CCC    | TZP               | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.006 | True | True  | None     | C4H10   | CCCC   | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.007 | True | True  | None     | C3H8O   | CCCO   | None              | SinglePoint          | False                         | None          | None                                                  | None                                                  |
| dummyjob.008 | True | True  | None     | C6H14   | CCCCCC | None              | GeometryOptimization | True                          | None          | None                                                  | None                                                  |
| dummyjob.009 | True | True  | None     | C4H10O  | CCCOC  | None              | SinglePoint          | True                          | None          | None                                                  | None                                                  |
| dummyjob.010 | True | True  | None     | None    | None   | None              | GeometryOptimization | False                         | None          | Ar 0.0000000000       0.0000000000       0.0000000000 | Ar 1.6050000000       0.9266471820       2.6050000000 |"""
        )


class TestJobAnalysisWithChemicalSystem(TestJobAnalysis):

    @pytest.fixture(scope="class")
    def dummy_single_jobs(self):
        # Generate dummy jobs for a selection of molecules and input settings
        skip_if_no_scm_libbase()
        from scm.utils.conversions import plams_molecule_to_chemsys

        smiles = ["CC", "C", "O", "CO", "CCC", "CCCC", "CCCO", "CCCCCC", "CCCOC", "Sys"]
        jobs = []
        for i, s in enumerate(smiles):
            sett = Settings()
            sett.input.ams.task = "GeometryOptimization" if i % 2 else "SinglePoint"
            sett.input.ams.Properties.NormalModes = "True" if i % 3 else "False"
            if i < 5:
                sett.input.ADF.Basis.Type = "TZP"
                if i % 2:
                    sett.input.ADF.xc.gga = "pbe"
            else:
                sett.input.DFTB
            if s == "Sys":
                sett.input.ams.System.Atoms = [
                    "Ar 0.0000000000       0.0000000000       0.0000000000",
                    "Ar 1.6050000000       0.9266471820       2.6050000000",
                ]
                mol = None
            else:
                mol = from_smiles(s)
                mol = plams_molecule_to_chemsys(mol)

            jobs.append(DummySingleJob(wait=i / 100, molecule=mol, settings=sett, name="dummyjob"))

        jm = JobManager(JobManagerSettings())
        for j in jobs:
            j.run(jobmanager=jm)
            j.ok()

        yield jobs

        jm._clean()
        shutil.rmtree(jm.workdir)
