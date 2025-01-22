from typing import Optional, Sequence, Union, Dict, List, Callable, Any, NamedTuple, Tuple, Hashable
import os
import csv
from pathlib import Path
import numpy as np
from numbers import Number

from scm.plams.core.basejob import Job, SingleJob
from scm.plams.core.settings import Settings
from scm.plams.interfaces.adfsuite.ams import AMSJob
from scm.plams.core.errors import PlamsError
from scm.plams.core.functions import requires_optional_package, load, config
from scm.plams.tools.table_formatter import format_in_table
from scm.plams.interfaces.molecule.rdkit import to_smiles
from scm.plams.mol.molecule import Molecule
from scm.plams.interfaces.adfsuite.inputparser import InputParserFacade

try:
    from scm.libbase import UnifiedChemicalSystem as ChemicalSystem
    from scm.utils.conversions import chemsys_to_plams_molecule

    _has_scm_libbase = True
except ImportError:
    _has_scm_libbase = False

try:
    from pandas import DataFrame

    _has_pandas = True
except ImportError:
    _has_pandas = False

try:
    from scm.pisa.block import DriverBlock
    from scm.pisa.input_def import DRIVER_BLOCK_FILES, ENGINE_BLOCK_FILES

    _has_scm_pisa = True
except ImportError:
    _has_scm_pisa = False


__all__ = ["JobAnalysis"]


class JobAnalysis:
    """
    Analysis tool for Jobs, which generates tables of data consisting of fields and their respective value for each job.

    The jobs and fields which are included in the analysis are customizable, to allow for flexible comparison.
    """

    class _Field(NamedTuple):
        name: str
        group: str
        value_extractor: Callable[[Job], Any]

    _job_info_fields = [
        _Field(name="Path", group="job_info", value_extractor=lambda j: j.path),
        _Field(name="Name", group="job_info", value_extractor=lambda j: j.name),
        _Field(name="OK", group="job_info", value_extractor=lambda j: j.ok()),
        _Field(name="Check", group="job_info", value_extractor=lambda j: j.check()),
        _Field(name="ErrorMsg", group="job_info", value_extractor=lambda j: j.get_errormsg()),
    ]

    _job_parent_fields = [
        _Field(name="ParentPath", group="job_parent", value_extractor=lambda j: j.parent.path if j.parent else None),
        _Field(name="ParentName", group="job_parent", value_extractor=lambda j: j.parent.name if j.parent else None),
    ]

    _molecule_fields = [
        _Field(
            name="Formula",
            group="mol",
            value_extractor=lambda j: (
                JobAnalysis._mol_formula_extractor(j.molecule) if isinstance(j, SingleJob) else None
            ),
        ),
        _Field(
            name="Smiles",
            group="mol",
            value_extractor=lambda j: (
                JobAnalysis._mol_smiles_extractor(j.molecule) if isinstance(j, SingleJob) else None
            ),
        ),
    ]

    _timing_fields = [
        _Field(
            name="CPUTime",
            group="timing",
            value_extractor=lambda j: (
                j.results.readrkf("General", "CPUTime") if isinstance(j, AMSJob) and j.results is not None else None
            ),
        ),
        _Field(
            name="SysTime",
            group="timing",
            value_extractor=lambda j: (
                j.results.readrkf("General", "SysTime") if isinstance(j, AMSJob) and j.results is not None else None
            ),
        ),
        _Field(
            name="ElapsedTime",
            group="timing",
            value_extractor=lambda j: (
                j.results.readrkf("General", "ElapsedTime") if isinstance(j, AMSJob) and j.results is not None else None
            ),
        ),
    ]

    @staticmethod
    def _mol_formula_extractor(
        mol: Optional[Union[Molecule, Dict[str, Molecule], "ChemicalSystem", Dict[str, "ChemicalSystem"]]]
    ) -> Optional[str]:
        if isinstance(mol, dict):
            return ", ".join([f"{n}: {JobAnalysis._mol_formula_extractor(m)}" for n, m in mol.items()])
        elif isinstance(mol, Molecule):
            return mol.get_formula()
        elif _has_scm_libbase and isinstance(mol, ChemicalSystem):
            return mol.formula()
        return None

    @staticmethod
    def _mol_smiles_extractor(
        mol: Optional[Union[Molecule, Dict[str, Molecule], "ChemicalSystem", Dict[str, "ChemicalSystem"]]]
    ):
        if isinstance(mol, dict):
            return ", ".join([f"{n}: {JobAnalysis._mol_smiles_extractor(m)}" for n, m in mol.items()])
        elif isinstance(mol, Molecule):
            return to_smiles(mol)
        elif _has_scm_libbase and isinstance(mol, ChemicalSystem):
            return JobAnalysis._mol_smiles_extractor(chemsys_to_plams_molecule(mol))
        else:
            print(f"**** {mol}")
        return None

    _reserved_names = ["_jobs", "_fields", "_pisa_programs"]

    def __init__(self, paths: Optional[Sequence[Union[str, os.PathLike]]] = None, jobs: Optional[Sequence[Job]] = None):
        self._jobs: Dict[str, Job] = {}
        self._fields: Dict[str, JobAnalysis._Field] = {}
        self.add_job_info_fields()

        if _has_scm_pisa:
            self._pisa_programs = {value: key for key, value in ENGINE_BLOCK_FILES.items()}
            self._pisa_programs.update({value: key for key, value in DRIVER_BLOCK_FILES.items()})

        if jobs:
            for j in jobs:
                self.add_job(j)

        if paths:
            for p in paths:
                self.load_job(p)

    @property
    def jobs(self) -> Dict[str, Job]:
        """
        Jobs currently included in analysis.

        :return: Dictionary of the job path and the |Job|
        """
        return {k: v for k, v in self._jobs.items()}

    @property
    def field_names(self) -> List[str]:
        """
        Names of current fields, as they appear in the analysis.

        :return: list of field names
        """
        return [k for k in self._fields]

    @property
    def field_groups(self) -> Dict[str, List[str]]:
        """
        Groups of current fields.

        :return: dictionary with group names as keys and field names as values
        """
        groups = {}
        for field in self._fields.values():
            if field.group in groups:
                groups[field.group].append(field.name)
            else:
                groups[field.group] = [field.name]
        return groups

    def get_analysis(self) -> Dict[str, List]:
        """
        Gets analysis data. This is effectively a table in the form of a dictionary,
        where the keys are the field names and the values are a list of data for each job.

        :return: analysis data as a dictionary of field names/lists of job values
        """

        def safe_value(job: Job, value_extractor: Callable[[Job], Any]):
            try:
                return value_extractor(job)
            except Exception as e:
                return f"ERROR: {str(e)}"

        log_stdout = config.log.stdout
        log_file = config.log.file
        try:
            # Disable logging while fetching results
            config.log.stdout = 0
            config.log.file = 0
            return {col_name: self._get_field_analysis(col_name) for col_name in self._fields}
        finally:
            config.log.stdout = log_stdout
            config.log.file = log_file

    def _get_field_analysis(self, name) -> List:
        """
        Gets analysis data for field with a given name. This gives a list of data for the given field, with a  value for each job.

        :param: name of the field
        :return: analysis data as list of job values
        """
        if name not in self._fields:
            raise KeyError(f"Field with name '{name}' is not part of the analysis.")

        value_extractor = self._fields[name].value_extractor

        def safe_value(job: Job):
            try:
                return value_extractor(job)
            except Exception as e:
                return f"ERROR: {str(e)}"

        log_stdout = config.log.stdout
        log_file = config.log.file
        try:
            # Disable logging while fetching results
            config.log.stdout = 0
            config.log.file = 0
            return [safe_value(j) for j in self._jobs.values()]
        finally:
            config.log.stdout = log_stdout
            config.log.file = log_file

    @requires_optional_package("pandas")
    def to_dataframe(self) -> "DataFrame":
        """
        Converts analysis data to a dataframe. The column names are the field names and the column values are the values for each job.

        :return: analysis data as a dataframe
        """
        return DataFrame(self.get_analysis())

    def to_table(self, max_col_width: int = -1, max_rows: int = 30) -> str:
        """
        Converts analysis data to a pretty-printed table.

        :param max_col_width: can be integer positive value or -1, defaults to -1 (no maximum width)
        :param max_rows: can be integer positive value or -1, defaults to 30
        :return: string representation of the table
        """
        return format_in_table(self.get_analysis(), max_col_width=max_col_width, max_rows=max_rows)

    def to_csv_file(self, path: Union[str, os.PathLike]) -> None:
        """
        Write the analysis to a csv file with the specified path.

        :param path: path to save the csv file
        """
        data = self.get_analysis()
        keys = list(data.keys())
        num_rows = len(data[keys[0]]) if len(keys) > 0 else 0

        with open(path, mode="w", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(keys)
            for i in range(num_rows):
                row = [data[k][i] for k in keys]
                writer.writerow(row)

    def add_job(self, job: Job) -> None:
        """
        Add a job to the analysis. This adds a row to the analysis data.

        :param job: |Job| to add to the analysis
        """
        if job.path in self._jobs:
            raise KeyError(f"Job with path '{job.path}' has already been added to the analysis.")

        self._jobs[job.path] = job

    def remove_job(self, job: Union[str, os.PathLike, Job]) -> None:
        """
        Remove a job from the analysis. This removes a row from the analysis data.

        :param job: |Job| or path to a job to remove from the analysis
        """
        path = job.path if isinstance(job, Job) else str(os.path.abspath(job))
        if path not in self._jobs:
            raise KeyError(f"Job with path '{path}' is not part of the analysis.")

        self._jobs.pop(path)

    def load_job(self, path: Union[str, os.PathLike], loaders: Sequence[Callable[[str], Job]] = None) -> None:
        """
        Add job to the analysis by loading from a given path to the job folder.
        If no dill file is present in that location, the loaders will be used to load the given job from the folder.

        :param path: path to folder from which to load the job
        :param loaders: functions to try and load jobs, defaults to :meth:`~scm.plams.interfaces.adfsuite.ams.AMSJob.load_external` followed by |load_external|
        """
        path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"Cannot find job file in location '{path}'")

        dill_file = (path / path.name).with_suffix(".dill")

        job = None
        loaders = (
            loaders
            if loaders
            else [
                lambda p: AMSJob.load_external(p),
                lambda p: SingleJob.load_external(p),
            ]
        )
        if dill_file.exists():
            job = load(dill_file)
        else:
            for loader in loaders:
                try:
                    job = loader(str(path))
                    break
                except Exception:
                    pass
        if not job:
            raise PlamsError(f"Could not load job from path '{path}'")

        self.add_job(job)

    def filter_jobs(self, predicate: Callable[[Dict[str, Any]], bool]) -> None:
        """
        Remove any jobs from the analysis where the given predicate for field values evaluates to ``True``.
        In other words, this removes rows(s) from the analysis data where the filter function evaluates to ``True`` given a dictionary of the row data.

        :param predicate: filter function which takes a dictionary of field names and their values and evaluates to ``True``/``False``
        """
        analysis = self.get_analysis()
        for i, j in enumerate(self.jobs):
            data = {k: v[i] for k, v in analysis.items()}
            if predicate(data):
                self.remove_job(j)

    def add_field(self, name: str, value_extractor: Callable[[Job], Any], group: Optional[str] = None) -> None:
        """
        Add a field to the analysis. This adds a column to the analysis data.

        :param name: name of the field
        :param value_extractor: callable to extract the value for the field from a job
        :param group: an optional group that this field belongs to
        """
        if name in self._fields:
            raise KeyError(f"Field with name '{name}' has already been added to the analysis.")

        self._fields[name] = self._Field(name=name, group=group, value_extractor=value_extractor)

    def rename_field(self, name: str, new_name: str) -> None:
        """
        Rename a field in the analysis. This is the header of the column in the analysis data.

        :param name: current name of the field
        :param new_name: new name of the field
        """

        if name not in self._fields:
            raise KeyError(f"Field with name '{name}' is not part of the analysis.")

        field = self._fields.pop(name)
        self.add_field(new_name, field.value_extractor, field.group)

    def reorder_fields(self, order: Sequence[str]) -> None:
        """
        Reorder fields based upon the given sequence of field names. This is the order the columns will appear in the analysis data.

        Any specified fields will be placed first, with remaining fields placed after with their order unchanged.

        :param order: sequence of fields to be placed at the start of the field ordering
        """

        def key(field_name):
            try:
                return order.index(field_name)
            except ValueError:
                return len(order)

        self.sort_fields(key=key)

    def sort_fields(self, key: Callable[[str], Any]) -> None:
        """
        Sort fields according to a sort key. This is the order the columns will appear in the analysis data.

        :param key: sorting function which accepts the field name
        """
        sorted_keys = sorted(self._fields.keys(), key=key)
        self._fields = {k: self._fields[k] for k in sorted_keys}

    def remove_field(self, name: str) -> None:
        """
        Remove a field from the analysis. This removes a column from the analysis data.

        :param name: name of the field
        """
        if name not in self._fields:
            raise KeyError(f"Field with name '{name}' is not part of the analysis.")

        self._fields.pop(name)

    def remove_field_group(self, group: str) -> None:
        """
        Remove a field group from the analysis. This removes column(s) from the analysis data.

        :param group: name of the group to remove
        """
        names = []
        for field in self._fields.values():
            if field.group == group:
                names.append(field.name)

        for name in names:
            self.remove_field(name)

    def filter_fields(self, predicate: Callable[[List[Any]], bool]) -> None:
        """
        Remove any fields from the analysis where the given predicate evaluates to ``True`` for all values.
        In other words, this removes column(s) from the analysis data where the filter function evaluates to ``True``
        given all the row values.

        :param predicate: filter function which takes values and evaluates to ``True``/``False``
        """
        for n, vals in self.get_analysis().items():
            if predicate(vals):
                self.remove_field(n)

    def remove_empty_fields(self) -> None:
        """
        Remove field(s) from the analysis which have ``None`` for all values. This removes column(s) from the analysis data,
        where all rows have empty values.
        """
        self.filter_fields(lambda vals: all([v is None for v in vals]))

    def remove_uniform_fields(self, tol: float = 1e-08, ignore_empty: bool = False) -> None:
        """
        Remove field(s) from the analysis which evaluate the same for all values. This removes column(s) from the analysis data,
        where all rows have the same value.

        :param tol: absolute tolerance for numeric value comparison, all values must fall within this range
        :param ignore_empty: when ``True`` ignore ``None`` values in comparison, defaults to ``False``.
        """

        def is_uniform(vals: List[Any]):
            # Skip over None values if set to be ignored
            vals = [v for v in vals if v is not None or not ignore_empty]
            # Empty list is uniform
            if not vals:
                return True
            # Check if all numeric values, and if so evaluate range within tolerance
            if all([isinstance(v, Number) and not isinstance(v, bool) for v in vals]):
                return np.ptp(vals) <= tol

            try:
                return all([v == vals[0] for v in vals])
            except ValueError:
                return False

        self.filter_fields(lambda vals: is_uniform(vals))

    def _add_preconfigured_fields(self, fields):
        for field in fields:
            if field.name not in self._fields:
                self._fields[field.name] = field

    def _remove_preconfigured_fields(self, fields):
        for field in fields:
            if field.name in self._fields:
                self._fields.pop(field.name)

    def add_job_info_fields(self) -> None:
        """
        Adds the set of default job information fields to the analysis.

        Namely:

        * Path: :attr:`~scm.plams.core.basejob.Job.path`
        * Name: :attr:`~scm.plams.core.basejob.Job.name`
        * OK: :meth:`~scm.plams.core.basejob.Job.ok`
        * Check: :meth:`~scm.plams.core.basejob.Job.check`
        * ErrorMsg: :meth:`~scm.plams.core.basejob.Job.get_erromsg`
        """
        self._add_preconfigured_fields(self._job_info_fields)

    def remove_job_info_fields(self) -> None:
        """
        Removes the set of default job information fields from the analysis.

        Namely:

        * Path: :attr:`~scm.plams.core.basejob.Job.path`
        * Name: :attr:`~scm.plams.core.basejob.Job.name`
        * OK: :meth:`~scm.plams.core.basejob.Job.ok`
        * Check: :meth:`~scm.plams.core.basejob.Job.check`
        * ErrorMsg: :meth:`~scm.plams.core.basejob.Job.get_erromsg`
        """
        self._remove_preconfigured_fields(self._job_info_fields)

    def add_job_parent_fields(self) -> None:
        """
        Adds the set of default job parent fields to the analysis.

        Namely:

        * ParentPath: :attr:`~scm.plams.core.basejob.Job.path` of :attr:`~scm.plams.core.basejob.Job.parent`
        * ParentName: :attr:`~scm.plams.core.basejob.Job.name` of :attr:`~scm.plams.core.basejob.Job.parent`
        """
        self._add_preconfigured_fields(self._job_parent_fields)

    def remove_job_parent_fields(self) -> None:
        """
        Removes the set of default job information fields from the analysis.

        Namely:

        * ParentPath: :attr:`~scm.plams.core.basejob.Job.path` of :attr:`~scm.plams.core.basejob.Job.parent`
        * ParentName: :attr:`~scm.plams.core.basejob.Job.name` of :attr:`~scm.plams.core.basejob.Job.parent`
        """
        self._add_preconfigured_fields(self._job_parent_fields)

    def add_molecule_fields(self) -> None:
        """
        Adds the set of default molecule fields to the analysis.

        Namely:

        * Formula: :method:`~scm.plams.mol.molecule.Molecule.get_formula` of :attr:`~scm.plams.core.basejob.SingleJob.molecule` and derived classes
        * Smiles: :func:`~scm.plams.interfaces.molecule.rdkit.to_smiles` with :attr:`~scm.plams.core.basejob.SingleJob.molecule` and derived classes

        Jobs without a molecule will not have these fields populated.
        Where jobs have multiple molecules, these will be returned as comma-separated values.
        """
        self._add_preconfigured_fields(self._molecule_fields)

    def remove_molecule_fields(self) -> None:
        """
        Removes the set of default molecule fields from the analysis.

        Namely:

        * Formula: :method:`~scm.plams.mol.molecule.Molecule.get_formula` of :attr:`~scm.plams.core.basejob.SingleJob.molecule` and derived classes
        * Smiles: :func:`~scm.plams.interfaces.molecule.rdkit.to_smiles` with :attr:`~scm.plams.core.basejob.SingleJob.molecule` and derived classes

        Jobs without a molecule will not have these fields populated.
        Where jobs have multiple molecules, these will be returned as comma-separated values.
        """
        self._remove_preconfigured_fields(self._molecule_fields)

    def add_timing_fields(self) -> None:
        """
        Adds the set of default timing fields to the analysis.

        Namely:

        * CPUTime: :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_readrkf` with General/CPUTime for :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results` and derived classes
        * SysTime: :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_readrkf` with General/SysTime for :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results` and derived classes
        * ElapsedTime: :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_readrkf` with General/ElapsedTime for :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results` and derived classes

        Jobs which do not derive from |AMSJob| or have no results will not have these fields populated.
        """
        self._add_preconfigured_fields(self._timing_fields)

    def remove_timing_fields(self) -> None:
        """
        Removes the set of default timing fields from the analysis.

        Namely:

        * CPUTime: :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_readrkf` with General/CPUTime for :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results` and derived classes
        * SysTime: :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_readrkf` with General/SysTime for :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results` and derived classes
        * ElapsedTime: :meth:`~scm.plams.interfaces.adfsuite.ams.AMSResults.get_readrkf` with General/ElapsedTime for :attr:`~scm.plams.interfaces.adfsuite.ams.AMSJob.results` and derived classes

        Jobs which do not derive from |AMSJob| or have no results will not have these fields populated.
        """
        self._remove_preconfigured_fields(self._timing_fields)

    def add_settings_field(self, key_tuple: Tuple[Hashable, ...]) -> None:
        """
        Add a field for a nested key from the job settings to the analysis.
        The name of the field will be a period-delimited string of the key path e.g. ("input", "ams", "task") will appear as field ``Input.Ams.Task``.

        :param key_tuple: nested tuple of keys in the settings object
        """
        self.add_field(
            ".".join([str(k).title() for k in key_tuple]),
            lambda j, k=key_tuple: self._get_job_settings(j).get_nested(k),
            group="settings",
        )

    def add_settings_fields(
        self,
        predicate: Optional[Callable[[Tuple[Hashable, ...]], bool]] = None,
        flatten_list: bool = True,
        group: Optional[str] = "settings",
    ) -> None:
        """
        Add a field for all nested keys which satisfy the predicate from the job settings to the analysis.
        The name of the fields will be a period-delimited string of the key path e.g. ("input", "ams", "task") will appear as field ``Input.Ams.Task``.

        :param predicate: optional predicate which evaluates to ``True`` or ``False`` given a nested key, by default will be ``True`` for every key
        :param flatten_list: whether to flatten lists in settings objects
        :param group: an optional group that this set of settings fields belongs to
        """

        all_blocks = set()
        all_keys = {}  # Use dict as a sorted set for keys
        for job in self._jobs.values():
            settings = self._get_job_settings(job)
            blocks = set(settings.block_keys(flatten_list=flatten_list))
            keys = {k: None for k in settings.nested_keys(flatten_list=flatten_list)}
            all_blocks = all_blocks.union(blocks)
            all_keys.update(keys)

        fields = []
        predicate = predicate if predicate else lambda _: True
        for key in all_keys:
            # Take only final nested keys i.e. those which are not block keys and satisfy the predicate
            if key not in all_blocks and predicate(key):
                name = ".".join([str(k).title() for k in key])
                field = self._Field(
                    name=name, value_extractor=lambda j, k=key: self._get_job_settings(j).get_nested(k), group=group
                )
                fields.append(field)

        self._add_preconfigured_fields(fields)

    def _get_job_settings(self, job: Job) -> Settings:
        """
        Get job settings converting any PISA input block to a standard settings object.
        """
        # Convert any PISA settings blocks to standard
        settings = Settings()
        if job.settings is not None:
            if isinstance(job.settings, Settings):
                settings = job.settings.copy()
            if _has_scm_pisa:
                if hasattr(job.settings, "input") and isinstance(job.settings.input, DriverBlock):
                    # Note use own input parser facade here to use caching
                    program = self._pisa_programs[job.settings.input.name].name.split(".")[0]
                    settings.input = InputParserFacade().to_settings(
                        program=program, text_input=job.settings.input.get_input_string()
                    )
        return settings

    def add_settings_input_fields(self, include_system_block: bool = False, flatten_list: bool = True) -> None:
        """
        Add a field for each input key in the :attr:`~scm.plams.core.basejob.Job.settings` object across all currently added jobs.

        :param include_system_block: whether to include keys for the system block, defaults to ``False``
        :param flatten_list: whether to flatten lists in settings objects
        """

        def predicate(key_tuple: Tuple[Hashable, ...]):
            if len(key_tuple) == 0 or str(key_tuple[0]).lower() != "input":
                return False

            return (
                len(key_tuple) < 3
                or str(key_tuple[1]).lower() != "ams"
                or str(key_tuple[2]).lower() != "system"
                or include_system_block
            )

        self.add_settings_fields(predicate, flatten_list, group="settings_input")

    def remove_settings_fields(self) -> None:
        """
        Remove all fields which contain ``settings`` in the field group name from the analysis.
        """
        for group in self.field_groups.keys():
            if "settings" in group.lower():
                self.remove_field_group(group)

    def __str__(self) -> str:
        return format_in_table(self.get_analysis(), max_col_width=12, max_rows=5)

    def __getitem__(self, name: str) -> List[Any]:
        """
        Get analysis data for a given field.

        :param name: name of the field
        :return: list of values for each job
        """
        return self._get_field_analysis(name)

    def __setitem__(self, name: str, value: Callable[[Job], Any]) -> None:
        """
        Set analysis for given field.

        :param name: name of the field
        :param value: callable to extract the value for the field from a job
        """
        if not callable(value):
            raise TypeError(f"To set a field, the value must be a callable which accepts a Job.")

        if name in self._fields:
            self.remove_field(name)
        self.add_field(name, value_extractor=value)

    def __delitem__(self, name: str) -> None:
        """
        Delete analysis for given field.

        :param name: name of the field
        """
        self.remove_field(name)

    def __getattr__(self, name: str) -> List[Any]:
        """
        Fallback to get analysis for given field when an attribute is not present.

        :param name: name of the field
        :return: list of values for each job
        """
        try:
            return self[name]
        except KeyError:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute or analysis field with name '{name}'"
            )

    def __setattr__(self, name, value) -> None:
        """
        Fallback to set analysis for given field.

        :param name: name of the field
        :param value: callable to extract the value for the field from a job
        """
        if name in self._reserved_names or hasattr(self.__class__, name):
            super().__setattr__(name, value)
        else:
            self[name] = value

    def __delattr__(self, name) -> None:
        """
        Fallback to set analysis for given field.

        :param name: name of the field
        """
        if name in self._reserved_names or hasattr(self.__class__, name):
            super().__delattr__(name)
        else:
            try:
                del self[name]
            except KeyError:
                raise AttributeError(
                    f"'{self.__class__.__name__}' object has no attribute or analysis field with name '{name}'"
                )

    def __dir__(self):
        """
        Return standard attributes, plus dynamically added field names which can be accessed via dot notation.
        """
        return [x for x in super().__dir__()] + [
            k for k in self._fields.keys() if isinstance(k, str) and k.isidentifier()
        ]
