# Changelog
Notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

The versioning for this project aligns with that of AMS. 
The format for this is the original year of the main release as a prefix, and an incremental postfix for each sub-release, starting at 101.
In general, subsequent sub-releases after the main release contain only bug fixes.
For example, 2024.101 is the major release of 2024, and 2024.102 is the first bugfix release.

This changelog is effective from the 2025 releases.

## [Unreleased]

### Added
* Methods `get_system`, `get_input_system` and `get_main_system` to `AMSResults`, which return an AMS `ChemicalSystem` instead of a PLAMS `Molecule` 
* `AMSJob` can accept an AMS `ChemicalSystem` instead of a PLAMS `Molecule` as an input system
* Specific `ConfigSettings` and related settings classes with explicitly defined fields
* Support for work functions: `AMSResults.get_work_function_results` and `plot_work_function`
* Example on `MoleculeFormats`
* Script `generate_example.sh` to generate documentation pages from notebook examples
* GitHub workflows for CI and publishing to PyPI

### Changed
* Functions for optional packages (e.g. RDKit, ASE) are available even when these packages are not installed, but will raise an `MissingOptionalPackageError` when called
* `AMSResults.get_main_ase_atoms` also includes atomic charges
* Global `config` is initialized with a `ConfigSettings` instead of loading from the standard `plams_defaults` file
* `init` and `finish` functions are now optional
* `Job.status` is a `JobStatus` string enum
* Supercell and RDKit properties are no longer serialized to AMS input
* Restructuring of examples and conversion of various examples to notebooks

### Fixed
* `Molecule.properties.charge` is a numeric instead of string type when loading molecule from a file
* `Molecule.delete_all_bonds` removes the reference molecule from the removed bond instances
* `SingleJob.load` returns the correctly loaded job
* `AMSJob.check` handles a `NoneType` status, returning `False`


### Deprecated
* `plams` launch script is deprecated in favour of simply running with `amspython`

### Removed
* Legacy `BANDJob`, `DFTBJob`, `UFFJob`, `MOPACJob`, `ReaxFFJob`, `CSHessianADFJob` and `ADFJob` have been removed
* Exception classes `AMSPipeDecodeError`, `AMSPipeError`, `AMSPipeInvalidArgumentError`, `AMSPipeLogicError`, `AMSPipeRuntimeError`, `AMSPipeUnknownArgumentError`, `AMSPipeUnknownMethodError`, `AMSPipeUnknownVersionError`, were moved from scm.plams to scm.amspipe.




