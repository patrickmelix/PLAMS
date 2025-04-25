from setuptools import find_packages, setup
from glob import glob
import os

# This minimal setup.py exists alongside the pyproject.toml for legacy reasons.
# Currently, the "artificial" prefix "scm." is added to the (sub)package names, which is not supported via the pyproject.toml.
# ToDo: the package should be restructured with a directory structure that reflects this, then the setuptools package finding used.
# See: https://setuptools.pypa.io/en/latest/userguide/pyproject_config.html
examples_files = [
    os.path.relpath(path)  # relative to package_dir
    for path in glob("examples/**", recursive=True)
    if os.path.isfile(path)
]
test_files = [
    os.path.relpath(path)  # relative to package_dir
    for path in glob("unit_tests/**", recursive=True)
    if os.path.isfile(path)
]

setup(
    packages=["scm.plams"] + ["scm.plams." + i for i in find_packages(".")],
    package_dir={"scm.plams": "."},
    package_data={"scm.plams": [".flake8"] + examples_files + test_files},
)
