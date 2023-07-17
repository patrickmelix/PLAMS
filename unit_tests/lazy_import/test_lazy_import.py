import re
import sys
import threading
import pytest

from scm.plams.lazy_import.lazy_import import LazyImporter
from scm.plams.core.functions import init
from scm.utils.runsubprocess import RunSubprocess


def test_lazy_import():
    # module defined in this test folder for testing purposes
    assert "mock_module" not in sys.modules
    # module raises an ImportError when executed, so we are sure it is lazily imported
    mock_module = LazyImporter.import_module("mock_module")
    assert "mock_module" in sys.modules
    try:
        # on accessing attribute, module is actually imported and it places a flag in
        # the sys.modules dictionary so we can test it is was actually executed
        mock_module.foo
    # module raises ImportError at the end
    except ImportError:
        assert sys.modules["test_lazy_import"] == 5


def test_star_import_no_slow_libraries():
    # the importtime commandline option prints an overview of all imported libraries and the time it took to import
    # in the following format:
    # import time: self [us] | cumulative | imported package
    # import time:       145 |        145 |   _io
    # import time:        25 |         25 |   marshal
    # import time:       199 |        199 |   posix
    # import time:       346 |        713 | _frozen_importlib_external
    success, _, error = RunSubprocess('sh -c ""$AMSBIN/amspython" -X importtime -c \'from scm.plams import *\'"')
    assert success
    # we have a list of 'slow' external libraries that we don't want to import from the root plams __init__.py
    # these libraries are either lazily imported using scm.plans.lazy_import, by inlining the imports in function
    # definitions instead of the top of the module, or by simply not importing certain modules from the root __init__.py
    slow_libraries = ["numpy", "scipy", "ase", "rdkit", "dill", "networkx", "natsort"]
    for lib in slow_libraries:
        # this regex is not guaranteed to never have false positives, and might have to be adjusted in the future
        assert not re.search(pattern=r"\|\s*" + lib, string=error)
    # if this tests fails, you can easily visualize the import flow with tuna https://github.com/nschloe/tuna
    # amspython -X importtime -c 'from scm.plams import *' 2> test_imports.log && amspython -m tuna test_imports.log


def test_warning():
    t = threading.Thread()
    # should raise the warning, since lazy modules are present in the module cache
    with pytest.warns(ImportWarning):
        t.start()
        t.join()


@pytest.mark.filterwarnings("error::pytest.PytestUnhandledThreadExceptionWarning")
def test_lazy_module_resolve_imports(tmp_path):
    from scm.plams.lazy_import import numpy as np
    from scm.plams.lazy_import import ase

    def access_module_attrs():
        for _ in range(1000):
            np.ndarray
            ase.Atom

    # Make sure to test with more threads than there are cores on the machine
    threads = [threading.Thread(target=access_module_attrs) for _ in range(256)]

    init(folder=tmp_path)  # calls LazyImporter.resolve_imports under the hood.
    for t in threads:
        t.start()
    for t in threads:
        t.join()
