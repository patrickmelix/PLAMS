from __future__ import annotations

from importlib import util, abc
import sys
from types import ModuleType
from typing import TYPE_CHECKING


def lazy_import(name: str) -> ModuleType:
    """Taken from https://docs.python.org/3.11/library/importlib.html#implementing-lazy-imports.
    Returns a ``ModuleType`` object that will only import the actual module when the first
    attribute on it is accessed. So for example:

    ..code:: python3

        np = lazy_import('numpy')

    will not actually import the module, but it will check if this module exists through
    the ``find_spec`` call. Note that this does not guarantee that the module can be successfully
    imported, but only that a module with that name exists.

    ..code:: python3

        np = lazy_import('numpy')

        def dot_product(vec1, vec2):
            return np.dot(vec1, vec2)

    In the example above the numpy module will only be imported when ``dot_product`` is first called.
    This may cause an import error at an unexpected location, so the unit tests should always cover
    the imports of all dependencies.

    Large advantage is that actual importing of slow external libraries can be postponed, so we can still
    expose all modules by importing them in the main __init__.py, without worrying about always importing
    slow external libraries that are only needed in specific parts of the code.

    :param name: name of the module to import
    :type name: str
    :return: Lazily loaded module
    :rtype: ModuleType
    """
    spec = util.find_spec(name)
    if spec is None:
        raise ImportError(f"Can not find spec for module {name}")
    assert isinstance(spec.loader, abc.Loader)
    loader = util.LazyLoader(spec.loader)
    spec.loader = loader
    module = util.module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


if TYPE_CHECKING:
    import numpy
    import networkx
    import dill
    import ase
else:
    numpy = lazy_import("numpy")
    networkx = lazy_import("networkx")
    dill = lazy_import("dill")
    ase = lazy_import("ase")
