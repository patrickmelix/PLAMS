from __future__ import annotations

import sys
from importlib import abc, util
from types import ModuleType
from typing import TYPE_CHECKING
import threading
import warnings

_old_threading_start = threading.Thread.start


def monkey_patch_thread_start(self: threading.Thread) -> None:
    lazy_modules = LazyImporter.get_modules()
    if len(lazy_modules) > 0:
        msg = f"""\


Thread started with lazily imported modules: {lazy_modules} present in sys.modules.
Likely caused by importing (indirectly) from scm.plams, which lazily imports some slow external libraries.
If you are planning to use the lazy modules inside of the thread, you should call:

from scm.plams.lazy_import import LazyImporter
LazyImporter.resolve_imports()

Which is normally done in the plams 'init' method. Not resolving the lazily imported modules before using them
in threads can result in an AttributeError.

If the modules will not be used by the thread, you can ingore this warning.
        """
        old_warning_filters = warnings.filters
        try:
            warnings.filters = []
            warnings.warn(msg, ImportWarning)
        finally:
            warnings.filters = old_warning_filters
    _old_threading_start(self)


threading.Thread.start = monkey_patch_thread_start  # type: ignore


class LazyImporter:
    _lazy_modules: list[str] = []

    @classmethod
    def resolve_imports(cls) -> None:
        for mod_name in cls.get_modules():
            sys.modules[mod_name].__name__  # getattr call triggers module loading
            cls._lazy_modules.remove(mod_name)

    @classmethod
    def get_modules(cls) -> list[str]:
        return [
            mod_name
            for mod_name in cls._lazy_modules
            if mod_name in sys.modules and isinstance(sys.modules[mod_name], util._LazyModule)  # type: ignore
        ]

    @classmethod
    def import_module(cls, name: str) -> ModuleType:
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
        cls._lazy_modules.append(name)
        return module


if TYPE_CHECKING:
    import ase
    import dill
    import networkx
    import numpy
else:
    # we have to make sure we don't import things twice
    if "numpy" in sys.modules:
        import numpy
    else:
        numpy = LazyImporter.import_module("numpy")
    if "networkx" in sys.modules:
        import networkx
    else:
        networkx = LazyImporter.import_module("networkx")
    if "dill" in sys.modules:
        import dill
    else:
        dill = LazyImporter.import_module("dill")
    if "ase" in sys.modules:
        import ase
    else:
        ase = LazyImporter.import_module("ase")
