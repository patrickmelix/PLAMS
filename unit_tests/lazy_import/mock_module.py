import sys

sys.modules["test_lazy_import"] = 5
raise ImportError
