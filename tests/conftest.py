"""The research session for the .sage cases: its preparser prelude binds Integer to the owned integers.

qc_sage_session lowers a .sage file after the conftests have run and binds the installed preparser's
runtime_namespace() in it; the preamble registers its literal constructors there when it is imported.
"""

import importlib

importlib.import_module("dzack_research.preamble.all")
