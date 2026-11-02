"""Import shim: the port needs ``sage_extensions.all``; this image has no ``sageparse``.

``sage_extensions.all`` imports ``sageparse.preparser.research``, which is absent from this VM
(not on PyPI, not in any site-packages here). Everything the port actually uses from that
namespace is reachable through the ordinary category modules, so install a stand-in module
``sage_extensions.all`` exposing exactly those names. Nothing in the repository is modified.
"""
import sys
import types

import sage.all  # noqa: F401  (initializes Sage's cyclic imports first)

from sage.rings.integer_ring import ZZ as ZZ  # noqa: F401
from sage_extensions.categories.abstract_categories.cat import Cat  # noqa: F401
from sage_extensions.categories.lattices import Lattices
from sage_extensions.categories.sets import Sets

stub = types.ModuleType("sage_extensions.all")
stub.Lattices = Lattices
stub.Sets = Sets
stub.Cat = Cat
stub.ZZ = ZZ
sys.modules.setdefault("sage_extensions.all", stub)
