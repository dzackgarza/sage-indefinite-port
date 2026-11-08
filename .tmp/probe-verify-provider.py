"""Confirm inside the container which sage-extensions tree is executing, and that it carries the fix."""
import hashlib
import inspect

import sage_extensions

base = inspect.getfile(sage_extensions).rsplit("/", 1)[0]
print("module:", inspect.getfile(sage_extensions))
for rel in ("categories/rational_integral_stabilizers.py",
            "categories/modules/orbit_spans.py",
            "categories/group/groups.py"):
    path = f"{base}/{rel}"
    print(f"  {rel:48s} sha256 {hashlib.sha256(open(path, 'rb').read()).hexdigest()[:16]}")

stabilizers = open(f"{base}/categories/rational_integral_stabilizers.py").read()
spans = open(f"{base}/categories/modules/orbit_spans.py").read()
print("  stabilizer(): composite disjunct first :",
      "and self.stable_containing_bound() is not None\n            )\n            or rational_group.is_finite() is True" in stabilizers)
print("  :378 containing is None first          :", "if containing is None and group.is_finite() is not True:" in stabilizers)
print("  orbit_spans containing is None first   :", "if containing is None and group.is_finite() is True:" in spans)
