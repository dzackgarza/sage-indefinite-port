"""Which recognition routes exist in the installed GAP, and is one merely not loaded?"""
from sage.all import libgap

def show(label, expr):
    try:
        print(f"{label}: {libgap.eval(expr)}", flush=True)
    except Exception as error:
        print(f"{label}: raised {type(error).__name__}: {str(error)[:120]}", flush=True)

print("GAP version:", libgap.eval("GAPInfo.Version"), flush=True)
show("IsBound(IsomorphismTypeOfMatrixGroup)", "IsBound(IsomorphismTypeOfMatrixGroup)")
show("IsBound(IsomorphismTypeOfMatrixGroupGeneric)", "IsBound(IsomorphismTypeOfMatrixGroupGeneric)")
show("IsBound(NiceMonomorphism)", "IsBound(NiceMonomorphism)")
show("IsBound(SizeOfMatrixGroup)", "IsBound(SizeOfMatrixGroup)")
show("IsBound(IsomorphismPermGroup)", "IsBound(IsomorphismPermGroup)")

print("--- installed packages (name only) ---", flush=True)
show("PackagesInfo names", "RecNames(GAPInfo.PackagesInfo)")
for package in ("recog", "atlasrep", "ctbllib", "forms", "genss", "io"):
    show(f'LoadPackage("{package}")', f'LoadPackage("{package}")')

print("--- after loading ---", flush=True)
show("IsBound(IsomorphismTypeOfMatrixGroup)", "IsBound(IsomorphismTypeOfMatrixGroup)")
show("IsBound(IsomorphismTypeOfMatrixGroupGeneric)", "IsBound(IsomorphismTypeOfMatrixGroupGeneric)")
show("IsBound(RecogniseMatrixGroup)", "IsBound(RecogniseMatrixGroup)")
