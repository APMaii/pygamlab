from typing import Any, Iterable
# Import safely (so missing libraries don’t break the module)
from typing import Any, Iterable

from collections.abc import Iterable

from collections.abc import Iterable
from typing import Any
from ..Primatom import GAM_Atom

try:
    from ase import Atoms as ASE_Atoms
except ImportError:
    ASE_Atoms = None

try:
    from pymatgen.core import Structure as PymatgenStructure,Molecule as PymatgenMolecule

except ImportError:
    PymatgenStructure = None

try:
    import MDAnalysis as mda
except ImportError:
    mda = None

try:
    #import openbabel
    from openbabel import pybel
except ImportError:
    #openbabel = None
    pybel = None

# Your custom class
try:
    from ..Primatom import GAM_Atom 

except ImportError:
    GAM_Atom = None

try:
    from rdkit import Chem
except ImportError:
    Chem= None


# ------------------------------
# Individual check functions
# ------------------------------
def is_ase(obj: Any) -> bool:
    """
    Check whether the object is an ASE Atoms instance.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if `obj` is an instance of ASE Atoms and ASE is available; False otherwise.
    """
    
    return ASE_Atoms is not None and isinstance(obj, ASE_Atoms)



def is_gam(obj: Any) -> bool:
    """
    Check if the object is a `GAM_Atom` or an iterable of `GAM_Atom` objects.

    This function supports single `GAM_Atom` instances, as well as lists, tuples,
    numpy arrays, or any iterable container of `GAM_Atom` objects.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if `obj` is a `GAM_Atom` or an iterable of `GAM_Atom`; False otherwise.
    """
    if GAM_Atom is None:
        return False

    # Single atom
    if isinstance(obj, GAM_Atom):
        return True

    # Check for iterable of GAM_Atom
    if isinstance(obj, Iterable) and not isinstance(obj, (str, bytes)):
        try:
            return all(isinstance(x, GAM_Atom) for x in obj)
        except Exception:
            return False

    return False




def is_pymatgen(obj: Any) -> bool:
    """
    Check if the object is a Pymatgen `Structure` or `Molecule` instance.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if `obj` is a `pymatgen.core.structure.Structure` or
        `pymatgen.core.structure.Molecule`; False otherwise.
    """
    return (
        (PymatgenStructure is not None and isinstance(obj, PymatgenStructure)) or
        (PymatgenMolecule is not None and isinstance(obj, PymatgenMolecule))
    )

def is_mdanalysis(obj: Any) -> bool:
    """
    Check if the object is an MDAnalysis `Universe` or `AtomGroup`.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if `obj` is an instance of `mda.Universe` or `mda.core.groups.AtomGroup`; False otherwise.
    """
    if mda is None:
        return False
    return isinstance(obj, (mda.Universe, mda.core.groups.AtomGroup))

def is_openbabel(obj: Any) -> bool:
    """
    Check if the object is an OpenBabel/Pybel molecule.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if `obj` is a `pybel.Molecule` (or OpenBabel OBMol if supported); False otherwise.
    """
    if pybel is None: #and openbabel is None:
        return False
    return (
        (pybel is not None and isinstance(obj, pybel.Molecule))
        #or (openbabel is not None and isinstance(obj, openbabel.OBMol))
    )

def is_rdkit(obj: Any) -> bool:
    """
    Check if the object is an RDKit molecule.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if `obj` is an instance of `rdkit.Chem.Mol`; False otherwise.
    """
    if Chem is None:
        return False
    return isinstance(obj, Chem.Mol)





# ------------------------------
# Global format checker
# ------------------------------
def detect_format(obj: Any) -> str:
    """
    Detect the molecular format/library of the given object.

    This function attempts to identify which molecular representation or
    computational chemistry library the object belongs to.

    Parameters
    ----------
    obj : Any
        Object representing molecular data.

    Returns
    -------
    str
        A string indicating the detected format:
            - 'ASE' : ASE Atoms object
            - 'GAM' : GAM_Atom or iterable of GAM_Atom
            - 'pymatgen' : Pymatgen Structure or Molecule
            - 'MDAnalysis' : MDAnalysis Universe or AtomGroup
            - 'OpenBabel' : OpenBabel/Pybel molecule
            - 'RDKit' : RDKit molecule
            - 'Unknown' : Format could not be determined
    """
    if is_ase(obj):
        return "ASE"
    elif is_gam(obj):
        return "GAM"
    elif is_pymatgen(obj):
        return "pymatgen"
    elif is_mdanalysis(obj):
        return "MDAnalysis"
    elif is_openbabel(obj):
        return "OpenBabel"
    elif is_rdkit(obj):
        return "RDKit"
    else:
        return "Unknown"
