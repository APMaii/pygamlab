from pymatgen.core import Element, Molecule, Structure, Lattice
from typing import List
from ..Primatom import GAM_Atom 
import numpy as np
from .checker import is_ase, is_gam, is_pymatgen, is_mdanalysis, is_openbabel , is_rdkit
from ase import Atoms

# Import required libraries for each format
try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available. Install with: pip install rdkit")

try:
    import MDAnalysis as mda
    from MDAnalysis import Universe
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False
    print("MDAnalysis not available. Install with: pip install MDAnalysis")

try:
    import mdtraj as md
    MDTRAJ_AVAILABLE = True
except ImportError:
    MDTRAJ_AVAILABLE = False
    print("MDTraj not available. Install with: pip install mdtraj")

try:
    from openbabel import pybel
    import openbabel as ob
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False
    print("OpenBabel not available. Install with: pip install openbabel-wheel")

try:
    import ase
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    print("ASE not available. Install with: pip install ase")





def gam_to_ase(atoms: List[GAM_Atom]) -> 'ase.Atoms':
        """
        Convert a list of GAM_Atom objects to an ASE Atoms object.

        Parameters
        ----------
        atoms : List[GAM_Atom]
            List of GAM_Atom objects to convert.

        Returns
        -------
        ase.Atoms
            ASE Atoms object representing the same atomic symbols and coordinates.

        Raises
        ------
        TypeError
            If the input is not a list of GAM_Atom objects.
        """
        if is_gam(atoms) is False:
            raise TypeError("Input is not a list of GAM_Atom objects.")
        symbols = [atom.element for atom in atoms]
        positions = [[atom.x, atom.y, atom.z] for atom in atoms]
        
        gam_atoms = Atoms(symbols=symbols, positions=positions)
        
        return gam_atoms

def ase_to_gam(atoms: 'ase.Atoms') -> List[GAM_Atom]:
    """
    Convert an ASE Atoms object to a list of GAM_Atom objects.

    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object to convert.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with the same atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not an ASE Atoms object.
    """
    if is_ase(atoms) is False:
        raise TypeError("Input is not an ASE Atoms object.")
    atoms = []
    for i, atom in enumerate(atoms, start=1):
        atoms.append(GAM_Atom(
            id=i,
            element=atom.symbol,
            x=atom.position[0],
            y=atom.position[1],
            z=atom.position[2]
        ))
    return atoms        



def gam_to_pymatgen(atoms: List[GAM_Atom]) -> Molecule:
    """
    Convert a list of GAM_Atom objects to a pymatgen Molecule object.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of GAM_Atom objects to convert.

    Returns
    -------
    pymatgen.core.structure.Molecule
        Molecule object with the same atomic symbols and coordinates.
    """
    #if is_gam(atoms) is False:
        #raise TypeError("Input is not a list of GAM_Atom objects.")
    symbols = [atom.element for atom in atoms]
    coords = [[atom.x, atom.y, atom.z] for atom in atoms]
    return Molecule(symbols, coords)



def pymatgen_to_gam(sites) -> List[GAM_Atom]:
    """
    Convert a pymatgen Structure or Molecule to a list of GAM_Atom objects.

    Parameters
    ----------
    sites : pymatgen Structure or Molecule
        Pymatgen object to convert.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with the same atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not a pymatgen Structure or Molecule object.
    """
    if is_pymatgen(sites) is False:
        raise TypeError("Input is not a pymatgen Structure or Molecule object.")
    atoms = []
    for i, site in enumerate(sites, start=1):
        atom = GAM_Atom(
            id=i,
            element=getattr(site, "species_string", site.species_string if hasattr(site, "species_string") else site.element),
            x=site.x,
            y=site.y,
            z=site.z
        )
        atoms.append(atom)
    return atom



def gam_to_rdkit(atoms: List[GAM_Atom]) -> 'Chem.Mol':
    """
    Convert a list of GAM_Atom objects to an RDKit Mol object.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of GAM_Atom objects to convert.

    Returns
    -------
    rdkit.Chem.Mol
        RDKit Mol object with atom types and 3D coordinates.

    Raises
    ------
    TypeError
        If the input is not a list of GAM_Atom objects.
    ImportError
        If RDKit is not installed.
    """
    if is_gam(atoms) is False:
        raise TypeError("Input is not a list of GAM_Atom objects.")
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is not available. Please install it first.")
    mol = Chem.RWMol()
    atom_indices = []

    # Add atoms
    for atom in atoms:
        a = Chem.Atom(atom.element)
        idx = mol.AddAtom(a)
        atom_indices.append(idx)

    # Create a conformer for coordinates
    conf = Chem.Conformer(len(atoms))
    for i, atom in enumerate(atoms):
        conf.SetAtomPosition(i, (atom.x, atom.y, atom.z))
    mol.AddConformer(conf)

    return mol.GetMol()


def rdkit_to_gam(mol: 'Chem.Mol') -> List[GAM_Atom]:
    """
    Convert an RDKit Mol object to a list of GAM_Atom objects.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit Mol object to convert.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with the same atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not an RDKit Mol object.
    ImportError
        If RDKit is not installed.
    """
    if is_rdkit(mol) is False:
        raise TypeError("Input is not an RDKit Mol object.")
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is not available. Please install it first.")
    atoms = []
    conf = mol.GetConformer()
    for i, atom in enumerate(mol.GetAtoms(), start=1):
        pos = conf.GetAtomPosition(i-1)
        atoms.append(GAM_Atom(
            id=i,
            element=atom.GetSymbol(),
            x=pos.x,
            y=pos.y,
            z=pos.z
        ))
    return atoms





def gam_to_mdanalysis(atoms: List[GAM_Atom]) -> 'Universe':
    """
    Convert a list of GAM_Atom objects to an MDAnalysis Universe object.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of GAM_Atom objects to convert.

    Returns
    -------
    MDAnalysis.core.universe.Universe
        Universe object containing the atoms with positions and topology.

    Raises
    ------
    TypeError
        If the input is not a list of GAM_Atom objects.
    ImportError
        If MDAnalysis is not installed.
    """
    if is_gam(atoms) is False:
        raise TypeError("Input is not a list of GAM_Atom objects.")
    if not MDANALYSIS_AVAILABLE:
        raise ImportError("MDAnalysis is not available. Please install it first.")

    n_atoms = len(atoms)

    # Create an empty Universe with 1 residue
    u = Universe.empty(n_atoms=n_atoms,
                       n_residues=1,
                       atom_resindex=[0]*n_atoms,
                       trajectory=True)

    # Atom/residue info
    names = [atom.element for atom in atoms]
    ids = [atom.id for atom in atoms]

    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resids", [1])      # one residue only
    u.add_TopologyAttr("resnames", ["MOL"])
    u.add_TopologyAttr("ids", ids)

    # Coordinates
    positions = np.array([[atom.x, atom.y, atom.z] for atom in atoms], dtype=np.float32)
    u.atoms.positions = positions

    return u



def mdanalysis_to_gam(u: 'Universe') -> List[GAM_Atom]:
    """
    Convert an MDAnalysis Universe or AtomGroup to a list of GAM_Atom objects.

    Parameters
    ----------
    u : MDAnalysis Universe or AtomGroup
        Universe or AtomGroup to convert.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with the same atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not an MDAnalysis Universe or AtomGroup.
    ImportError
        If MDAnalysis is not installed.
    """
    if is_mdanalysis(u) is False:
        raise TypeError("Input is not an MDAnalysis Universe or AtomGroup object.")
    if not MDANALYSIS_AVAILABLE:
        raise ImportError("MDAnalysis is not available. Please install it first.")
    atoms = []
    for i, a in enumerate(u.atoms, start=1):
        x, y, z = a.position
        atoms.append(GAM_Atom(
            id=i,
            element=a.name,
            x=x,
            y=y,
            z=z
        ))
    return atoms







def gam_to_openbabel(atoms: List[GAM_Atom]) -> 'pybel.Molecule':
    """
    Convert a list of GAM_Atom objects to an OpenBabel pybel.Molecule object.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of GAM_Atom objects to convert.

    Returns
    -------
    pybel.Molecule
        OpenBabel molecule with atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not a list of GAM_Atom objects.
    ImportError
        If OpenBabel/pybel is not installed.
    """
    if is_gam(atoms) is False:
        raise TypeError("Input is not a list of GAM_Atom objects.")
    if not OPENBABEL_AVAILABLE:
        raise ImportError("OpenBabel is not available. Please install it first.")
    mol = ob.OBMol()
    for atom in atoms:
        a = mol.NewAtom()
        a.SetAtomicNum(ob.OBElementTable().GetAtomicNum(atom.element))
        a.SetVector(atom.x, atom.y, atom.z)
    return pybel.Molecule(mol)    

def openbabel_to_gam(mol: 'pybel.Molecule') -> List[GAM_Atom]:
    """
    Convert an OpenBabel pybel.Molecule object to a list of GAM_Atom objects.

    Parameters
    ----------
    mol : pybel.Molecule
        OpenBabel molecule to convert.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with the same atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not a pybel.Molecule object.
    ImportError
        If OpenBabel/pybel is not installed.
    """
    if is_openbabel(mol) is False:
        raise TypeError("Input is not an OpenBabel pybel.Molecule object.")
    if not OPENBABEL_AVAILABLE:
        raise ImportError("OpenBabel is not available. Please install it first.")
    atoms = []
    for atom in mol.atoms:
        atoms.append(GAM_Atom(
            id=atom.id,
            element=atom.element,
            x=atom.x,
            y=atom.y,
            z=atom.z
        ))
    return atoms        





def gam_to_pymatgen_structure(atoms: List[GAM_Atom]) -> Structure:
    """
    Convert a list of GAM_Atom objects to a pymatgen Structure object.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of GAM_Atom objects to convert.

    Returns
    -------
    pymatgen.core.structure.Structure
        Structure object representing the atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not a list of GAM_Atom objects.
    """
    if is_gam(atoms) is False:
        raise TypeError("Input is not a list of GAM_Atom objects.")
    symbols = [atom.element for atom in atoms]
    coords = [[atom.x, atom.y, atom.z] for atom in atoms]
    return Structure(symbols, coords)

def pymatgen_structure_to_gam(structure: Structure) -> List[GAM_Atom]:
    """
    Convert a pymatgen Structure object to a list of GAM_Atom objects.

    Parameters
    ----------
    structure : pymatgen Structure
        Pymatgen Structure object to convert.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with the same atomic symbols and coordinates.

    Raises
    ------
    TypeError
        If the input is not a pymatgen Structure object.
    """
    if is_pymatgen(structure) is False:
        raise TypeError("Input is not a pymatgen Structure object.")
    atoms = []
    for i, site in enumerate(structure.sites, start=1):
        atoms.append(GAM_Atom(
            id=i,
            element=site.species_string,
            x=site.x,
            y=site.y,
            z=site.z
        ))
    return atoms    



'''
def gam_to_mdtraj(atoms: List[GAM_Atom]) -> 'md.Trajectory':
    if is_gam(atoms) is False:
        raise TypeError("Input is not a list of GAM_Atom objects.")
    if not MDTRAJ_AVAILABLE:
        raise ImportError("MDTraj is not available. Please install it first.")
    return md.Trajectory(atoms)
    
    
def mdtraj_to_gam(traj: 'md.Trajectory') -> List[GAM_Atom]:
    if is_mdtraj(traj) is False:
        raise TypeError("Input is not an MDTraj Trajectory object.")
    if not MDTRAJ_AVAILABLE:
        raise ImportError("MDTraj is not available. Please install it first.")
    atoms = []
    for frame in traj:
        for atom in frame:
            atoms.append(GAM_Atom(
                id=atom.id,
                element=atom.element,
                x=atom.x,
                y=atom.y,
                z=atom.z
            ))
    return atoms    
'''






