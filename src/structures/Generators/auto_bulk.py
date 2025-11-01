from ..Primatom import GAM_Atom , GAM_Bond
import random
import numpy as np
from ase.build import bulk
from ase import Atoms
from typing import List, Optional, Dict, Tuple, Union, Any
import copy
from ase.visualize import view
import os
import ase
try:
    from ase import Atoms
    HAS_ASE = True
except ImportError:
    HAS_ASE = False
    
     

DEFAULT_LATTICE_DB = {
    'Cu': {'crystal_structure': 'fcc', 'a': 3.615},
    'Ni': {'crystal_structure': 'fcc', 'a': 3.52},
    'Co': {'crystal_structure': 'hcp', 'a': 2.51, 'c': 4.07},
    'Fe': {'crystal_structure': 'bcc', 'a': 2.866},
    'Al': {'crystal_structure': 'fcc', 'a': 4.05},
    'Ag': {'crystal_structure': 'fcc', 'a': 4.09},
    'Au': {'crystal_structure': 'fcc', 'a': 4.08},
    'Pb': {'crystal_structure': 'fcc', 'a': 4.95},
    'Mo': {'crystal_structure': 'bcc', 'a': 3.15},
    'W': {'crystal_structure': 'bcc', 'a': 3.165},
    'Ta': {'crystal_structure': 'bcc', 'a': 3.30},
    'Ti': {'crystal_structure': 'hcp', 'a': 2.95, 'c': 4.68},
    'Zn': {'crystal_structure': 'hcp', 'a': 2.66, 'c': 4.95},
    'Mg': {'crystal_structure': 'hcp', 'a': 3.21, 'c': 5.21},
    'Sn': {'crystal_structure': 'tetragonal', 'a': 5.83, 'c': 3.18},
    'Nb': {'crystal_structure': 'bcc', 'a': 3.30},
    'Pd': {'crystal_structure': 'fcc', 'a': 3.89},
    'Pt': {'crystal_structure': 'fcc', 'a': 3.92},
}



class AdvancedAlloys:

    """
    Generator for multi-element random alloys with customizable crystal structures and lattice constants.

    The `AdvancedAlloys` class allows the creation of *randomized substitutional alloys* 
    from multiple metallic or intermetallic elements, with user-defined 
    stoichiometric fractions, crystal types (e.g., fcc, bcc, hcp), and 
    lattice parameters. It can handle both simple and mixed lattice systems, 
    automatically constructing ASE-compatible supercells and mapping them 
    into GAM-compatible `GAM_Atom` objects.

    The class provides full control over alloy composition, randomization, 
    and reproducibility via random seeds. It also supports both isotropic 
    and anisotropic (hexagonal) lattice constants and offers 
    structure export and transformation operations (translation, rotation, etc.).

    Parameters
    ----------
    elements : list of str
        Chemical symbols of alloy components (e.g., ``['Cu', 'Ni', 'Co']``).
    fractions : list of float
        Fractional composition of each element. Must sum to 1.0.
    crystal_structures : list of str, optional
        Crystal structure type for each element 
        (e.g., ``['fcc', 'fcc', 'hcp']``).  
        If omitted, defaults are taken from the internal lattice database.
    lattice_constants : list of float or list of dict, optional
        Lattice constants for each element.  
        Each entry can be:
            - A float (for cubic lattices), e.g. ``3.615``  
            - A dict with keys like ``{'a': 2.95, 'c': 4.68}`` for hexagonal lattices  
        If omitted, values are retrieved from the internal lattice database.
    supercell_size : tuple of int, default=(3, 3, 3)
        Supercell replication factors along the (x, y, z) directions.
    seed : int, optional
        Random seed for reproducible alloy configurations.
    verbose : bool, default=True
        If True, prints diagnostic information during structure generation.
    metadata : dict, optional
        Additional metadata for traceability or simulation setup.

    Attributes
    ----------
    atoms : list[GAM_Atom]
        List of atomic objects in GAM-compatible format.
    ASE_atoms : ase.Atoms
        The underlying ASE `Atoms` object containing the generated structure.
    bonds : list[GAM_Bond]
        List of bonds (optional; currently not populated by default).
    elements : list[str]
        Element symbols defining the alloy system.
    fractions : list[float]
        Fractional composition of elements.
    crystal_structures : list[str]
        Crystal structure assigned to each element.
    lattice_constants : list[float or dict]
        Lattice constants for each component.
    supercell_size : tuple[int, int, int]
        Supercell replication dimensions.
    meta : dict
        Dictionary storing all structure metadata.
    origin : tuple[float, float]
        Origin coordinates for translation and rotation.
    rotation : float
        Rotation angle of the system (degrees).
    
    Methods
    -------
    _set_defaults()
        Set default crystal structures and lattice constants based on a predefined database.
    _build_advanced_structure() -> ase.Atoms
        Build the random alloy supercell using ASE’s bulk crystal generator.
    translate(dx, dy, dz=0.0)
        Translate all atoms by a displacement vector.
    rotate(angle_deg, about_center=True)
        Rotate the alloy structure around the z-axis.
    copy()
        Return a deep copy of the structure.
    get_atoms() -> list[GAM_Atom]
        Retrieve a deep copy of all atoms in GAM format.
    get_positions() -> list[tuple[float, float, float]]
        Return Cartesian coordinates of all atoms.
    get_elements() -> list[str]
        Return element symbols of all atoms.
    to_xyz(path: str)
        Save the generated structure to an XYZ file.
    to_ase() -> ase.Atoms
        Export the structure as an ASE `Atoms` object.
    
    Notes
    -----
    - If no crystal structure or lattice constant is provided for an element, 
      the class retrieves default values from an internal lattice database 
      (`DEFAULT_LATTICE_DB`).
    - All elements are randomly assigned to atomic sites in the supercell 
      according to the specified composition fractions.
    - Hexagonal lattices must define both 'a' and 'c' lattice constants.
    - The randomization is deterministic if a `seed` is provided.
    - The final composition may vary slightly due to integer rounding 
      of atom counts in the finite supercell.

    Examples
    --------
    >>> from pygamlab import AdvancedAlloys

    >>> # 1. Create a ternary fcc random alloy
    >>> alloy = AdvancedAlloys(
    ...     elements=['Cu', 'Ni', 'Co'],
    ...     fractions=[0.4, 0.3, 0.3],
    ...     crystal_structures=['fcc', 'fcc', 'fcc'],
    ...     lattice_constants=[3.61, 3.52, 3.55],
    ...     supercell_size=(4, 4, 4),
    ...     seed=42
    ... )
    >>> alloy.to_xyz("CuNiCo.xyz")

    >>> # 2. Create an hcp binary alloy with anisotropic lattice constants
    >>> alloy = AdvancedAlloys(
    ...     elements=['Ti', 'Al'],
    ...     fractions=[0.7, 0.3],
    ...     crystal_structures=['hcp', 'hcp'],
    ...     lattice_constants=[{'a': 2.95, 'c': 4.68}, {'a': 2.86, 'c': 4.57}],
    ...     supercell_size=(3, 3, 2)
    ... )
    >>> ase_atoms = alloy.to_ase()

    >>> # 3. Randomized quaternary alloy (defaults to fcc)
    >>> alloy = AdvancedAlloys(['Fe', 'Ni', 'Cr', 'Mn'], [0.25, 0.25, 0.25, 0.25])
    >>> alloy.to_xyz("FeNiCrMn.xyz")
    >>> print(len(alloy.atoms), "atoms generated.")
    """
    
    def __init__(
        self,
        elements: List[str],
        fractions: List[float],
        crystal_structures: Optional[List[str]] = None,
        lattice_constants: Optional[List[Union[float, Dict[str, float]]]] = None,
        supercell_size: Tuple[int, int, int] = (3, 3, 3),
        seed: Optional[int] = None,
        verbose: bool = True,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize advanced alloy structure with given parameters.

        Args:
            elements: List of chemical symbols, e.g. ['Cu', 'Ni', 'Co']
            fractions: List of fractions, must sum to 1.0
            crystal_structures: Optional list matching elements, e.g. ['fcc', 'fcc', 'hcp'].
                If None, use default from database or 'fcc'.
            lattice_constants: Optional list matching elements.
                Each entry can be a float (a), or dict with keys like {'a':.., 'c':..}.
                If None, use defaults.
            supercell_size: Tuple of ints for supercell replication.
            seed: Random seed for reproducibility.
            verbose: Print info.
            metadata: Optional metadata dictionary
        """
        # Required attributes
        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []
        self.metadata = metadata or {}
        
        # Store parameters
        self.elements = elements
        self.fractions = fractions
        self.crystal_structures = crystal_structures
        self.lattice_constants = lattice_constants
        self.supercell_size = supercell_size
        self.seed = seed
        self.verbose = verbose
        self.origin = (0.0, 0.0)   # default origin for translations/rotations
        self.rotation = 0.0  
        
        self.meta={
            'elements':elements ,
            'fractions':fractions ,
            'crystal_structures':crystal_structures ,
            'lattice_constants':lattice_constants ,
            'supercell_size':supercell_size ,
            'seed':seed ,
            'verbose':verbose }

        # Validation
        if len(elements) != len(fractions):
            raise ValueError("Length of elements and fractions must be equal.")
        if abs(sum(fractions) - 1.0) > 1e-6:
            raise ValueError("Fractions must sum to 1.0")

        if seed is not None:
            random.seed(seed)

        # Set defaults for crystal_structures and lattice_constants if not given
        self._set_defaults()
        
        # Build the structure
        self.ASE_atoms = self._build_advanced_structure()
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        
        

    def _set_defaults(self):
        """Set default crystal structures and lattice constants."""
        if self.crystal_structures is None:
            self.crystal_structures = []
            for el in self.elements:
                cs = DEFAULT_LATTICE_DB.get(el, {}).get('crystal_structure', 'fcc')
                self.crystal_structures.append(cs)
        if self.lattice_constants is None:
            self.lattice_constants = []
            for el in self.elements:
                lc = DEFAULT_LATTICE_DB.get(el, {}).get('a', None)
                self.lattice_constants.append(lc)
    
    def _build_advanced_structure(self) -> Atoms:
        """Build the advanced ASE atoms structure."""
        # For simplicity, build base structure from first element
        base_el = self.elements[0]
        base_cs = self.crystal_structures[0]
        base_lc = self.lattice_constants[0]

        if base_lc is None:
            raise ValueError(f"Lattice constant for base element {base_el} not specified and not found in database.")

        # Build base bulk
        if isinstance(base_lc, dict):
            # handle dict lattice constants (e.g. hexagonal)
            if base_cs.lower() == 'hcp' or base_cs.lower() == 'hexagonal':
                a = base_lc.get('a')
                c = base_lc.get('c')
                if a is None or c is None:
                    raise ValueError(f"Hexagonal lattice constants require both 'a' and 'c'")
                # hexagonal unit cell with ASE bulk
                atoms = bulk(base_el, 'hcp', a=a, c=c)
            else:
                raise NotImplementedError(f"Dict lattice constants not supported for crystal structure {base_cs}")
        else:
            atoms = bulk(base_el, base_cs, a=base_lc)

        # Create supercell
        atoms = atoms.repeat(self.supercell_size)

        total_atoms = len(atoms)
        num_atoms_per_element = [int(round(frac * total_atoms)) for frac in self.fractions]

        # Fix rounding issues
        diff = total_atoms - sum(num_atoms_per_element)
        if diff != 0:
            max_idx = num_atoms_per_element.index(max(num_atoms_per_element))
            num_atoms_per_element[max_idx] += diff

        # Assign elements randomly
        indices = list(range(total_atoms))
        random.shuffle(indices)

        start = 0
        for el, count in zip(self.elements, num_atoms_per_element):
            for i in indices[start:start+count]:
                atoms[i].symbol = el
            start += count
            
        if self.verbose:
            print(f"Generated alloy with composition:")
            for el, num in zip(self.elements, num_atoms_per_element):
                print(f"  {el}: {num} atoms ({num/total_atoms:.2%})")
            print(f"Total atoms: {total_atoms}")
            print(f"Crystal structure based on {base_el}: {base_cs} with lattice constant {base_lc}")
            print(f"Supercell size: {self.supercell_size}")

        return atoms
    
    
    def translate(self, dx: float, dy: float, dz: float = 0.0):
        """Translate all atoms by (dx, dy, dz)."""
        for atom in self.atoms:
            atom.x += dx
            atom.y += dy
            atom.z += dz
        
        # Update ASE atoms positions as well
        for i, atom in enumerate(self.ASE_atoms):
            atom.position += np.array([dx, dy, dz])

        # Update origin in metadata
        ox, oy = self.origin
        self.origin = (ox + dx, oy + dy)
        self.meta["origin"] = self.origin
        
        #return self.atoms
        
    
    def rotate(self, angle_deg: float, about_center: bool = True):
        """Rotate structure about z-axis."""
        angle_rad = np.radians(angle_deg)
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
        
        # Get center of mass
        positions = np.array([(atom.x, atom.y, atom.z) for atom in self.atoms])
        center = np.mean(positions, axis=0)
        
        if about_center:
            # Rotate about center of mass
            for atom in self.atoms:
                x_rel = atom.x - center[0]
                y_rel = atom.y - center[1]
                atom.x = center[0] + x_rel * cos_a - y_rel * sin_a
                atom.y = center[1] + x_rel * sin_a + y_rel * cos_a
        else:
            # Rotate about (0, 0)
            for atom in self.atoms:
                x_old, y_old = atom.x, atom.y
                atom.x = x_old * cos_a - y_old * sin_a
                atom.y = x_old * sin_a + y_old * cos_a
        
        # Update ASE atoms positions as well
        for i, atom in enumerate(self.ASE_atoms):
            atom.position[0] = self.atoms[i].x
            atom.position[1] = self.atoms[i].y
            atom.position[2] = self.atoms[i].z

        # Update rotation in metadata
        self.rotation = (self.rotation + angle_deg) % 360
        self.meta["rotation"] = self.rotation
        #return self.atoms
    
    def copy(self):
        """Create a deep copy of the structure."""
        return copy.deepcopy(self.atoms)
    
    
    def get_atoms(self):
        """
        Return a copy of the atoms in this molecule.

        Returns
        -------
        List[str]
            A list of atom symbols.
        """
        return copy.deepcopy(self.atoms)
    
    def get_positions(self) -> List[Tuple[float, float, float]]:
        """Return positions of all atoms as a list of (x, y, z) tuples."""
        return [atom.get_position() for atom in self.atoms]

    def get_elements(self) -> List[str]:
        """Return list of element symbols in the molecule."""
        return [atom.element for atom in self.atoms]
    

            
            
    def _to_GAM_Atoms(self,atoms):
        return [
            GAM_Atom(
                id=atom.index,
                element=atom.symbol,
                x=atom.position[0],
                y=atom.position[1],
                z=atom.position[2]
            )
            for atom in atoms
        ]
    
   
    
    def to_xyz(self, path: str) -> None:
        """Build (if needed) and save structure directly to XYZ file."""
        
        # Build only when saving
        if not self.atoms:
            raise ValueError("No atoms available to save.")
        
        # Prepare XYZ string directly
        lines = [str(len(self.atoms))]
        '''
        comment_parts = []
        for key, value in self.meta.items():
            comment_parts.append(f"{key}={value}")
        lines.append("Silicene structure: " + ", ".join(comment_parts))
        '''
        lines.append('Generated by PyGamlab')
        for atom in self.atoms:
            lines.append(f"{atom.element:2s} {atom.x:12.6f} {atom.y:12.6f} {atom.z:12.6f}")
        
        # Save directly
        with open(path, 'w') as f:
            f.write("\n".join(lines))

        if self.verbose:
            print(f"Structure saved to {path} (XYZ format)")
        
        

    def to_ase(self) -> Optional['ase.Atoms']:
        """
        Export structure to ASE Atoms object (if ASE is available).
        
        Returns:
            ASE Atoms object or None if ASE not available
        """
       
        if not HAS_ASE:
            print("ASE not available. Install with: pip install ase")
            return None
        
        if not self.atoms:
            return Atoms()
        
        symbols = [atom.element for atom in self.atoms]
        positions = [atom.position for atom in self.atoms]
        
        atoms_obj = Atoms(symbols=symbols, positions=positions)
        
        # Add metadata
        atoms_obj.info.update(self.meta)
        
        #or 
        #Just 
        #return self.ASE_atoms
        
        return atoms_obj
    
    
    


