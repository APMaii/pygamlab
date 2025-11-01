from ..Primatom import GAM_Atom , GAM_Bond 
#from Generators.gam_atom import GAM_Atom , GAM_Bond
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import copy
from typing import List, Tuple, Union, Optional, Dict, Any
import io


# Optional imports
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

try:
    import ase
    from ase import Atoms
    HAS_ASE = True
except ImportError:
    HAS_ASE = False







def graphene_unit_cell(lattice_constant: float = 2.46) -> Tuple[np.ndarray, np.ndarray]:

    """
    Flexible graphene nanostructure generator.
    
    Supports both periodic unit cells and finite flakes with various edge types,
    rotations, and coordinate transformations. Designed for extensibility to
    multilayers, defects, and strain.
    """
    


    
    a = lattice_constant
    
    # Two-atom basis in graphene unit cell
    basis = np.array([
        [0.0, 0.0, 0.0],                    # First carbon
        [a / 2.0, a * np.sqrt(3) / 6.0, 0.0]  # Second carbon
    ])
    
    # Lattice vectors for hexagonal lattice
    lattice_vectors = np.array([
        [a, 0.0, 0.0],                      # a1
        [a / 2.0, a * np.sqrt(3) / 2.0, 0.0]  # a2
    ])
    
    return basis, lattice_vectors






class Graphene:
    """
    Class for constructing and manipulating monolayer graphene nanostructures.

    This class provides a flexible framework for generating both periodic graphene
    unit cells and finite graphene flakes with customizable geometries, edge
    terminations, and orientations. The implementation is suitable for atomistic
    modeling, visualization, and export to external simulation packages (e.g.,
    ASE, LAMMPS, Quantum ESPRESSO).

    The `Graphene` class supports:
        - Generation of primitive unit cells or finite flakes with armchair or
          zigzag edges.
        - Parameterized lattice constant, width, and length.
        - Structural transformations such as translation and rotation.
        - Automatic bond construction via nearest-neighbor detection.
        - Export of atomic configurations to `.xyz` files or ASE objects.
    
    This class is designed to be extendable for modeling:
        - Multilayer graphene and heterostructures
        - Point or line defects (vacancies, dopants, grain boundaries)
        - Strained or rotated configurations (e.g., twisted bilayers)
        - Functionalized graphene or 2D carbon allotropes

    Parameters
    ----------
    lattice_constant : float, optional
        Lattice constant of graphene in Ångströms (default: 2.46 Å).
    width : float, optional
        Width of the graphene flake (in Å) along the x-direction (default: 10.0 Å).
    length : float, optional
        Length of the graphene flake (in Å) along the y-direction (default: 10.0 Å).
    edge_type : str, optional
        Edge termination type, either ``'armchair'`` or ``'zigzag'`` (default: 'armchair').
    rotation : float, optional
        Rotation angle of the entire structure (in degrees) about the z-axis (default: 0.0°).
    vacuum : float, optional
        Extra spacing in the z-direction (useful for periodic boundary simulations) (default: 10.0 Å).
    origin : tuple of float, optional
        (x, y) coordinates of the origin shift (default: (0.0, 0.0)).

    Attributes
    ----------
    atoms : list of GAM_Atom
        List of carbon atoms in the graphene structure.
    bonds : list of GAM_Bond
        List of covalent C–C bonds identified via nearest-neighbor search.
    lattice_vectors : ndarray of shape (2, 3)
        Lattice vectors of the graphene unit cell.
    meta : dict
        Metadata dictionary storing structure parameters and build status.
    origin : tuple of float
        Cartesian coordinates of the structure’s reference origin.

    Methods
    -------
    build() -> Graphene
        Constructs the graphene lattice based on current parameters.
    translate(dx, dy, dz=0.0)
        Translates all atoms by the specified displacement vector.
    rotate(angle_deg, about_center=True)
        Rotates the structure around the z-axis.
    nearest_neighbors(r_cut=None) -> Graphene
        Generates C–C bonds using a cutoff-based nearest-neighbor criterion.
    get_atoms() -> list[GAM_Atom]
        Returns a deep copy of the atomic list.
    get_positions() -> list[tuple[float, float, float]]
        Returns a list of atomic coordinates.
    get_elements() -> list[str]
        Returns a list of element symbols (typically all 'C').
    copy() -> Graphene
        Returns a deep copy of the current graphene structure.
    to_xyz(path: str)
        Writes the atomic configuration to an XYZ file.
    to_ase() -> ase.Atoms or None
        Exports the structure as an ASE `Atoms` object (if ASE is installed).

    Notes
    -----
    - The graphene lattice constant (2.46 Å) corresponds to a C–C bond length
      of approximately 1.42 Å.
    - The structure is generated in the xy-plane with z = 0 for all atoms.
    - Rotation and translation are applied after flake construction.
    - For periodic systems, the vacuum spacing ensures minimal z-overlap between layers.

    Examples
    --------
    >>> g = Graphene(width=20.0, length=15.0, edge_type='zigzag', rotation=15)
    >>> g.build()
    >>> g.to_xyz('graphene_flake.xyz')
    >>> atoms = g.to_ase()
    >>> print(len(atoms), "atoms in graphene flake")
    """
    
    def __init__(self, lattice_constant: float = 2.46, width: float = 10.0, 
                 length: float = 10.0, edge_type: str = 'armchair', 
                 rotation: float = 0.0, vacuum: float = 10.0, 
                 origin: Tuple[float, float] = (0.0, 0.0)):
        """
        Initialize graphene structure parameters.
        
        Args:
            lattice_constant: Lattice constant in Å (default: 2.46 Å)
            width: Width of bounding box in Å
            length: Length of bounding box in Å
            edge_type: 'armchair' or 'zigzag' edge termination
            rotation: Rotation angle in degrees about z-axis
            vacuum: Extra spacing in z-direction in Å
            origin: (x, y) coordinate shift
            periodic: If True, create primitive cell; if False, create finite flake
        """
        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []
        self.lattice_vectors: Optional[np.ndarray] = None
        self.origin=origin
        
        self.width=width
        self.lattice_constant=lattice_constant
        self.edge_type=edge_type
        self.rotation=rotation
        self.vacuum=vacuum
        self.origin=origin
        
        
        
        
        # Store parameters in metadata
        self.meta = {
            'lattice_constant': lattice_constant,
            'width': width,
            'length': length,
            'edge_type': edge_type,
            'rotation': rotation,
            'vacuum': vacuum,
            'origin': origin,
            'built': False
        }
        
        self.build()


    
    def build(self) -> 'Graphene':
        """
        Construct the graphene structure based on parameters.
        
        Returns:
            Self for method chaining
        """
        # Clear existing structure
        self.atoms.clear()
        self.bonds.clear()
        
        # Get unit cell
        basis, lattice_vectors = graphene_unit_cell(self.meta['lattice_constant'])
        
        #if self.meta['periodic']:
            # Create primitive unit cell
        #    self.lattice_vectors = lattice_vectors
        #    for i, pos in enumerate(basis):
        #        atom = GAM_Atom(i, "C", pos[0], pos[1], pos[2])
        #        self.atoms.append(atom)
        #else:
            # Create finite flake
        self._build_finite_flake(basis, lattice_vectors)
        
        # Apply rotation if specified
        if self.meta['rotation'] != 0.0:
            self.rotate(self.meta['rotation'], about_origin=True)
        
        # Apply origin shift
        if self.meta['origin'] != (0.0, 0.0):
            self.translate(self.meta['origin'][0], self.meta['origin'][1], 0.0)
        
        # Generate bonds based on nearest neighbors
        self.nearest_neighbors()
        
        self.meta['built'] = True
        
    
    def _build_finite_flake(self, basis: np.ndarray, lattice_vectors: np.ndarray):
        """Build finite graphene flake by tiling and trimming."""
        a1, a2 = lattice_vectors
        width, length = self.meta['width'], self.meta['length']
        edge_type = self.meta['edge_type']
        
        # Determine tiling parameters based on edge type
        if edge_type == 'armchair':
            # For armchair edges, align with x-axis
            n1 = int(np.ceil(width / np.linalg.norm(a1))) + 2
            n2 = int(np.ceil(length / np.linalg.norm(a2))) + 2
        else:  # zigzag
            # For zigzag edges, rotate lattice by 30 degrees conceptually
            n1 = int(np.ceil(width / np.linalg.norm(a1))) + 2
            n2 = int(np.ceil(length / np.linalg.norm(a2))) + 2
        
        # Generate atoms by tiling
        atom_id = 0
        positions = []
        
        for i in range(-n1//2, n1//2 + 1):
            for j in range(-n2//2, n2//2 + 1):
                for k, base_pos in enumerate(basis):
                    pos = base_pos + i * a1 + j * a2
                    positions.append(pos)
        
        positions = np.array(positions)
        
        # Apply edge-type specific orientation
        if edge_type == 'zigzag':
            # Rotate by 30 degrees for zigzag edges
            angle = np.radians(30)
            rotation_matrix = np.array([
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1]
            ])
            positions = positions @ rotation_matrix.T
        
        # Trim to desired size
        x_min, x_max = -width/2, width/2
        y_min, y_max = -length/2, length/2
        
        mask = ((positions[:, 0] >= x_min) & (positions[:, 0] <= x_max) &
                (positions[:, 1] >= y_min) & (positions[:, 1] <= y_max))
        
        valid_positions = positions[mask]
        
        # Create atoms
        for i, pos in enumerate(valid_positions):
            atom = GAM_Atom(i, "C", pos[0], pos[1], pos[2])
            self.atoms.append(atom)
    
    def nearest_neighbors(self, r_cut: Optional[float] = None) -> 'Graphene':
        """
        Build bonds based on nearest neighbor distances.
        
        Args:
            r_cut: Distance cutoff for bonding (default: 1.6 * C-C bond length)
        
        Returns:
            Self for method chaining
        """
        if r_cut is None:
            r_cut = 1.6 * 1.42  # ~1.6 times C-C bond length
        
        self.bonds.clear()
        bond_id = 0
        
        for i, atom1 in enumerate(self.atoms):
            for j, atom2 in enumerate(self.atoms[i+1:], i+1):
                distance = np.linalg.norm(atom1.position - atom2.position)
                if distance <= r_cut:
                    bond = GAM_Bond(bond_id, atom1.id, atom2.id)
                    self.bonds.append(bond)
                    bond_id += 1
        
        return self
    
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

    def copy(self) -> 'Graphene':
        """Create deep copy of the structure."""
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
    

    def to_xyz(self, path: str) -> None:
        """Build (if needed) and save structure directly to XYZ file."""
        
        # Build only when saving
        if not self.meta.get('built', False):
            self.build()
        
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
        
        return atoms_obj








    
    
    
    
    
    
    
 
    
    