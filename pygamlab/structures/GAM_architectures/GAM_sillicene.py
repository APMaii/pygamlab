from ..Primatom import GAM_Atom , GAM_Bond , GAM_Molecule
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import copy
import io
from typing import List, Tuple, Optional, Union, Dict, Any


# Optional imports
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

try:
    import ase
    from ase import Atoms
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False




def silicene_unit_cell(lattice_constant: float = 3.87) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """
    Generate silicene unit cell with 2-atom basis and lattice vectors.
    
    Args:
        lattice_constant: Lattice constant in Angstroms
    
    Returns:
        Tuple of (basis_positions, lattice_vectors)
    """
    # Silicene has hexagonal lattice like graphene
    a = lattice_constant
    
    # Lattice vectors for hexagonal cell
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([a/2, a*np.sqrt(3)/2, 0.0])
    lattice_vectors = [a1, a2]
    
    # Two-atom basis (A and B sublattices)
    basis = [
        np.array([0.0, 0.0, 0.0]),  # A sublattice
        np.array([a/2, a/(2*np.sqrt(3)), 0.0])  # B sublattice
    ]
    
    return basis, lattice_vectors


class Silicene:
    """
    Class for constructing and manipulating monolayer silicene nanostructures.

    This class provides a flexible generator for 2D silicene sheets and flakes,
    incorporating its characteristic buckled honeycomb lattice. Silicene, the
    silicon analogue of graphene, exhibits a low-buckled geometry due to the
    larger atomic radius of silicon, leading to mixed sp²–sp³ hybridization.
    Its electronic structure, mechanical flexibility, and tunable bandgap make
    it a key candidate for post-silicon electronics and 2D semiconductor devices.

    The `Silicene` class supports systematic construction of finite or extended
    structures with customizable lattice parameters, edge orientations, and
    buckling amplitude. It also provides methods for geometric transformations
    and bond detection, making it suitable for atomistic simulations, DFT input
    generation, and molecular dynamics preprocessing.

    Supported features include:
        - Finite nanoflake or extended monolayer generation
        - Armchair or zigzag edge terminations
        - Buckled atomic geometry (adjustable buckling amplitude)
        - Rotation and translation transformations
        - Distance-based bond detection
        - Export to XYZ and ASE-compatible structures

    The implementation can be extended to:
        - Multilayer silicene and silicene heterostructures
        - Defect engineering (vacancies, dopants)
        - Strain-dependent bandgap analysis
        - Functionalization or oxidation studies

    Parameters
    ----------
    lattice_constant : float, optional
        Lattice constant of silicene (Å). Default is 3.87 Å.
    width : float, optional
        Width of the structure in the armchair direction (Å). Default is 10.0 Å.
    length : float, optional
        Length of the structure in the zigzag direction (Å). Default is 10.0 Å.
    edge_type : {'zigzag', 'armchair'}, optional
        Edge orientation of the finite flake. Default is 'armchair'.
    rotation : float, optional
        Rotation angle of the entire structure about the z-axis (degrees). Default is 0.0°.
    vacuum : float, optional
        Out-of-plane vacuum spacing for periodic systems (Å). Default is 15.0 Å.
    buckling_height : float, optional
        Out-of-plane displacement between A and B sublattices (Å). Default is 0.44 Å.
    origin : tuple of float, optional
        (x, y) coordinate shift applied after construction. Default is (0.0, 0.0).

    Attributes
    ----------
    atoms : list of GAM_Atom
        List of Si atoms in the structure.
    bonds : list of GAM_Bond
        List of Si–Si bonds identified using a distance cutoff.
    lattice_vectors : list of numpy.ndarray
        Lattice vectors defining the silicene unit cell geometry.
    meta : dict
        Metadata containing structure parameters and descriptors.
    origin : tuple of float
        Current origin of the flake in Cartesian coordinates.

    Methods
    -------
    build() -> None
        Construct the silicene lattice and populate atom and bond lists.
    nearest_neighbors(r_cut: float | None = None) -> None
        Identify nearest neighbors based on a distance cutoff and store bonds.
    translate(dx: float, dy: float, dz: float = 0.0) -> None
        Translate all atomic coordinates by the specified vector.
    rotate(angle_deg: float, about_center: bool = True) -> None
        Rotate the silicene structure about the z-axis.
    copy() -> list[GAM_Atom]
        Return a deep copy of all atoms in the structure.
    get_atoms() -> list[GAM_Atom]
        Retrieve a copy of all atomic objects in the structure.
    get_positions() -> list[tuple[float, float, float]]
        Return atomic positions as a list of Cartesian coordinates.
    get_elements() -> list[str]
        Return the list of element symbols (typically all 'Si').
    to_xyz(path: str) -> None
        Export the silicene structure to an `.xyz` file.
    to_ase() -> ase.Atoms
        Convert the structure to an ASE `Atoms` object for use in simulation workflows.

    Notes
    -----
    - Silicene adopts a buckled honeycomb lattice with two atoms per unit cell.
    - Typical Si–Si bond length is ~2.28 Å, and the buckling height is ~0.44 Å.
    - The default cutoff for nearest-neighbor detection corresponds to roughly
      1.6× the Si–Si bond length.
    - The armchair and zigzag directions are rotated by 30° relative to each
      other, similar to graphene but with buckled atomic planes.
    - The `vacuum` spacing should be increased (≥15 Å) when simulating isolated layers
      to prevent interlayer interactions in periodic boundary conditions.

    Examples
    --------
    >>> s = Silicene(width=15.0, length=12.0, edge_type='zigzag', buckling_height=0.44)
    >>> s.build()
    >>> s.to_xyz("silicene_flake.xyz")
    >>> atoms = s.to_ase()
    >>> print(f"Silicene structure with {len(atoms)} atoms successfully generated.")
    """
    
    def __init__(self, lattice_constant: float = 3.87, width: float = 10.0, 
                 length: float = 10.0, edge_type: str = 'armchair', 
                 rotation: float = 0.0, vacuum: float = 15.0, 
                 buckling_height: float = 0.44, origin: Tuple[float, float] = (0.0, 0.0)):
        """
        Initialize silicene structure parameters.
        
        Args:
            lattice_constant: Lattice constant in Angstroms (~3.87 Å)
            width: Width of bounding box in Angstroms
            length: Length of bounding box in Angstroms
            edge_type: 'armchair' or 'zigzag' edge termination
            rotation: Rotation angle in degrees about z-axis
            vacuum: Z-spacing for periodic cells in Angstroms
            buckling_height: Intrinsic buckling displacement in Angstroms (~0.44 Å)
            origin: (x0, y0) coordinate shift
            periodic: If True, create periodic cell; if False, finite flake
        """
        self.lattice_constant = lattice_constant
        self.width = width
        self.length = length
        self.edge_type = edge_type.lower()
        self.rotation = rotation
        self.vacuum = vacuum
        self.buckling_height = buckling_height
        self.origin = origin
        
        # Initialize structure containers
        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []
        self.lattice_vectors: Optional[List[np.ndarray]] = None
        
        # Metadata dictionary
        self.meta: Dict[str, Any] = {
            'lattice_constant': lattice_constant,
            'width': width,
            'length': length,
            'edge_type': edge_type,
            'rotation': rotation,
            'vacuum': vacuum,
            'buckling_height': buckling_height,
            'origin': origin
        }
        
        # Build structure
        self.build()
    
    def build(self) -> None:
        """Construct the silicene flake with buckling."""
        self.atoms.clear()
        self.bonds.clear()
        
        # Get unit cell
        basis, lattice_vectors = silicene_unit_cell(self.lattice_constant)
        
        #if self.periodic:
         #   self.lattice_vectors = [vec.copy() for vec in lattice_vectors]
        
        # Determine tiling dimensions
        a1, a2 = lattice_vectors
        n1 = max(1, int(np.ceil(self.width / np.linalg.norm(a1))))
        n2 = max(1, int(np.ceil(self.length / np.linalg.norm(a2))))
        
        atom_id = 0
        positions = []
        
        # Generate atoms by tiling unit cells
        for i in range(-n1, n1 + 1):
            for j in range(-n2, n2 + 1):
                cell_origin = i * a1 + j * a2
                
                for k, basis_pos in enumerate(basis):
                    pos = cell_origin + basis_pos
                    
                    # Check if position is within bounding box (for finite flakes)
                    #if not self.periodic:
                     #   if (abs(pos[0]) > self.width/2 or abs(pos[1]) > self.length/2):
                     #       continue
                    
                    positions.append((pos, k))  # k indicates sublattice (0=A, 1=B)
        
        # Apply rotation if specified
        if self.rotation != 0.0:
            angle_rad = np.radians(self.rotation)
            cos_theta = np.cos(angle_rad)
            sin_theta = np.sin(angle_rad)
            rotation_matrix = np.array([
                [cos_theta, -sin_theta, 0],
                [sin_theta, cos_theta, 0],
                [0, 0, 1]
            ])
        
        # Create atoms with buckling
        for pos, sublattice in positions:
            if self.rotation != 0.0:
                pos = rotation_matrix @ pos
            
            # Apply buckling: A sublattice up, B sublattice down
            z_offset = self.buckling_height/2 if sublattice == 0 else -self.buckling_height/2
            
            # Apply origin shift
            x = pos[0] + self.origin[0]
            y = pos[1] + self.origin[1]
            z = pos[2] + z_offset
            
            atom = GAM_Atom(atom_id, "Si", x, y, z, {"sublattice": sublattice})
            self.atoms.append(atom)
            atom_id += 1
        
        # Generate bonds
        self.nearest_neighbors()
    
    def nearest_neighbors(self, r_cut: Optional[float] = None) -> None:
        """
        Build bonds using distance cutoff.
        
        Args:
            r_cut: Cutoff radius for bond detection (default: 2.28 * 1.6 Å)
        """
        if r_cut is None:
            r_cut = 2.28 * 1.6  # Typical Si-Si bond length × safety factor
        
        self.bonds.clear()
        bond_id = 0
        
        for i, atom1 in enumerate(self.atoms):
            for j, atom2 in enumerate(self.atoms[i+1:], i+1):
                pos1 = atom1.position
                pos2 = atom2.position
                distance = np.linalg.norm(pos2 - pos1)
                
                if distance <= r_cut:
                    bond = GAM_Bond(bond_id, atom1.id, atom2.id)
                    self.bonds.append(bond)
                    bond_id += 1
    
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


        
    
    def copy(self) -> 'Silicene':
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
        

        

    
    def to_ase(self) -> 'ase.Atoms':
        """Export to ASE Atoms object if available."""
        if not ASE_AVAILABLE:
            raise ImportError("ASE not available. Install with: pip install ase")
        
        positions = []
        symbols = []
        
        for atom in self.atoms:
            positions.append([atom.x, atom.y, atom.z])
            symbols.append(atom.element)
        
        atoms = Atoms(symbols=symbols, positions=positions)
        
        #if self.lattice_vectors and self.periodic:
            # Set up periodic cell
        #    cell = np.array(self.lattice_vectors + [[0, 0, self.vacuum]])
        #    atoms.set_cell(cell)
        #    atoms.set_pbc([True, True, False])
        
        return atoms








