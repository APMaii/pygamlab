from ..Primatom import GAM_Atom , GAM_Bond 
#from Generators.gam_atom import GAM_Atom , GAM_Bond
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import copy
import io
from typing import List, Tuple, Union, Optional, Dict, Any

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




def phosphorene_unit_cell(lattice_a: float = 3.313, lattice_b: float = 4.376, 
                         puckering_height: float = 2.14) -> Tuple[List[Tuple], np.ndarray]:
    """
    Generate phosphorene unit cell with 4 atoms in puckered arrangement.
    
    Args:
        lattice_a: Lattice constant along armchair direction (Å)
        lattice_b: Lattice constant along zigzag direction (Å)  
        puckering_height: Out-of-plane puckering distance (Å)
    
    Returns:
        basis: List of (x, y, z, sublayer) tuples for 4 atoms
        lattice_vectors: 2x2 array of lattice vectors
    """
    # Orthorhombic unit cell with 4 atoms
    # Two sublayers with alternating puckering
    basis = [
        (0.0, 0.0, puckering_height/2, 0),           # Upper sublayer
        (lattice_a/2, 0.0, -puckering_height/2, 1), # Lower sublayer
        (0.0, lattice_b/2, -puckering_height/2, 1), # Lower sublayer
        (lattice_a/2, lattice_b/2, puckering_height/2, 0)  # Upper sublayer
    ]
    
    lattice_vectors = np.array([
        [lattice_a, 0.0],
        [0.0, lattice_b]
    ])
    
    return basis, lattice_vectors


class Phosphorene:
    """
    Class for constructing and manipulating monolayer phosphorene nanostructures.

    This class provides a general framework for generating monolayer black
    phosphorus (phosphorene) in both periodic and finite geometries. Phosphorene,
    unlike graphene, exhibits an anisotropic puckered structure and a
    direction-dependent electronic bandgap, making it essential to correctly
    handle its orthorhombic lattice and out-of-plane buckling.

    The `Phosphorene` class allows for systematic construction of nanoribbons,
    flakes, or extended sheets with user-defined edge types, rotation angles,
    and lattice dimensions. It is designed for integration with atomistic
    modeling workflows (e.g., DFT, MD, or TB simulations).

    Supported features include:
        - Generation of finite monolayer flakes or extended sheets
        - Armchair or zigzag edge terminations
        - Rotation and translation transformations
        - Puckering height adjustment for structural realism
        - Bond generation based on nearest-neighbor distance criteria
        - Export to XYZ and ASE-compatible structures

    The implementation is easily extendable to:
        - Multilayer phosphorene and van der Waals heterostructures
        - Vacancy or substitutional defects
        - Strained configurations and edge reconstructions
        - Functionalization or doping studies

    Parameters
    ----------
    lattice_a : float, optional
        Lattice constant along the armchair direction (Å). Default is 3.313 Å.
    lattice_b : float, optional
        Lattice constant along the zigzag direction (Å). Default is 4.376 Å.
    width : float, optional
        Width of the generated flake (Å). Default is 10.0 Å.
    length : float, optional
        Length of the generated flake (Å). Default is 10.0 Å.
    edge_type : {'zigzag', 'armchair'}, optional
        Edge orientation of the flake. Default is 'zigzag'.
    rotation : float, optional
        Rotation angle (in degrees) about the z-axis. Default is 0.0°.
    vacuum : float, optional
        Additional spacing in the z-direction for isolated layer simulations (Å). Default is 15.0 Å.
    puckering_height : float, optional
        Vertical displacement between upper and lower sublayers in the puckered structure (Å). Default is 2.14 Å.
    origin : tuple of float, optional
        (x, y) coordinate shift applied after construction. Default is (0.0, 0.0).

    Attributes
    ----------
    atoms : list of GAM_Atom
        List of phosphorus atoms in the structure.
    bonds : list of GAM_Bond
        List of covalent P–P bonds computed via nearest-neighbor search.
    lattice_vectors : ndarray of shape (2, 3)
        Lattice vectors defining the orthorhombic cell geometry.
    meta : dict
        Metadata dictionary containing structural parameters and statistics.
    origin : tuple of float
        Current origin of the structure in Cartesian coordinates.

    Methods
    -------
    build() -> None
        Construct the phosphorene lattice from the specified parameters.
    nearest_neighbors(r_cut: float | None = None) -> None
        Identify and store P–P bonds using a cutoff-based distance search.
    translate(dx: float, dy: float, dz: float = 0.0) -> None
        Translate all atoms by the given displacement vector.
    rotate(angle_deg: float, about_center: bool = True) -> None
        Rotate the structure around the z-axis.
    copy() -> list[GAM_Atom]
        Return a deep copy of the current atomic structure.
    get_atoms() -> list[GAM_Atom]
        Retrieve a copy of all atoms in the structure.
    get_positions() -> list[tuple[float, float, float]]
        Get all atomic positions as a list of (x, y, z) coordinates.
    get_elements() -> list[str]
        Get a list of all element symbols (typically all 'P').
    to_xyz(path: str) -> None
        Export the phosphorene structure to an `.xyz` file.
    to_ase() -> ase.Atoms
        Convert the current structure into an ASE `Atoms` object for use with
        simulation packages.

    Notes
    -----
    - Phosphorene exhibits a puckered honeycomb lattice with four atoms per
      unit cell and anisotropic in-plane lattice constants.
    - The armchair direction corresponds to the short lattice vector (a-axis),
      and the zigzag direction to the longer one (b-axis).
    - The default cutoff for nearest neighbors corresponds to approximately
      1.6 × the typical P–P bond length (~2.2 Å).
    - For periodic systems, ensure that the `vacuum` spacing is sufficiently
      large to avoid interlayer interactions.

    Examples
    --------
    >>> p = Phosphorene(width=15.0, length=12.0, edge_type='armchair', rotation=10)
    >>> p.build()
    >>> p.to_xyz("phosphorene_flake.xyz")
    >>> ase_atoms = p.to_ase()
    >>> print(f"Phosphorene flake with {len(ase_atoms)} atoms constructed successfully.")
    """
    
    def __init__(self, lattice_a: float = 3.313, lattice_b: float = 4.376,
                 width: float = 10.0, length: float = 10.0, edge_type: str = 'zigzag',
                 rotation: float = 0.0, vacuum: float = 15.0, puckering_height: float = 2.14,
                 origin: Tuple[float, float] = (0.0, 0.0)):
        """
        Initialize phosphorene structure parameters.
        
        Args:
            lattice_a: Lattice constant along armchair direction (Å)
            lattice_b: Lattice constant along zigzag direction (Å)
            width: Structure width (Å)
            length: Structure length (Å)
            edge_type: 'zigzag' or 'armchair' edge orientation
            rotation: Rotation angle about z-axis (degrees)
            vacuum: Vacuum spacing in z-direction (Å)
            puckering_height: Out-of-plane puckering distance (Å)
            origin: (x0, y0) coordinate shift after construction
            periodic: If True, use periodic boundary conditions
        """
        self.lattice_a = lattice_a
        self.lattice_b = lattice_b
        self.width = width
        self.length = length
        self.edge_type = edge_type.lower()
        self.rotation = rotation
        self.vacuum = vacuum
        self.puckering_height = puckering_height
        self.origin = origin
        
        # Structure data
        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []
        self.lattice_vectors: Optional[np.ndarray] = None
        self.meta: Dict[str, Any] = {
            'lattice_a':lattice_a,
            'lattice_b':lattice_b,
            'width':width,
            'length':length,
            'edge_type':edge_type,
            'rotation': rotation,
            'vacuum':vacuum,
            'puckering_height':puckering_height,
            'origin':origin}
        
        # Build structure
        self.build()
    
    def build(self):
        """Construct the phosphorene structure."""
        # Get unit cell basis and lattice vectors
        basis, lattice_vecs = phosphorene_unit_cell(
            self.lattice_a, self.lattice_b, self.puckering_height
        )
        
        #if self.periodic:
        #    self.lattice_vectors = lattice_vecs
        
        # Calculate number of unit cells needed
        n_a = int(np.ceil(self.width / self.lattice_a)) + 1
        n_b = int(np.ceil(self.length / self.lattice_b)) + 1
        
        # Generate atoms by tiling unit cell
        atom_id = 0
        temp_atoms = []
        
        for i in range(n_a):
            for j in range(n_b):
                for x_rel, y_rel, z_rel, sublayer in basis:
                    x = i * self.lattice_a + x_rel
                    y = j * self.lattice_b + y_rel
                    z = z_rel
                    
                    # Apply edge type orientation
                    if self.edge_type == 'armchair':
                        x, y = y, x  # Swap coordinates for armchair edges
                    
                    temp_atoms.append((x, y, z, atom_id, sublayer))
                    atom_id += 1
        
        # Filter atoms within specified dimensions
        filtered_atoms = []
        for x, y, z, aid, sublayer in temp_atoms:
            if 0 <= x <= self.width and 0 <= y <= self.length:
                filtered_atoms.append((x, y, z, aid, sublayer))
        
        # Apply rotation if specified
        if self.rotation != 0.0:
            angle_rad = np.radians(self.rotation)
            cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
            
            rotated_atoms = []
            for x, y, z, aid, sublayer in filtered_atoms:
                x_rot = x * cos_a - y * sin_a
                y_rot = x * sin_a + y * cos_a
                rotated_atoms.append((x_rot, y_rot, z, aid, sublayer))
            filtered_atoms = rotated_atoms
        
        # Apply origin shift
        ox, oy = self.origin
        final_atoms = [(x + ox, y + oy, z, aid, sublayer) 
                      for x, y, z, aid, sublayer in filtered_atoms]
        
        # Create Atom objects
        self.atoms = []
        for x, y, z, aid, sublayer in final_atoms:
            atom = GAM_Atom(aid, "P", x, y, z, {"sublayer": sublayer})
            self.atoms.append(atom)
        
        # Generate bonds
        self.nearest_neighbors()
        
        # Store metadata
        self.meta = {
            "lattice_a": self.lattice_a,
            "lattice_b": self.lattice_b,
            "width": self.width,
            "length": self.length,
            "edge_type": self.edge_type,
            "rotation": self.rotation,
            "vacuum": self.vacuum,
            "puckering_height": self.puckering_height,
            "origin": self.origin,
            "n_atoms": len(self.atoms),
            "n_bonds": len(self.bonds)
        }
    
    def nearest_neighbors(self, r_cut: Optional[float] = None):
        """Determine bonds based on distance cutoff."""
        if r_cut is None:
            r_cut = 2.2 * 1.6  # Typical P-P bond length × safety factor
        
        self.bonds = []
        bond_id = 0
        
        for i, atom1 in enumerate(self.atoms):
            for j, atom2 in enumerate(self.atoms[i+1:], i+1):
                pos1 = atom1.position
                pos2 = atom2.position
                dist = np.linalg.norm(pos2 - pos1)
                
                if dist <= r_cut:
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
            
            

    

    def to_ase(self):
        """Convert to ASE Atoms object (if ASE is available)."""
        if not ASE_AVAILABLE:
            raise ImportError("ASE not available. Install with: pip install ase")
        
        symbols = [atom.element for atom in self.atoms]
        positions = [[atom.x, atom.y, atom.z] for atom in self.atoms]
        
        atoms = Atoms(symbols=symbols, positions=positions)
        
        #if self.periodic and self.lattice_vectors is not None:
        #    # Set up periodic boundary conditions
        #    cell = np.zeros((3, 3))
        #    cell[:2, :2] = self.lattice_vectors
        #    cell[2, 2] = self.vacuum
        #    atoms.set_cell(cell)
        #    atoms.set_pbc([True, True, False])
        
        return atoms















