from ..Primatom import GAM_Atom , GAM_Bond 
#from gam_atom import GAM_Atom , GAM_Bond
#from gam_atom import REFERENCE_PERIODIC_TABLE_CONFIG
import copy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
from typing import Dict, List, Tuple, Optional, Any, Union
import json
from dataclasses import dataclass, asdict
from ase import Atoms
#from ..Primatom.gam_atom import REFERENCE_PERIODIC_TABLE_CONFIG

REFERENCE_PERIODIC_TABLE_CONFIG = {
    'H':  {'number': 1,  'mass': 1.008,   'config': '1s1',                 'radius': 0.53, 'lattice_constant': None},
    'He': {'number': 2,  'mass': 4.003,   'config': '1s2',                 'radius': 0.31, 'lattice_constant': None},
    'Li': {'number': 3,  'mass': 6.941,   'config': '[He] 2s1',            'radius': 1.67, 'lattice_constant': 3.49},  # bcc
    'Be': {'number': 4,  'mass': 9.012,   'config': '[He] 2s2',            'radius': 1.12, 'lattice_constant': 2.29},  # hcp
    'B':  {'number': 5,  'mass': 10.811,  'config': '[He] 2s2 2p1',        'radius': 0.85, 'lattice_constant': 8.73},
    'C':  {'number': 6,  'mass': 12.011,  'config': '[He] 2s2 2p2',        'radius': 0.67, 'lattice_constant': 3.57},  # diamond
    'N':  {'number': 7,  'mass': 14.007,  'config': '[He] 2s2 2p3',        'radius': 0.56, 'lattice_constant': None},
    'O':  {'number': 8,  'mass': 15.999,  'config': '[He] 2s2 2p4',        'radius': 0.48, 'lattice_constant': None},
    'F':  {'number': 9,  'mass': 18.998,  'config': '[He] 2s2 2p5',        'radius': 0.42, 'lattice_constant': None},
    'Ne': {'number': 10, 'mass': 20.180,  'config': '[He] 2s2 2p6',        'radius': 0.38, 'lattice_constant': None},
    'Na': {'number': 11, 'mass': 22.990,  'config': '[Ne] 3s1',            'radius': 1.90, 'lattice_constant': 4.23},  # bcc
    'Mg': {'number': 12, 'mass': 24.305,  'config': '[Ne] 3s2',            'radius': 1.45, 'lattice_constant': 3.21},  # hcp
    'Al': {'number': 13, 'mass': 26.982,  'config': '[Ne] 3s2 3p1',        'radius': 1.18, 'lattice_constant': 4.05},  # fcc
    'Si': {'number': 14, 'mass': 28.085,  'config': '[Ne] 3s2 3p2',        'radius': 1.11, 'lattice_constant': 5.43},  # diamond
    'P':  {'number': 15, 'mass': 30.974,  'config': '[Ne] 3s2 3p3',        'radius': 1.06, 'lattice_constant': 7.17},
    'S':  {'number': 16, 'mass': 32.06,   'config': '[Ne] 3s2 3p4',        'radius': 1.02, 'lattice_constant': None},
    'Cl': {'number': 17, 'mass': 35.453,  'config': '[Ne] 3s2 3p5',        'radius': 0.99, 'lattice_constant': None},
    'Ar': {'number': 18, 'mass': 39.948,  'config': '[Ne] 3s2 3p6',        'radius': 0.95, 'lattice_constant': None},
    'K':  {'number': 19, 'mass': 39.098,  'config': '[Ar] 4s1',            'radius': 2.43, 'lattice_constant': 5.23},  # bcc
    'Ca': {'number': 20, 'mass': 40.078,  'config': '[Ar] 4s2',            'radius': 1.94, 'lattice_constant': 5.58},  # fcc
    'Fe': {'number': 26, 'mass': 55.845,  'config': '[Ar] 3d6 4s2',        'radius': 1.56, 'lattice_constant': 2.87},  # bcc
    'Cu': {'number': 29, 'mass': 63.546,  'config': '[Ar] 3d10 4s1',       'radius': 1.45, 'lattice_constant': 3.61},  # fcc
    'Zn': {'number': 30, 'mass': 65.38,   'config': '[Ar] 3d10 4s2',       'radius': 1.42, 'lattice_constant': 2.66},  # hcp
    'Ag': {'number': 47, 'mass': 107.868, 'config': '[Kr] 4d10 5s1',       'radius': 1.65, 'lattice_constant': 4.09},  # fcc
    'Cd': {'number': 48, 'mass': 112.411, 'config': '[Kr] 4d10 5s2',       'radius': 1.55, 'lattice_constant': 2.98},  # hcp
    'W':  {'number': 74, 'mass': 183.84,  'config': '[Xe] 4f14 5d4 6s2',   'radius': 1.35, 'lattice_constant': 3.16},  # bcc
    'Pt': {'number': 78, 'mass': 195.084, 'config': '[Xe] 4f14 5d9 6s1',   'radius': 1.39, 'lattice_constant': 3.92},  # fcc
    'Au': {'number': 79, 'mass': 196.967, 'config': '[Xe] 4f14 5d10 6s1',  'radius': 1.44, 'lattice_constant': 4.08},  # fcc
    'Ta': {'number': 73, 'mass': 180.948, 'config': '[Xe] 4f14 5d3 6s2',   'radius': 1.43, 'lattice_constant': 3.30},  # bcc
    'Nb': {'number': 41, 'mass': 92.906,  'config': '[Kr] 4d4 5s1',        'radius': 1.46, 'lattice_constant': 3.30},  # bcc
    'Mo': {'number': 42, 'mass': 95.94,   'config': '[Kr] 4d5 5s1',        'radius': 1.39, 'lattice_constant': 3.15},  # bcc
}


@dataclass
class AtomicData:
    """Store atomic data for different elements."""
    
    # Atomic radii in Angstroms (covalent radii)
    ATOMIC_RADII={element : properties['radius'] for element,properties in REFERENCE_PERIODIC_TABLE_CONFIG.items()}

    
    # Lattice parameters in Angstroms
    LATTICE_CONSTANTS = {element : properties['lattice_constant'] for element,properties in REFERENCE_PERIODIC_TABLE_CONFIG.items()}



class Nanoparticle_Generator:
    """
    General-purpose nanoparticle structure generator for atomistic simulations.

    The `Nanoparticle_Generator` class provides a flexible and modular tool for
    constructing atomic-scale models of nanoparticles with arbitrary geometries,
    materials, and orientations. It supports the creation of finite clusters,
    spherical or polyhedral nanoparticles, and core–shell architectures with
    user-defined crystallography and dimensions.

    This class can be used to prepare atomic coordinates for density functional
    theory (DFT), molecular dynamics (MD), or machine learning (ML) workflows.
    It is particularly suited for generating input structures for nanomaterials
    research involving metals, semiconductors, and hybrid systems.

    Supported features include:
        - Generation of finite nanoparticles with specified radius, shape, and lattice
        - Core–shell or multi-layered structures
        - Selection of crystalline structure (FCC, BCC, HCP, etc.)
        - Randomized orientation and rotation of particles
        - Atom-based filtering for specific elements or sublattices
        - Export to standard formats such as XYZ or ASE-compatible structures

    The implementation can be extended to support:
        - Nanorods, nanowires, and anisotropic morphologies
        - Surface relaxation or reconstruction
        - Doping, vacancies, and alloy generation
        - Integration with LAMMPS, VASP, or ASE simulation interfaces

    Parameters
    ----------
    element : str
        Chemical symbol of the atomic species (e.g., 'Au', 'Ag', 'Si').
    lattice_constant : float
        Lattice constant of the crystal in Ångströms.
    structure_type : {'FCC', 'BCC', 'HCP', 'diamond', ...}
        Crystal structure used for atomic packing.
    radius : float
        Radius of the nanoparticle in Ångströms.
    center : tuple of float, optional
        Cartesian coordinates of the nanoparticle center. Default is (0.0, 0.0, 0.0).
    core_element : str, optional
        If specified, the core region of the nanoparticle will use this element.
    core_radius : float, optional
        Radius of the core region (Å) when building core–shell particles.
    random_orientation : bool, optional
        If True, applies a random rotation to the generated particle. Default is False.
    vacuum : float, optional
        Optional vacuum spacing in Ångströms for periodic export. Default is 10.0 Å.
    save_path : str, optional
        Path for saving generated structures in `.xyz` or `.cif` format.
    seed : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    atoms : list of GAM_Atom
        List of atoms with positions, element type, and optional metadata.
    bonds : list of GAM_Bond
        List of identified bonds (if neighbor detection is used).
    lattice_vectors : numpy.ndarray
        Lattice vectors defining the crystal unit cell.
    meta : dict
        Metadata dictionary containing generation parameters and history.
    center : tuple of float
        Center of the nanoparticle in Cartesian coordinates.
    radius : float
        Radius of the nanoparticle in Ångströms.

    Methods
    -------
    build() -> None
        Construct the nanoparticle structure from the given parameters.
    generate_core_shell(core_element: str, core_radius: float) -> None
        Create a core–shell nanoparticle by replacing atoms within the core radius.
    apply_rotation(angle: float, axis: str = 'z') -> None
        Rotate the nanoparticle by a given angle around a specified axis.
    translate(dx: float, dy: float, dz: float) -> None
        Translate all atoms by the specified vector.
    nearest_neighbors(r_cut: float | None = None) -> None
        Identify bonds based on a distance cutoff.
    filter_atoms(element: str) -> list[GAM_Atom]
        Return atoms belonging to a specific element.
    get_positions() -> list[tuple[float, float, float]]
        Return atomic positions as (x, y, z) tuples.
    get_elements() -> list[str]
        Return a list of all atomic species in the nanoparticle.
    to_xyz(path: str) -> None
        Export the nanoparticle to an XYZ file.
    to_ase() -> ase.Atoms
        Convert the nanoparticle to an ASE-compatible `Atoms` object.

    Notes
    -----
    - The atomic cutoff radius for neighbor detection should be set according
      to the typical bond length of the chosen element (e.g., 1.2× nearest-neighbor distance).
    - The nanoparticle surface morphology depends on both lattice type and radius.
    - Randomized orientation can be used to eliminate orientation bias in ML datasets.
    - Large nanoparticles (>10,000 atoms) may require optimized neighbor-search algorithms.
    - The class is designed to be modular, allowing integration with simulation workflows
      such as VASP, Quantum ESPRESSO, LAMMPS, or CP2K.

    Examples
    --------
    >>> # Generate a 5 nm gold nanoparticle
    >>> from pygamlab import Nanoparticle_Generator
    >>> npg = Nanoparticle_Generator(element='Au', lattice_constant=4.08,
    ...                              structure_type='FCC', radius=25.0)
    >>> npg.build()
    >>> npg.to_xyz('Au_nanoparticle.xyz')

    >>> # Generate a core–shell Ag@Au nanoparticle
    >>> cs = Nanoparticle_Generator(element='Au', lattice_constant=4.08,
    ...                             structure_type='FCC', radius=30.0)
    >>> cs.generate_core_shell(core_element='Ag', core_radius=15.0)
    >>> cs.to_xyz('AgAu_core_shell.xyz')
    """
    
    def __init__(self, element: str, size_nm: float, 
                 shape: str = 'sphere', doping: Optional[Dict[str, float]] = None,
                 coating: Optional[Tuple[str, float]] = None, 
                 crystal_structure: str = 'FCC', **kwargs):
        """
        Initialize a nanoparticle with specified parameters.
        
        Args:
            element: Primary element for the nanoparticle (e.g., 'Au')
            size_nm: Particle size in nanometers
            shape: Shape type ('sphere', 'rod', 'cube', 'octahedron')
            doping: Optional dictionary of doping elements with percentages
            coating: Optional tuple of (material, thickness_nm)
            crystal_structure: Crystal structure ('FCC', 'BCC', 'hexagonal', 'diamond')
            **kwargs: Additional parameters (surface_roughness, lattice_constant, etc.)
        """
        # Validate size_nm to ensure it's a positive value
        if size_nm <= 0:
            raise ValueError(f"Particle size must be a positive value. Provided size: {size_nm} nm.")
        
        
        # Validate doping percentages
        if doping:
            for element, percentage in doping.items():
                if percentage > 100:
                    raise ValueError(f"Doping percentage for {element} exceeds 100%. Provided: {percentage}%.")
        
        
        
        # Store basic parameters
        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []
        self.element = element
        self.elements = [element]  # Keep for backward compatibility
        self.size_nm = size_nm
        self.shape = shape.lower()
        self.doping = doping or {}
        self.coating = coating
        self.crystal_structure = crystal_structure.upper()
        
        # Store additional parameters
        self.surface_roughness = kwargs.get('surface_roughness', 0.0)
        self.lattice_constant = kwargs.get('lattice_constant', None)
        self.custom_parameters = kwargs
        
        # Initialize atomic data
        self.atomic_data = AtomicData()
        
        # Store all parameters in a dictionary
        self.meta = {
            'element': self.element,
            'elements': self.elements,
            'size_nm': self.size_nm,
            'shape': self.shape,
            'doping': self.doping,
            'coating': self.coating,
            'crystal_structure': self.crystal_structure,
            'surface_roughness': self.surface_roughness,
            'lattice_constant': self.lattice_constant,
            **kwargs
        }
        
        # Generate the nanoparticle structure
        self.positions = None
        self.atom_types = None
        self._build()
        
        
    
    def _get_lattice_constant(self) -> float:
        """
        Get the lattice constant for the primary element.
        
        Returns:
            float: Lattice constant in Angstroms
        """
        if self.lattice_constant is not None:
            return self.lattice_constant
        
        if self.element in self.atomic_data.LATTICE_CONSTANTS:
            return self.atomic_data.LATTICE_CONSTANTS[self.element]
        else:
            # Default lattice constant
            return 4.0
    
    def _generate_unit_cell_positions(self) -> Tuple[np.ndarray, List[str]]:
        """
        Generate unit cell positions based on crystal structure.
        
        Returns:
            Tuple of (positions array, atom types list)
        """
        a = self._get_lattice_constant()
        
        if self.crystal_structure == 'FCC':
            # Face-centered cubic
            positions = np.array([
                [0, 0, 0],
                [0.5, 0.5, 0],
                [0.5, 0, 0.5],
                [0, 0.5, 0.5]
            ]) * a
            atom_types = [self.element] * 4
            
        elif self.crystal_structure == 'BCC':
            # Body-centered cubic
            positions = np.array([
                [0, 0, 0],
                [0.5, 0.5, 0.5]
            ]) * a
            atom_types = [self.element] * 2
            
        elif self.crystal_structure == 'HEXAGONAL':
            # Hexagonal close-packed
            c_a_ratio = self.custom_parameters.get('c_a_ratio', 1.633)
            c = a * c_a_ratio
            positions = np.array([
                [0, 0, 0],
                [2/3, 1/3, 0.5]
            ]) * np.array([a, a, c])
            atom_types = [self.element] * 2
            
        elif self.crystal_structure == 'DIAMOND':
            # Diamond cubic
            positions = np.array([
                [0, 0, 0],
                [0.25, 0.25, 0.25],
                [0.5, 0.5, 0],
                [0.75, 0.75, 0.25],
                [0.5, 0, 0.5],
                [0.75, 0.25, 0.75],
                [0, 0.5, 0.5],
                [0.25, 0.75, 0.75]
            ]) * a
            atom_types = [self.element] * 8
            
        else:
            raise ValueError(f"Unsupported crystal structure: {self.crystal_structure}")
        
        return positions, atom_types
    
    def _generate_lattice_positions(self) -> Tuple[np.ndarray, List[str]]:
        """
        Generate lattice positions by replicating unit cells.
        
        Returns:
            Tuple of (positions array, atom types list)
        """
        unit_positions, unit_atom_types = self._generate_unit_cell_positions()
        a = self._get_lattice_constant()
        
        # Determine number of unit cells needed
        size_angstrom = self.size_nm * 10  # Convert nm to Angstrom
        n_cells = int(np.ceil(size_angstrom / a)) + 1
        
        all_positions = []
        all_atom_types = []
        
        # Replicate unit cell
        for i in range(-n_cells, n_cells + 1):
            for j in range(-n_cells, n_cells + 1):
                for k in range(-n_cells, n_cells + 1):
                    offset = np.array([i, j, k]) * a
                    shifted_positions = unit_positions + offset
                    
                    all_positions.append(shifted_positions)
                    all_atom_types.extend(unit_atom_types)
        
        positions = np.vstack(all_positions)
        
        return positions, all_atom_types
    
    def _apply_shape_constraint(self, positions: np.ndarray) -> np.ndarray:
        """
        Apply shape constraints to filter atomic positions.
        
        Args:
            positions: Array of atomic positions
            
        Returns:
            Boolean mask for positions within the shape
        """
        size_angstrom = self.size_nm * 10 / 2  # Half-size in Angstrom
        
        if self.shape == 'sphere':
            distances = np.linalg.norm(positions, axis=1)
            mask = distances <= size_angstrom
            
        elif self.shape == 'cube':
            mask = np.all(np.abs(positions) <= size_angstrom, axis=1)
            
        elif self.shape == 'rod':
            # Rod along z-axis
            aspect_ratio = self.custom_parameters.get('aspect_ratio', 3.0)
            rod_length = size_angstrom * aspect_ratio
            rod_radius = size_angstrom / aspect_ratio
            
            xy_distances = np.linalg.norm(positions[:, :2], axis=1)
            z_mask = np.abs(positions[:, 2]) <= rod_length
            xy_mask = xy_distances <= rod_radius
            mask = z_mask & xy_mask
            
        elif self.shape == 'octahedron':
            # Regular octahedron
            mask = (np.abs(positions[:, 0]) + np.abs(positions[:, 1]) + 
                   np.abs(positions[:, 2])) <= size_angstrom
            
        else:
            raise ValueError(f"Unsupported shape: {self.shape}")
        
        return mask
    
    def _apply_doping(self, atom_types: List[str]) -> List[str]:
        """
        Apply doping by randomly replacing atoms.
        
        Args:
            atom_types: Original atom types
            
        Returns:
            Modified atom types with doping
        """
        if not self.doping:
            return atom_types
        
        doped_types = atom_types.copy()
        n_atoms = len(atom_types)
        
        for dopant, percentage in self.doping.items():
            n_dopant = int(n_atoms * percentage / 100)
            indices = random.sample(range(n_atoms), n_dopant)
            
            for idx in indices:
                doped_types[idx] = dopant
        
        return doped_types
    
    def _add_coating(self, positions: np.ndarray, atom_types: List[str]) -> Tuple[np.ndarray, List[str]]:
        """
        Add coating layer to the nanoparticle.
        
        Args:
            positions: Core positions
            atom_types: Core atom types
            
        Returns:
            Tuple of (all positions, all atom types) including coating
        """
        if not self.coating:
            return positions, atom_types
        
        coating_material, coating_thickness_nm = self.coating
        coating_thickness_angstrom = coating_thickness_nm * 10
        
        # Find surface atoms (simplified approach)
        center = np.mean(positions, axis=0)
        distances_from_center = np.linalg.norm(positions - center, axis=1)
        surface_threshold = np.percentile(distances_from_center, 90)
        surface_mask = distances_from_center >= surface_threshold
        surface_positions = positions[surface_mask]
        
        # Generate coating positions
        coating_positions = []
        coating_types = []
        
        for pos in surface_positions:
            direction = (pos - center) / np.linalg.norm(pos - center)
            # Add multiple layers
            n_layers = max(1, int(coating_thickness_angstrom / 2.0))
            for layer in range(1, n_layers + 1):
                coat_pos = pos + direction * layer * 2.0
                coating_positions.append(coat_pos)
                coating_types.append(coating_material)
        
        if coating_positions:
            all_positions = np.vstack([positions, np.array(coating_positions)])
            all_types = atom_types + coating_types
        else:
            all_positions = positions
            all_types = atom_types
        
        return all_positions, all_types
    
    def _apply_surface_roughness(self, positions: np.ndarray) -> np.ndarray:
        """
        Apply surface roughness by adding random displacement.
        
        Args:
            positions: Original positions
            
        Returns:
            Positions with surface roughness applied
        """
        if self.surface_roughness <= 0:
            return positions
        
        # Find surface atoms
        center = np.mean(positions, axis=0)
        distances_from_center = np.linalg.norm(positions - center, axis=1)
        surface_threshold = np.percentile(distances_from_center, 80)
        surface_mask = distances_from_center >= surface_threshold
        
        # Apply random displacement to surface atoms
        roughness_positions = positions.copy()
        surface_displacement = np.random.normal(0, self.surface_roughness, 
                                              (np.sum(surface_mask), 3))
        roughness_positions[surface_mask] += surface_displacement
        
        return roughness_positions
    
    def _build(self):
        """Generate the complete nanoparticle structure."""
        # Generate initial lattice
        positions, atom_types = self._generate_lattice_positions()
        
        # Apply shape constraint
        shape_mask = self._apply_shape_constraint(positions)
        positions = positions[shape_mask]
        atom_types = [atom_types[i] for i, mask in enumerate(shape_mask) if mask]
        
        # Apply doping
        atom_types = self._apply_doping(atom_types)
        
        # Add coating
        positions, atom_types = self._add_coating(positions, atom_types)
        
        # Apply surface roughness
        positions = self._apply_surface_roughness(positions)
        
        # Center the nanoparticle
        center = np.mean(positions, axis=0)
        positions -= center
        
        self.positions = positions
        self.atom_types = atom_types
        
        for i in range(len(self.atom_types)):
            x, y, z = self.positions[i]  # unpack coordinates
            each_atom = GAM_Atom(id=i, element=self.atom_types[i], x=x, y=y, z=z)
            self.atoms.append(each_atom)
            
            
        
        
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
    
    def add_element(self, element: str, percentage: float, distribution: str = 'random') -> None:
        """
        Add an element to create an alloy nanoparticle.
        
        Args:
            element: Element symbol to add (e.g., 'Ag', 'Cu')
            percentage: Percentage of the new element (0-100)
            distribution: Distribution method ('random', 'surface', 'core', 'layered')
        """
        if not 0 <= percentage <= 100:
            raise ValueError(f"Percentage must be between 0 and 100. Provided: {percentage}%")
        
        if element not in self.atomic_data.ATOMIC_RADII:
            raise ValueError(f"Element {element} not supported. Available elements: {list(self.atomic_data.ATOMIC_RADII.keys())}")
        
        # Add to elements list
        if element not in self.elements:
            self.elements.append(element)
        
        # Calculate number of atoms to replace
        n_atoms = len(self.atom_types)
        n_to_replace = int(n_atoms * percentage / 100)
        
        if distribution == 'random':
            # Random distribution throughout the nanoparticle
            indices = random.sample(range(n_atoms), n_to_replace)
            for idx in indices:
                self.atom_types[idx] = element
                
        elif distribution == 'surface':
            # Replace atoms on the surface
            center = np.mean(self.positions, axis=0)
            distances_from_center = np.linalg.norm(self.positions - center, axis=1)
            surface_threshold = np.percentile(distances_from_center, 80)
            surface_indices = np.where(distances_from_center >= surface_threshold)[0]
            
            if len(surface_indices) >= n_to_replace:
                selected_indices = random.sample(list(surface_indices), n_to_replace)
            else:
                selected_indices = surface_indices
                # Fill remaining with random selection
                remaining = n_to_replace - len(surface_indices)
                other_indices = [i for i in range(n_atoms) if i not in surface_indices]
                if other_indices:
                    selected_indices.extend(random.sample(other_indices, min(remaining, len(other_indices))))
            
            for idx in selected_indices:
                self.atom_types[idx] = element
                
        elif distribution == 'core':
            # Replace atoms in the core
            center = np.mean(self.positions, axis=0)
            distances_from_center = np.linalg.norm(self.positions - center, axis=1)
            core_threshold = np.percentile(distances_from_center, 20)
            core_indices = np.where(distances_from_center <= core_threshold)[0]
            
            if len(core_indices) >= n_to_replace:
                selected_indices = random.sample(list(core_indices), n_to_replace)
            else:
                selected_indices = core_indices
                # Fill remaining with random selection
                remaining = n_to_replace - len(core_indices)
                other_indices = [i for i in range(n_atoms) if i not in core_indices]
                if other_indices:
                    selected_indices.extend(random.sample(other_indices, min(remaining, len(other_indices))))
            
            for idx in selected_indices:
                self.atom_types[idx] = element
                
        elif distribution == 'layered':
            # Create layered structure (core-shell or alternating layers)
            center = np.mean(self.positions, axis=0)
            distances_from_center = np.linalg.norm(self.positions - center, axis=1)
            
            # Sort atoms by distance from center
            sorted_indices = np.argsort(distances_from_center)
            
            # Replace atoms in specific layers
            layer_size = n_atoms // 4  # Divide into 4 layers
            for i in range(min(n_to_replace, n_atoms)):
                layer_idx = (i // layer_size) % 4
                if layer_idx % 2 == 1:  # Replace in odd-numbered layers
                    self.atom_types[sorted_indices[i]] = element
                    
        else:
            raise ValueError(f"Unsupported distribution method: {distribution}. Use 'random', 'surface', 'core', or 'layered'")
        
        # Update atoms list
        for i, atom_type in enumerate(self.atom_types):
            if i < len(self.atoms):
                self.atoms[i].element = atom_type
            else:
                # Create new atom if needed
                x, y, z = self.positions[i]
                new_atom = GAM_Atom(id=i, element=atom_type, x=x, y=y, z=z)
                self.atoms.append(new_atom)
        
        # Update metadata
        self.meta['elements'] = self.elements
        print(f"Added {element} ({percentage}%) with {distribution} distribution")
    
    def get_composition(self) -> Dict[str, float]:
        """
        Get the current composition of the nanoparticle as percentages.
        
        Returns:
            Dictionary mapping element symbols to their percentages
        """
        if not self.atom_types:
            return {}
        
        unique_elements, counts = np.unique(self.atom_types, return_counts=True)
        total_atoms = len(self.atom_types)
        
        composition = {}
        for element, count in zip(unique_elements, counts):
            composition[element] = (count / total_atoms) * 100
        
        return composition
    
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
        """
        Export structure to ASE Atoms object (if ASE is available).
        
        Returns:
            ASE Atoms object or None if ASE not available
        """

        symbols = [atom.element for atom in self.atoms]
        positions = [atom.position for atom in self.atoms]
        
        atoms_obj = Atoms(symbols=symbols, positions=positions)
        
        # Add metadata
        atoms_obj.info.update(self.meta)
        
        return atoms_obj
    
    
    
    def _to_xyz(self, file_path: str) -> None:
        """
        Save nanoparticle atomic positions to an XYZ file.
        
        Args:
            file_path: Path to save the XYZ file
        """
        if self.positions is None:
            raise ValueError("Nanoparticle structure not generated")
        
        with open(file_path, 'w') as f:
            # Write number of atoms
            f.write(f"{len(self.positions)}\n")
            
            # Write comment line with nanoparticle info
            comment = (f"Nanoparticle: {', '.join(self.elements)}, "
                      f"Size: {self.size_nm} nm, Shape: {self.shape}, "
                      f"Crystal: {self.crystal_structure}")
            f.write(f"{comment}\n")
            
            # Write atomic coordinates
            for i, (pos, atom_type) in enumerate(zip(self.positions, self.atom_types)):
                f.write(f"{atom_type:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        
        print(f"Nanoparticle structure saved to {file_path}")
    
    def _to_ase(self):
        """
        Return an ASE Atoms object (requires ASE installation).
        
        Returns:
            ASE Atoms object containing the nanoparticle structure
        """
        if self.positions is None:
            raise ValueError("Nanoparticle structure not generated")
        
        try:
            from ase import Atoms
            atoms = Atoms(symbols=self.atom_types, positions=self.positions)
            return atoms
        except ImportError:
            print("ASE not available. Install with: pip install ase")
            # Return a simple dictionary instead
            return {
                'symbols': self.atom_types,
                'positions': self.positions,
                'cell': None,
                'pbc': [False, False, False]
            }
    
    
    
    
    def get_info(self) -> Dict[str, Any]:
        """
        Get comprehensive information about the nanoparticle.
        
        Returns:
            Dictionary containing nanoparticle information
        """
        if self.positions is None:
            return self.meta
        
        # Calculate additional properties
        unique_elements, counts = np.unique(self.atom_types, return_counts=True)
        composition = dict(zip(unique_elements, counts))
        
        center_of_mass = np.mean(self.positions, axis=0)
        max_distance = np.max(np.linalg.norm(self.positions - center_of_mass, axis=1))
        actual_size_nm = max_distance * 2 / 10  # Convert to nm
        
        info = self.meta.copy()
        info.update({
            'n_atoms': len(self.positions),
            'composition': composition,
            'actual_size_nm': actual_size_nm,
            'center_of_mass': center_of_mass.tolist(),
            'bounding_box': {
                'x_range': [float(self.positions[:, 0].min()), float(self.positions[:, 0].max())],
                'y_range': [float(self.positions[:, 1].min()), float(self.positions[:, 1].max())],
                'z_range': [float(self.positions[:, 2].min()), float(self.positions[:, 2].max())]
            }
        })
        
        return info
    

    def save_parameters(self, file_path: str) -> None:
        """
        Save nanoparticle parameters to a JSON file.
        
        Args:
            file_path: Path to save the parameters file
        """
        info = self.get_info()
        with open(file_path, 'w') as f:
            json.dump(info, f, indent=2, default=str)
        print(f"Parameters saved to {file_path}")





    


