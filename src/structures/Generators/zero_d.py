from ..Primatom import GAM_Atom , GAM_Bond
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.cluster.icosahedron import Icosahedron
from ase.cluster.wulff import wulff_construction
from ase.visualize import view
import matplotlib.pyplot as plt
from ase.io import write
from ase.io.pov import get_bondpairs
import random
from typing import List, Tuple, Dict, Optional, Union , Any
from ase.atom import Atom
from ase.data import atomic_numbers, covalent_radii
from ase.constraints import FixAtoms
import itertools
import copy





class Nano_ZeroD_Builder:
    """
    High-level atomistic builder for 0D nanostructures.

    The `Nano_ZeroD_Builder` class provides a unified interface for constructing
    atomistic models of zero-dimensional (0D) nanostructures, including:

        - **Quantum dots** (semiconducting nanocrystals with quantum confinement)
        - **Nanoclusters** (small metallic or covalent atomic clusters)
        - **Nanoparticles** (larger crystalline particles, optionally faceted using Wulff construction)

    This class supports both geometric and crystallographic control, allowing
    for the creation of realistic nanoscale systems based on their lattice type,
    surface facets, and atomic arrangement. The resulting structures can be used
    as input for DFT, MD, or ML simulations.

    Features
    --------
    - Build 0D nanostructures from bulk crystals or predefined motifs
    - Support for FCC, BCC, HCP, diamond, zincblende, and related lattices
    - Integration with ASE (Atomic Simulation Environment)
    - Alloy, defect, and strain engineering
    - Ligand surface functionalization
    - Core–shell nanoparticle generation
    - Randomized variation generation for dataset expansion

    Parameters
    ----------
    material : str
        Chemical symbol of the primary element (e.g., 'Au', 'CdSe').
    structure_type : {'quantum_dot', 'nanocluster', 'nanoparticle'}
        Type of nanostructure to build.
    size : float
        Size parameter — interpreted as:
            - Radius (Å) for quantum dots or nanoparticles
            - Number of atomic shells for nanoclusters
    crystal_structure : str, optional
        Underlying lattice type. Defaults to `'fcc'`.
    lattice_constant : float, optional
        Lattice constant in Ångströms. Required for bulk-based generation.
    vacuum : float, optional
        Vacuum spacing around the nanostructure (Å). Default is 10.0 Å.
    surfaces : list of tuple[int], optional
        Miller indices of the surfaces for Wulff construction (only for nanoparticles).
    surface_energies : list of float, optional
        Surface energies (J/m²) corresponding to `surfaces`.
    noshells : int, optional
        Number of shells for nanocluster structures (if not using `size`).
    radius_type : {'sphere'}, optional
        Defines the nanoparticle geometry. Currently only `'sphere'` is supported.
    **kwargs : dict
        Additional keyword arguments passed to ASE structure constructors.

    Attributes
    ----------
    ASE_atoms : ase.Atoms
        The underlying ASE `Atoms` object representing the structure.
    atoms : list[GAM_Atom]
        A list of GAM_Atom objects for GAM-Lab compatibility.
    bonds : list[GAM_Bond]
        List of identified bonds (if applicable).
    meta : dict
        Metadata dictionary containing build parameters and history.
    structure_type : str
        Type of structure ('quantum_dot', 'nanocluster', 'nanoparticle').
    radius_type : str
        Shape type used for spherical generation.

    Methods
    -------
    create_alloy(elements, compositions)
        Randomly substitutes atomic species to generate alloyed structures.
    add_defects(defect_type, concentration, elements=None)
        Introduce vacancies, interstitials, or substitutional defects.
    create_variations(base_structure, size_range, n_structures, **kwargs)
        Generate multiple structural variations within a size range.
    apply_strain(strain, direction='xyz')
        Apply uniaxial or isotropic strain to the structure.
    modify_surface(ligands, coverage=1.0)
        Add ligand atoms to surface sites for surface modification.
    create_core_shell(core_material, shell_material, shell_thickness, shell_lattice_constant=None)
        Create a core–shell nanostructure by wrapping the core with another material.
    translate(dx, dy, dz=0.0)
        Translate all atoms by a specified vector.
    rotate(angle_deg, about_center=True)
        Rotate the structure about the z-axis or its center of mass.
    copy()
        Create a deep copy of the structure.
    get_atoms()
        Return a copy of the GAM_Atom list.
    get_positions()
        Return atomic positions as a list of (x, y, z) tuples.
    get_elements()
        Return a list of atomic symbols.
    _to_GAM_Atoms(atoms)
        Internal conversion from ASE atoms to GAM_Atom objects.

    Notes
    -----
    - The generated nanostructures can serve as initial geometries for
      relaxation, energy minimization, or electronic structure calculations.
    - Surface atoms are automatically centered with optional vacuum padding.
    - The Wulff construction minimizes total surface energy for given facets.
    - Alloying, defect creation, and surface modification are stochastic; use a fixed random seed for reproducibility.
    - The core–shell builder assumes isotropic growth; anisotropic extension can be implemented by modifying the shell generator.

    Examples
    --------
    >>> from pygamlab import Nano_ZeroD_Builder
    >>> # Generate a gold nanocluster with 5 atomic shells
    >>> builder = Nano_ZeroD_Builder(material='Au', structure_type='nanocluster', size=5)
    >>> atoms = builder.ASE_atoms
    >>> builder.to_xyz('Au_cluster.xyz')

    >>> # Build a 3 nm CdSe quantum dot
    >>> qd = Nano_ZeroD_Builder(material='CdSe', structure_type='quantum_dot',
    ...                         crystal_structure='zincblende',
    ...                         lattice_constant=6.08, size=15.0)
    >>> qd.apply_strain(0.02)
    >>> qd.to_xyz('CdSe_QD_strained.xyz')

    >>> # Create a core–shell Ag@Au nanoparticle
    >>> cs = Nano_ZeroD_Builder(material='Ag', structure_type='nanoparticle',
    ...                         size=25.0, surfaces=[(1, 1, 1), (1, 0, 0)],
    ...                         surface_energies=[1.2, 1.5])
    >>> cs.create_core_shell(core_material='Ag', shell_material='Au', shell_thickness=5.0)
    >>> cs.to_xyz('AgAu_core_shell.xyz')
    """

    def __init__(self, material: str, structure_type: str,
               size: float,
               crystal_structure: str = "fcc",
               lattice_constant: float = None,
               vacuum: float = 10.0,
               surfaces=None,
               surface_energies=None,
               noshells: int = None,
               radius_type: str = "sphere",
               **kwargs):
        """
        Initialize the builder.

        :param material: Chemical symbol, e.g., "Au", "CdSe"
        :param lattice_constant: Lattice constant in Å (optional)
        :param crystal_structure: 'fcc', 'bcc', 'hcp', 'diamond', 'zincblende', etc.
        :param structure_type: 'quantum_dot', 'nanocluster', or 'nanoparticle'
        :param size: Radius (Å) for quantum_dot/nanoparticle, or number of shells for nanocluster
        :param vacuum: Vacuum spacing in Å
        :param surfaces: List of Miller indices for Wulff construction (nanoparticles)
        :param surface_energies: Corresponding surface energies (J/m²)
        :param noshells: Number of atomic shells for nanocluster (if not using size)
        :param radius_type: Shape for quantum_dot/nanoparticle: 'sphere' (currently supported)
        :param kwargs: Additional ASE builder arguments
        :return: ASE Atoms object

        """
        self.material = material
        self.lattice_constant = lattice_constant
        self.crystal_structure = crystal_structure

        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []


        self.structure_type=structure_type
        self.size=size
        self.surfaces=surfaces
        self.surface_energies=surface_energies
        self.noshells=noshells
        self.radius_type=radius_type




        self.meta={'material':material ,
        'lattice_constant':lattice_constant ,
        'crystal_structure': crystal_structure,
        'structure_type':structure_type ,
        'size': size,
        'surfaces':surfaces ,
        'surface_energies': surface_energies,
        'noshells':noshells ,
        'radius_type': radius_type}





        structure_type = structure_type.lower()

        if structure_type == "nanocluster":
            shells = noshells if noshells is not None else int(size)
            atoms = Icosahedron(self.material, noshells=shells)
            self.atoms=self._to_GAM_Atoms(atoms)

        elif structure_type == "quantum_dot":
            # Build bulk lattice first
            atoms_bulk = bulk(self.material,
                              crystalstructure=self.crystal_structure,
                              a=self.lattice_constant)
            positions = atoms_bulk.get_positions()
            center = positions.mean(axis=0)
            mask = np.linalg.norm(positions - center, axis=1) <= size
            atoms = atoms_bulk[mask]
            self.atoms=atoms

        elif structure_type == "nanoparticle":
            if surfaces is None or surface_energies is None:
                raise ValueError("For nanoparticle Wulff construction, provide surfaces and surface_energies.")
            atoms = wulff_construction(self.material,
                                       surfaces, surface_energies,
                                       size,  # number of atoms approx
                                       self.crystal_structure,
                                       rounding="above")
            
            #here changed
            #self.atoms=atoms

        else:
            raise ValueError(f"Unknown 0D structure type: {structure_type}")

        atoms.center(vacuum=vacuum)
        #here changed
        #self.atoms.center(vacuum=vacuum)
    

        self.ASE_atoms = atoms
        self.atoms=self._to_GAM_Atoms(atoms)

    def create_alloy(self,elements: List[str], 
                    compositions: List[float]) -> GAM_Atom:
        """
        Create an alloy by randomly substituting atoms.
        
        :param base_atoms: Base structure
        :param elements: List of elements to include
        :param compositions: Corresponding atomic fractions
        :return: Alloyed structure
        """
        if len(elements) != len(compositions) or abs(sum(compositions) - 1.0) > 1e-6:
            raise ValueError("Invalid composition specification")
            
        atoms = self.ASE_atoms
        n_atoms = len(self.atoms)
        indices = list(range(n_atoms))
        random.shuffle(indices)
        
        current_idx = 0
        for element, comp in zip(elements, compositions):
            n_element = int(comp * n_atoms)
            for idx in indices[current_idx:current_idx + n_element]:
                atoms[idx].symbol = element
            current_idx += n_element
        
        self.ASE_atoms=atoms
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return atoms
        #return self.atoms


    def add_defects(self, defect_type: str,
                    concentration: float, elements: List[str] = None) -> GAM_Atom:
        """
        Add defects to the structure.
        
        :param defect_type: 'vacancy', 'interstitial', or 'substitution'
        :param concentration: Defect concentration (0-1)
        :param elements: Elements for substitution/interstitial
        :return: Structure with defects
        """
        if concentration < 0 or concentration > 1:
            raise ValueError("Concentration must be between 0 and 1")
            
        if defect_type not in ['vacancy', 'interstitial', 'substitution']:
            raise ValueError("Defect type must be 'vacancy', 'interstitial', or 'substitution'")
            
        if defect_type in ['substitution', 'interstitial'] and not elements:
            raise ValueError(f"Elements must be provided for {defect_type} defects")
            
        defected_atoms = self.ASE_atoms.copy()
        n_atoms = len(self.ASE_atoms)
        n_defects = int(concentration * n_atoms)
        
        # Ensure we don't try to remove more atoms than we have
        if defect_type == 'vacancy' and n_defects >= n_atoms:
            raise ValueError(f"Cannot remove {n_defects} atoms from structure with only {n_atoms} atoms")
        
        if defect_type == 'vacancy':
            # Randomly select atoms to remove
            indices_to_remove = random.sample(range(n_atoms), n_defects)
            # Create new atoms object without the selected indices
            remaining_indices = [i for i in range(n_atoms) if i not in indices_to_remove]
            defected_atoms = self.ASE_atoms[remaining_indices]
            
        elif defect_type == 'substitution':
            # Randomly select atoms to substitute
            indices = random.sample(range(n_atoms), n_defects)
            for idx in indices:
                defected_atoms[idx].symbol = random.choice(elements)
                
        elif defect_type == 'interstitial':
            # Add interstitial atoms at random positions
            positions = self.ASE_atoms.get_positions()
            
            # Get the bounding box of existing atoms
            min_pos = np.min(positions, axis=0)
            max_pos = np.max(positions, axis=0)
            
            for _ in range(n_defects):
                # Generate random position within the bounding box
                pos = min_pos + np.random.rand(3) * (max_pos - min_pos)
                # Add some random offset to avoid overlapping with existing atoms
                pos += (np.random.rand(3) - 0.5) * 2.0  # ±1 Å random offset
                
                new_atom = Atom(random.choice(elements), pos)
                defected_atoms.append(new_atom)
               
        self.ASE_atoms = defected_atoms
        self.atoms = self._to_GAM_Atoms(self.ASE_atoms)



    #remove.....
    def create_variations(self, base_structure: str, size_range: Tuple[float, float],
                         n_structures: int, **kwargs) -> List[GAM_Atom]:
        """
        Create multiple variations of a structure.
        
        :param base_structure: Structure type
        :param size_range: (min_size, max_size)
        :param n_structures: Number of structures to generate
        :return: List of structures
        """
        structures = []
        sizes = np.linspace(size_range[0], size_range[1], n_structures)
        
        for size in sizes:
            atoms = self.create(base_structure, size=size, **kwargs)
            structures.append(atoms)
            
        return structures

    def apply_strain(self, strain: float, 
                    direction: str = 'xyz') -> GAM_Atom:
        """
        Apply strain to the structure.
        
        :param atoms: Input structure
        :param strain: Strain factor (positive for tension, negative for compression)
        :param direction: Direction(s) to apply strain
        :return: Strained structure
        """
        strained_atoms = self.ASE_atoms.copy()
        cell = self.ASE_atoms.get_cell()
        
        for i, axis in enumerate(['x', 'y', 'z']):
            if axis in direction:
                cell[i] *= (1 + strain)
                
        strained_atoms.set_cell(cell, scale_atoms=True)

        self.ASE_atoms=strained_atoms
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)


        #return strained_atoms
        #return self.atoms


    def modify_surface(self, ligands: List[str],
                      coverage: float = 1.0) -> GAM_Atom:
        """
        Modify surface with ligands.
        
        :param atoms: Input structure
        :param ligands: List of ligand elements
        :param coverage: Surface coverage (0-1)
        :return: Surface modified structure
        """
        modified = self.ASE_atoms.copy()
        positions = self.ASE_atoms.get_positions()
        center = positions.mean(axis=0)
        
        # Find surface atoms (simple distance-based criterion)
        distances = np.linalg.norm(positions - center, axis=1)
        max_dist = np.max(distances) * 0.9
        surface_indices = np.where(distances >= max_dist)[0]
        
        # Add ligands
        n_ligands = int(len(surface_indices) * coverage)
        selected_indices = random.sample(list(surface_indices), n_ligands)
        
        for idx in selected_indices:
            pos = positions[idx]
            direction = pos - center
            direction = direction / np.linalg.norm(direction)
            ligand_pos = pos + direction * 2.0  # 2 Å bond length
            modified.append(Atom(random.choice(ligands), ligand_pos))
            
            
        self.ASE_atoms=modified
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)


        #return modified
        #return self.atoms


    def create_core_shell(self, core_material: str, shell_material: str,
                         shell_thickness: float, shell_lattice_constant: float = None) -> GAM_Atom:
        """
        Create a core-shell nanostructure.
        
        :param core_material: Material for the core (e.g., "Au", "CdSe")
        :param shell_material: Material for the shell (e.g., "Ag", "CdS")
        :param shell_thickness: Thickness of the shell in Å
        :param shell_lattice_constant: Lattice constant for shell material (optional)
        :return: Core-shell structure
        """
        # Get current structure as reference
        current_positions = self.ASE_atoms.get_positions()
        center = np.mean(current_positions, axis=0)
        
        # Calculate distances from center
        distances = np.linalg.norm(current_positions - center, axis=1)
        max_core_distance = np.max(distances)
        
        # Create new atoms object starting with core
        core_shell_atoms = self.ASE_atoms.copy()
        
        # Determine shell lattice constant
        if shell_lattice_constant is None:
            # Use default values for common materials
            default_lattice_constants = {
                "Ag": 4.09, "Au": 4.08, "Cu": 3.61, "Pt": 3.92,
                "CdS": 5.83, "CdSe": 6.08, "ZnS": 5.41, "ZnSe": 5.67
            }
            shell_lattice_constant = default_lattice_constants.get(shell_material, 4.0)
        
        # Create shell by adding atoms around the core
        shell_radius = max_core_distance + shell_thickness
        
        # Generate shell atoms using a simple cubic grid approach
        # This is a simplified approach - for more accurate structures, 
        # you might want to use proper crystal structure generation
        
        # Calculate grid spacing based on shell material
        grid_spacing = shell_lattice_constant / 2.0
        
        # Generate grid points for shell
        x_range = np.arange(center[0] - shell_radius, center[0] + shell_radius + grid_spacing, grid_spacing)
        y_range = np.arange(center[1] - shell_radius, center[1] + shell_radius + grid_spacing, grid_spacing)
        z_range = np.arange(center[2] - shell_radius, center[2] + shell_radius + grid_spacing, grid_spacing)
        
        shell_atoms_added = 0
        for x in x_range:
            for y in y_range:
                for z in z_range:
                    pos = np.array([x, y, z])
                    distance_from_center = np.linalg.norm(pos - center)
                    
                    # Check if position is in shell region (between core surface and shell surface)
                    if (distance_from_center > max_core_distance and 
                        distance_from_center <= shell_radius):
                        
                        # Check if position is not too close to existing atoms
                        min_distance_to_existing = np.min(np.linalg.norm(current_positions - pos, axis=1))
                        
                        if min_distance_to_existing > 1.5:  # Minimum 1.5 Å separation
                            new_atom = Atom(shell_material, pos)
                            core_shell_atoms.append(new_atom)
                            shell_atoms_added += 1
        
        # Update the structure
        self.ASE_atoms = core_shell_atoms
        self.atoms = self._to_GAM_Atoms(self.ASE_atoms)
        
        # Update metadata
        self.meta['core_material'] = core_material
        self.meta['shell_material'] = shell_material
        self.meta['shell_thickness'] = shell_thickness
        self.meta['shell_lattice_constant'] = shell_lattice_constant
        
        #return self.atoms


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
    








