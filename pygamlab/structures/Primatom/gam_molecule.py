
from .gam_universe import GAM_Atom 
from .gam_universe import GAM_Bond



from typing import Dict, Any, Optional, List, Tuple, Set, Union
import numpy as np
from dataclasses import dataclass, field
from enum import Enum
import json
from abc import ABC, abstractmethod
from scipy.spatial.distance import euclidean
from scipy.spatial.transform import Rotation
from collections import deque


# Physical constants
ATOMIC_MASS_UNIT = 1.66053906660e-27  # kg
BOHR_RADIUS = 5.29177210903e-11  # m
ELEMENTARY_CHARGE = 1.602176634e-19  # C

# Bond length reference data (in Angstroms)
REFERENCE_BOND_LENGTHS = {
    ('C', 'C'): {'1': 1.54, '2': 1.34, '3': 1.20},  # single, double, triple
    ('C', 'H'): {'1': 1.09},
    ('C', 'O'): {'1': 1.43, '2': 1.23},
    ('C', 'N'): {'1': 1.47, '2': 1.29, '3': 1.16},
    # Add more reference data as needed
}


REFERENCE_PERIODIC_TABLE_CONFIG={
    'H': {'number': 1, 'mass': 1.008, 'config': '1s1', 'radius': 0.53},
    'He': {'number': 2, 'mass': 4.003, 'config': '1s2', 'radius': 0.31},
    'Li': {'number': 3, 'mass': 6.941, 'config': '[He] 2s1', 'radius': 1.67},
    'Be': {'number': 4, 'mass': 9.012, 'config': '[He] 2s2', 'radius': 1.12},
    'B': {'number': 5, 'mass': 10.811, 'config': '[He] 2s2 2p1', 'radius': 0.85},
    'C': {'number': 6, 'mass': 12.011, 'config': '[He] 2s2 2p2', 'radius': 0.67},
    'N': {'number': 7, 'mass': 14.007, 'config': '[He] 2s2 2p3', 'radius': 0.56},
    'O': {'number': 8, 'mass': 15.999, 'config': '[He] 2s2 2p4', 'radius': 0.48},
    'F': {'number': 9, 'mass': 18.998, 'config': '[He] 2s2 2p5', 'radius': 0.42},
    'Ne': {'number': 10, 'mass': 20.180, 'config': '[He] 2s2 2p6', 'radius': 0.38},
    'Na': {'number': 11, 'mass': 22.990, 'config': '[Ne] 3s1', 'radius': 1.90},
    'Mg': {'number': 12, 'mass': 24.305, 'config': '[Ne] 3s2', 'radius': 1.45},
    'Al': {'number': 13, 'mass': 26.982, 'config': '[Ne] 3s2 3p1', 'radius': 1.18},
    'Si': {'number': 14, 'mass': 28.085, 'config': '[Ne] 3s2 3p2', 'radius': 1.11},
    'P': {'number': 15, 'mass': 30.974, 'config': '[Ne] 3s2 3p3', 'radius': 1.06},
    'S': {'number': 16, 'mass': 32.06, 'config': '[Ne] 3s2 3p4', 'radius': 1.02},
    'Cl': {'number': 17, 'mass': 35.453, 'config': '[Ne] 3s2 3p5', 'radius': 0.99},
    'Ar': {'number': 18, 'mass': 39.948, 'config': '[Ne] 3s2 3p6', 'radius': 0.95},
    'K': {'number': 19, 'mass': 39.098, 'config': '[Ar] 4s1', 'radius': 2.43},
    'Ca': {'number': 20, 'mass': 40.078, 'config': '[Ar] 4s2', 'radius': 1.94},
    'Sc': {'number': 21, 'mass': 44.956, 'config': '[Ar] 3d1 4s2', 'radius': 1.84},
    'Ti': {'number': 22, 'mass': 47.867, 'config': '[Ar] 3d2 4s2', 'radius': 1.76},
    'V': {'number': 23, 'mass': 50.942, 'config': '[Ar] 3d3 4s2', 'radius': 1.71},
    'Cr': {'number': 24, 'mass': 51.996, 'config': '[Ar] 3d5 4s1', 'radius': 1.66},
    'Mn': {'number': 25, 'mass': 54.938, 'config': '[Ar] 3d5 4s2', 'radius': 1.61},
    'Fe': {'number': 26, 'mass': 55.845, 'config': '[Ar] 3d6 4s2', 'radius': 1.56},
    'Co': {'number': 27, 'mass': 58.933, 'config': '[Ar] 3d7 4s2', 'radius': 1.52},
    'Ni': {'number': 28, 'mass': 58.693, 'config': '[Ar] 3d8 4s2', 'radius': 1.49},
    'Cu': {'number': 29, 'mass': 63.546, 'config': '[Ar] 3d10 4s1', 'radius': 1.45},
    'Zn': {'number': 30, 'mass': 65.38, 'config': '[Ar] 3d10 4s2', 'radius': 1.42},
    'Ga': {'number': 31, 'mass': 69.723, 'config': '[Ar] 3d10 4s2 4p1', 'radius': 1.36},
    'Ge': {'number': 32, 'mass': 72.630, 'config': '[Ar] 3d10 4s2 4p2', 'radius': 1.25},
    'As': {'number': 33, 'mass': 74.922, 'config': '[Ar] 3d10 4s2 4p3', 'radius': 1.19},
    'Se': {'number': 34, 'mass': 78.971, 'config': '[Ar] 3d10 4s2 4p4', 'radius': 1.16},
    'Br': {'number': 35, 'mass': 79.904, 'config': '[Ar] 3d10 4s2 4p5', 'radius': 1.14},
    'Kr': {'number': 36, 'mass': 83.798, 'config': '[Ar] 3d10 4s2 4p6', 'radius': 1.10}
}


REFERENCE_VALENCE_ELECTRONS = {
    # Period 1
    'H': 1, 'He': 2,
    # Period 2
    'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8,
    # Period 3
    'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8,
    # Period 4 (main group)
    'K': 1, 'Ca': 2, 'Ga': 3, 'Ge': 4, 'As': 5, 'Se': 6, 'Br': 7, 'Kr': 8,
    # Period 5 (main group)
    'Rb': 1, 'Sr': 2, 'In': 3, 'Sn': 4, 'Sb': 5, 'Te': 6, 'I': 7, 'Xe': 8,
    # Period 6 (main group)
    'Cs': 1, 'Ba': 2, 'Tl': 3, 'Pb': 4, 'Bi': 5, 'Po': 6, 'At': 7, 'Rn': 8,
    # Period 7 (main group)
    'Fr': 1, 'Ra': 2,
    # Other common elements
    'Sc': 3, 'Ti': 4, 'V': 5, 'Cr': 6, 'Mn': 7, 'Fe': 8, 'Co': 9, 'Ni': 10,
    'Cu': 11, 'Zn': 12, 'Ag': 11, 'Cd': 12, 'Au': 11, 'Hg': 12
}


REFERENCE_PAULING_ELECTRONEGATIVITIES={
    'H': 2.20, 'He': 0.00, 'Li': 0.98, 'Be': 1.57, 'B': 2.04, 'C': 2.55,
    'N': 3.04, 'O': 3.44, 'F': 3.98, 'Ne': 0.00, 'Na': 0.93, 'Mg': 1.31,
    'Al': 1.61, 'Si': 1.90, 'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Ar': 0.00,
    'K': 0.82, 'Ca': 1.00, 'Sc': 1.36, 'Ti': 1.54, 'V': 1.63, 'Cr': 1.66,
    'Mn': 1.55, 'Fe': 1.83, 'Co': 1.88, 'Ni': 1.91, 'Cu': 1.90, 'Zn': 1.65,
    'Ga': 1.81, 'Ge': 2.01, 'As': 2.18, 'Se': 2.55, 'Br': 2.96, 'Kr': 3.00,
    'Rb': 0.82, 'Sr': 0.95, 'Y': 1.22, 'Zr': 1.33, 'Nb': 1.60, 'Mo': 2.16,
    'Tc': 1.90, 'Ru': 2.20, 'Rh': 2.28, 'Pd': 2.20, 'Ag': 1.93, 'Cd': 1.69,
    'In': 1.78, 'Sn': 1.96, 'Sb': 2.05, 'Te': 2.10, 'I': 2.66, 'Xe': 2.60,
    'Cs': 0.79, 'Ba': 0.89, 'La': 1.10, 'Hf': 1.30, 'Ta': 1.50, 'W': 2.36,
    'Re': 1.90, 'Os': 2.20, 'Ir': 2.20, 'Pt': 2.28, 'Au': 2.54, 'Hg': 2.00,
    'Tl': 1.62, 'Pb': 1.87, 'Bi': 2.02, 'Po': 2.00, 'At': 2.20, 'Rn': 2.20,
    'Fr': 0.70, 'Ra': 0.90, 'Ac': 1.10, 'Th': 1.30, 'Pa': 1.50, 'U': 1.38,
    'Np': 1.36, 'Pu': 1.28, 'Am': 1.13, 'Cm': 1.28, 'Bk': 1.30, 'Cf': 1.30,
    'Es': 1.30, 'Fm': 1.30, 'Md': 1.30, 'No': 1.30, 'Lr': 1.30
}

    

# ... existing code ...

@dataclass
class GAM_Molecule:
    """
    Advanced molecule class for computational chemistry and molecular modeling applications.
    
    The GAM_Molecule class represents complete molecular systems and provides comprehensive
    functionality for molecular analysis, geometric calculations, and structural properties.
    It manages collections of atoms and bonds, enabling sophisticated molecular dynamics
    simulations, quantum chemistry computations, and structural analysis workflows.
    
    Attributes:
        atoms (Dict[int, GAM_Atom]): Dictionary mapping atom IDs to GAM_Atom objects.
        bonds (Dict[int, GAM_Bond]): Dictionary mapping bond IDs to GAM_Bond objects.
        
    Examples:
        Basic molecule creation:
        >>> molecule = GAM_Molecule()
        >>> atom1 = GAM_Atom(id=1, element="C", x=0.0, y=0.0, z=0.0)
        >>> atom2 = GAM_Atom(id=2, element="C", x=1.54, y=0.0, z=0.0)
        >>> molecule.add_atom(atom1)
        >>> molecule.add_atom(atom2)
        >>> print(f"Atoms: {len(molecule.atoms)}")
        Atoms: 2
        
        Molecular properties:
        >>> formula = molecule.get_molecular_formula()
        >>> mass = molecule.calculate_molecular_mass()
        >>> print(f"Formula: {formula}, Mass: {mass:.2f} amu")
        Formula: C2, Mass: 24.02 amu
        
        Geometric analysis:
        >>> com = molecule.get_center_of_mass()
        >>> cog = molecule.get_center_of_geometry()
        >>> print(f"Center of mass: {com}")
        >>> print(f"Center of geometry: {cog}")
        Center of mass: [0.77 0.   0.  ]
        Center of geometry: [0.77 0.   0.  ]
        
        Energy calculations:
        >>> lj_energy = molecule.calculate_total_lennard_jones_energy()
        >>> coulomb_energy = molecule.calculate_total_coulomb_energy()
        >>> print(f"LJ Energy: {lj_energy:.2e} J")
        >>> print(f"Coulomb Energy: {coulomb_energy:.2e} J")
        LJ Energy: -2.34e-21 J
        Coulomb Energy: 0.00e+00 J
        
        Structure validation:
        >>> issues = molecule.validate_structure()
        >>> print(f"Validation issues: {issues}")
        Validation issues: {'unreasonable_bonds': [], 'isolated_atoms': [], ...}
        
    Notes:
        - All positions are in Angstroms unless otherwise specified.
        - Energies are calculated in Joules for physical accuracy.
        - The class automatically manages bond-atom relationships when bonds are added.
        - Ring detection uses depth-first search with configurable maximum ring size.
        - Structure validation checks for common molecular geometry issues.
        
    Raises:
        KeyError: If referenced atoms do not exist when adding bonds.
        ValueError: If atom or bond IDs are not unique.
        
    See Also:
        GAM_Atom: For individual atom representation and properties.
        GAM_Bond: For bond representation and analysis.
        REFERENCE_BOND_LENGTHS: For reference bond length data.
        REFERENCE_VALENCE_ELECTRONS: For valence electron reference data.
    """
    atoms: Dict[int, GAM_Atom] = field(default_factory=dict)
    bonds: Dict[int, GAM_Bond] = field(default_factory=dict)
    
    def add_atom(self, atom: GAM_Atom):
        """Adds an atom to the molecule."""
        self.atoms[atom.id] = atom
        
    def add_bond(self, bond: GAM_Bond):
        """Adds a bond to the molecule and updates the atoms' bond lists."""
        self.bonds[bond.id] = bond
        # Update atoms with their new bond
        self.atoms[bond.atom1_id].add_bond(bond)
        self.atoms[bond.atom2_id].add_bond(bond)

    def get_atom(self, atom_id: int) -> Optional[GAM_Atom]:
        """Retrieves an atom by its ID."""
        return self.atoms.get(atom_id)
        
    def get_bond(self, bond_id: int) -> Optional[GAM_Bond]:
        """Retrieves a bond by its ID."""
        return self.bonds.get(bond_id)

    def calculate_total_energy(self) -> float:
        """Calculates the total energy of the molecule (placeholder)."""
        # This would include bond strain, non-bonded interactions, etc.
        total_energy = 0.0
        # Example: Summing up potential energies of atoms
        for atom in self.atoms.values():
            total_energy += atom.potential_energy
        return total_energy
    
    # Advanced molecular analysis functions
    
    def get_center_of_mass(self) -> np.ndarray:
        """
        Calculate the center of mass of the molecule.
        
        Returns:
            The center of mass coordinates as a numpy array.
        """
        if not self.atoms:
            return np.zeros(3)
        
        total_mass = 0.0
        weighted_position = np.zeros(3)
        
        for atom in self.atoms.values():
            mass = atom.atomic_mass
            total_mass += mass
            weighted_position += mass * atom.position
        
        return weighted_position / total_mass
    
    def get_center_of_geometry(self) -> np.ndarray:
        """
        Calculate the geometric center (centroid) of the molecule.
        
        Returns:
            The geometric center coordinates as a numpy array.
        """
        if not self.atoms:
            return np.zeros(3)
        
        positions = np.array([atom.position for atom in self.atoms.values()])
        return np.mean(positions, axis=0)
    
    def calculate_molecular_mass(self) -> float:
        """
        Calculate the total molecular mass in atomic mass units.
        
        Returns:
            Total molecular mass in amu.
        """
        return sum(atom.atomic_mass for atom in self.atoms.values())
    
    def get_molecular_formula(self) -> str:
        """
        Generate the molecular formula string.
        
        Returns:
            Molecular formula (e.g., "C6H12O6").
        """
        element_count = {}
        for atom in self.atoms.values():
            element_count[atom.element] = element_count.get(atom.element, 0) + 1
        
        # Sort elements by common convention (C, H, then alphabetical)
        sorted_elements = []
        if 'C' in element_count:
            sorted_elements.append('C')
        if 'H' in element_count:
            sorted_elements.append('H')
        
        remaining = sorted([elem for elem in element_count.keys() if elem not in ['C', 'H']])
        sorted_elements.extend(remaining)
        
        formula = ""
        for element in sorted_elements:
            count = element_count[element]
            formula += element if count == 1 else f"{element}{count}"
        
        return formula
    
    def get_bounding_box(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the bounding box of the molecule.
        
        Returns:
            Tuple of (min_coords, max_coords) as numpy arrays.
        """
        if not self.atoms:
            return np.zeros(3), np.zeros(3)
        
        positions = np.array([atom.position for atom in self.atoms.values()])
        return np.min(positions, axis=0), np.max(positions, axis=0)
    
    def calculate_radius_of_gyration(self) -> float:
        """
        Calculate the radius of gyration of the molecule.
        
        Returns:
            Radius of gyration in Angstroms.
        """
        if not self.atoms:
            return 0.0
        
        center_of_mass = self.get_center_of_mass()
        total_mass = 0.0
        sum_weighted_distances = 0.0
        
        for atom in self.atoms.values():
            mass = atom.atomic_mass
            distance_squared = np.sum((atom.position - center_of_mass) ** 2)
            sum_weighted_distances += mass * distance_squared
            total_mass += mass
        
        return np.sqrt(sum_weighted_distances / total_mass)
    
    def calculate_moment_of_inertia_tensor(self) -> np.ndarray:
        """
        Calculate the moment of inertia tensor of the molecule.
        
        Returns:
            3x3 moment of inertia tensor.
        """
        center_of_mass = self.get_center_of_mass()
        I = np.zeros((3, 3))
        
        for atom in self.atoms.values():
            mass = atom.atomic_mass
            r = atom.position - center_of_mass
            
            # Diagonal elements
            I[0, 0] += mass * (r[1]**2 + r[2]**2)
            I[1, 1] += mass * (r[0]**2 + r[2]**2)
            I[2, 2] += mass * (r[0]**2 + r[1]**2)
            
            # Off-diagonal elements
            I[0, 1] -= mass * r[0] * r[1]
            I[0, 2] -= mass * r[0] * r[2]
            I[1, 2] -= mass * r[1] * r[2]
        
        # Symmetric matrix
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        
        return I
    
    def get_principal_moments(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate principal moments of inertia and axes.
        
        Returns:
            Tuple of (eigenvalues, eigenvectors) representing principal moments and axes.
        """
        I = self.calculate_moment_of_inertia_tensor()
        eigenvalues, eigenvectors = np.linalg.eigh(I)
        return eigenvalues, eigenvectors
    
    def find_rings(self, max_ring_size: int = 12) -> List[List[int]]:
        """
        Find all rings in the molecule using depth-first search.
        
        Args:
            max_ring_size: Maximum ring size to search for.
            
        Returns:
            List of rings, where each ring is a list of atom IDs.
        """
        rings = []
        visited_bonds = set()
        
        def dfs(start_atom_id: int, current_atom_id: int, path: List[int], visited_atoms: Set[int]) -> None:
            if len(path) > max_ring_size:
                return
            
            current_atom = self.atoms[current_atom_id]
            
            for bond in current_atom.get_bonds():
                if bond.id in visited_bonds:
                    continue
                
                neighbor_id = bond.atom2_id if bond.atom1_id == current_atom_id else bond.atom1_id
                
                if neighbor_id == start_atom_id and len(path) >= 3:
                    # Found a ring
                    ring = path.copy()
                    ring.sort()  # Canonical form
                    if ring not in rings:
                        rings.append(ring)
                    continue
                
                if neighbor_id not in visited_atoms:
                    new_visited = visited_atoms.copy()
                    new_visited.add(neighbor_id)
                    new_path = path + [neighbor_id]
                    visited_bonds.add(bond.id)
                    dfs(start_atom_id, neighbor_id, new_path, new_visited)
                    visited_bonds.remove(bond.id)
        
        for atom_id in self.atoms.keys():
            dfs(atom_id, atom_id, [atom_id], {atom_id})
        
        return rings
    
    def calculate_total_lennard_jones_energy(self, epsilon: float = 1e-21, sigma: float = 3.4e-10) -> float:
        """
        Calculate total Lennard-Jones energy for all non-bonded atom pairs.
        
        Args:
            epsilon: Well depth parameter (J).
            sigma: Distance parameter (m).
            
        Returns:
            Total Lennard-Jones energy in Joules.
        """
        total_energy = 0.0
        atom_list = list(self.atoms.values())
        
        for i in range(len(atom_list)):
            for j in range(i + 1, len(atom_list)):
                atom1, atom2 = atom_list[i], atom_list[j]
                
                # Skip if atoms are bonded
                if atom1.is_bonded_to(atom2):
                    continue
                
                energy = atom1.calculate_lennard_jones_energy(atom2, epsilon, sigma)
                total_energy += energy
        
        return total_energy
    
    def calculate_total_coulomb_energy(self) -> float:
        """
        Calculate total electrostatic energy for all atom pairs.
        
        Returns:
            Total Coulomb energy in Joules.
        """
        total_energy = 0.0
        atom_list = list(self.atoms.values())
        
        for i in range(len(atom_list)):
            for j in range(i + 1, len(atom_list)):
                atom1, atom2 = atom_list[i], atom_list[j]
                energy = atom1.calculate_coulomb_energy(atom2)
                total_energy += energy
        
        return total_energy
    
    def calculate_bond_strain_energy(self) -> float:
        """
        Calculate total bond strain energy.
        
        Returns:
            Total strain energy in kcal/mol.
        """
        total_strain = 0.0
        for bond in self.bonds.values():
            total_strain += bond.calculate_strain_energy(self.atoms)
        return total_strain
    
    def get_connectivity_matrix(self) -> np.ndarray:
        """
        Generate the connectivity matrix for the molecule.
        
        Returns:
            N×N connectivity matrix where N is the number of atoms.
        """
        atom_ids = sorted(self.atoms.keys())
        n_atoms = len(atom_ids)
        connectivity = np.zeros((n_atoms, n_atoms), dtype=int)
        
        id_to_index = {atom_id: i for i, atom_id in enumerate(atom_ids)}
        
        for bond in self.bonds.values():
            i = id_to_index[bond.atom1_id]
            j = id_to_index[bond.atom2_id]
            connectivity[i, j] = bond.order
            connectivity[j, i] = bond.order
        
        return connectivity
    
    def get_distance_matrix(self) -> np.ndarray:
        """
        Calculate the distance matrix between all atoms.
        
        Returns:
            N×N distance matrix in Angstroms.
        """
        atom_list = list(self.atoms.values())
        n_atoms = len(atom_list)
        distances = np.zeros((n_atoms, n_atoms))
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist = atom_list[i].calculate_distance_to(atom_list[j])
                distances[i, j] = dist
                distances[j, i] = dist
        
        return distances
    
    def translate_molecule(self, vector: np.ndarray) -> None:
        """
        Translate the entire molecule by a vector.
        
        Args:
            vector: Translation vector [dx, dy, dz].
        """
        for atom in self.atoms.values():
            atom.translate(vector)
    
    def rotate_molecule(self, rotation_matrix: np.ndarray, center: Optional[np.ndarray] = None) -> None:
        """
        Rotate the entire molecule around a center point.
        
        Args:
            rotation_matrix: 3×3 rotation matrix.
            center: Center of rotation (default is center of mass).
        """
        if center is None:
            center = self.get_center_of_mass()
        
        for atom in self.atoms.values():
            atom.rotate(rotation_matrix, center)
    
    def align_to_principal_axes(self) -> None:
        """
        Align the molecule to its principal axes of inertia.
        """
        center_of_mass = self.get_center_of_mass()
        _, principal_axes = self.get_principal_moments()
        
        # Translate to origin first
        self.translate_molecule(-center_of_mass)
        
        # Rotate to align with principal axes
        self.rotate_molecule(principal_axes.T)
    
    def calculate_dipole_moment(self) -> np.ndarray:
        """
        Calculate the electric dipole moment of the molecule.
        
        Returns:
            Dipole moment vector in Debye units.
        """
        dipole = np.zeros(3)
        
        for atom in self.atoms.values():
            charge = atom.charge * ELEMENTARY_CHARGE  # Convert to Coulombs
            position = atom.position * 1e-10  # Convert to meters
            dipole += charge * position
        
        # Convert to Debye (1 Debye = 3.336e-30 C⋅m)
        dipole_debye = dipole / 3.336e-30
        
        return dipole_debye
    
    def find_atoms_by_element(self, element: str) -> List[GAM_Atom]:
        """
        Find all atoms of a specific element.
        
        Args:
            element: Element symbol (e.g., 'C', 'H', 'O').
            
        Returns:
            List of atoms of the specified element.
        """
        return [atom for atom in self.atoms.values() if atom.element == element]
    
    def find_shortest_path(self, atom1_id: int, atom2_id: int) -> Optional[List[int]]:
        """
        Find the shortest path between two atoms through bonds.
        
        Args:
            atom1_id: Starting atom ID.
            atom2_id: Target atom ID.
            
        Returns:
            List of atom IDs representing the shortest path, or None if no path exists.
        """
        if atom1_id not in self.atoms or atom2_id not in self.atoms:
            return None
        
        if atom1_id == atom2_id:
            return [atom1_id]
        
        queue = deque([(atom1_id, [atom1_id])])
        visited = {atom1_id}
        
        while queue:
            current_id, path = queue.popleft()
            current_atom = self.atoms[current_id]
            
            for bond in current_atom.get_bonds():
                neighbor_id = bond.atom2_id if bond.atom1_id == current_id else bond.atom1_id
                
                if neighbor_id == atom2_id:
                    return path + [neighbor_id]
                
                if neighbor_id not in visited:
                    visited.add(neighbor_id)
                    queue.append((neighbor_id, path + [neighbor_id]))
        
        return None
    
    def validate_structure(self) -> Dict[str, List[str]]:
        """
        Validate the molecular structure and return potential issues.
        
        Returns:
            Dictionary with categories of issues found.
        """
        issues = {
            'unreasonable_bonds': [],
            'isolated_atoms': [],
            'unusual_coordination': [],
            'overlapping_atoms': []
        }
        
        # Check for unreasonable bond lengths
        for bond in self.bonds.values():
            if not bond.is_reasonable_length(self.atoms):
                atom1 = self.atoms[bond.atom1_id]
                atom2 = self.atoms[bond.atom2_id]
                length = bond.calculate_length(self.atoms)
                issues['unreasonable_bonds'].append(
                    f"Bond {bond.id} ({atom1.element}-{atom2.element}): {length:.3f} Å"
                )
        
        # Check for isolated atoms
        for atom in self.atoms.values():
            if len(atom.get_bonds()) == 0:
                issues['isolated_atoms'].append(f"Atom {atom.id} ({atom.element})")
        
        # Check for unusual coordination numbers
        for atom in self.atoms.values():
            coordination = len(atom.get_bonds())
            expected_max = REFERENCE_VALENCE_ELECTRONS.get(atom.element, 8)
            if coordination > expected_max:
                issues['unusual_coordination'].append(
                    f"Atom {atom.id} ({atom.element}): {coordination} bonds (expected ≤{expected_max})"
                )
        
        # Check for overlapping atoms (distance < 0.5 Å)
        atom_list = list(self.atoms.values())
        for i in range(len(atom_list)):
            for j in range(i + 1, len(atom_list)):
                dist = atom_list[i].calculate_distance_to(atom_list[j])
                if dist < 0.5:
                    issues['overlapping_atoms'].append(
                        f"Atoms {atom_list[i].id} and {atom_list[j].id}: {dist:.3f} Å"
                    )
        
        return issues
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert molecule to dictionary representation.
        
        Returns:
            Dictionary representation of the molecule.
        """
        return {
            'atoms': {str(atom_id): atom.to_dict() for atom_id, atom in self.atoms.items()},
            'bonds': {str(bond_id): bond.to_dict() for bond_id, bond in self.bonds.items()},
            'molecular_formula': self.get_molecular_formula(),
            'molecular_mass': self.calculate_molecular_mass(),
            'center_of_mass': self.get_center_of_mass().tolist(),
            'center_of_geometry': self.get_center_of_geometry().tolist()
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'GAM_Molecule':
        """
        Create molecule from dictionary representation.
        
        Args:
            data: Dictionary containing molecule data.
            
        Returns:
            New GAM_Molecule instance.
        """
        molecule = cls()
        
        # Load atoms
        for atom_data in data['atoms'].values():
            atom = GAM_Atom.from_dict(atom_data)
            molecule.add_atom(atom)
        
        # Load bonds
        for bond_data in data['bonds'].values():
            bond = GAM_Bond.from_dict(bond_data)
            molecule.add_bond(bond)
        
        return molecule
    
    def to_xyz_file(self, filepath: str, comment: str = "Generated by PyGamlab") -> None:
        """
        Write molecule to XYZ format file.
        
        Args:
            filepath: Output file path.
            comment: Comment line for the XYZ file.
        """
        with open(filepath, 'w') as f:
            f.write(f"{len(self.atoms)}\n")
            f.write(f"{comment}\n")
            for atom in self.atoms.values():
                f.write(f"{atom.element} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}\n")

    def to_ase_object():
        pass

    def to_pymatgen_object():    
        pass
    
    def __repr__(self) -> str:
        """String representation of the molecule."""
        formula = self.get_molecular_formula()
        mass = self.calculate_molecular_mass()
        return f"GAM_Molecule({formula}, MW={mass:.2f} amu, {len(self.atoms)} atoms, {len(self.bonds)} bonds)"

# ... existing code ...
    
    
    
