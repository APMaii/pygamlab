
#from .gam_bond import GAM_Bond


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
    'Au': {'number': 2, 'mass': 4.003, 'config': '1s2', 'radius': 0.31},
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
    'Kr': {'number': 36, 'mass': 83.798, 'config': '[Ar] 3d10 4s2 4p6', 'radius': 1.10},
    'Cd': {'number': 48, 'mass': 112.411, 'config': '[Kr] 4d10 5s2', 'radius': 1.55},
    'Pt': {'number': 78, 'mass': 195.084, 'config': '[Xe] 4f14 5d9 6s1', 'radius': 1.77},
    'Ag': {'number': 47, 'mass': 107.868, 'config': '[Kr] 4d10 5s1', 'radius': 1.65},
    'W': {'number': 74, 'mass': 183.84, 'config': '[Xe] 4f14 5d4 6s2', 'radius': 0.53},
    'Nb': {'number': 41,'mass': 92.906, 'config': '[Kr] 4d4 5s1', 'radius': 1.46},
    'Mo': {'number': 42,'mass': 95.94, 'config': '[Kr] 4d5 5s1', 'radius': 1.39},
    'Ta': {'number': 73,'mass': 180.948, 'config': '[Xe] 4f14 5d3 6s2', 'radius': 1.43}

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
    'Cu': 11, 'Zn': 12, 'Ag': 11, 'Cd': 12, 'Au': 11, 'Hg': 12, 'Pt': 10, 'Nb': 5,
    'Mo': 6, 'Ta': 5
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

class ElectronicState(Enum):
    """Enumeration for electronic states of atoms."""
    GROUND_STATE = "ground"
    EXCITED_STATE = "excited"
    IONIZED = "ionized"

@dataclass
class AtomicOrbital:
    """Represents an atomic orbital with quantum numbers."""
    n: int  # Principal quantum number
    l: int  # Angular momentum quantum number
    m: int  # Magnetic quantum number
    s: float  # Spin quantum number
    occupation: float = 0.0  # Electron occupation number

class AtomicData:
    """Static class containing fundamental atomic data."""
    _periodic_table: Dict[str, Dict[str, Any]] =REFERENCE_PERIODIC_TABLE_CONFIG  

    @classmethod
    def get_atomic_data(cls, element: str) -> Dict[str, Any]:
        """Get atomic data for an element."""
        return cls._periodic_table.get(element, {})

class ForceField(ABC):
    """Abstract base class for force field implementations."""
    @abstractmethod
    def calculate_energy(self, atom: 'GAM_Atom', neighbors: List['GAM_Atom']) -> float:
        """Calculate the energy contribution of an atom."""
        pass




class GAM_Atom:
    """
    Advanced atom class for computational chemistry and molecular modeling applications.
    
    The GAM_Atom class represents individual atoms within molecular systems and provides
    comprehensive functionality for scientific calculations including molecular dynamics
    simulations, quantum chemistry computations, and structural analysis. It maintains
    backward compatibility with simple atomic representations while offering advanced
    features for sophisticated scientific workflows.
    
    Attributes:
        id (int): Unique identifier for the atom within a molecular system.
        element (str): Chemical element symbol (e.g., 'C', 'H', 'O', 'N').
        x (float): X-coordinate position in Angstroms.
        y (float): Y-coordinate position in Angstroms.
        z (float): Z-coordinate position in Angstroms.
        charge (float): Formal charge of the atom in elementary charge units.
        electronic_state (ElectronicState): Electronic state enumeration (ground, excited, ionized).
        atomic_number (int): Atomic number from periodic table.
        atomic_mass (float): Atomic mass in atomic mass units (AMU).
        electron_configuration (str): Electron configuration string.
        covalent_radius (float): Covalent radius in Angstroms.
        velocity (np.ndarray): 3D velocity vector in Angstroms/ps.
        force (np.ndarray): 3D force vector for molecular dynamics.
        potential_energy (float): Potential energy contribution in Joules.
        kinetic_energy (float): Kinetic energy in Joules.
        metadata (Dict[str, Any]): Additional metadata for custom properties.
        orbitals (List[AtomicOrbital]): List of atomic orbitals with quantum numbers.
        bonds (List[GAM_Bond]): List of bonds connected to this atom.
        
    Examples:
        Basic atom creation:
        >>> atom = GAM_Atom(id=1, element="C", x=0.0, y=0.0, z=0.0)
        >>> print(atom.element)
        C
        
        Advanced atom with charge and electronic state:
        >>> from ElectronicState import ElectronicState
        >>> charged_atom = GAM_Atom(
        ...     id=2, element="O", x=1.5, y=0.0, z=0.0, 
        ...     charge=-0.5, electronic_state=ElectronicState.EXCITED_STATE
        ... )
        
        Distance calculation between atoms:
        >>> atom1 = GAM_Atom(id=1, element="C", x=0.0, y=0.0, z=0.0)
        >>> atom2 = GAM_Atom(id=2, element="C", x=1.54, y=0.0, z=0.0)
        >>> distance = atom1.calculate_distance_to(atom2)
        >>> print(f"Distance: {distance:.2f} Å")
        Distance: 1.54 Å
        
        Molecular dynamics velocity initialization:
        >>> atom = GAM_Atom(id=1, element="H", x=0.0, y=0.0, z=0.0)
        >>> atom.set_velocity_from_temperature(300.0)  # 300 K
        >>> print(f"Velocity magnitude: {np.linalg.norm(atom.velocity):.3f} Å/ps")
        
    Notes:
        - All positions are in Angstroms unless otherwise specified.
        - Charges are in elementary charge units (e).
        - Energies are calculated in Joules for physical accuracy.
        - The class automatically initializes atomic properties from reference data.
        - Bond management requires integration with GAM_Bond class.
        
    Raises:
        ValueError: If element symbol is not found in periodic table reference data.
        ValueError: If position vector has insufficient dimensions.
        ValueError: If velocity array has incorrect shape or type.
        
    See Also:
        GAM_Bond: For bond representation and management.
        AtomicData: For periodic table reference data.
        ElectronicState: For electronic state enumeration.
    """
    
    def __init__(self, 
                 id: int,
                 element: str = "C",
                 x: float = 0.0,
                 y: float = 0.0,
                 z: float = 0.0,
                 metadata: Optional[Dict[str, Any]] = None,
                 charge: float = 0.0,
                 electronic_state: ElectronicState = ElectronicState.GROUND_STATE):
        """Initialize an Atom instance with both original and advanced properties."""
        # Original properties
        self.id = id
        self.element = element.capitalize()
        self.x = x
        self.y = y
        self.z = z

        #----------
        self.metadata = metadata or {}

        # Advanced properties
        self.charge = charge
        self.electronic_state = electronic_state
        self._velocity = np.zeros(3)
        
        # Initialize atomic properties from periodic table
        atomic_data = AtomicData.get_atomic_data(self.element)
        if not atomic_data:
            raise ValueError(f"Element '{self.element}' not found in periodic table data.")
    

        self.atomic_type: str = self.element  # Default to element name
        self.atomic_number = atomic_data.get('number', 0)
        self.atomic_mass = atomic_data.get('mass', 0.0)
        self.electron_configuration = atomic_data.get('config', '')
        self.covalent_radius = atomic_data.get('radius', 0.0)
        
        # Initialize quantum mechanical properties
        self.orbitals: List[AtomicOrbital] = []
        self._initialize_orbitals()
        
        # Dynamic properties
        self.force = np.zeros(3)
        self.potential_energy = 0.0
        self.kinetic_energy = 0.0

        # New optional properties
        self._bonds: List['GAM_Bond'] = []
        self.bond_partners: List['GAM_Atom'] = []

    def _initialize_orbitals(self):
        """Initialize atomic orbitals based on electron configuration."""
        pass

    @property
    def position(self) -> np.ndarray:
        """Get position as numpy array."""
        return np.array([self.x, self.y, self.z])

    @position.setter
    def position(self, pos: Union[List[float], Tuple[float, ...], np.ndarray]):
        """Set position from array-like."""
        if len(pos) < 3:
            raise ValueError("Position must be a 3D vector")
        self.x, self.y, self.z = float(pos[0]), float(pos[1]), float(pos[2])

    @property
    def velocity(self) -> np.ndarray:
        """Get velocity vector."""
        return self._velocity.copy()

    @velocity.setter
    def velocity(self, vel: np.ndarray):
        """Set velocity vector with bounds checking."""
        if not isinstance(vel, np.ndarray) or vel.shape != (3,):
            raise ValueError("Velocity must be a 3D numpy array")
        self._velocity = vel.copy()

    def rotate(self, rotation_matrix: np.ndarray, center: Optional[np.ndarray] = None) -> None:
        """
        Rotate the atom around a center point using a rotation matrix.
        
        Args:
            rotation_matrix: 3x3 rotation matrix
            center: Center of rotation (default is origin)
        """
        if center is None:
            center = np.zeros(3)
        
        pos = self.position - center
        rotated_pos = np.dot(rotation_matrix, pos)
        self.position = rotated_pos + center

    def rotate_around_axis(self, axis: np.ndarray, angle_degrees: float, center: Optional[np.ndarray] = None) -> None:
        """
        Rotate the atom around an arbitrary axis by a given angle.
        
        Args:
            axis: The axis of rotation (will be normalized)
            angle_degrees: Rotation angle in degrees
            center: Center of rotation (default is origin)
        """
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        rotation = Rotation.from_rotvec(axis * np.radians(angle_degrees))
        self.rotate(rotation.as_matrix(), center)



    #------------------------------------------------

    def reflect(self, plane_normal: np.ndarray) -> None:
        """
        Reflect the atom through a plane defined by its normal vector.
        
        Args:
            plane_normal: Normal vector of the reflection plane
        """
        plane_normal = np.array(plane_normal)
        plane_normal = plane_normal / np.linalg.norm(plane_normal)
        reflection_matrix = np.eye(3) - 2 * np.outer(plane_normal, plane_normal)
        self.position = np.dot(reflection_matrix, self.position)

    
    def mirror(self, plane_normal: np.ndarray, plane_point: Optional[np.ndarray] = None) -> None:
        """
        Reflects the atom's position across a plane.
        
        Args:
            plane_normal: The normal vector of the reflection plane.
            plane_point: A point on the reflection plane. If None, the origin is used.
        """
        plane_normal = np.array(plane_normal)
        plane_normal = plane_normal / np.linalg.norm(plane_normal) # Normalize
        
        if plane_point is None:
            plane_point = np.zeros(3)
        
        # Vector from a point on the plane to the atom
        vec_to_atom = self.position - plane_point
        
        # Calculate the distance from the atom to the plane
        distance_to_plane = np.dot(vec_to_atom, plane_normal)
        
        # The new position is the old position minus twice the projection onto the normal vector
        reflection_vector = 2 * distance_to_plane * plane_normal
        self.position = self.position - reflection_vector



    #---instead of mirror and reflect
    def reflect_across_plane(self, plane_normal: np.ndarray, plane_point: Optional[np.ndarray] = None) -> None:
        """
        Reflects the atom's position across a plane defined by its normal and a point.

        Args:
            plane_normal: The normal vector of the reflection plane.
            plane_point: A point on the reflection plane. If None, the origin is used.
        """
        plane_normal = np.array(plane_normal)
        plane_normal = plane_normal / np.linalg.norm(plane_normal)

        if plane_point is None:
            plane_point = np.zeros(3)

        vec_to_atom = self.position - plane_point
        distance_to_plane = np.dot(vec_to_atom, plane_normal)

        reflection_vector = 2 * distance_to_plane * plane_normal
        self.position = self.position - reflection_vector


    def inversion_center(self, center: Optional[np.ndarray] = None) -> None:
        """
        Inverts the atom's position through a central point.
        
        Args:
            center: The center of inversion. If None, the origin is used.
        """
        if center is None:
            center = np.zeros(3)
        
        # Inversion is equivalent to reflecting through the center
        self.position = center - (self.position - center)


    def translate(self, vector: np.ndarray) -> None:
        """
        Translate the atom by a vector.
        
        Args:
            vector: Translation vector [dx, dy, dz]
        """
        self.position = self.position + np.array(vector)

    def calculate_distance_to(self, other: 'GAM_Atom') -> float:
        """Calculate distance to another atom in Angstroms."""
        return np.linalg.norm(self.position - other.position)

    def calculate_angle(self, atom2: 'GAM_Atom', atom3: 'GAM_Atom') -> float:
        """Calculate angle between three atoms (self-atom2-atom3) in degrees."""
        v1 = self.position - atom2.position
        v2 = atom3.position - atom2.position
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))

    def calculate_dihedral(self, atom2: 'GAM_Atom', atom3: 'GAM_Atom', atom4: 'GAM_Atom') -> float:
        """
        Calculate dihedral angle between four atoms.
        
        Returns:
            Dihedral angle in degrees
        """
        b1 = self.position - atom2.position
        b2 = atom3.position - atom2.position
        b3 = atom4.position - atom3.position

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        
        return np.degrees(np.arctan2(y, x))

    #I think this is correct
    def calculate_improper(self, atom2: 'GAM_Atom', atom3: 'GAM_Atom', atom4: 'GAM_Atom') -> float:
        """
        Calculate improper dihedral angle (atom1 is central).

        Returns:
            Improper angle in degrees
        """
        v1 = atom2.position - self.position
        v2 = atom3.position - self.position
        v3 = atom4.position - self.position

        # Normal vector to the plane defined by atoms self, 2, and 3
        plane_normal = np.cross(v1, v2)

        # Angle between the normal and the vector to atom 4
        cos_angle = np.dot(plane_normal, v3) / (np.linalg.norm(plane_normal) * np.linalg.norm(v3))

        # The improper angle is the arccosine of this value
        return np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))


    def is_bonded_to(self, other: 'GAM_Atom') -> bool:
        """Check if this atom is bonded to another atom."""
        return other.id in self._bonds

    def calculate_kinetic_energy(self) -> float:
        """Calculate kinetic energy in Joules."""
        mass_kg = self.atomic_mass * ATOMIC_MASS_UNIT
        return 0.5 * mass_kg * np.sum(self._velocity ** 2)


    def get_valence_electrons(self) -> int:
        """
        Get the number of valence electrons for the atom.
        
        Returns:
            The number of valence electrons.
        """
        return REFERENCE_VALENCE_ELECTRONS.get(self.element, 0)
    

    def get_bonding_partners(self) -> List['GAM_Atom']:
        """
        Returns a list of GAM_Atom objects this atom is bonded to.
        
        NOTE: This requires a parent Molecule or System class to manage the full list of atoms.
              The implementation here is a placeholder.
        """
        return self._bonds
    

    def apply_lennard_jones_potential(self, neighbor: 'GAM_Atom', epsilon: float, sigma: float) -> float:
        """
        Calculates the Lennard-Jones potential energy with a neighbor atom.
        
        Args:
            neighbor: The neighboring GAM_Atom.
            epsilon: The well depth of the potential.
            sigma: The distance at which the potential is zero.
            
        Returns:
            The Lennard-Jones potential energy contribution.
        """
        r = self.calculate_distance_to(neighbor)
        if r == 0:
            return float('inf')
            
        term1 = (sigma / r) ** 12
        term2 = (sigma / r) ** 6
        
        return 4 * epsilon * (term1 - term2)
    
    def calculate_vdw_energy(self, neighbor: 'GAM_Atom', epsilon: float, sigma: float) -> float:
        """
        Calculates the Van der Waals interaction energy using the Lennard-Jones potential.

        Args:
            neighbor: The neighboring GAM_Atom.
            epsilon: The well depth of the potential, representing the strength of the interaction.
            sigma: The distance at which the potential is zero.

        Returns:
            The Van der Waals potential energy contribution in Joules.
        """
        r = self.calculate_distance_to(neighbor) * 1e-10  # Convert to meters
        if r == 0:
            return float('inf')  # Prevent division by zero
        
        term1 = (sigma / r) ** 12
        term2 = (sigma / r) ** 6

        # The Lennard-Jones equation: E_LJ = 4*epsilon * [ (sigma/r)^12 - (sigma/r)^6 ]
        return 4 * epsilon * (term1 - term2)
    
    #combine calculate vdw and leonard jones in one function belwo :

    def calculate_lennard_jones_energy(self, neighbor: 'GAM_Atom', epsilon: float, sigma: float) -> float:
        """
        Calculates the Lennard-Jones potential energy with a neighbor atom.

        Args:
            neighbor: The neighboring GAM_Atom.
            epsilon: The well depth of the potential (e.g., in Joules).
            sigma: The distance at which the potential is zero (e.g., in meters).

        Returns:
            The Lennard-Jones potential energy contribution in Joules.
        """
        r = self.calculate_distance_to(neighbor) * 1e-10  # Convert to meters
        if r == 0:
            return float('inf')

        term1 = (sigma / r) ** 12
        term2 = (sigma / r) ** 6

        return 4 * epsilon * (term1 - term2)

    
    def calculate_electronegativity(self, scale: str = 'Pauling') -> float:
        """
        Calculates or retrieves the atom's electronegativity from a lookup table.
        
        Args:
            scale: The electronegativity scale to use ('Pauling' or other).
            
        Returns:
            The electronegativity value.
        """
        # A simple lookup table for the Pauling scale
        PAULING_ELECTRONEGATIVITIES = REFERENCE_PAULING_ELECTRONEGATIVITIES
        
        
        if scale == 'Pauling':
            return PAULING_ELECTRONEGATIVITIES.get(self.element, 0.0)
        else:
            raise ValueError(f"Electronegativity scale '{scale}' not supported.")
        
    def set_velocity_from_temperature(self, temperature_k: float) -> None:
        """
        Sets the atom's velocity based on a Maxwell-Boltzmann distribution for a given temperature.

        Args:
            temperature_k: The desired temperature in Kelvin.
        """
        # kB is the Boltzmann constant in J/K
        BOLTZMANN_CONSTANT = 1.380649e-23
        mass_kg = self.atomic_mass * ATOMIC_MASS_UNIT

        # Standard deviation for the Maxwell-Boltzmann distribution
        sigma = np.sqrt(BOLTZMANN_CONSTANT * temperature_k / mass_kg)

        # Generate random velocity components from a normal distribution
        self._velocity = np.random.normal(loc=0.0, scale=sigma, size=3)


    
    
    def calculate_coulomb_energy(self, neighbor: 'GAM_Atom') -> float:
        """
        Calculates the electrostatic Coulomb energy between two atoms.

        Args:
            neighbor: The neighboring GAM_Atom.

        Returns:
            The electrostatic potential energy in Joules.
        """
        # k_e is Coulomb's constant in N·m²/C²
        COULOMB_CONSTANT = 8.9875517923e9
        
        r = self.calculate_distance_to(neighbor) * 1e-10  # Convert to meters
        if r == 0:
            return float('inf')  # Prevent division by zero
        
        # Charges are in elementary charge units, convert to Coulombs
        q1 = self.charge * ELEMENTARY_CHARGE
        q2 = neighbor.charge * ELEMENTARY_CHARGE

        # Coulomb's Law: E_coulomb = k_e * (q1 * q2) / r
        return COULOMB_CONSTANT * (q1 * q2) / r

    def apply_force(self, force: np.ndarray, dt: float):
        """Apply force for time step dt using Velocity Verlet integration."""
        self._velocity += 0.5 * (force / self.atomic_mass) * dt
        new_pos = self.position + self._velocity * dt
        self.x, self.y, self.z = new_pos
        self.force = force

    def calculate_electron_density(self, point: np.ndarray) -> float:
        """Calculate electron density at a point in space."""
        r = np.linalg.norm(point - self.position)
        return self.atomic_number * np.exp(-r / BOHR_RADIUS)

    def to_dict(self) -> Dict[str, Any]:
        """Convert atom to dictionary representation."""
        return {
            'id': self.id,
            'element': self.element,
            'x': self.x,
            'y': self.y,
            'z': self.z,
            'charge': self.charge,
            'electronic_state': self.electronic_state.value,
            'metadata': self.metadata,
            'atomic_number': self.atomic_number,
            'atomic_mass': self.atomic_mass,
            'electron_configuration': self.electron_configuration,
            'velocity': self._velocity.tolist(),
            'covalent_radius': self.covalent_radius
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'GAM_Atom':
        """Create atom instance from dictionary representation."""
        atom = cls(
            id=data['id'],
            element=data['element'],
            x=data['x'],
            y=data['y'],
            z=data['z'],
            metadata=data['metadata'],
            charge=data.get('charge', 0.0),
            electronic_state=ElectronicState(data.get('electronic_state', 'ground'))
        )
        if 'velocity' in data:
            atom.velocity = np.array(data['velocity'])
        return atom
    
    def to_json(self, filepath: str):
        """
        Saves the atom's data to a JSON file.
        
        Args:
            filepath: The path to the output JSON file.
        """
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)

    @classmethod
    def from_json(cls, filepath: str) -> 'GAM_Atom':
        """
        Loads atom data from a JSON file and reconstructs the GAM_Atom object.
        
        Args:
            filepath: The path to the input JSON file.
            
        Returns:
            A new GAM_Atom object with the loaded data.
        """
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)
    

    #-----for BONDS
    def add_bond(self, bond: 'GAM_Bond'):
        """Adds a GAM_Bond object to the atom's list of bonds."""
        if bond not in self._bonds:
            self._bonds.append(bond)

    def remove_bond(self, bond: 'GAM_Bond'):
        """Removes a GAM_Bond object from the atom's list of bonds."""
        if bond in self._bonds:
            self._bonds.remove(bond)

    def get_bonds(self) -> List['GAM_Bond']:
        """Returns a list of all GAM_Bond objects connected to this atom."""
        return self._bonds

    def get_neighbor_atoms(self, atoms: Dict[int, 'GAM_Atom']) -> List['GAM_Atom']:
        """
        Returns a list of GAM_Atom objects directly connected via a bond.

        Args:
            atoms: Dictionary mapping atom IDs to Atom objects.
        """
        neighbors = []
        for bond in self._bonds:
            if bond.atom1_id == self.id:
                neighbors.append(atoms[bond.atom2_id])
            else:
                neighbors.append(atoms[bond.atom1_id])
        return neighbors
    


    def __repr__(self) -> str:
        """String representation maintaining original format with additional info."""
        return (f"Atom(id={self.id}, element='{self.element}', "
                f"x={self.x:.3f}, y={self.y:.3f}, z={self.z:.3f}, "
                f"charge={self.charge:.3f}, state={self.electronic_state.value})")





class GAM_Bond:
    """
    Advanced bond class for computational chemistry and molecular modeling applications.
    
    The GAM_Bond class represents covalent bonds between atoms within molecular systems
    and provides comprehensive functionality for structural analysis, geometric calculations,
    and bond property evaluation. It supports various bond orders, strain energy calculations,
    and advanced geometric transformations essential for molecular dynamics simulations
    and quantum chemistry computations.
    
    Attributes:
        id (int): Unique identifier for the bond within a molecular system.
        atom1_id (int): ID of the first atom in the bond.
        atom2_id (int): ID of the second atom in the bond.
        order (int): Bond order (1=single, 2=double, 3=triple).
        _length (Optional[float]): Cached bond length in Angstroms.
        _energy (Optional[float]): Cached bond energy in Joules.
        
    Examples:
        Basic bond creation:
        >>> bond = GAM_Bond(id=1, atom1_id=1, atom2_id=2, order=1)
        >>> print(f"Bond order: {bond.order}")
        Bond order: 1
        
        Bond length calculation:
        >>> atoms = {1: GAM_Atom(1, "C", 0.0, 0.0, 0.0), 
        ...           2: GAM_Atom(2, "C", 1.54, 0.0, 0.0)}
        >>> bond = GAM_Bond(1, 1, 2, 1)
        >>> length = bond.calculate_length(atoms)
        >>> print(f"Bond length: {length:.2f} Å")
        Bond length: 1.54 Å
        
        Bond validation:
        >>> is_reasonable = bond.is_reasonable_length(atoms)
        >>> print(f"Reasonable bond length: {is_reasonable}")
        Reasonable bond length: True
        
        Strain energy calculation:
        >>> strain_energy = bond.calculate_strain_energy(atoms)
        >>> print(f"Strain energy: {strain_energy:.3f} kcal/mol")
        Strain energy: 0.000 kcal/mol
        
    Notes:
        - Bond lengths are calculated in Angstroms.
        - Strain energies are calculated in kcal/mol using harmonic approximation.
        - Bond validation uses reference bond lengths with 20% tolerance.
        - All geometric calculations require access to the complete atom dictionary.
        
    Raises:
        ValueError: If fixed atom ID in rotation operations is not part of the bond.
        ValueError: If atoms associated with the bond do not exist in the provided dictionary.
        
    See Also:
        GAM_Atom: For individual atom representation and properties.
        GAM_Molecule: For complete molecular system management.
        REFERENCE_BOND_LENGTHS: For reference bond length data.
    """
    
    def __init__(self, id: int, atom1_id: int, atom2_id: int, order: int = 1):
        self.id = id
        self.atom1_id = atom1_id
        self.atom2_id = atom2_id
        self.order = order
        self._length: Optional[float] = None
        self._energy: Optional[float] = None

    @property
    def bond_vector(self, atoms: Dict[int, GAM_Atom]) -> np.ndarray:
        """Calculate the bond vector from atom1 to atom2."""
        atom1 = atoms[self.atom1_id]
        atom2 = atoms[self.atom2_id]
        return atom2.position - atom1.position

    def calculate_length(self, atoms: Dict[int, GAM_Atom]) -> float:
        """Calculate the bond length in Angstroms."""
        return np.linalg.norm(self.bond_vector(atoms))

    def is_reasonable_length(self, atoms: Dict[int, GAM_Atom]) -> bool:
        """
        Check if the bond length is reasonable based on reference data.
        
        Args:
            atoms: Dictionary mapping atom IDs to Atom objects
            
        Returns:
            bool: True if the bond length is within 20% of reference value
        """
        atom1 = atoms[self.atom1_id]
        atom2 = atoms[self.atom2_id]
        
        # Sort elements to match reference data format
        elements = tuple(sorted([atom1.element, atom2.element]))
        ref_lengths = REFERENCE_BOND_LENGTHS.get(elements, {})
        ref_length = ref_lengths.get(str(self.order))
        
        if ref_length is None:
            return True  # No reference data available
            
        actual_length = self.calculate_length(atoms)
        return abs(actual_length - ref_length) / ref_length < 0.2
    
    #maybe mroe reliable:
    def is_reasonable_length(self, atoms: Dict[int, GAM_Atom]) -> bool:
        """
        Check if the bond length is reasonable based on reference data.
        ...
        """
        atom1 = atoms[self.atom1_id]
        atom2 = atoms[self.atom2_id]

        elements = tuple(sorted([atom1.element, atom2.element]))
        ref_lengths = REFERENCE_BOND_LENGTHS.get(elements, {})
        # Corrected lookup to use integer key
        ref_length = ref_lengths.get(self.order)

        if ref_length is None:
            return True

        actual_length = self.calculate_length(atoms)
        return abs(actual_length - ref_length) / ref_length < 0.2




    def calculate_angle_with(self, other: 'GAM_Bond', atoms: Dict[int, GAM_Atom]) -> float:
        """
        Calculate angle between this bond and another bond.
        
        Args:
            other: Another Bond object
            atoms: Dictionary mapping atom IDs to Atom objects
            
        Returns:
            float: Angle in degrees between the bonds
        """
        v1 = self.bond_vector(atoms)
        v2 = other.bond_vector(atoms)
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))

    def rotate_around_axis(self, 
                         axis: np.ndarray, 
                         angle_degrees: float,
                         atoms: Dict[int, GAM_Atom],
                         fixed_atom_id: Optional[int] = None) -> None:
        """
        Rotate the bond around an axis while keeping one atom fixed.
        
        Args:
            axis: Rotation axis
            angle_degrees: Rotation angle in degrees
            atoms: Dictionary mapping atom IDs to Atom objects
            fixed_atom_id: ID of the atom to keep fixed (defaults to atom1_id)
        """
        if fixed_atom_id is None:
            fixed_atom_id = self.atom1_id

        if fixed_atom_id not in {self.atom1_id, self.atom2_id}:
            raise ValueError(f"Fixed atom ID must be either {self.atom1_id} or {self.atom2_id}.")
            
        moving_atom_id = self.atom2_id if fixed_atom_id == self.atom1_id else self.atom1_id
        
        # Get the fixed atom's position as rotation center
        center = atoms[fixed_atom_id].position
        
        # Rotate the moving atom
        atoms[moving_atom_id].rotate_around_axis(axis, angle_degrees, center)

    def calculate_strain_energy(self, atoms: Dict[int, GAM_Atom]) -> float:
        """
        Calculate the strain energy of the bond based on deviation from ideal length.
        
        Args:
            atoms: Dictionary mapping atom IDs to Atom objects
            
        Returns:
            float: Strain energy in kcal/mol
        """
        atom1 = atoms[self.atom1_id]
        atom2 = atoms[self.atom2_id]
        
        elements = tuple(sorted([atom1.element, atom2.element]))
        ref_lengths = REFERENCE_BOND_LENGTHS.get(elements, {})
        ref_length = ref_lengths.get(str(self.order))
        
        if ref_length is None:
            return 0.0
            
        actual_length = self.calculate_length(atoms)
        # Simple harmonic approximation
        force_constant = 100.0  # kcal/mol/Å² (approximate)
        return 0.5 * force_constant * (actual_length - ref_length) ** 2
    
    def get_elements(self, atoms: Dict[int, 'GAM_Atom']) -> Tuple[str, str]:
        """
        Returns a tuple of the elements involved in the bond.

        Args:
            atoms: Dictionary mapping atom IDs to Atom objects.

        Returns:
            A sorted tuple of the elements (e.g., ('C', 'H')).
        """
        atom1 = atoms[self.atom1_id]
        atom2 = atoms[self.atom2_id]
        return tuple(sorted([atom1.element, atom2.element]))
    
    def get_bond_type(self, atoms: Dict[int, 'GAM_Atom']) -> str:
        """
        Returns a descriptive string for the bond type.

        Args:
            atoms: Dictionary mapping atom IDs to Atom objects.

        Returns:
            A string describing the bond (e.g., 'C-C single bond').
        """
        elements = self.get_elements(atoms)
        bond_order_map = {1: 'single', 2: 'double', 3: 'triple'}
        bond_type = bond_order_map.get(self.order, 'unspecified')
        return f"{elements[0]}-{elements[1]} {bond_type} bond"
    

    
    def is_aromatic(self, atoms: Dict[int, 'GAM_Atom']) -> bool:
        """
        Checks if the bond is likely part of an aromatic system.
        This is a simple heuristic based on bond order and element types.

        Args:
            atoms: Dictionary mapping atom IDs to Atom objects.

        Returns:
            True if the bond is aromatic, False otherwise.
        """
        if self.order == 1:
            elements = self.get_elements(atoms)
            # Aromatic bonds are often C-C bonds with an order of 1.5,
            # but in some representations, they are labeled as order 1.
            if elements == ('C', 'C'):
                # Check for a more subtle indicator of aromaticity, such as sp2 hybridization
                atom1 = atoms[self.atom1_id]
                atom2 = atoms[self.atom2_id]
                # A common heuristic is checking if both atoms have 3 neighbors
                if len(atom1._bonds) == 3 and len(atom2._bonds) == 3:
                    return True
        return False
    

    def calculate_dihedral_angle(self, atom3_id: int, atom4_id: int, atoms: Dict[int, 'GAM_Atom']) -> float:
        """
        Calculates the dihedral angle between this bond and a plane defined by
        the third and fourth atoms.

        Args:
            atom3_id: The ID of the third atom.
            atom4_id: The ID of the fourth atom.
            atoms: Dictionary mapping atom IDs to Atom objects.

        Returns:
            The dihedral angle in degrees.
        """
        atom1 = atoms[self.atom1_id]
        atom2 = atoms[self.atom2_id]
        atom3 = atoms[atom3_id]
        atom4 = atoms[atom4_id]

        b1 = atom1.position - atom2.position
        b2 = atom3.position - atom2.position
        b3 = atom4.position - atom3.position

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        m1 = np.cross(n1, b2 / np.linalg.norm(b2))

        x = np.dot(n1, n2)
        y = np.dot(m1, n2)

        return np.degrees(np.arctan2(y, x))
    
    def calculate_angle_with_plane(self, plane_normal: np.ndarray, atoms: Dict[int, 'GAM_Atom']) -> float:
        """
        Calculates the angle between the bond vector and a user-defined plane.

        Args:
            plane_normal: The normal vector of the plane.
            atoms: Dictionary mapping atom IDs to Atom objects.

        Returns:
            The angle in degrees.
        """
        bond_vector = self.bond_vector(atoms)
        normal_vec = plane_normal / np.linalg.norm(plane_normal)
        
        # The angle between the vector and the normal of the plane
        cos_angle_to_normal = np.dot(bond_vector, normal_vec) / (np.linalg.norm(bond_vector))
        
        # The angle with the plane is 90 degrees minus the angle with the normal
        angle_with_plane_rad = np.pi / 2 - np.arccos(np.clip(cos_angle_to_normal, -1.0, 1.0))
        
        return np.degrees(angle_with_plane_rad)
    


    def to_dict(self) -> Dict[str, Any]:
        """
        Serializes the bond's data into a dictionary.
        
        Returns:
            A dictionary representation of the bond.
        """
        return {
            'id': self.id,
            'atom1_id': self.atom1_id,
            'atom2_id': self.atom2_id,
            'order': self.order
        }
    
    def get_atoms(self, atoms: Dict[int, 'GAM_Atom']) -> Tuple['GAM_Atom', 'GAM_Atom']:
        """
        Retrieves the two GAM_Atom objects connected by this bond.

        Args:
            atoms: Dictionary mapping atom IDs to Atom objects.

        Returns:
            A tuple containing the two GAM_Atom objects.
        """
        atom1 = atoms.get(self.atom1_id)
        atom2 = atoms.get(self.atom2_id)
        if not all([atom1, atom2]):
            raise ValueError("One or more atoms associated with this bond do not exist.")
        return atom1, atom2

    def get_neighboring_bonds(self, atoms: Dict[int, 'GAM_Atom']) -> List['GAM_Bond']:
        """
        Finds all other bonds that share an atom with this bond.

        Args:
            atoms: Dictionary mapping atom IDs to Atom objects.

        Returns:
            A list of GAM_Bond objects that are neighbors to this bond.
        """
        atom1, atom2 = self.get_atoms(atoms)
        neighboring_bonds = []
        for bond in atom1.get_bonds():
            if bond.id != self.id:
                neighboring_bonds.append(bond)
        for bond in atom2.get_bonds():
            if bond.id != self.id and bond not in neighboring_bonds:
                neighboring_bonds.append(bond)
        return neighboring_bonds

    

    def is_part_of_ring(self, atoms: Dict[int, 'GAM_Atom'], max_ring_size: int = 7) -> bool:
        """
        Checks if the bond is part of a closed ring structure using Breadth-First Search (BFS).
        ...
        """
        atom1, atom2 = self.get_atoms(atoms)
    
        # Use a more efficient deque for the queue
        queue = deque([(atom2, 1, {atom2.id})])  # (current_atom, path_length, visited_set)
    
        while queue:
            current_atom, path_length, visited = queue.popleft()
    
            if path_length >= max_ring_size:
                continue
    
            for bond in current_atom.get_bonds():
                if bond.id == self.id:
                    continue
    
                neighbor = atoms[bond.atom2_id if bond.atom1_id == current_atom.id else bond.atom1_id]
    
                if neighbor.id == atom1.id:
                    return True  # Found a path back to atom1, forming a ring
    
                if neighbor.id not in visited:
                    new_visited = visited.copy()
                    new_visited.add(neighbor.id)
                    queue.append((neighbor, path_length + 1, new_visited))
        return False
    


    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'GAM_Bond':
        """
        Reconstructs a GAM_Bond object from a dictionary.

        Args:
            data: A dictionary containing the bond's data.

        Returns:
            A new GAM_Bond object.
        """
        return cls(
            id=data['id'],
            atom1_id=data['atom1_id'],
            atom2_id=data['atom2_id'],
            order=data['order']
        )
    
    def to_json(self, filepath: str):
        """
        Saves the bond's data to a JSON file.
        
        Args:
            filepath: The path to the output JSON file.
        """
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)  
    


    def __repr__(self) -> str:
        return f"Bond(id={self.id}, atom1_id={self.atom1_id}, atom2_id={self.atom2_id}, order={self.order})"
    
    

    
 





    
