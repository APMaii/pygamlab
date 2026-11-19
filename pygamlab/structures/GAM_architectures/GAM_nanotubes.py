from ..Primatom import GAM_Atom , GAM_Bond , GAM_Molecule
#from gam_atom import GAM_Atom , GAM_Bond
import numpy as np
from math import gcd, sqrt, sin, cos, pi, floor, atan, atan2
from typing import Tuple, List, Dict, Optional, Union, Any
from dataclasses import dataclass
from enum import Enum
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
from ase import Atoms



class NanotubeType(Enum):
    """Nanotube chirality classification"""
    ZIGZAG = "zigzag"
    ARMCHAIR = "armchair" 
    CHIRAL = "chiral"

@dataclass
class AtomProperties:
    """Properties for different atom types"""
    symbol: str
    radius: float  # Atomic radius in Angstroms
    color: Tuple[int, int, int]  # RGB color for visualization
    mass: float  # Atomic mass in amu
    bond_length: float  # Typical bond length in Angstroms








class Nanotube_Generator:
    """
    Advanced chiral nanotube structure generator with extensive customization options.

    The `Nanotube_Generator` class provides a robust and modular framework for constructing
    single-walled (SWCNT), double-walled (DWCNT), or multi-walled (MWCNT) nanotube structures
    with arbitrary chirality, composition, and geometry. It is designed for computational
    materials science applications such as density functional theory (DFT), molecular dynamics (MD),
    and machine learning (ML) datasets requiring atomically resolved nanotube structures.

    This generator supports a wide range of materials beyond carbon, including silicon, 
    germanium, boron nitride (BN), silicon carbide (SiC), transition metal dichalcogenides (MoS₂, WS₂),
    and III–V compounds such as GaN and InN. It enables detailed control over structural,
    chemical, and defect configurations to simulate realistic experimental conditions.

    Key Features
    ------------
    - Generation of chiral (n, m) nanotubes with precise geometric parameters
    - Calculation of diameter, chiral angle, translation vector, and lattice constants
    - Support for multi-walled nanotubes via concentric wall generation
    - Introduction of various structural defects (vacancies, Stone–Wales, substitutional)
    - Doping and alloying with customizable dopant types and concentrations
    - Application of uniaxial, compressive, or radial strain
    - Surface functionalization with chemical groups (e.g., –OH, –COOH)
    - Export to XYZ and ASE-compatible formats

    Parameters
    ----------
    n, m : int
        Chiral indices defining the nanotube (n, m). Determines the nanotube’s 
        chirality, electronic properties, and geometry.
    length : float, default=10.0
        Length of the nanotube in Ångströms.
    atom_type : str, default='C'
        Atomic species of the nanotube (e.g., 'C', 'BN', 'MoS2').
    custom_atom : AtomProperties, optional
        Custom atom definition, if not included in predefined atom data.
    defects : dict, optional
        Defect configuration. Example:
        ``{'vacancies': 0.02, 'stone_wales': 0.01, 'substitutions': {'atom': 'B', 'ratio': 0.03}}``
    doping : dict, optional
        Doping configuration. Example:
        ``{'dopant': 'N', 'concentration': 0.02, 'pattern': 'periodic'}``
    strain : dict, optional
        Strain parameters defining deformation mode and magnitude. Example:
        ``{'type': 'tensile', 'magnitude': 0.05}``
    multi_wall : list of tuple(int, int), optional
        List of chiral indices defining additional walls for multi-walled nanotubes.
        Example: ``[(10,10), (15,15), (20,20)]``
    functionalization : dict, optional
        Parameters for surface functionalization.
        Example: ``{'groups': ['COOH', 'OH'], 'density': 0.1}``
    verbose : bool, default=False
        Print detailed generation progress and statistics.

    Attributes
    ----------
    atoms : list of GAM_Atom
        List of atomic objects representing all atoms in the nanotube.
    bonds : list of GAM_Bond
        List of interatomic bonds with computed lengths and connectivity.
    meta : dict
        Metadata containing nanotube properties, lattice parameters, and configuration.
    last_generated : dict
        Cached data of the most recently generated nanotube.
    generation_stats : dict
        Summary statistics such as atom count, bond lengths, and composition.
    origin : tuple of float
        Origin of the nanotube (x, y) coordinates.
    rotation : float
        Cumulative rotation angle applied to the nanotube (in degrees).

    Methods
    -------
    classify_nanotube(n, m) -> NanotubeType
        Classify the nanotube as zigzag, armchair, or chiral.
    calculate_properties(n, m, bond_length) -> Dict[str, float]
        Compute geometric properties including diameter, chiral angle, and unit cell parameters.
    _generate_base_structure(n, m, length, atom_props, sign) -> Tuple[List, List]
        Construct the base atomic structure and bonds for a single nanotube wall.
    _generate_bonds(atoms, max_distance) -> List[Dict]
        Identify and create bonds between atoms based on cutoff distance.
    _apply_strain(atoms, strain_params) -> List[Dict]
        Apply tensile, compressive, or radial strain transformations.
    _introduce_defects(atoms, bonds, defect_params, atom_props) -> Tuple[List, List]
        Introduce vacancy, substitutional, or Stone–Wales defects into the structure.
    _apply_doping(atoms, doping_params, atom_props) -> List[Dict]
        Substitute atoms with dopants according to given pattern and concentration.
    _create_multi_wall(atoms, bonds, wall_specs, length, base_atom_props) -> Tuple[List, List]
        Build multi-walled nanotube structures with appropriate inter-wall spacing.
    _add_functional_groups(atoms, bonds, functionalization) -> Tuple[List, List]
        Add surface functional groups to outer atoms.
    _calculate_stats(nanotube_data) -> Dict
        Calculate statistical summaries of the generated structure.
    _print_stats() -> None
        Display generation statistics in human-readable format.
    translate(dx, dy, dz=0.0) -> None
        Translate the nanotube by a specified displacement vector.
    rotate(angle_deg, about_center=True) -> None
        Rotate the nanotube around the z-axis by a given angle.
    copy() -> List[GAM_Atom]
        Return a deep copy of the nanotube’s atomic configuration.
    get_atoms() -> List[GAM_Atom]
        Return a list of atom objects representing the nanotube.
    get_positions() -> List[Tuple[float, float, float]]
        Retrieve all atomic positions.
    get_elements() -> List[str]
        Retrieve the list of atomic element symbols.
    to_xyz(filename, nanotube_data=None) -> None
        Export the nanotube structure to XYZ format.
    to_ase() -> ase.Atoms
        Convert the structure to an ASE `Atoms` object for simulation workflows.
    get_supported_atoms() -> List[str]
        Return all supported predefined atomic types.

    Notes
    -----
    - Chirality determines the nanotube’s electronic properties: armchair nanotubes (n = m)
      are metallic, zigzag (m = 0) can be semiconducting or metallic depending on n, and
      chiral tubes exhibit varying electronic characteristics.
    - Defects and doping can significantly modify mechanical, optical, and electrical behavior.
    - Strain and functionalization options are useful for simulating real-world experimental conditions.
    - The geometry generation algorithm is based on the unrolled graphene lattice method
      with periodic closure along the chiral vector.
    - Interlayer spacing for multi-walled nanotubes defaults to 3.4 Å (graphitic spacing).

    Examples
    --------
    >>> # Generate a pristine (10,10) armchair carbon nanotube
    >>> from pygamlab import Nanotube_Generator
    >>> cnt = Nanotube_Generator(n=10, m=10, length=20.0, atom_type='C', verbose=True)

    >>> # Generate a nitrogen-doped (12,0) zigzag BN nanotube
    >>> doped_bn = Nanotube_Generator(n=12, m=0, atom_type='BN',
    ...     doping={'dopant': 'N', 'concentration': 0.03, 'pattern': 'random'})

    >>> # Create a triple-walled CNT
    >>> mwcnt = Nanotube_Generator(n=5, m=5, multi_wall=[(5,5), (10,10), (15,15)], length=30.0)
    >>> mwcnt.to_xyz("multiwalled_CNT.xyz")
    """
    
    # Predefined atom properties
    ATOM_DATA = {
        'C': AtomProperties('C', 0.70, (64, 64, 64), 12.011, 1.42),
        'Si': AtomProperties('Si', 1.11, (240, 200, 160), 28.085, 2.35),
        'Ge': AtomProperties('Ge', 1.20, (102, 143, 143), 72.630, 2.44),
        'BN': AtomProperties('BN', 0.83, (255, 181, 181), 12.51, 1.45),  # Boron Nitride
        'SiC': AtomProperties('SiC', 0.91, (255, 215, 0), 20.048, 1.89),
        'GaN': AtomProperties('GaN', 0.87, (138, 43, 226), 41.865, 1.94),
        'InN': AtomProperties('InN', 0.92, (166, 166, 171), 128.825, 2.15),
        'AlN': AtomProperties('AlN', 0.85, (191, 166, 166), 20.491, 1.89),
        'WS2': AtomProperties('WS2', 1.35, (255, 215, 0), 247.97, 2.41),  # Tungsten disulfide
        'MoS2': AtomProperties('MoS2', 1.30, (138, 43, 226), 160.07, 2.38), # Molybdenum disulfide
    }
    
    
    def __init__(self,
                         n: int,
                         m: int,
                         length: float = 10.0,
                         atom_type: str = 'C',
                         custom_atom: Optional[AtomProperties] = None,
                         defects: Optional[Dict] = None,
                         doping: Optional[Dict] = None,
                         strain: Optional[Dict] = None,
                         multi_wall: Optional[List[Tuple[int, int]]] = None,
                         functionalization: Optional[Dict] = None,
                         verbose: bool = False) -> Dict:
        """
        Generate advanced nanotube structure with multiple customization options
        
        Parameters:
        -----------
        n, m : int
            Chiral indices (n,m)
        length : float
            Nanotube length in Angstroms
        atom_type : str
            Type of atoms ('C', 'Si', 'BN', etc.)
        custom_atom : AtomProperties, optional
            Custom atom properties if not in predefined list
        defects : dict, optional
            Defect parameters: {'vacancies': ratio, 'stone_wales': ratio, 'substitutions': {...}}
        doping : dict, optional
            Doping parameters: {'dopant': 'B', 'concentration': 0.01, 'pattern': 'random'}
        strain : dict, optional
            Strain parameters: {'type': 'tensile/compressive', 'magnitude': 0.05}
        
        multi_wall : List[Tuple[int, int]], optional
            List of (n,m) for multi-walled nanotubes: [(10,10), (15,15), (20,20)]
        functionalization : dict, optional
            Surface functionalization: {'groups': ['COOH', 'OH'], 'density': 0.1}
        verbose : bool
            Print detailed information
            
        Returns:
        --------
        dict: Complete nanotube data structure for pygamlab
        """
        
        self.last_generated = None
        self.generation_stats = {}
        
        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []
        
        # Initialize transformation attributes
        self.origin = (0.0, 0.0)
        self.rotation = 0.0
        
        
        # Get atom properties
        if custom_atom:
            atom_props = custom_atom
        else:
            atom_props = self.ATOM_DATA.get(atom_type, self.ATOM_DATA['C'])
        
        # Ensure n >= m
        if n < m:
            n, m = m, n
            sign = -1
        else:
            sign = 1
            
        # Calculate properties
        properties = self.calculate_properties(n, m, atom_props.bond_length)
        nanotube_type = self.classify_nanotube(n, m)
        
        if verbose:
            print(f"Generating {nanotube_type.value} nanotube ({n},{m})")
            print(f"Diameter: {properties['diameter']:.2f} Å")
            print(f"Chiral angle: {properties['chiral_angle']:.2f}°")
        
        # Generate base structure
        self.atoms, self.bonds = self._generate_base_structure(n, m, length, atom_props, sign)
        
        
        
        # Apply modifications
        if strain:
            self.atoms = self._apply_strain(self.atoms, strain)
            
        if defects:
            self.atoms, self.bonds = self._introduce_defects(self.atoms, self.bonds, defects, atom_props)
            
        if doping:
            print('doping')
            self.atoms = self._apply_doping(self.atoms, doping, atom_props)

        if functionalization:
            self.atoms, self.bonds = self._add_functional_groups(self.atoms, self.bonds, functionalization)
            
        if multi_wall:
            self.atoms, self.bonds = self._create_multi_wall(self.atoms, self.bonds, multi_wall, length, atom_props)
        
        
        
        dict_atoms=self.atoms
        
        self.atoms=[ GAM_Atom(
                id=atom['id'],
                element=atom['symbol'],
                x=atom['position'][0],
                y=atom['position'][1],
                z=atom['position'][2],
            )
            for atom in dict_atoms
        ]
        
        
        
        # Generate final structure
        self.meta = {
            'atoms': self.atoms,
            'dict_atoms':dict_atoms,
            'bonds': self.bonds,
            'properties': properties,
            'nanotube_type': nanotube_type.value,
            'indices': (n, m),
            'atom_type': atom_type,
            'atom_properties': atom_props,
            'length': length,
            'num_atoms': len(self.atoms),
            'metadata': {
                'defects': defects,
                'doping': doping,
                'strain': strain,
                'multi_wall': multi_wall,
                'functionalization': functionalization
            }
        }
        
        self.last_generated = self.meta
        self.generation_stats = self._calculate_stats(self.meta)
        

        
        if verbose:
            self._print_stats()
            
        #return nanotube_data
    
    
    
    
    def classify_nanotube(self, n: int, m: int) -> NanotubeType:
        """Classify nanotube type based on (n,m) indices"""
        if m == 0:
            return NanotubeType.ZIGZAG
        elif n == m:
            return NanotubeType.ARMCHAIR
        else:
            return NanotubeType.CHIRAL
    
    def calculate_properties(self, n: int, m: int, bond_length: float) -> Dict[str, float]:
        """Calculate geometric properties of the nanotube"""
        # Ensure n >= m for consistency
        if n < m:
            n, m = m, n
            
        sq3 = sqrt(3.0)
        a = sq3 * bond_length  # Lattice parameter
        
        # Chiral vector length
        l_chiral = sqrt(n*n + m*m + n*m)
        
        # Diameter and circumference
        diameter = a * l_chiral / pi
        circumference = pi * diameter
        
        # Chiral angle in degrees
        if n == m:
            chiral_angle = 30.0  # Armchair
        elif m == 0:
            chiral_angle = 0.0   # Zigzag
        else:
            chiral_angle = atan(sqrt(3) * m / (2*n + m)) * 180.0 / pi
            
        # Translation vector length
        d_gcd = gcd(n, m)
        if (n - m) % (3 * d_gcd) == 0:
            d_r = 3 * d_gcd
        else:
            d_r = d_gcd
            
        translation_length = sq3 * a * sqrt(n*n + m*m + n*m) / d_r
        
        # Number of atoms per unit cell
        atoms_per_cell = 2 * (n*n + m*m + n*m) / d_gcd
        
        return {
            'diameter': diameter,
            'circumference': circumference,
            'chiral_angle': chiral_angle,
            'translation_length': translation_length,
            'atoms_per_cell': int(atoms_per_cell),
            'lattice_parameter': a,
            'chiral_vector_length': l_chiral * a
        }
    
    
    
    def _generate_base_structure(self, n, m, length, atom_props, sign):
        """Generate the base nanotube atomic structure"""
        sq3 = sqrt(3.0)
        a = sq3 * atom_props.bond_length
        l2 = n * n + m * m + n * m
        l1 = sqrt(l2)
        
        nd = gcd(n, m)
        if (n - m) % (3 * nd) == 0:
            ndr = 3 * nd
        else:
            ndr = nd
            
        nr = (2 * m + n) // ndr
        ns = -(2 * n + m) // ndr
        nn = 2 * l2 // ndr
        
        # Find translation vector
        ichk = 0
        if nr == 0:
            n60 = 1
        else:
            n60 = nr * 4
            
        absn = abs(n60)
        nnp = []
        nnq = []
        for i in range(-absn, absn + 1):
            for j in range(-absn, absn + 1):
                j2 = nr * j - ns * i
                if j2 == 1:
                    j1 = m * i - n * j
                    if j1 > 0 and j1 < nn:
                        ichk += 1
                        nnp.append(i)
                        nnq.append(j)
                        
        if ichk == 0:
            raise RuntimeError('Translation vector not found!')
        if ichk >= 2:
            raise RuntimeError('Multiple translation vectors found!')
            
        nnnp, nnnq = nnp[0], nnq[0]
        
        lp = nnnp * nnnp + nnnq * nnnq + nnnp * nnnq
        r = a * sqrt(lp)
        c = a * l1
        t = sq3 * c / ndr
        rs = c / (2.0 * pi)
        
        # Angular parameters
        q1 = atan((sq3 * m) / (2 * n + m))
        q2 = atan((sq3 * nnnq) / (2 * nnnp + nnnq))
        q3 = q1 - q2
        q4 = 2.0 * pi / nn
        q5 = atom_props.bond_length * cos((pi / 6.0) - q1) / c * 2.0 * pi
        
        h1 = abs(t) / abs(sin(q3))
        h2 = atom_props.bond_length * sin((pi / 6.0) - q1)
        
        # Generate atomic positions
        atoms = []
        atom_id = 0
        
        for i in range(nn):
            k = floor(i * abs(r) / h1)
            
            # First atom
            x1 = rs * cos(i * q4)
            y1 = rs * sin(i * q4)
            z1 = (i * abs(r) - k * h1) * sin(q3)
            
            kk2 = abs(floor((z1 + 0.0001) / t))
            if z1 >= t - 0.0001:
                z1 -= t * kk2
            elif z1 < 0:
                z1 += t * kk2
                
            atoms.append({
                'id': atom_id,
                'symbol': atom_props.symbol,
                'position': [x1, y1, sign * z1],
                'color': atom_props.color,
                'radius': atom_props.radius,
                'mass': atom_props.mass
            })
            atom_id += 1
            
            # Second atom
            z3 = (i * abs(r) - k * h1) * sin(q3) - h2
            x2 = rs * cos(i * q4 + q5)
            y2 = rs * sin(i * q4 + q5)
            
            if z3 >= 0 and z3 < t:
                z2 = z3
            else:
                z2 = (i * abs(r) - (k + 1) * h1) * sin(q3) - h2
                kk = abs(floor(z2 / t))
                if z2 >= t - 0.0001:
                    z2 -= t * kk
                elif z2 < 0:
                    z2 += t * kk
                    
            atoms.append({
                'id': atom_id,
                'symbol': atom_props.symbol,
                'position': [x2, y2, sign * z2],
                'color': atom_props.color,
                'radius': atom_props.radius,
                'mass': atom_props.mass
            })
            atom_id += 1
            
        
        # Replicate for desired length
        num_cells = max(1, int(length / t))
        if num_cells > 1:
            base_atoms = atoms[:]
            for cell in range(1, num_cells):
                for base_atom in base_atoms:
                    new_atom = base_atom.copy()
                    new_atom['id'] = atom_id
                    new_pos = new_atom['position'].copy()
                    new_pos[2] += cell * t * sign
                    new_atom['position'] = new_pos
                    atoms.append(new_atom)
                    atom_id += 1
                    
                    
        
        
        # Generate bonds
        bonds = self._generate_bonds(atoms, atom_props.bond_length * 1.2)  # 20% tolerance
        
        
        return atoms, bonds
    

    
    def _generate_bonds(self, atoms, max_distance):
        """Generate bonds between atoms based on distance"""
        bonds = []
        bond_id = 0
        
        for i, atom1 in enumerate(atoms):
            for j, atom2 in enumerate(atoms[i+1:], i+1):
                pos1 = np.array(atom1['position'])
                pos2 = np.array(atom2['position'])
                distance = np.linalg.norm(pos1 - pos2)
                
                if distance <= max_distance:
                    bonds.append({
                        'id': bond_id,
                        'atom1_id': atom1['id'],
                        'atom2_id': atom2['id'],
                        'length': distance,
                        'type': 'single'
                    })
                    bond_id += 1
                    
        return bonds
    
    def _apply_strain(self, atoms, strain_params):
        """Apply strain to the nanotube"""
        strain_type = strain_params.get('type', 'tensile')
        magnitude = strain_params.get('magnitude', 0.0)
        
        if strain_type == 'tensile':
            # Stretch along z-axis
            for atom in atoms:
                atom['position'][2] *= (1 + magnitude)
        elif strain_type == 'compressive':
            # Compress along z-axis
            for atom in atoms:
                atom['position'][2] *= (1 - magnitude)
        elif strain_type == 'radial':
            # Radial expansion/compression
            for atom in atoms:
                x, y = atom['position'][0], atom['position'][1]
                r = sqrt(x*x + y*y)
                if r > 0:
                    factor = (1 + magnitude)
                    atom['position'][0] = x * factor
                    atom['position'][1] = y * factor
                    
        return atoms
    
    def _introduce_defects(self, atoms, bonds, defect_params, atom_props):
        """Introduce various defects into the structure"""
        modified_atoms = atoms[:]
        modified_bonds = bonds[:]
        
        # Vacancy defects
        if 'vacancies' in defect_params:
            vacancy_ratio = defect_params['vacancies']
            num_vacancies = int(len(atoms) * vacancy_ratio)
            vacancy_indices = np.random.choice(len(atoms), num_vacancies, replace=False)
            
            # Remove atoms and associated bonds
            modified_atoms = [atom for i, atom in enumerate(atoms) if i not in vacancy_indices]
            vacancy_atom_ids = {atoms[i]['id'] for i in vacancy_indices}
            modified_bonds = [bond for bond in bonds 
                            if bond['atom1_id'] not in vacancy_atom_ids 
                            and bond['atom2_id'] not in vacancy_atom_ids]
        
        # Stone-Wales defects (bond rotations)
        if 'stone_wales' in defect_params:
            sw_ratio = defect_params['stone_wales']
            # Implementation would involve identifying adjacent hexagons and rotating bonds
            # This is complex and would require topological analysis
            pass
        
        # Substitutional defects
        if 'substitutions' in defect_params:
            sub_params = defect_params['substitutions']
            sub_atom_type = sub_params.get('atom', 'B')
            sub_ratio = sub_params.get('ratio', 0.01)
            
            if sub_atom_type in self.ATOM_DATA:
                sub_props = self.ATOM_DATA[sub_atom_type]
                num_substitutions = int(len(modified_atoms) * sub_ratio)
                sub_indices = np.random.choice(len(modified_atoms), num_substitutions, replace=False)
                
                for idx in sub_indices:
                    modified_atoms[idx]['symbol'] = sub_props.symbol
                    modified_atoms[idx]['color'] = sub_props.color
                    modified_atoms[idx]['radius'] = sub_props.radius
                    modified_atoms[idx]['mass'] = sub_props.mass
        
        return modified_atoms, modified_bonds
    
    def _apply_doping(self, atoms, doping_params, atom_props):
        """Apply doping to the nanotube"""
        dopant = doping_params.get('dopant', 'B')
        concentration = doping_params.get('concentration', 0.01)
        pattern = doping_params.get('pattern', 'random')
        

        
        if dopant not in self.ATOM_DATA:
            return atoms
            
        dopant_props = self.ATOM_DATA[dopant]
        num_dopants = int(len(atoms) * concentration)
        
        if pattern == 'random':
            doping_indices = np.random.choice(len(atoms), num_dopants, replace=False)
        elif pattern == 'periodic':
            step = len(atoms) // num_dopants
            doping_indices = list(range(0, len(atoms), step))[:num_dopants]
        else:
            doping_indices = np.random.choice(len(atoms), num_dopants, replace=False)
        
        modified_atoms = atoms[:]
        for idx in doping_indices:
            modified_atoms[idx]['symbol'] = dopant_props.symbol
            modified_atoms[idx]['color'] = dopant_props.color
            modified_atoms[idx]['radius'] = dopant_props.radius
            modified_atoms[idx]['mass'] = dopant_props.mass
            
        return modified_atoms
    

    
    
    
    def _create_multi_wall(self, atoms, bonds, wall_specs, length, base_atom_props):
        """Create multi-walled nanotube"""
        all_atoms = atoms[:]
        all_bonds = bonds[:]
        
        # Start with the innermost wall (already generated)
        current_atom_count = len(atoms)
        current_bond_count = len(bonds)
        
        for i, (n, m) in enumerate(wall_specs):
            # Generate additional wall with appropriate spacing
            wall_atoms, wall_bonds = self._generate_base_structure(
                n, m, length, base_atom_props, 1
            )
            
            # Adjust radial position for proper interlayer spacing (typically ~3.4 Å)
            spacing = 3.4 * (i + 1)
            for atom in wall_atoms:
                x, y = atom['position'][0], atom['position'][1]
                r = sqrt(x*x + y*y)
                if r > 0:
                    factor = (r + spacing) / r
                    atom['position'][0] = x * factor
                    atom['position'][1] = y * factor
                
                # Update atom ID to avoid conflicts
                atom['id'] = current_atom_count + atom['id']
            
            # Update bond atom IDs
            for bond in wall_bonds:
                bond['id'] = current_bond_count + bond['id']
                bond['atom1_id'] += current_atom_count
                bond['atom2_id'] += current_atom_count
            
            all_atoms.extend(wall_atoms)
            all_bonds.extend(wall_bonds)
            
            # Update counters for next iteration
            current_atom_count += len(wall_atoms)
            current_bond_count += len(wall_bonds)
        
        return all_atoms, all_bonds
    
    def _add_functional_groups(self, atoms, bonds, functionalization):
        """Add functional groups to the nanotube surface"""
        # This is a placeholder implementation
        # In a real implementation, you would add specific functional groups
        # like -COOH, -OH, -NH2, etc. to surface atoms
        
        groups = functionalization.get('groups', [])
        density = functionalization.get('density', 0.1)
        
        if not groups or density <= 0:
            return atoms, bonds
            
        # For now, just return the original structure
        # A full implementation would involve:
        # 1. Identifying surface atoms
        # 2. Adding functional group atoms
        # 3. Creating new bonds
        # 4. Updating atom IDs and positions
        
        return atoms, bonds
    
    def _calculate_stats(self, nanotube_data):
        """Calculate statistics for the generated nanotube"""
        atoms = nanotube_data['dict_atoms']
        bonds = nanotube_data['bonds']
        
        # Atom type distribution
        atom_counts = {}
        for atom in atoms:
            symbol = atom['symbol']
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        # Bond length statistics
        bond_lengths = [bond['length'] for bond in bonds]
        
        stats = {
            'total_atoms': len(atoms),
            'total_bonds': len(bonds),
            'atom_distribution': atom_counts,
            'average_bond_length': np.mean(bond_lengths) if bond_lengths else 0,
            'min_bond_length': np.min(bond_lengths) if bond_lengths else 0,
            'max_bond_length': np.max(bond_lengths) if bond_lengths else 0,
            'bond_length_std': np.std(bond_lengths) if bond_lengths else 0
        }
        
        return stats
    
    def _print_stats(self):
        """Print generation statistics"""
        if not self.generation_stats:
            return
            
        stats = self.generation_stats
        print("\n=== Nanotube Generation Statistics ===")
        print(f"Total atoms: {stats['total_atoms']}")
        print(f"Total bonds: {stats['total_bonds']}")
        print(f"Atom distribution: {stats['atom_distribution']}")
        print(f"Average bond length: {stats['average_bond_length']:.3f} Å")
        print(f"Bond length range: {stats['min_bond_length']:.3f} - {stats['max_bond_length']:.3f} Å")
        print("=" * 40)
        
        
        
        
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


    
    def to_xyz(self, filename: str, nanotube_data: Optional[Dict] = None):
        """Export nanotube to XYZ format"""
        if nanotube_data is None:
            nanotube_data = self.last_generated
            
        if not nanotube_data:
            raise ValueError("No nanotube data to export")
            
        atoms = nanotube_data['atoms']
        
        with open(filename, 'w') as f:
            f.write(f"{len(atoms)}\n")
            #f.write(f"Nanotube ({nanotube_data['indices'][0]},{nanotube_data['indices'][1]}) - {nanotube_data['atom_type']}\n")
            f.write('Generated ')

            for atom in atoms:
                x, y, z = atom['position']
                f.write(f"{atom['symbol']} {x:.6f} {y:.6f} {z:.6f}\n")
                
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
    
    def get_supported_atoms(self) -> List[str]:
        """Get list of supported atom types"""
        return list(self.ATOM_DATA.keys())
        



