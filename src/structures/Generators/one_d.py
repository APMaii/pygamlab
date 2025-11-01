from ..Primatom import GAM_Atom , GAM_Bond

from math import gcd, sqrt
import numpy as np
from ase.atoms import Atoms
import copy
from ase.visualize import view   # optional, to visualize
from ase.io import write 
#from one_d import ase_nanotube
from typing import List, Tuple, Dict, Optional, Union , Any
import numpy as np
from ase.build import bulk
from ase import Atoms
from ase.atom import Atom
from ase.geometry import cell_to_cellpar
from ase.visualize import view
from ase.io import write, read
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io.trajectory import Trajectory
import random



from ase.atoms import Atoms
import copy

try:
    import ase
    from ase import Atoms
    HAS_ASE = True
except ImportError:
    HAS_ASE = False




def ase_nanotube(n, m, length=1.0, bond=1.42, symbol='C', verbose=False,
             vacuum=None):
    """Create an atomic structure.

    Creates a single-walled nanotube whose structure is specified using the
    standardized (n, m) notation.

    Parameters
    ----------
    n : int
        n in the (n, m) notation.
    m : int
        m in the (n, m) notation.
    length : float, optional
        Length of the nanotube in Angstroms.
    bond : float, optional
        Bond length between neighboring atoms.
    symbol : str, optional
        Chemical element to construct the nanotube from.
    verbose : bool, optional
        If True, will display key geometric parameters.

    Returns
    -------
    ase.atoms.Atoms
        An ASE Atoms object corresponding to the specified molecule.

    Examples
    --------
    >>> from ase.build import nanotube
    >>> atoms1 = nanotube(6, 0, length=4)
    >>> atoms2 = nanotube(3, 3, length=6, bond=1.4, symbol='Si')
    """
    if n < m:
        m, n = n, m
        sign = -1
    else:
        sign = 1

    nk = 6000
    sq3 = sqrt(3.0)
    a = sq3 * bond
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
        raise RuntimeError('not found p, q strange!!')
    if ichk >= 2:
        raise RuntimeError('more than 1 pair p, q strange!!')

    nnnp = nnp[0]
    nnnq = nnq[0]

    if verbose:
        print('the symmetry vector is', nnnp, nnnq)

    lp = nnnp * nnnp + nnnq * nnnq + nnnp * nnnq
    r = a * sqrt(lp)
    c = a * l1
    t = sq3 * c / ndr

    if 2 * nn > nk:
        raise RuntimeError('parameter nk is too small!')

    rs = c / (2.0 * np.pi)

    if verbose:
        print('radius=', rs, t)

    q1 = np.arctan((sq3 * m) / (2 * n + m))
    q2 = np.arctan((sq3 * nnnq) / (2 * nnnp + nnnq))
    q3 = q1 - q2

    q4 = 2.0 * np.pi / nn
    q5 = bond * np.cos((np.pi / 6.0) - q1) / c * 2.0 * np.pi

    h1 = abs(t) / abs(np.sin(q3))
    h2 = bond * np.sin((np.pi / 6.0) - q1)

    ii = 0
    x, y, z = [], [], []
    for i in range(nn):
        x1, y1, z1 = 0, 0, 0

        k = np.floor(i * abs(r) / h1)
        x1 = rs * np.cos(i * q4)
        y1 = rs * np.sin(i * q4)
        z1 = (i * abs(r) - k * h1) * np.sin(q3)
        kk2 = abs(np.floor((z1 + 0.0001) / t))
        if z1 >= t - 0.0001:
            z1 -= t * kk2
        elif z1 < 0:
            z1 += t * kk2
        ii += 1

        x.append(x1)
        y.append(y1)
        z.append(z1)
        z3 = (i * abs(r) - k * h1) * np.sin(q3) - h2
        ii += 1

        if z3 >= 0 and z3 < t:
            x2 = rs * np.cos(i * q4 + q5)
            y2 = rs * np.sin(i * q4 + q5)
            z2 = (i * abs(r) - k * h1) * np.sin(q3) - h2
            x.append(x2)
            y.append(y2)
            z.append(z2)
        else:
            x2 = rs * np.cos(i * q4 + q5)
            y2 = rs * np.sin(i * q4 + q5)
            z2 = (i * abs(r) - (k + 1) * h1) * np.sin(q3) - h2
            kk = abs(np.floor(z2 / t))
            if z2 >= t - 0.0001:
                z2 -= t * kk
            elif z2 < 0:
                z2 += t * kk
            x.append(x2)
            y.append(y2)
            z.append(z2)

    ntotal = 2 * nn
    X = []
    for i in range(ntotal):
        X.append([x[i], y[i], sign * z[i]])

    # Calculate number of unit cells needed for desired length
    num_cells = max(1, int(length / t))
    if num_cells > 1:
        xx = X[:]
        for mnp in range(2, num_cells + 1):
            for i in range(len(xx)):
                X.append(xx[i][:2] + [xx[i][2] + (mnp - 1) * t])

    transvec = t
    numatom = ntotal * num_cells
    diameter = rs * 2
    chiralangle = np.arctan((sq3 * n) / (2 * m + n)) / np.pi * 180

    cell = [[0, 0, 0], [0, 0, 0], [0, 0, length * t]]
    atoms = Atoms(symbol + str(numatom),
                  positions=X,
                  cell=cell,
                  pbc=[False, False, True])
    if vacuum:
        atoms.center(vacuum, axis=(0, 1))
    if verbose:
        print('translation vector =', transvec)
        print('diameter = ', diameter)
        print('chiral angle = ', chiralangle)
    return atoms








#========================
#========================
#========================
#========================
#========================
#========================



class Nano_OneD_Builder:
    """
    Builder class for 1D nanomaterials such as nanowires, nanorods, and nanotubes.

    This class provides flexible tools to construct, modify, and export 
    one-dimensional nanostructures based on atomic-scale models. 
    It supports metallic, semiconductor, and carbon-based nanostructures, 
    including multi-walled carbon nanotubes (MWCNTs), nanowire arrays, 
    and nanotube bundles.

    Structures can be generated from ASE's bulk and nanotube generators,
    and then modified using defect insertion, alloying, or strain application.
    The generated geometries can be exported in multiple formats (XYZ, CIF, PDB, etc.)
    or directly returned as ASE `Atoms` objects or GAM-compatible atom lists.

    Parameters
    ----------
    material : str
        Chemical symbol of the element (e.g., "Au", "C", "Si").
    structure_type : str
        Type of nanostructure to generate: one of 
        ``"nanowire"``, ``"nanorod"``, or ``"nanotube"``.
    lattice_constant : float, optional
        Lattice constant of the bulk crystal in Ångströms. 
        Required for crystalline materials other than carbon nanotubes.
    crystal_structure : str, default="fcc"
        ASE bulk structure type (e.g., ``"fcc"``, ``"bcc"``, ``"hcp"``, ``"diamond"``).
    direction : tuple of int, default=(1, 0, 0)
        Miller indices defining the nanowire growth direction.
    radius : float, default=1.0
        Radius of the 1D structure in nanometers.
    length : float, default=5.0
        Length of the structure in nanometers.
    vacuum : float, default=10.0
        Vacuum spacing in Ångströms.
    **kwargs :
        Additional parameters reserved for future extensions.

    Attributes
    ----------
    ASE_atoms : ase.Atoms
        The atomic structure generated using ASE.
    atoms : list of GAM_Atom
        List of atoms in GAM-compatible format.
    meta : dict
        Metadata describing structure parameters.
    structure_type : str
        The chosen type of nanostructure ("nanowire", "nanorod", or "nanotube").
    direction : tuple of int
        Growth direction (Miller indices).
    radius, length, vacuum : float
        Structural dimensions and vacuum spacing.
    
    Notes
    -----
    - Nanotubes are generated using ASE's `ase.build.nanotube()` function.
    - Nanowires/nanorods are cut from bulk crystals and shaped cylindrically.
    - All distances are internally converted to Ångströms.
    - The builder can introduce defects, create alloys, or apply strain
      to simulate realistic material behavior.

    Examples
    --------
    >>> builder = Nano_OneD_Builder(material="Au", structure_type="nanowire",
    ...                             lattice_constant=4.08, crystal_structure="fcc",
    ...                             radius=1.5, length=10.0)
    >>> builder.add_defects(defect_type="vacancy", concentration=0.02)
    >>> builder.apply_strain(0.05, direction='z')
    >>> builder.export("gold_nanowire", formats=['xyz', 'cif'])
    >>> ase_atoms = builder.to_ase()
    """

    def __init__(self, material: str,structure_type: str,
                 lattice_constant: float = None, crystal_structure: str = "fcc",
                 direction=(1, 0, 0),
               radius=1.0, length=5.0, vacuum=10.0, **kwargs):
        """
        Create a 1D nanomaterial.

        :param material: Chemical symbol, e.g., "Au", "C", "Si"
        :param lattice_constant: Lattice constant in Å (optional)
        :param crystal_structure: ASE bulk structure type ("fcc", "bcc", "hcp", "diamond", etc.)
        :param structure_type: "nanowire", "nanorod", or "nanotube"
        :param direction: Miller indices tuple for crystal growth direction
        :param radius: Radius in nanometers
        :param length: Length in nanometers
        :param vacuum: Vacuum spacing in Å
        :return: ASE Atoms object
        """
        self.material = material
        self.lattice_constant = lattice_constant
        self.crystal_structure = crystal_structure



        self.structure_type=structure_type
        self.direction=direction
        self.radius=radius
        self.length=length
        self.vacuum=vacuum


        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []


        self.meta={'material':material,
        'lattice_constant':lattice_constant ,
        'crystal_structure':crystal_structure ,
        'structure_type': structure_type,
        'direction': direction,
        'radius': radius,
        'length': length,
        'vacuum':vacuum }


        radius_ang = radius * 10.0  # nm → Å
        length_ang = length * 10.0  # nm → Å

        if structure_type.lower() == "nanotube":
            atoms = ase_nanotube(6, 6, length=length_ang, vacuum=vacuum, symbol=self.material)
            
            self.ASE_atoms=atoms
            
            self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        elif structure_type.lower() in ["nanowire", "nanorod"]:
            # Create bulk structure
            atoms = bulk(self.material, crystalstructure=self.crystal_structure,
                         a=self.lattice_constant)
            # Make it long in the chosen direction
            repetitions = (int(length_ang / cell_to_cellpar(atoms.cell)[0]) + 1,
                           int(length_ang / cell_to_cellpar(atoms.cell)[1]) + 1,
                           int(length_ang / cell_to_cellpar(atoms.cell)[2]) + 1)
            atoms = atoms.repeat(repetitions)

            # Select atoms within cylindrical radius
            center = atoms.get_center_of_mass()
            positions = atoms.get_positions()
            distances = np.linalg.norm(positions[:, :2] - center[:2], axis=1)
            mask = distances <= radius_ang
            atoms = atoms[mask]

            # Add vacuum
            atoms.center(vacuum=vacuum, axis=(0, 1))

            self.ASE_atoms=atoms
            self.atoms=self._to_GAM_Atoms(self.ASE_atoms)


        else:
            raise ValueError(f"Unknown 1D structure type: {structure_type}")

        #return atoms

    def _create(self, material: str,
           structure_type: str,
           lattice_constant: float = None,
           crystal_structure: str = "fcc",
           direction=(1, 0, 0),
           radius: float = 1.0,
           length: float = 5.0,
           vacuum: float = 10.0):
        radius_ang = radius * 10.0  # nm → Å
        length_ang = length * 10.0  # nm → Å

        if structure_type.lower() == "nanotube":
            atoms = ase_nanotube(6, 6, length=length_ang, vacuum=vacuum, symbol=material)
        elif structure_type.lower() in ["nanowire", "nanorod"]:
            atoms = bulk(material, crystalstructure=crystal_structure,a = lattice_constant)

            reps = (int(length_ang / cell_to_cellpar(atoms.cell)[0]) + 1,
                    int(length_ang / cell_to_cellpar(atoms.cell)[1]) + 1,
                    int(length_ang / cell_to_cellpar(atoms.cell)[2]) + 1)
            atoms = atoms.repeat(reps)

            center = atoms.get_center_of_mass()
            distance = np.linalg.norm(atoms.get_positions()[:, :2] - center[:2], axis=1)
            atoms = atoms[distance <= radius_ang]

            atoms.center(vacuum=vacuum, axis=(0, 1))

        else:
            raise ValueError(f"Unknown 1D structure type: {structure_type}")
        return atoms
    
    

    
    
    

    def add_defects(self, defect_type='vacancy', 
                    concentration=0.05, seed=None, min_distance=2.0):
        """
        Add defects to the structure.
        
        Parameters:
        -----------
        defect_type : str
            Type of defect ('vacancy', 'interstitial')
        concentration : float
            Concentration of defects (0-1)
        seed : int
            Random seed for reproducibility
        min_distance : float
            Minimum distance between atoms for interstitial defects (in Angstroms)
        """
        if seed is not None:
            random.seed(seed)
            
        num_atoms = len(self.ASE_atoms)
        num_defects = int(num_atoms * concentration)
        
        if defect_type == 'vacancy':
            # Create vacancies by removing random atoms
            indices = list(range(num_atoms))
            to_remove = random.sample(indices, num_defects)
            del self.ASE_atoms[to_remove]
            
        elif defect_type == 'interstitial':
            # Add interstitial atoms at appropriate positions
            existing_positions = self.ASE_atoms.get_positions()
            cell = self.ASE_atoms.get_cell()
            
            # Get the actual bounds of the structure (not the full cell)
            if self.structure_type.lower() == "nanotube":
                # For nanotubes, find the cylindrical region
                center = np.mean(existing_positions, axis=0)
                # Use x,y coordinates to find radius
                distances_from_center = np.sqrt((existing_positions[:, 0] - center[0])**2 + 
                                             (existing_positions[:, 1] - center[1])**2)
                max_radius = np.max(distances_from_center) + 2.0  # Add buffer
                
                # Find z bounds
                z_min, z_max = np.min(existing_positions[:, 2]), np.max(existing_positions[:, 2])
                z_range = z_max - z_min
                
                attempts = 0
                max_attempts = num_defects * 100  # Limit attempts to avoid infinite loop
                
                for _ in range(num_defects):
                    if attempts >= max_attempts:
                        print(f"Warning: Could only place {_}/{num_defects} interstitial atoms due to space constraints")
                        break
                        
                    while attempts < max_attempts:
                        attempts += 1
                        
                        # Generate random position within the nanotube region
                        angle = random.uniform(0, 2 * np.pi)
                        radius = random.uniform(0, max_radius)
                        z_pos = random.uniform(z_min - 1.0, z_max + 1.0)  # Slight extension beyond bounds
                        
                        x_pos = center[0] + radius * np.cos(angle)
                        y_pos = center[1] + radius * np.sin(angle)
                        new_pos = np.array([x_pos, y_pos, z_pos])
                        
                        # Check distance from existing atoms
                        distances = np.linalg.norm(existing_positions - new_pos, axis=1)
                        if np.min(distances) >= min_distance:
                            # Position is good, add the atom
                            self.ASE_atoms.append(Atom(self.material, position=new_pos))
                            existing_positions = np.vstack([existing_positions, new_pos])
                            break
                            
            else:  # nanowire/nanorod
                # For nanowires/nanorods, find the cylindrical region
                center = np.mean(existing_positions, axis=0)
                # Use x,y coordinates to find radius
                distances_from_center = np.sqrt((existing_positions[:, 0] - center[0])**2 + 
                                             (existing_positions[:, 1] - center[1])**2)
                max_radius = np.max(distances_from_center) + 2.0  # Add buffer
                
                # Find z bounds
                z_min, z_max = np.min(existing_positions[:, 2]), np.max(existing_positions[:, 2])
                
                attempts = 0
                max_attempts = num_defects * 100
                
                for _ in range(num_defects):
                    if attempts >= max_attempts:
                        print(f"Warning: Could only place {_}/{num_defects} interstitial atoms due to space constraints")
                        break
                        
                    while attempts < max_attempts:
                        attempts += 1
                        
                        # Generate random position within the nanowire region
                        angle = random.uniform(0, 2 * np.pi)
                        radius = random.uniform(0, max_radius)
                        z_pos = random.uniform(z_min - 1.0, z_max + 1.0)
                        
                        x_pos = center[0] + radius * np.cos(angle)
                        y_pos = center[1] + radius * np.sin(angle)
                        new_pos = np.array([x_pos, y_pos, z_pos])
                        
                        # Check distance from existing atoms
                        distances = np.linalg.norm(existing_positions - new_pos, axis=1)
                        if np.min(distances) >= min_distance:
                            # Position is good, add the atom
                            self.ASE_atoms.append(Atom(self.material, position=new_pos))
                            existing_positions = np.vstack([existing_positions, new_pos])
                            break
        
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)




    def create_alloy(self, elements: list, 
                     concentrations: list, seed=None):
        """
        Create an alloy by replacing atoms with different elements.
        
        Parameters:
        -----------
        elements : list
            List of element symbols to add
        concentrations : list
            List of concentrations for each element (should sum to <= 1)
        seed : int
            Random seed for reproducibility
        """
        if seed is not None:
            random.seed(seed)
            
        if len(elements) != len(concentrations):
            raise ValueError("Number of elements must match number of concentrations")
            
        if sum(concentrations) > 1:
            raise ValueError("Sum of concentrations must be <= 1")
            
        num_atoms = len(self.ASE_atoms)
        original_symbols = self.ASE_atoms.get_chemical_symbols()
        
        # Calculate number of atoms for each element
        num_each = [int(num_atoms * conc) for conc in concentrations]
        
        # Randomly select atoms to replace
        all_indices = list(range(num_atoms))
        random.shuffle(all_indices)
        
        current_idx = 0
        for element, num in zip(elements, num_each):
            indices_to_change = all_indices[current_idx:current_idx + num]
            for idx in indices_to_change:
                self.ASE_atoms[idx].symbol = element
            current_idx += num
            
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

    def apply_strain(self,strain: float, direction: str = 'z'):
        """
        Apply strain to the structure.
        
        Parameters:
        -----------
        strain : float
            Strain value (positive for tension, negative for compression)
        direction : str
            Direction to apply strain ('x', 'y', or 'z')
        """
        cell = self.ASE_atoms.get_cell()
        positions = self.ASE_atoms.get_positions()
        
        # Apply strain in the specified direction
        if direction == 'x':
            cell[0] *= (1 + strain)
            positions[:, 0] *= (1 + strain)
        elif direction == 'y':
            cell[1] *= (1 + strain)
            positions[:, 1] *= (1 + strain)
        else:  # 'z'
            cell[2] *= (1 + strain)
            positions[:, 2] *= (1 + strain)
            
        self.ASE_atoms.set_cell(cell)
        self.ASE_atoms.set_positions(positions)
        
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        
    def create_multiple_nanowires(self, num_wires: int, spacing: float = 2.0, 
                                arrangement: str = 'linear', **wire_params):
        """
        Create multiple parallel nanowires/nanorods arranged in a specific pattern.
        
        Parameters:
        -----------
        num_wires : int
            Number of nanowires to create
        spacing : float
            Spacing between wires in nanometers
        arrangement : str
            Arrangement pattern ('linear', 'hexagonal', 'square')
        wire_params : dict
            Parameters to pass to create() method for each wire
        
        Returns:
        --------
        ase.atoms.Atoms
            Combined structure of multiple nanowires
        """
        spacing_ang = spacing * 10.0  # nm → Å
        all_atoms = []
        
        # Create positions based on arrangement
        if arrangement == 'linear':
            positions = [(i * spacing_ang, 0, 0) for i in range(num_wires)]
        
        elif arrangement == 'hexagonal':
            positions = []
            layers = int(np.ceil(np.sqrt(num_wires / 2)))
            current_wire = 0
            for layer in range(layers):
                y = layer * spacing_ang * np.sqrt(3) / 2
                offset = (layer % 2) * spacing_ang / 2
                x_positions = np.arange(layers) * spacing_ang + offset
                for x in x_positions:
                    if current_wire < num_wires:
                        positions.append((x, y, 0))
                        current_wire += 1
        
        elif arrangement == 'square':
            side_length = int(np.ceil(np.sqrt(num_wires)))
            positions = []
            for i in range(side_length):
                for j in range(side_length):
                    if len(positions) < num_wires:
                        positions.append((i * spacing_ang, j * spacing_ang, 0))
        
        else:
            raise ValueError(f"Unknown arrangement pattern: {arrangement}")
        
        # Create individual wires and combine them
        combined_atoms = None
        for pos in positions:
            # Create single wire
            wire = self._create(material = self.material,structure_type = "nanowire", lattice_constant = self.lattice_constant,crystal_structure = self.crystal_structure,direction = self.direction, radius = self.radius, length = self.length, vacuum = self.vacuum)

            # Translate wire to its position
            wire.translate(pos)
            
            if combined_atoms is None:
                combined_atoms = wire
            else:
                combined_atoms.extend(wire)


        self.ASE_atoms=combined_atoms
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
    
    
    
    def create_mwcnt(self, num_walls: int, initial_n: int = 6, initial_m: int = 6, 
                    wall_spacing: float = 0.34, length: float = 10.0, vacuum: float = 10.0):
        """
        Create a multi-walled carbon nanotube (MWCNT) with concentric tubes.
        
        Parameters:
        -----------
        num_walls : int
            Number of walls in the MWCNT
        initial_n, initial_m : int
            Chiral indices for the innermost tube
        wall_spacing : float
            Spacing between walls in nm (default 0.34 nm, typical for graphite)
        length : float
            Length of the MWCNT in nm
        vacuum : float
            Vacuum spacing in Angstroms
        
        Returns:
        --------
        ase.atoms.Atoms
            Multi-walled carbon nanotube structure with concentric tubes
        """
        if self.material != "C":
            raise ValueError("Multi-walled nanotubes are only supported for carbon")
            
        wall_spacing_ang = wall_spacing * 10.0  # nm → Å
        length_ang = length * 10.0  # nm → Å
        
        # Calculate initial radius
        initial_radius = self._calculate_tube_radius(initial_n, initial_m)
        
        # Create innermost tube
        combined_atoms = None
        current_n, current_m = initial_n, initial_m
        
        for wall in range(num_walls):
            # Calculate indices for next tube to maintain proper spacing
            if wall > 0:
                target_radius = initial_radius + wall * wall_spacing_ang
                # Find appropriate n,m that gives radius close to target
                best_diff = float('inf')
                best_n = current_n
                best_m = current_m
                
                # Search for best n,m combination
                for dn in range(1, 6):
                    for dm in range(1, 6):
                        test_n = current_n + dn
                        test_m = current_m + dm
                        test_radius = self._calculate_tube_radius(test_n, test_m)
                        diff = abs(test_radius - target_radius)
                        if diff < best_diff:
                            best_diff = diff
                            best_n = test_n
                            best_m = test_m
                
                current_n, current_m = best_n, best_m
            
            # Create current wall
            current_tube = ase_nanotube(current_n, current_m, 
                                      length=length_ang, 
                                      vacuum=vacuum, 
                                      symbol="C")
            
            # Center the tube
            cell_center = current_tube.cell.diagonal() / 2
            positions = current_tube.get_positions()
            com = positions.mean(axis=0)
            translation = cell_center - com
            current_tube.translate(translation)
            
            if combined_atoms is None:
                combined_atoms = current_tube
            else:
                # Ensure new tube is centered relative to existing structure
                existing_com = combined_atoms.get_positions().mean(axis=0)
                new_com = current_tube.get_positions().mean(axis=0)
                translation = existing_com - new_com
                current_tube.translate(translation)
                combined_atoms.extend(current_tube)
        
        # Final centering of entire structure
        cell_center = combined_atoms.cell.diagonal() / 2
        positions = combined_atoms.get_positions()
        com = positions.mean(axis=0)
        translation = cell_center - com
        combined_atoms.translate(translation)



        self.ASE_atoms=combined_atoms
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        
        
    def _create_mwcnt(self, num_walls: int, initial_n: int = 6, initial_m: int = 6, 
                    wall_spacing: float = 0.34, length: float = 10.0, vacuum: float = 10.0):
        """
        Create a multi-walled carbon nanotube (MWCNT) with concentric tubes.
        
        Parameters:
        -----------
        num_walls : int
            Number of walls in the MWCNT
        initial_n, initial_m : int
            Chiral indices for the innermost tube
        wall_spacing : float
            Spacing between walls in nm (default 0.34 nm, typical for graphite)
        length : float
            Length of the MWCNT in nm
        vacuum : float
            Vacuum spacing in Angstroms
        
        Returns:
        --------
        ase.atoms.Atoms
            Multi-walled carbon nanotube structure with concentric tubes
        """
        if self.material != "C":
            raise ValueError("Multi-walled nanotubes are only supported for carbon")
            
        wall_spacing_ang = wall_spacing * 10.0  # nm → Å
        length_ang = length * 10.0  # nm → Å
        
        # Calculate initial radius
        initial_radius = self._calculate_tube_radius(initial_n, initial_m)
        
        # Create innermost tube
        combined_atoms = None
        current_n, current_m = initial_n, initial_m
        
        for wall in range(num_walls):
            # Calculate indices for next tube to maintain proper spacing
            if wall > 0:
                target_radius = initial_radius + wall * wall_spacing_ang
                # Find appropriate n,m that gives radius close to target
                best_diff = float('inf')
                best_n = current_n
                best_m = current_m
                
                # Search for best n,m combination
                for dn in range(1, 6):
                    for dm in range(1, 6):
                        test_n = current_n + dn
                        test_m = current_m + dm
                        test_radius = self._calculate_tube_radius(test_n, test_m)
                        diff = abs(test_radius - target_radius)
                        if diff < best_diff:
                            best_diff = diff
                            best_n = test_n
                            best_m = test_m
                
                current_n, current_m = best_n, best_m
            
            # Create current wall
            current_tube = ase_nanotube(current_n, current_m, 
                                      length=length_ang, 
                                      vacuum=vacuum, 
                                      symbol="C")
            
            # Center the tube
            cell_center = current_tube.cell.diagonal() / 2
            positions = current_tube.get_positions()
            com = positions.mean(axis=0)
            translation = cell_center - com
            current_tube.translate(translation)
            
            if combined_atoms is None:
                combined_atoms = current_tube
            else:
                # Ensure new tube is centered relative to existing structure
                existing_com = combined_atoms.get_positions().mean(axis=0)
                new_com = current_tube.get_positions().mean(axis=0)
                translation = existing_com - new_com
                current_tube.translate(translation)
                combined_atoms.extend(current_tube)
        
        # Final centering of entire structure
        cell_center = combined_atoms.cell.diagonal() / 2
        positions = combined_atoms.get_positions()
        com = positions.mean(axis=0)
        translation = cell_center - com
        combined_atoms.translate(translation)



        #self.ASE_atoms=combined_atoms
        #self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        
        #return self.ASE_atoms
        return combined_atoms

    def _calculate_tube_radius(self, n: int, m: int, bond_length: float = 1.42):
        """
        Calculate the radius of a carbon nanotube with given (n,m) indices.
        
        Parameters:
        -----------
        n, m : int
            Chiral indices
        bond_length : float
            C-C bond length in Angstroms
            
        Returns:
        --------
        float
            Tube radius in Angstroms
        """
        a = bond_length * np.sqrt(3)  # graphene lattice constant
        return a * np.sqrt(n*n + m*m + n*m) / (2 * np.pi)
    
    
    def create_nanotube_bundle(self, num_tubes: int, tube_type: str = 'SWCNT',
                             arrangement: str = 'hexagonal', spacing: float = 0.34,
                             **tube_params):
        """
        Create a bundle of carbon nanotubes arranged in a specific pattern.
        
        Parameters:
        -----------
        num_tubes : int
            Number of tubes in the bundle
        tube_type : str
            Type of tubes ('SWCNT' or 'MWCNT')
        arrangement : str
            Arrangement pattern ('hexagonal', 'square', 'circular')
        spacing : float
            Spacing between tubes in nm
        tube_params : dict
            Parameters for tube creation (passed to create() or create_mwcnt())
        
        Returns:
        --------
        ase.atoms.Atoms
            Bundle of nanotubes
        """
        if self.material != "C":
            raise ValueError("Nanotube bundles are only supported for carbon")
            
        spacing_ang = spacing * 10.0  # nm → Å
        
        # Create a single tube first to get its radius
        if tube_type == 'SWCNT':
            template_tube = self._create(structure_type="nanotube", **tube_params)
        else:  # MWCNT
            template_tube = self._create_mwcnt(**tube_params)
        
        positions = template_tube.get_positions()
        print("Template tube positions:", positions)
        tube_radius = np.max(np.sqrt(positions[:,0]**2 + positions[:,1]**2))
        effective_spacing = max(spacing_ang, 2 * tube_radius)
        
        # Calculate positions based on arrangement
        if arrangement == 'hexagonal':
            # For hexagonal packing, calculate required layers
            positions = []
            if num_tubes == 1:
                positions = [(0, 0, 0)]
            else:
                # Calculate number of shells needed
                shell = 1
                tubes_in_shells = 1
                while tubes_in_shells < num_tubes:
                    shell += 1
                    tubes_in_shells += 6 * (shell - 1)
                
                # Generate hexagonal grid
                for n in range(shell):
                    if n == 0:
                        if num_tubes > 0:
                            positions.append((0, 0, 0))
                    else:
                        # Add tubes in current shell
                        for i in range(6):  # 6 sides
                            for j in range(n):  # points per side
                                angle = i * np.pi / 3
                                next_angle = (i + 1) * np.pi / 3
                                
                                # Interpolate between corners
                                t = j / n
                                r1 = n * effective_spacing
                                x = r1 * (np.cos(angle) * (1-t) + np.cos(next_angle) * t)
                                y = r1 * (np.sin(angle) * (1-t) + np.sin(next_angle) * t)
                                
                                if len(positions) < num_tubes:
                                    positions.append((x, y, 0))
                                    
        elif arrangement == 'square':
            side_length = int(np.ceil(np.sqrt(num_tubes)))
            center_offset = (side_length - 1) * effective_spacing / 2
            positions = []
            for i in range(side_length):
                for j in range(side_length):
                    if len(positions) < num_tubes:
                        x = i * effective_spacing - center_offset
                        y = j * effective_spacing - center_offset
                        positions.append((x, y, 0))
                        
        elif arrangement == 'circular':
            positions = []
            if num_tubes == 1:
                positions = [(0, 0, 0)]
            else:
                # Place tubes in concentric circles
                tubes_placed = 0
                shell = 0
                while tubes_placed < num_tubes:
                    if shell == 0:
                        positions.append((0, 0, 0))
                        tubes_placed += 1
                    else:
                        # Calculate number of tubes that fit in this shell
                        circumference = 2 * np.pi * shell * effective_spacing
                        tubes_in_shell = min(
                            int(circumference / effective_spacing),
                            num_tubes - tubes_placed
                        )
                        
                        if tubes_in_shell > 0:
                            for i in range(tubes_in_shell):
                                angle = 2 * np.pi * i / tubes_in_shell
                                r = shell * effective_spacing
                                x = r * np.cos(angle)
                                y = r * np.sin(angle)
                                positions.append((x, y, 0))
                                tubes_placed += 1
                    shell += 1
        
        else:
            raise ValueError(f"Unknown arrangement pattern: {arrangement}")
        
        # Create individual tubes and combine them
        combined_atoms = None
        for pos in positions:
            # Create tube
            if tube_type == 'SWCNT':
                tube = self._create(structure_type="nanotube", **tube_params)
            else:  # MWCNT
                tube = self._create_mwcnt(**tube_params)
            
            # Center each tube before translation
            tube_com = tube.get_positions().mean(axis=0)
            tube.translate(-tube_com)
            
            # Translate to final position
            tube.translate(pos)
            
            if combined_atoms is None:
                combined_atoms = tube
            else:
                combined_atoms.extend(tube)
        
        # Center the entire bundle
        if combined_atoms is not None:
            cell_center = combined_atoms.cell.diagonal() / 2
            bundle_com = combined_atoms.get_positions().mean(axis=0)
            translation = cell_center - bundle_com
            combined_atoms.translate(translation)
        

        self.ASE_atoms=combined_atoms
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        
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


    
    

    def export(self,  filename: str, formats=None, 
               trajectory=False, frame_interval=1):
        """
        Export structure in various formats with additional options.
        
        Parameters:
        -----------
        filename : str
            Base filename without extension
        formats : list
            List of format strings (default: ['xyz', 'cif', 'vasp', 'pdb'])
        trajectory : bool
            Whether to save as trajectory (useful for dynamics)
        frame_interval : int
            Interval between frames for trajectory
        """
        if formats is None:
            formats = ['xyz', 'cif', 'vasp', 'pdb']
            
        for fmt in formats:
            output_file = f"{filename}.{fmt}"
            try:
                write(output_file, self.ASE_atoms)
                print(f"Successfully exported to {output_file}")
            except Exception as e:
                print(f"Failed to export to {fmt} format: {str(e)}")
                
        if trajectory:
            traj_file = f"{filename}.traj"
            traj = Trajectory(traj_file, 'w')
            traj.write(self.ASE_atoms)
            print(f"Saved trajectory to {traj_file}")


    def to_xyz(self, filename: str) -> None:
        """
        Save the structure to a file.
        
        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object to save
        filename : str
            Output filename
        format : str
            Output format 'xyz'
        """
        write(filename, self.ASE_atoms, format='xyz')
        
        
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
        
        #or 
        #Just 
        #return self.ASE_atoms
        
        return atoms_obj
        
        
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














