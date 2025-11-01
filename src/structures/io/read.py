import json
import re
from typing import List, Dict, Optional, Tuple, Any
from ..Primatom import GAM_Atom
import json
from typing import List, Dict, Optional, Tuple
from ..Primatom import GAM_Atom , GAM_Bond , GAM_Molecule




def read_structure(file_path: str, fmt: Optional[str] = None) -> List[GAM_Atom]:
    """
    Universal import function to read molecular structure files into GAM_Atom objects.

    Parameters
    ----------
    file_path : str
        Path to the structure file.
    fmt : str, optional
        File format to read ('xyz', 'cif', 'poscar', 'pdb', 'lammps', 'gaussian', 'pygam').
        If None, the format is auto-detected from the file extension.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects representing the atomic structure.

    Raises
    ------
    ValueError
        If the format cannot be determined or is unsupported.
    """
    if fmt is None:
        # Auto-detect format from file extension
        file_ext = file_path.lower().split('.')[-1]
        if file_ext == 'xyz':
            fmt = 'xyz'
        elif file_ext == 'cif':
            fmt = 'cif'
        elif file_path.endswith('POSCAR') or file_ext == 'poscar':
            fmt = 'poscar'
        elif file_ext in ['pdb']:
            fmt = 'pdb'
        elif file_ext in ['lammps', 'data']:
            fmt = 'lammps'
        elif file_ext in ['gjf', 'com', 'gaussian']:
            fmt = 'gaussian'
        elif file_ext in ['pygam', 'pygamlab']:
            fmt = 'pygam'
        else:
            raise ValueError(f"Cannot auto-detect format for file: {file_path}")
    
    fmt = fmt.lower()
    
    if fmt == "xyz":
        return read_xyz(file_path)
    elif fmt == "cif":
        return read_cif(file_path)
    elif fmt == "poscar":
        return read_poscar(file_path)
    elif fmt == "pdb":
        return read_pdb(file_path)
    elif fmt == "lammps":
        return read_lammps(file_path)
    elif fmt in ("gjf", "gaussian", "com"):
        return read_gaussian(file_path)
    elif fmt in ('pygamlab','pygam','gamlab'):
        return read_pygam(file_path)
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def read_xyz(file_path: str) -> List[GAM_Atom]:
    """
    Read an XYZ file and return a list of GAM_Atom objects.

    XYZ format:
    Line 1: number of atoms
    Line 2: comment
    Following lines: element X Y Z

    Parameters
    ----------
    file_path : str
        Path to the XYZ file.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with positions read from the file.

    Raises
    ------
    ValueError
        If the file is invalid or coordinates cannot be parsed.
    """
    atoms = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Remove empty lines and strip whitespace
    lines = [line.strip() for line in lines if line.strip()]
    
    if len(lines) < 2:
        raise ValueError("Invalid XYZ file: too few lines")
    
    # Read number of atoms
    try:
        n_atoms = int(lines[0])
    except ValueError:
        raise ValueError("Invalid XYZ file: first line must be number of atoms")
    
    # Skip comment line and read atoms
    atom_lines = lines[2:2+n_atoms]
    
    for i, line in enumerate(atom_lines):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid atom line {i+3}: {line}")
        
        element = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            raise ValueError(f"Invalid coordinates on line {i+3}: {line}")
        
        # Create GAM_Atom with basic properties
        atom = GAM_Atom(
            id=i+1,
            element=element,
            x=x, y=y, z=z
        )
        atoms.append(atom)
    
    print(f"Successfully read {len(atoms)} atoms from XYZ file: {file_path}")
    return atoms


def read_cif(file_path: str) -> List[GAM_Atom]:
    """
    Read a CIF (Crystallographic Information File) and return a list of GAM_Atom objects.

    Notes
    -----
    This is a basic CIF reader focusing only on atomic coordinates.

    Parameters
    ----------
    file_path : str
        Path to the CIF file.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with fractional or Cartesian coordinates.

    Raises
    ------
    ValueError
        If coordinates cannot be parsed.
    """
    atoms = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Find the atom site loop
    in_atom_loop = False
    atom_data_started = False
    label_idx, x_idx, y_idx, z_idx = None, None, None, None
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('loop_'):
            in_atom_loop = False
            atom_data_started = False
            continue
        
        if in_atom_loop and line.startswith('_atom_site_'):
            if 'label' in line:
                label_idx = len([l for l in lines[:lines.index(line+'\n')] if l.strip().startswith('_atom_site_')])
            elif 'fract_x' in line:
                x_idx = len([l for l in lines[:lines.index(line+'\n')] if l.strip().startswith('_atom_site_')])
            elif 'fract_y' in line:
                y_idx = len([l for l in lines[:lines.index(line+'\n')] if l.strip().startswith('_atom_site_')])
            elif 'fract_z' in line:
                z_idx = len([l for l in lines[:lines.index(line+'\n')] if l.strip().startswith('_atom_site_')])
            continue
        
        if line.startswith('_atom_site_label'):
            in_atom_loop = True
            label_idx = 0
            continue
        
        if in_atom_loop and not line.startswith('_') and line:
            if all(idx is not None for idx in [label_idx, x_idx, y_idx, z_idx]):
                parts = line.split()
                if len(parts) > max(label_idx, x_idx, y_idx, z_idx):
                    element = parts[label_idx].strip('0123456789')
                    try:
                        x = float(parts[x_idx])
                        y = float(parts[y_idx]) 
                        z = float(parts[z_idx])
                        
                        atom = GAM_Atom(
                            id=len(atoms)+1,
                            element=element,
                            x=x, y=y, z=z
                        )
                        atoms.append(atom)
                    except (ValueError, IndexError):
                        continue
    
    print(f"Successfully read {len(atoms)} atoms from CIF file: {file_path}")
    return atoms


def read_poscar(file_path: str) -> List[GAM_Atom]:
    """
    Read a POSCAR (VASP) file and return a list of GAM_Atom objects.

    Parameters
    ----------
    file_path : str
        Path to the POSCAR file.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with coordinates converted to Cartesian if necessary.

    Raises
    ------
    ValueError
        If the file format is invalid or missing critical lines.
    """
    atoms = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    lines = [line.strip() for line in lines if line.strip()]
    
    if len(lines) < 8:
        raise ValueError("Invalid POSCAR file: too few lines")
    
    # Skip comment line and scaling factor
    # Read lattice vectors (lines 2-4)
    lattice = []
    for i in range(2, 5):
        vec = [float(x) for x in lines[i].split()]
        lattice.append(vec)
    
    # Read element names and counts
    elements = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    
    # Check if selective dynamics is present
    coord_start = 8
    if lines[7].lower().startswith('s'):  # Selective dynamics
        coord_start = 9
    
    coord_type = lines[coord_start-1].lower()
    is_cartesian = coord_type.startswith('c') or coord_type.startswith('k')
    
    # Read atomic coordinates
    atom_idx = 0
    for elem_idx, (element, count) in enumerate(zip(elements, counts)):
        for i in range(count):
            if coord_start + atom_idx >= len(lines):
                break
            
            coord_line = lines[coord_start + atom_idx].split()
            x, y, z = float(coord_line[0]), float(coord_line[1]), float(coord_line[2])
            
            # Convert fractional to Cartesian if needed
            if not is_cartesian:
                # Manual matrix multiplication for fractional to Cartesian conversion
                cart_x = x * lattice[0][0] + y * lattice[1][0] + z * lattice[2][0]
                cart_y = x * lattice[0][1] + y * lattice[1][1] + z * lattice[2][1]
                cart_z = x * lattice[0][2] + y * lattice[1][2] + z * lattice[2][2]
                x, y, z = cart_x, cart_y, cart_z
            
            atom = GAM_Atom(
                id=atom_idx+1,
                element=element,
                x=x, y=y, z=z
            )
            atoms.append(atom)
            atom_idx += 1
    
    print(f"Successfully read {len(atoms)} atoms from POSCAR file: {file_path}")
    return atoms


def read_pdb(file_path: str) -> List[GAM_Atom]:
    """
    Read a PDB (Protein Data Bank) file and return a list of GAM_Atom objects.

    Parameters
    ----------
    file_path : str
        Path to the PDB file.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with atomic coordinates.

    Notes
    -----
    The function reads 'ATOM' and 'HETATM' records and attempts to extract element symbols.

    """
    atoms = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            # PDB format: columns are fixed-width
            try:
                atom_id = int(line[6:11].strip())
                element = line[76:78].strip()
                if not element:  # fallback to atom name if element field is empty
                    element = line[12:16].strip().strip('0123456789')
                
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                atom = GAM_Atom(
                    id=atom_id,
                    element=element,
                    x=x, y=y, z=z
                )
                atoms.append(atom)
                
            except (ValueError, IndexError):
                continue
    
    print(f"Successfully read {len(atoms)} atoms from PDB file: {file_path}")
    return atoms


def read_lammps(file_path: str) -> List[GAM_Atom]:
    """
    Read a LAMMPS data file and return a list of GAM_Atom objects.

    Parameters
    ----------
    file_path : str
        Path to the LAMMPS data file.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with coordinates and charges.

    Notes
    -----
    Atom types are currently used as placeholders for element symbols.
    A proper type-to-element mapping may be needed for accurate representation.

    Raises
    ------
    ValueError
        If the 'Atoms' section cannot be found.
    """
    atoms = []
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Find Atoms section
    atoms_section_start = None
    for i, line in enumerate(lines):
        if line.strip().lower() == 'atoms':
            atoms_section_start = i + 1
            break
    
    if atoms_section_start is None:
        raise ValueError("No 'Atoms' section found in LAMMPS file")
    
    # Read atom data
    for i in range(atoms_section_start, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) < 7:
            continue
        
        try:
            atom_id = int(parts[0])
            # parts[1] is molecule ID, parts[2] is atom type
            charge = float(parts[3])
            x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
            
            # For simplicity, use atom type as element (would need type mapping in real use)
            element = f"T{parts[2]}"  # Placeholder - would need atom type to element mapping
            
            atom = GAM_Atom(
                id=atom_id,
                element=element,
                x=x, y=y, z=z,
                charge=charge
            )
            atoms.append(atom)
            
        except (ValueError, IndexError):
            continue
    
    print(f"Successfully read {len(atoms)} atoms from LAMMPS file: {file_path}")
    return atoms


def read_gaussian(file_path: str) -> List[GAM_Atom]:
    """
    Read a Gaussian input or output file and return a list of GAM_Atom objects.

    Parameters
    ----------
    file_path : str
        Path to the Gaussian file.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects with atomic coordinates.

    Raises
    ------
    ValueError
        If a geometry section cannot be detected in the file.
    """
    atoms = []
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    
    # Find geometry section (after charge and multiplicity line)
    geom_start = None
    for i, line in enumerate(lines):
        line = line.strip()
        # Look for charge multiplicity line (e.g., "0 1")
        if re.match(r'^-?\d+\s+\d+\s*$', line):
            geom_start = i + 1
            break
    
    if geom_start is None:
        raise ValueError("No geometry section found in Gaussian file")
    
    # Read atoms until blank line
    atom_id = 1
    for i in range(geom_start, len(lines)):
        line = lines[i].strip()
        if not line:  # Empty line ends geometry section
            break
        
        parts = line.split()
        if len(parts) >= 4:
            try:
                element = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                
                atom = GAM_Atom(
                    id=atom_id,
                    element=element,
                    x=x, y=y, z=z
                )
                atoms.append(atom)
                atom_id += 1
                
            except (ValueError, IndexError):
                continue
    
    print(f"Successfully read {len(atoms)} atoms from Gaussian file: {file_path}")
    return atoms


def read_pygam(file_path: str) -> Tuple[List[GAM_Atom], Dict[str, Any]]:
    """
    Read a PYGAM format file and return atoms and associated metadata.

    Parameters
    ----------
    file_path : str
        Path to the PYGAM file.

    Returns
    -------
    Tuple[List[GAM_Atom], Dict[str, Any]]
        - List[GAM_Atom]: Atomic positions and properties.
        - Dict[str, Any]: Metadata dictionary containing system, cell, timestep, bonds, dynamics, and additional metadata.

    Notes
    -----
    PYGAM format includes sections for atoms, bonds, dynamics parameters, and JSON metadata.
    """
    atoms = []
    metadata = {
        'system': '',
        'cell': None,
        'timestep': 0,
        'bonds': [],
        'dynamics': {},
        'metadata': {}
    }
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    current_section = None
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('#') or not line:
            continue
        
        if line.startswith('SYSTEM'):
            metadata['system'] = line.split(None, 1)[1] if len(line.split()) > 1 else ''
        elif line.startswith('NATOMS'):
            # Number of atoms - we'll count them as we read
            continue
        elif line.startswith('CELL'):
            parts = line.split()[1:]
            if len(parts) >= 3:
                metadata['cell'] = [float(x) for x in parts[:3]]
        elif line.startswith('TIMESTEP'):
            metadata['timestep'] = int(line.split()[1])
        elif line.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 11:
                try:
                    atom_id = int(parts[1])
                    element = parts[2]
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    vx, vy, vz = float(parts[6]), float(parts[7]), float(parts[8])
                    charge = float(parts[9])
                    electronic_state = parts[10]
                    mass = float(parts[11]) if len(parts) > 11 else None
                    
                    atom = GAM_Atom(
                        id=atom_id,
                        element=element,
                        x=x, y=y, z=z,
                        charge=charge
                    )
                    
                    # Set velocity if GAM_Atom supports it
                    if hasattr(atom, '_velocity'):
                        atom._velocity = [vx, vy, vz]
                    
                    # Set electronic state and mass if supported
                    if hasattr(atom, 'electronic_state'):
                        atom.electronic_state.name = electronic_state
                    if hasattr(atom, 'atomic_mass') and mass is not None:
                        atom.atomic_mass = mass
                    
                    atoms.append(atom)
                    
                except (ValueError, IndexError):
                    continue
        elif line.startswith('BOND'):
            parts = line.split()
            if len(parts) >= 5:
                try:
                    atom1, atom2 = int(parts[1]), int(parts[2])
                    order = int(parts[3])
                    length = float(parts[4])
                    metadata['bonds'].append((atom1, atom2, order, length))
                except (ValueError, IndexError):
                    continue
        elif line.startswith('METADATA'):
            # Parse JSON metadata
            json_str = line[8:].strip()  # Remove "METADATA" prefix
            try:
                metadata['metadata'] = json.loads(json_str)
            except json.JSONDecodeError:
                continue
        elif ' ' in line and not line.startswith('BOND') and not line.startswith('ATOM'):
            # Dynamics parameters
            parts = line.split(None, 1)
            if len(parts) == 2:
                try:
                    key, value = parts[0].lower(), float(parts[1])
                    metadata['dynamics'][key] = value
                except ValueError:
                    metadata['dynamics'][parts[0].lower()] = parts[1]
    
    print(f"Successfully read {len(atoms)} atoms from PYGAM file: {file_path}")
    return atoms, metadata


# Convenience function that returns just atoms from PYGAM
def read_pygam_atoms(file_path: str) -> List[GAM_Atom]:
    """
    Read a PYGAM file and return only the list of atoms.

    Parameters
    ----------
    file_path : str
        Path to the PYGAM file.

    Returns
    -------
    List[GAM_Atom]
        List of GAM_Atom objects from the PYGAM file.
    """
    """Read PYGAM file and return only the atoms list."""
    atoms, _ = read_pygam(file_path)
    return atoms








def from_pygam(pygam_str: str) -> Tuple[List[GAM_Atom], Dict[str, any]]:
    """
    Parse a PYGAM-formatted string and return atoms and metadata.

    Parameters
    ----------
    pygam_str : str
        Content of a PYGAM file as a string.

    Returns
    -------
    Tuple[List[GAM_Atom], Dict[str, any]]
        - List[GAM_Atom]: List of atoms with positions, charges, velocities, and states.
        - Dict[str, any]: Dictionary containing metadata, dynamics, and bonds information.
    """
    atoms: List[GAM_Atom] = []
    metadata: Dict[str, any] = {}
    dynamics: Dict[str, float] = {}
    bonds: List[Tuple[int, int, int, float]] = []

    lines = pygam_str.strip().splitlines()
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith("ATOM"):
            parts = line.split()
            id = int(parts[1])
            element = parts[2]
            x, y, z = map(float, parts[3:6])
            vx, vy, vz = map(float, parts[6:9])
            charge = float(parts[9])
            state = parts[10]
            mass = float(parts[11])
            atom = GAM_Atom(id=id, element=element, x=x, y=y, z=z,
                            charge=charge, electronic_state=state)
            atom._velocity[:] = [vx, vy, vz]
            atom.atomic_mass = mass
            atoms.append(atom)

        elif line.startswith("BOND"):
            _, a1, a2, order, length = line.split()
            bonds.append((int(a1), int(a2), int(order), float(length)))

        elif line.startswith("METADATA"):
            meta_str = line[len("METADATA"):].strip()
            metadata = json.loads(meta_str)

        else:
            parts = line.split()
            if len(parts) == 2 and parts[0].isupper():
                try:
                    dynamics[parts[0]] = float(parts[1])
                except ValueError:
                    pass

    return atoms, {"metadata": metadata, "dynamics": dynamics, "bonds": bonds}




def import_pygam(file_path: str) -> Tuple[List[GAM_Atom], Dict[str, any]]:
    """
    Import atoms and metadata from a PYGAM file on disk.

    Parameters
    ----------
    file_path : str
        Path to the PYGAM file.

    Returns
    -------
    Tuple[List[GAM_Atom], Dict[str, any]]
        - List[GAM_Atom]: Atoms read from the file.
        - Dict[str, any]: Metadata dictionary including dynamics and bonds.

    Notes
    -----
    This function reads the file and uses `from_pygam` internally for parsing.
    """
    with open(file_path, "r") as f:
        pygam_str = f.read()
    atoms, info = from_pygam(pygam_str)
    print(f".pygam file successfully imported: {file_path}")
    return atoms, info



