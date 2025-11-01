
from typing import List, Optional
from collections import Counter
import json
from typing import List, Dict, Optional, Tuple
from ..Primatom import GAM_Atom 





def export(atoms: list[GAM_Atom], fmt: str, file_path: Optional[str] = None, **kwargs) -> None:
    """
    Export a list of GAM_Atom objects to a file in a specified molecular format.

    This is a universal export function supporting multiple formats, including
    XYZ, CIF, POSCAR (VASP), PDB, LAMMPS, Gaussian, and the custom .pygam format.

    Parameters
    ----------
    atoms : list[GAM_Atom]
        List of atoms to export.
    fmt : str
        Target format for export. Supported values include:
        'xyz', 'cif', 'poscar', 'pdb', 'lammps', 'gaussian', 'gjf', 
        'pygamlab', 'pygam', 'gamlab'.
    file_path : str, optional
        File path to save the output. If None, a default path will be generated
        based on the chosen format.
    **kwargs : dict
        Additional keyword arguments passed to the corresponding format-specific exporter.

    Raises
    ------
    ValueError
        If an unsupported format is specified.
    """
    fmt = fmt.lower()
    
    # Generate default filename if not provided
    if file_path is None:
        if fmt in ("gjf", "gaussian"):
            file_path = "structure.gjf"
        elif fmt == "poscar":
            file_path = "POSCAR"
        elif fmt == "lammps":
            file_path = "data.lammps"
        elif fmt in ('pygamlab','pygam','gamlab'):
            file_path = "structure.pygam"
        else:
            file_path = f"structure.{fmt}"
    
    if fmt == "xyz":
        to_xyz(atoms, file_path=file_path, **kwargs)
    elif fmt == "cif":
        to_cif(atoms, file_path=file_path, **kwargs)
    elif fmt == "poscar":
        to_poscar(atoms, file_path=file_path, **kwargs)
    elif fmt == "pdb":
        to_pdb(atoms, file_path=file_path, **kwargs)
    elif fmt == "lammps":
        to_lammps(atoms, file_path=file_path, **kwargs)
    elif fmt in ("gjf", "gaussian"):
        to_gaussian(atoms, file_path=file_path, **kwargs)
    elif fmt in ('pygamlab','pygam','gamlab'):
        export_pygam(atoms, file_path=file_path, **kwargs)
    else:
        raise ValueError(f"Unsupported format: {fmt}")





def to_pygam(atoms: List[GAM_Atom],
             system: str = "GAM system",
             cell: Optional[List[float]] = None,
             timestep: int = 0,
             bonds: Optional[List[Tuple[int, int, int, float]]] = None,
             dynamics: Optional[Dict[str, float]] = None,
             metadata: Optional[dict] = None) -> str:
    """
    Convert a list of GAM_Atom objects and optional simulation info to a .pygam formatted string.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of atoms to include in the output.
    system : str
        Name or description of the system.
    cell : list of float, optional
        Unit cell dimensions [a, b, c] in Angstroms.
    timestep : int
        Simulation timestep (default is 0).
    bonds : list of tuples, optional
        Bond information as (atom1_id, atom2_id, order, length).
    dynamics : dict, optional
        Dynamics-related metadata such as temperature, pressure, etc.
    metadata : dict, optional
        Arbitrary additional metadata stored in JSON format.

    Returns
    -------
    str
        Formatted string in .pygam file format.
    """

    lines = [
        "#========================================================",
        "# PYGAM STRUCTURE FILE (version 1.0)",
        "#========================================================",
        f"SYSTEM  {system}",
        f"NATOMS  {len(atoms)}",
        f"CELL    {' '.join(f'{c:.6f}' for c in (cell or [0.0, 0.0, 0.0]))}",
        "UNITS   Angstrom eV fs",
        f"TIMESTEP  {timestep}",
        "",
        "#--------------------------------------------------------",
        "# Atoms section",
        "# Format: ATOM id element x y z vx vy vz charge state mass",
        "#--------------------------------------------------------"
    ]

    for atom in atoms:
        vx, vy, vz = atom._velocity.tolist()
        lines.append(
            f"ATOM {atom.id} {atom.element} "
            f"{atom.x:.6f} {atom.y:.6f} {atom.z:.6f} "
            f"{vx:.6f} {vy:.6f} {vz:.6f} "
            f"{atom.charge:.3f} {atom.electronic_state.name} {atom.atomic_mass:.3f}"
        )

    if bonds:
        lines.extend([
            "",
            "#--------------------------------------------------------",
            "# Bonds section",
            "# Format: BOND atom1 atom2 order length",
            "#--------------------------------------------------------"
        ])
        for b1, b2, order, length in bonds:
            lines.append(f"BOND {b1} {b2} {order} {length:.3f}")

    if dynamics:
        lines.extend([
            "",
            "#--------------------------------------------------------",
            "# Dynamics section",
            "#--------------------------------------------------------"
        ])
        for key, value in dynamics.items():
            lines.append(f"{key.upper()}   {value}")

    if metadata:
        lines.extend([
            "",
            "#--------------------------------------------------------",
            "# Metadata block (JSON)",
            "#--------------------------------------------------------",
            "METADATA " + json.dumps(metadata, indent=2)
        ])

    lines.append("")
    lines.append("#========================================================")
    lines.append("# END OF STEP")
    lines.append("#========================================================")

    return "\n".join(lines)






def export_pygam(atoms: List[GAM_Atom],
                 file_path: str = "structure.pygam",
                 system: str = "GAM system",
                 cell: Optional[List[float]] = None,
                 timestep: int = 0,
                 bonds: Optional[List[Tuple[int, int, int, float]]] = None,
                 dynamics: Optional[Dict[str, float]] = None,
                 metadata: Optional[dict] = None):
    """
    Export a list of GAM_Atom objects to a .pygam file on disk.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of atoms to export.
    file_path : str
        Destination file path for the .pygam file.
    system : str
        Name or description of the system.
    cell : list of float, optional
        Unit cell dimensions [a, b, c] in Angstroms.
    timestep : int
        Simulation timestep (default is 0).
    bonds : list of tuples, optional
        Bond information as (atom1_id, atom2_id, order, length).
    dynamics : dict, optional
        Dynamics-related metadata such as temperature, pressure, etc.
    metadata : dict, optional
        Arbitrary additional metadata stored in JSON format.

    Notes
    -----
    This function writes the .pygam file directly to disk.
    """

    pygam_str = to_pygam(atoms, system, cell, timestep, bonds, dynamics, metadata)
    with open(file_path, "w") as f:
        f.write(pygam_str)
    print(f".pygam file successfully exported: {file_path}")




# ---------- XYZ ----------
def to_xyz(atoms: List[GAM_Atom], file_path: str = "structure.xyz", comment: str = "Generated by GAM_Atoms") -> None:
    """
    Export a list of GAM_Atom objects to the XYZ file format.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of atoms to export.
    file_path : str
        Destination XYZ file path.
    comment : str
        Optional comment line added after the atom count.
    """
    lines = [str(len(atoms))]
    lines.append(comment)
    for atom in atoms:
        lines.append(f"{atom.element} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}")
    
    with open(file_path, "w") as f:
        f.write("\n".join(lines))
    print(f"XYZ file successfully exported: {file_path}")


# ---------- CIF (very minimal) ----------
def to_cif(atoms: List[GAM_Atom], 
           file_path: str = "structure.cif",
           cell: Optional[List[List[float]]] = None, 
           name: str = "GAM_Structure") -> None:
    """
    Export a list of GAM_Atom objects to a CIF (Crystallographic Information File).

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of atoms to export.
    file_path : str
        Destination CIF file path.
    cell : list of lists, optional
        3x3 unit cell matrix. Required for proper fractional coordinates.
    name : str
        CIF data block name.
    
    Notes
    -----
    Coordinates are converted to fractional units if a cell is provided.
    """
    lines = [f"data_{name}"]
    if cell:
        a, b, c = cell
        lines += [
            f"_cell_length_a    {a[0]:.6f}",
            f"_cell_length_b    {b[1]:.6f}",
            f"_cell_length_c    {c[2]:.6f}",
            "_cell_angle_alpha 90.000",
            "_cell_angle_beta  90.000",
            "_cell_angle_gamma 90.000",
        ]
    lines += [
        "loop_",
        "_atom_site_label",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z"
    ]
    for atom in atoms:
        # assuming Cartesian coords and unit cell provided
        if cell:
            x = atom.x / cell[0][0]
            y = atom.y / cell[1][1]
            z = atom.z / cell[2][2]
        else:
            x, y, z = atom.x, atom.y, atom.z
        lines.append(f"{atom.element} {x:.6f} {y:.6f} {z:.6f}")
    
    with open(file_path, "w") as f:
        f.write("\n".join(lines))
    print(f"CIF file successfully exported: {file_path}")


# ---------- POSCAR (VASP) ----------
def to_poscar(atoms: List[GAM_Atom], 
              file_path: str = "POSCAR",
              cell: Optional[List[List[float]]] = None, 
              comment: str = "Generated by GAM_Atoms") -> None:
    """
    Export a list of GAM_Atom objects to VASP POSCAR format.

    Parameters
    ----------
    atoms : List[GAM_Atom]
        List of atoms to export.
    file_path : str
        Destination POSCAR file path.
    cell : list of lists
        3x3 lattice vectors (required).
    comment : str
        Optional comment line in the POSCAR file.

    Raises
    ------
    ValueError
        If `cell` is not provided.
    """
    if cell is None:
        raise ValueError("POSCAR export requires a cell (3x3 matrix).")

    # count elements
    counts = Counter(atom.element for atom in atoms)
    elements = list(counts.keys())
    numbers = list(counts.values())

    lines = [comment, "1.0"]
    # lattice
    for vec in cell:
        lines.append(f"{vec[0]:.6f} {vec[1]:.6f} {vec[2]:.6f}")
    # element line
    lines.append(" ".join(elements))
    # counts line
    lines.append(" ".join(str(n) for n in numbers))
    lines.append("Cartesian")

    for atom in atoms:
        lines.append(f"{atom.x:.6f} {atom.y:.6f} {atom.z:.6f}")
    
    with open(file_path, "w") as f:
        f.write("\n".join(lines))
    print(f"POSCAR file successfully exported: {file_path}")




# --------------------------------
def to_pdb(atoms: list[GAM_Atom], file_path: str = "structure.pdb", title: str="Generated by GAM_Atoms") -> None:
    """
    Export a list of GAM_Atom objects to the PDB file format.

    Parameters
    ----------
    atoms : list[GAM_Atom]
        List of atoms to export.
    file_path : str
        Destination PDB file path.
    title : str
        TITLE line in the PDB file.
    """
    lines = [f"TITLE     {title}"]
    for i, atom in enumerate(atoms, start=1):
        lines.append(
            f"ATOM  {i:5d} {atom.element:>2}   MOL     1    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00"
        )
    lines.append("END")
    
    with open(file_path, "w") as f:
        f.write("\n".join(lines))
    print(f"PDB file successfully exported: {file_path}")



def to_lammps(atoms: list[GAM_Atom], file_path: str = "data.lammps", cell: list[list[float]] = None) -> None:
    """
    Export a list of GAM_Atom objects to LAMMPS data format.

    Parameters
    ----------
    atoms : list[GAM_Atom]
        List of atoms to export.
    file_path : str
        Destination LAMMPS data file path.
    cell : list of lists, optional
        Simulation box dimensions (currently not used in example).
    """
    counts = len(atoms)
    lines = [f"{counts} atoms", "", "Atoms"]
    for i, atom in enumerate(atoms, start=1):
        # example: atom-ID molecule-ID atom-type q x y z
        lines.append(f"{i} 1 {i} {atom.charge:.3f} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}")
    
    with open(file_path, "w") as f:
        f.write("\n".join(lines))
    print(f"LAMMPS data file successfully exported: {file_path}")



def to_gaussian(atoms: list[GAM_Atom], 
                file_path: str = "structure.gjf",
                charge: int=0, 
                multiplicity: int=1,
                route: str="#P HF/6-31G(d) OPT") -> None:
    """
    Export a list of GAM_Atom objects to Gaussian input format (.gjf).

    Parameters
    ----------
    atoms : list[GAM_Atom]
        List of atoms to export.
    file_path : str
        Destination Gaussian input file path.
    charge : int
        Total charge of the system (default: 0).
    multiplicity : int
        Spin multiplicity of the system (default: 1).
    route : str
        Gaussian route section specifying method and keywords.
    """
    lines = [route, "", "Generated by GAM_Atoms", ""]
    lines.append(f"{charge} {multiplicity}")
    for atom in atoms:
        lines.append(f"{atom.element} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}")
    lines.append("")
    
    with open(file_path, "w") as f:
        f.write("\n".join(lines))
    print(f"Gaussian input file successfully exported: {file_path}")









