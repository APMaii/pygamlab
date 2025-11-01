from ..Primatom import GAM_Atom , GAM_Bond
import numpy as np
from ase import Atoms
from typing import Tuple, Optional, Dict, Union, Literal , List
import warnings
from ase.visualize import view
from ase.io import write
import matplotlib.pyplot as plt
from ase.io.pov import get_bondpairs
import copy
from ase.build import graphene
from ase.build import mx2


from ase.atoms import Atoms
import copy

try:
    import ase
    from ase import Atoms
    HAS_ASE = True
except ImportError:
    HAS_ASE = False





LAYERED_MATERIALS = {
    # Group 6 TMDs (most stable in 2H phase)
    'mos2': {
        'formula': 'MoS2',
        'lattice_constant': 3.160,
        'thickness': 6.15,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Molybdenum disulfide - semiconductor'
    },
    'mose2': {
        'formula': 'MoSe2',
        'lattice_constant': 3.289,
        'thickness': 6.46,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Molybdenum diselenide - semiconductor'
    },
    'mote2': {
        'formula': 'MoTe2',
        'lattice_constant': 3.518,
        'thickness': 6.97,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Molybdenum ditelluride - semiconductor'
    },
    'ws2': {
        'formula': 'WS2',
        'lattice_constant': 3.153,
        'thickness': 6.18,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Tungsten disulfide - semiconductor'
    },
    'wse2': {
        'formula': 'WSe2',
        'lattice_constant': 3.282,
        'thickness': 6.49,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Tungsten diselenide - semiconductor'
    },
    'wte2': {
        'formula': 'WTe2',
        'lattice_constant': 3.496,
        'thickness': 7.04,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Tungsten ditelluride - semimetal'
    },
    
    # Group 5 TMDs (can be metallic in 1T phase)
    'vs2': {
        'formula': 'VS2',
        'lattice_constant': 3.160,
        'thickness': 5.76,
        'preferred_phase': '1T',
        'stable_phases': ['1T', '2H'],
        'description': 'Vanadium disulfide - metallic'
    },
    'vse2': {
        'formula': 'VSe2',
        'lattice_constant': 3.355,
        'thickness': 6.10,
        'preferred_phase': '1T',
        'stable_phases': ['1T', '2H'],
        'description': 'Vanadium diselenide - metallic'
    },
    'nbs2': {
        'formula': 'NbS2',
        'lattice_constant': 3.313,
        'thickness': 5.985,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Niobium disulfide - metallic'
    },
    'nbse2': {
        'formula': 'NbSe2',
        'lattice_constant': 3.445,
        'thickness': 6.26,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Niobium diselenide - superconductor'
    },
    'tas2': {
        'formula': 'TaS2',
        'lattice_constant': 3.314,
        'thickness': 5.895,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Tantalum disulfide - metallic/CDW'
    },
    'tase2': {
        'formula': 'TaSe2',
        'lattice_constant': 3.436,
        'thickness': 6.25,
        'preferred_phase': '2H',
        'stable_phases': ['2H', '1T'],
        'description': 'Tantalum diselenide - superconductor'
    },
    
    # Group 4 TMDs
    'tis2': {
        'formula': 'TiS2',
        'lattice_constant': 3.407,
        'thickness': 5.695,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Titanium disulfide - metallic'
    },
    'tise2': {
        'formula': 'TiSe2',
        'lattice_constant': 3.535,
        'thickness': 6.008,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Titanium diselenide - semimetal/CDW'
    },
    'zrs2': {
        'formula': 'ZrS2',
        'lattice_constant': 3.658,
        'thickness': 5.83,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Zirconium disulfide - semiconductor'
    },
    'zrse2': {
        'formula': 'ZrSe2',
        'lattice_constant': 3.77,
        'thickness': 6.13,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Zirconium diselenide - semiconductor'
    },
    'hfs2': {
        'formula': 'HfS2',
        'lattice_constant': 3.615,
        'thickness': 5.85,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Hafnium disulfide - semiconductor'
    },
    'hfse2': {
        'formula': 'HfSe2',
        'lattice_constant': 3.743,
        'thickness': 6.15,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Hafnium diselenide - semiconductor'
    },
    
    # Group 7 TMDs
    'res2': {
        'formula': 'ReS2',
        'lattice_constant': 3.144,
        'thickness': 6.05,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Rhenium disulfide - anisotropic semiconductor'
    },
    'rese2': {
        'formula': 'ReSe2',
        'lattice_constant': 3.300,
        'thickness': 6.34,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Rhenium diselenide - anisotropic semiconductor'
    },
    
    # Group 10 TMDs
    'pts2': {
        'formula': 'PtS2',
        'lattice_constant': 3.540,
        'thickness': 5.036,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Platinum disulfide - semiconductor'
    },
    'ptse2': {
        'formula': 'PtSe2',
        'lattice_constant': 3.728,
        'thickness': 5.081,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Platinum diselenide - semimetal'
    },
    
    # Metal halides
    'cri3': {
        'formula': 'CrI3',
        'lattice_constant': 6.867,
        'thickness': 6.60,
        'preferred_phase': '2H',
        'stable_phases': ['2H'],
        'description': 'Chromium triiodide - 2D ferromagnet'
    },
    'crbr3': {
        'formula': 'CrBr3',
        'lattice_constant': 6.30,
        'thickness': 6.12,
        'preferred_phase': '2H',
        'stable_phases': ['2H'],
        'description': 'Chromium tribromide - 2D ferromagnet'
    },
    'vi3': {
        'formula': 'VI3',
        'lattice_constant': 6.86,
        'thickness': 6.74,
        'preferred_phase': '2H',
        'stable_phases': ['2H'],
        'description': 'Vanadium triiodide - 2D ferromagnet'
    },
    
    # Other layered materials
    'bi2se3': {
        'formula': 'Bi2Se3',  # Simplified for mx2 function
        'lattice_constant': 4.138,
        'thickness': 9.54,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Bismuth selenide - topological insulator'
    },
    'bi2te3': {
        'formula': 'Bi2Te3',  # Simplified for mx2 function
        'lattice_constant': 4.386,
        'thickness': 10.15,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Bismuth telluride - topological insulator'
    },
    'sb2te3': {
        'formula': 'Sb2Te3',  # Simplified for mx2 function
        'lattice_constant': 4.264,
        'thickness': 10.33,
        'preferred_phase': '1T',
        'stable_phases': ['1T'],
        'description': 'Antimony telluride - topological insulator'
    }
}



# Material parameters database - all designed to work with the graphene function
NANOSHEET_MATERIALS = {
    'graphene': {
        'formula': 'C2',
        'lattice_constant': 2.460,
        'thickness': 0.0,
        'description': 'Pure graphene - carbon honeycomb lattice'
    },
    'silicene': {
        'formula': 'Si2', 
        'lattice_constant': 3.86,
        'thickness': 0.44,
        'description': 'Silicon analogue of graphene with buckling'
    },
    'germanene': {
        'formula': 'Ge2',
        'lattice_constant': 4.02, 
        'thickness': 0.64,
        'description': 'Germanium analogue of graphene with buckling'
    },
    'stanene': {
        'formula': 'Sn2',
        'lattice_constant': 4.67,
        'thickness': 0.85,
        'description': 'Tin analogue of graphene with buckling'
    },
    'plumbene': {
        'formula': 'Pb2',
        'lattice_constant': 4.97,
        'thickness': 0.95,
        'description': 'Lead analogue of graphene with buckling'
    },
    'boron_nitride': {
        'formula': 'BN',
        'lattice_constant': 2.504,
        'thickness': 0.0,
        'description': 'Hexagonal boron nitride (white graphene)'
    },
    'aluminum_nitride': {
        'formula': 'AlN',
        'lattice_constant': 3.111,
        'thickness': 0.0,
        'description': 'Aluminum nitride honeycomb structure'
    },
    'gallium_nitride': {
        'formula': 'GaN', 
        'lattice_constant': 3.189,
        'thickness': 0.0,
        'description': 'Gallium nitride honeycomb structure'
    },
    'phosphorene': {
        'formula': 'P2',  # Simplified to work with graphene function
        'lattice_constant': 3.314,
        'thickness': 1.05,  # Half the actual buckling for the graphene function
        'description': 'Black phosphorus (simplified honeycomb approximation)'
    },
    'arsenene': {
        'formula': 'As2',
        'lattice_constant': 3.608,
        'thickness': 1.4,
        'description': 'Arsenic analogue of phosphorene'
    },
    'antimonene': {
        'formula': 'Sb2', 
        'lattice_constant': 4.120,
        'thickness': 1.65,
        'description': 'Antimony analogue of phosphorene'
    },
    'bismuthene': {
        'formula': 'Bi2',
        'lattice_constant': 4.54,
        'thickness': 1.86,
        'description': 'Bismuth analogue of phosphorene'
    },
    'beryllium_oxide': {
        'formula': 'BeO',
        'lattice_constant': 2.698,
        'thickness': 0.0,
        'description': 'Beryllium oxide honeycomb structure'
    },
    'zinc_oxide': {
        'formula': 'ZnO',
        'lattice_constant': 3.249,
        'thickness': 0.6,
        'description': 'Zinc oxide honeycomb structure with buckling'
    },
    'cadmium_sulfide': {
        'formula': 'CdS',
        'lattice_constant': 4.136,
        'thickness': 0.8,
        'description': 'Cadmium sulfide honeycomb structure'
    }
}



class Nano_TwoD_Builder:

    """
    Comprehensive builder for 2D nanomaterials such as nanosheets, monolayers, and van der Waals heterostructures.

    The `Nano_TwoD_Builder` class provides an advanced atomistic modeling framework 
    for constructing two-dimensional (2D) nanostructures, including:

        - **Nanosheets:** elemental or binary 2D lattices (e.g., graphene, silicene, phosphorene)
        - **Transition Metal Dichalcogenides (TMDs):** MX₂ layered compounds with tunable 1T/2H phases
        - **Heterostructures:** stacked multilayer systems with controlled interlayer spacing

    The builder supports both *honeycomb-like* (graphene family) and *layered* (TMD) topologies.
    It integrates with ASE’s `graphene()` and `mx2()` generators and enables
    direct conversion to GAM-compatible formats for downstream simulation pipelines.

    Features
    --------
    - Generate 2D nanosheets and TMD monolayers
    - Build van der Waals heterostructures from arbitrary combinations of layers
    - Support for custom formulas, lattice constants, and layer thicknesses
    - Predefined constructors for common 2D materials (graphene, MoS₂, WS₂, etc.)
    - ASE-based export and visualization (XYZ, CIF, PDB)
    - Structure translation and rotation operations
    - GAM-compatible atom and bond representation

    Parameters
    ----------
    material : str
        Name or symbol of the 2D material (e.g., "graphene", "MoS2", "BN").
    structure_type : {'nanosheet', 'tmd'}
        Type of 2D nanostructure to build.

    Attributes
    ----------
    ASE_atoms : ase.Atoms
        Underlying ASE structure representing the generated material.
    atoms : list[GAM_Atom]
        List of atoms in GAM-compatible format.
    bonds : list[GAM_Bond]
        Optional list of identified bonds.
    layered_materials : dict
        Database of transition metal dichalcogenides (TMDs) and related 2D materials.
    nanosheet_materials : dict
        Database of single-element or binary 2D sheet materials.
    meta : dict
        Metadata dictionary containing build parameters and structure details.
    structure_type : str
        Selected structure type ('nanosheet' or 'tmd').
    material : str
        Material name used for generation.

    Methods
    -------
    create_layered_material(phase='2H', size=(10,10,1), vacuum=10.0, ...)
        Build a TMD or layered 2D material (e.g., MoS₂, WSe₂) using ASE’s `mx2()` function.
    create_nanosheet(size=(10,10,1), vacuum=10.0, ...)
        Build a honeycomb-like nanosheet (e.g., graphene, silicene) using ASE’s `graphene()` function.
    create_layered_heterostructure(materials_and_phases, layer_spacing=6.5, ...)
        Stack multiple TMD layers with custom interlayer spacing.
    create_nanosheet_heterostructure(materials, layer_spacing=3.4, ...)
        Build van der Waals heterostructures combining 2D sheets (e.g., graphene/hBN).
    create_custom_sheet(formula, lattice_constant, thickness=0.0, ...)
        Define a custom 2D lattice from scratch.
    list_layered_materials(category=None)
        Display all available layered (TMD) materials in the internal database.
    list_nanosheet_materials()
        Display all available 2D sheet materials and their physical properties.
    translate(dx, dy, dz=0.0)
        Translate atomic coordinates in 3D space.
    rotate(angle_deg, about_center=True)
        Rotate the sheet around the z-axis (in-plane rotation).
    copy()
        Deep copy the atomic structure and metadata.
    get_atoms()
        Return a deep copy of all atoms as GAM_Atom objects.
    get_positions()
        Retrieve atomic positions as (x, y, z) tuples.
    get_elements()
        Return list of atomic symbols.
    to_xyz(filename)
        Export structure to `.xyz` file format.
    to_ase()
        Convert to an ASE `Atoms` object for external use.

    Notes
    -----
    - The builder automatically selects between `graphene()` and `mx2()` constructors 
      based on the specified `structure_type`.
    - Both nanosheet and layered materials can be stacked into heterostructures
      with adjustable vacuum spacing and interlayer distances.
    - Built-in materials include metallic, semiconducting, and superconducting 2D systems:
      MoS₂, WS₂, MoSe₂, WSe₂, NbSe₂, TaSe₂, PtS₂, CrI₃, graphene, silicene, BN, phosphorene, etc.
    - The `meta` dictionary records all parameters used in construction, 
      enabling reproducibility and dataset traceability.
    - All units are expressed in Ångströms unless otherwise noted.

    Examples
    --------
    >>> from pygamlab import Nano_TwoD_Builder

    >>> # 1. Generate a graphene sheet
    >>> sheet = Nano_TwoD_Builder(material="graphene", structure_type="nanosheet")
    >>> sheet.to_xyz("graphene.xyz")

    >>> # 2. Build a MoS2 monolayer in the 2H phase
    >>> tmd = Nano_TwoD_Builder(material="MoS2", structure_type="tmd")
    >>> tmd.create_layered_material(phase='2H', size=(5,5,1), vacuum=15.0)
    >>> tmd.to_xyz("MoS2_2H.xyz")

    >>> # 3. Construct a van der Waals heterostructure: Graphene/hBN
    >>> het = Nano_TwoD_Builder(material="graphene", structure_type="nanosheet")
    >>> het.create_nanosheet_heterostructure(["graphene", "boron_nitride"], layer_spacing=3.35)
    >>> het.to_xyz("graphene_hbn_heterostructure.xyz")

    >>> # 4. Create a custom 2D material
    >>> custom = Nano_TwoD_Builder(material="CustomX", structure_type="nanosheet")
    >>> custom.create_custom_sheet(formula="X2", lattice_constant=3.6, thickness=0.2)
    >>> custom.to_xyz("CustomX_sheet.xyz")
    """

    def __init__(self,material: str, structure_type:str):

        
        self.atoms: List[GAM_Atom] = []
        self.bonds: List[GAM_Bond] = []
        

        self.structure_type = structure_type
        self.material = material

        self.meta={'material':material,
        'structure_type': structure_type}




        self.layered_materials = LAYERED_MATERIALS.copy()
        self.nanosheet_materials = NANOSHEET_MATERIALS.copy()

        if structure_type.lower() == "nanosheet":
            self.ASE_atoms=self.create_nanosheet()
            self.atoms=self._to_GAM_Atoms(self.ASE_atoms) 


        elif structure_type.lower() == "tmd":
            self.ASE_atoms=self.create_layered_material()
            self.atoms=self._to_GAM_Atoms(self.ASE_atoms) 

        else:
            raise ValueError("Invalid structure type. Must be 'nanosheet', 'nanoribbon', or 'nanotube'.")



    def create_layered_material(self,phase: Optional[str] = None,
                              size: Tuple[int, int, int] = (10,10, 1),
                              vacuum: Optional[float] = 10.0,
                              custom_formula: Optional[str] = None,
                              custom_lattice_constant: Optional[float] = None,
                              custom_thickness: Optional[float] = None) -> Atoms:
        
        


        
        """
        Create a layered material using the ASE mx2 function.
        
        Parameters
        ----------
        material : str
            Material name from database
        phase : str, optional
            Crystal phase ('2H' or '1T'). Uses preferred phase if None.
        size : tuple
            Supercell size (nx, ny, nz)
        vacuum : float, optional
            Vacuum spacing in z-direction
        custom_formula : str, optional
            Override default formula
        custom_lattice_constant : float, optional
            Override default lattice constant
        custom_thickness : float, optional
            Override default thickness
            
        Returns
        -------
        Atoms
            Layered material structure
        """
        material_key = self.material.lower().replace(' ', '').replace('-', '')
        
        if material_key not in self.layered_materials:
            available = ', '.join(self.layered_materials.keys())
            raise ValueError(f"Material '{self.material}' not found. Available: {available}")
        
        props = self.layered_materials[material_key]
        
        # Determine phase
        if phase is None:
            phase = props['preferred_phase']
        elif phase not in props['stable_phases']:
            stable = ', '.join(props['stable_phases'])
            warnings.warn(f"Phase '{phase}' may not be stable for {self.material}. "
                         f"Stable phases: {stable}")
        
        # Get parameters (with custom overrides)
        formula = custom_formula or props['formula']
        a = custom_lattice_constant or props['lattice_constant']
        thickness = custom_thickness or props['thickness']
        
        # Create structure using mx2 function
        atoms = mx2(
            formula=formula,
            kind=phase,
            a=a,
            thickness=thickness,
            size=size,
            vacuum=vacuum
        )

        self.meta={
            'material':self.material,
            'structure_type': self.structure_type,
            'phase':phase,
            'size':size,
            'vacuum': vacuum,
            'formula': custom_formula,
            'lattice_constant': custom_lattice_constant,
            'thickness' : custom_thickness,
        }

        return atoms
    
    # Convenience methods for common TMDs
    def mos2_sheet(self, phase='2H', size=(1, 1, 1), vacuum=10.0):
        """Create MoS2 sheet."""
        
        if self.material.lower()!='mos2':
            raise ValueError(f"Your material is {self.material} and you want to create mos2, build this class again.")
            
        #self.ASE_atoms = self.create_layered_material('mos2', phase=phase, size=size, vacuum=vacuum)
        self.ASE_atoms = self.create_layered_material(phase=phase, size=size, vacuum=vacuum)
        
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('mos2', phase=phase, size=size, vacuum=vacuum)
    
    def ws2_sheet(self, phase='2H', size=(1, 1, 1), vacuum=10.0):
        if self.material.lower()!='ws2':
            raise ValueError(f"Your material is {self.material} and you want to create ws2, build this class again.")
            
        """Create WS2 sheet."""
        self.ASE_atoms = self.create_layered_material( phase=phase, size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('ws2', phase=phase, size=size, vacuum=vacuum)
    
    def mose2_sheet(self, phase='2H', size=(1, 1, 1), vacuum=10.0):
        if self.material.lower()!='mose2':
            raise ValueError(f"Your material is {self.material} and you want to create mose2, build this class again.")
        """Create MoSe2 sheet."""
        self.ASE_atoms = self.create_layered_material(phase=phase, size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('mose2', phase=phase, size=size, vacuum=vacuum)
    
    def wse2_sheet(self, phase='2H', size=(1, 1, 1), vacuum=10.0):
        if self.material.lower()!='wse2':
            raise ValueError(f"Your material is {self.material} and you want to create wse2, build this class again.")
        """Create WSe2 sheet."""
        self.ASE_atoms = self.create_layered_material(phase=phase, size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('wse2', phase=phase, size=size, vacuum=vacuum)
    
    def nbse2_sheet(self, phase='2H', size=(1, 1, 1), vacuum=10.0):
        if self.material.lower()!='nbse2':
            raise ValueError(f"Your material is {self.material} and you want to create nbse2, build this class again.")
        """Create NbSe2 sheet (superconductor)."""
        self.ASE_atoms = self.create_layered_material(phase=phase, size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('nbse2', phase=phase, size=size, vacuum=vacuum)
    
    def tase2_sheet(self, phase='2H', size=(1, 1, 1), vacuum=10.0):
        if self.material.lower()!='tase2':
            raise ValueError(f"Your material is {self.material} and you want to create tase2, build this class again.")
        """Create TaSe2 sheet (superconductor)."""
        self.ASE_atoms = self.create_layered_material(phase=phase, size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('tase2', phase=phase, size=size, vacuum=vacuum)
    
    def pts2_sheet(self, phase='1T', size=(1, 1, 1), vacuum=10.0):
        if self.material.lower()!='pts2':
            raise ValueError(f"Your material is {self.material} and you want to create pts2, build this class again.")
        """Create PtS2 sheet (semiconductor)."""
        self.ASE_atoms = self.create_layered_material(phase, size, vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('pts2', phase=phase, size=size, vacuum=vacuum)
    
    def cri3_sheet(self, phase='2H', size=(1, 1, 1), vacuum=10.0):
        if self.material.lower()!='cris3':
            raise ValueError(f"Your material is {self.material} and you want to create cris3, build this class again.")
        """Create CrI3 sheet (2D ferromagnet)."""
        self.ASE_atoms = self.create_layered_material(phase=phase, size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)
        #return self.atoms
        #return self.create_layered_material('cri3', phase=phase, size=size, vacuum=vacuum)
    


    
    def create_layered_heterostructure(self,
                             materials_and_phases: list,
                             layer_spacing: float = 6.5,
                             size: Tuple[int, int, int] = (1, 1, 1),
                             vacuum: Optional[float] = 10.0) -> GAM_Atom:
        """
        Create a heterostructure from multiple layered materials.
        
        Parameters
        ----------
        materials_and_phases : list
            List of tuples: [(material, phase), ...] or list of material names
        layer_spacing : float
            Distance between layers
        size : tuple
            Supercell size
        vacuum : float, optional
            Total vacuum spacing
            
        Returns
        -------
        Atoms
            Heterostructure
        """

        

        if len(materials_and_phases) < 2:
            raise ValueError("Need at least 2 materials for heterostructure")
        
        layers = []
        z_offsets = []
        current_z = 0
        
        for i, item in enumerate(materials_and_phases):
            # Handle both (material, phase) tuples and simple material names
            if isinstance(item, (list, tuple)) and len(item) == 2:
                material, phase = item
            else:
                material, phase = item, None
            
            # Create layer
            self.material=material
            layer = self.create_layered_material(phase=phase, size=size, vacuum=None)
            layers.append(layer)
            z_offsets.append(current_z)
            
            # Calculate next layer position
            if i < len(materials_and_phases) - 1:
                material_key = material.lower().replace(' ', '').replace('-', '')
                layer_thickness = self.layered_materials[material_key]['thickness']
                current_z += layer_thickness + layer_spacing
        
        # Combine layers
        all_positions = []
        all_symbols = []
        
        for layer, z_offset in zip(layers, z_offsets):
            positions = layer.get_positions()
            positions[:, 2] += z_offset
            all_positions.extend(positions)
            all_symbols.extend(layer.get_chemical_symbols())
        
        # Create combined structure
        cell = layers[0].cell.copy()
        total_height = current_z + (vacuum if vacuum else 15.0)
        cell[2, 2] = total_height
        
        heterostructure = Atoms(
            symbols=all_symbols,
            positions=all_positions,
            cell=cell,
            pbc=(1, 1, 0)
        )
        
        if vacuum:
            heterostructure.center(vacuum, axis=2)


        self.ASE_atoms=heterostructure
        
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)




        self.meta={'materials_and_phases':materials_and_phases,
                   'structure_type':self.structure_type,
                   'layer_spacing':layer_spacing,
                   'size':size,
                   'vacuum':vacuum}
        
        
        #return self.atoms
    
    
    
    def list_layered_materials(self, category: Optional[str] = None) -> None:
        """
        Display available materials.
        
        Parameters
        ----------
        category : str, optional
            Filter by category ('tmd', 'halide', 'topological', etc.)
        """
        print("Available Layered Materials (using ASE mx2 function):")
        print("=" * 65)
        
        for name, props in self.layered_materials.items():
            if category and category.lower() not in props['description'].lower():
                continue
                
            print(f"{name.upper()}:")
            print(f"  Formula: {props['formula']}")
            print(f"  Lattice constant: {props['lattice_constant']:.3f} Å")
            print(f"  Thickness: {props['thickness']:.2f} Å")
            print(f"  Preferred phase: {props['preferred_phase']}")
            print(f"  Stable phases: {', '.join(props['stable_phases'])}")
            print(f"  Description: {props['description']}")
            print()
    


    



   #=====================      ====================   
   #=====================     ====================  
   #=====================     ====================  
   #===================== Nnaosheet 
   #=====================     ====================  
   #=====================     ====================  
   #======================    ====================  
    def create_nanosheet(self, 
                        size: Tuple[int, int, int] = (10, 10, 1),
                        vacuum: Optional[float] = 10.0,
                        custom_formula: Optional[str] = None,
                        custom_lattice_constant: Optional[float] = None,
                        custom_thickness: Optional[float] = None) -> Atoms:
        

        """
        Create a nanosheet using the ASE graphene function with modified parameters.
        
        Parameters
        ----------
        material : str
            Material name from the database
        size : tuple
            Supercell size (nx, ny, nz)
        vacuum : float, optional
            Vacuum spacing in z-direction
        custom_formula : str, optional
            Override the default formula
        custom_lattice_constant : float, optional
            Override the default lattice constant
        custom_thickness : float, optional
            Override the default thickness
            
        Returns
        -------
        Atoms
            ASE Atoms object of the nanosheet
        """
        material_key = self.material.lower().replace(' ', '_').replace('-', '_')
        
        if material_key not in self.nanosheet_materials:
            available = ', '.join(self.nanosheet_materials.keys())
            raise ValueError(f"Material '{self.material}' not found. Available: {available}")
        
        # Get material properties
        props = self.nanosheet_materials[material_key]
        
        # Use custom parameters if provided
        formula = custom_formula or props['formula']
        a = custom_lattice_constant or props['lattice_constant'] 
        thickness = custom_thickness or props['thickness']
        
        # Call the original graphene function with modified parameters
        atoms = graphene(
            formula=formula,
            a=a,
            thickness=thickness,
            size=size,
            vacuum=vacuum
        )

        self.meta={'material':self.material,
                   'structure_type':self.structure_type,
                   'size':size,
                   'vacuum': vacuum,
                   'formula':custom_formula,
                   'lattice_constant':custom_lattice_constant,
                   'thickness':custom_thickness}
        
        #self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        return atoms
    
    def graphene_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create graphene nanosheet."""
        self.ASE_atoms=self.create_nanosheet('graphene', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        
        #return self.atoms
    
    def silicene_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create silicene nanosheet.""" 
        self.ASE_atoms=self.create_nanosheet('silicene', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        #return self.atoms   
    
    def germanene_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create germanene nanosheet."""
        self.ASE_atoms= self.create_nanosheet('germanene', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        #return self.atoms 
    
    def stanene_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create stanene nanosheet."""
        self.ASE_atoms= self.create_nanosheet('stanene', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        #return self.atoms 
    
    def boron_nitride_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create hexagonal boron nitride nanosheet."""
        self.ASE_atoms= self.create_nanosheet('boron_nitride', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        #return self.atoms 
    
    def phosphorene_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create phosphorene nanosheet (honeycomb approximation)."""
        self.ASE_atoms= self.create_nanosheet('phosphorene', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        #return self.atoms 
    
    def arsenene_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create arsenene nanosheet."""
        self.ASE_atoms= self.create_nanosheet('arsenene', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        #return self.atoms 
    
    def antimonene_sheet(self, size=(1, 1, 1), vacuum=10.0):
        """Create antimonene nanosheet."""
        self.ASE_atoms= self.create_nanosheet('antimonene', size=size, vacuum=vacuum)
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        #return self.atoms 
    
    def create_custom_sheet(self,
                          formula: str,
                          lattice_constant: float,
                          thickness: float = 0.0,
                          size: Tuple[int, int, int] = (1, 1, 1),
                          vacuum: Optional[float] = 10.0) -> GAM_Atom:
        """
        Create a custom nanosheet with specified parameters.
        
        Parameters
        ----------
        formula : str
            Chemical formula (e.g., 'C2', 'Si2', 'BN')
        lattice_constant : float
            Lattice constant in Angstrom
        thickness : float
            Thickness/buckling parameter
        size : tuple
            Supercell size
        vacuum : float, optional
            Vacuum spacing
            
        Returns
        -------
        Atoms
            Custom nanosheet structure
        """
        self.ASE_atoms=graphene(
            formula=formula,
            a=lattice_constant,
            thickness=thickness,
            size=size,
            vacuum=vacuum
        )
        
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)



        #return self.atoms
        
    
    
    def create_nanosheet_heterostructure(self,
                             materials: list,
                             layer_spacing: float = 3.4,
                             size: Tuple[int, int, int] = (1, 1, 1),
                             vacuum: Optional[float] = 10.0) -> GAM_Atom:
        """
        Create a heterostructure by stacking different materials.
        
        Parameters
        ----------
        materials : list
            List of material names to stack
        layer_spacing : float
            Distance between layers in Angstrom
        size : tuple
            Supercell size for each layer
        vacuum : float, optional
            Total vacuum spacing
            
        Returns
        -------
        Atoms
            Heterostructure combining multiple materials
        """
        if len(materials) < 2:
            raise ValueError("Need at least 2 materials for heterostructure")
        
        layers = []
        z_offsets = []
        current_z = 0
        
        for i, material in enumerate(materials):
            # Create layer without vacuum
            layer = self.create_nanosheet(material, size=size, vacuum=None)
            layers.append(layer)
            z_offsets.append(current_z)
            
            # Calculate next layer position
            if i < len(materials) - 1:
                layer_thickness = self.nanosheet_materials[material.lower().replace(' ', '_').replace('-', '_')]['thickness']
                current_z += max(layer_thickness, 0.1) + layer_spacing
        
        # Combine all layers
        all_positions = []
        all_symbols = []
        
        for layer, z_offset in zip(layers, z_offsets):
            positions = layer.get_positions()
            positions[:, 2] += z_offset
            all_positions.extend(positions)
            all_symbols.extend(layer.get_chemical_symbols())
        
        # Create new cell with appropriate z-dimension
        cell = layers[0].cell.copy()
        total_height = current_z + (vacuum if vacuum else 10.0)
        cell[2, 2] = total_height
        
        heterostructure = Atoms(
            symbols=all_symbols,
            positions=all_positions,
            cell=cell,
            pbc=(1, 1, 0)
        )
        
        if vacuum:
            heterostructure.center(vacuum, axis=2)

        self.ASE_atoms = heterostructure
        
        self.atoms=self._to_GAM_Atoms(self.ASE_atoms)

        self.meta = {
            'materials': materials,
            'structure_type': self.structure_type,
            'layer_spacing': layer_spacing,
            'size': size,
            'vacuum': vacuum
        }
        
        #return self.atoms
    

    def list_nanosheet_materials(self) -> None:
        """Display all available materials and their properties."""
        print("Available 2D Materials (using ASE graphene function):")
        print("=" * 60)
        
        for name, props in self.nanosheet_materials.items():
            print(f"{name.replace('_', ' ').title()}:")
            print(f"  Formula: {props['formula']}")
            print(f"  Lattice constant: {props['lattice_constant']:.3f} Å")
            print(f"  Thickness/Buckling: {props['thickness']:.3f} Å")
            print(f"  Description: {props['description']}")
            print()
            
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
    
    
    
   









