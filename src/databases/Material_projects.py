
import json
from typing import List, Optional
from mp_api.client import MPRester



    
class MaterialsProject_Explorer:
    """
    A unified interface to search and fetch various material properties from the Materials Project.

    This class allows searching by elements, chemical formula, or specific MP-IDs, and fetching
    electronic, mechanical, thermodynamic, and structural properties.

    Attributes
    ----------
    api_key : str, optional
        API key for accessing Materials Project (if required).
    mpr : MPRester
        MPRester instance for querying the API.
    results : dict
        Stores the latest fetched data for a material.

    Examples
    --------
    >>> explorer = MaterialsProject_Explorer(api_key="YOUR_API_KEY")
    >>> results = explorer.search_materials(formula="MoS2")
    >>> explorer.fetch_electronic_properties("mp-123")
    {'mp-id': 'mp-123', 'band_gap': 1.23, ...}
    >>> all_data = explorer.fetch_all_data("mp-123")
    >>> all_data['electronic_prop']['band_gap']
    1.23
    >>> explorer.close()
    """

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize the MaterialsProject_Explorer.

        Parameters
        ----------
        api_key : str, optional
            API key for Materials Project access. If None, assumes local access or environment key.
        """
        self.api_key = api_key
        self.mpr = MPRester(api_key) if api_key else MPRester()
        self.results = {
            'electronic_prop': {},
            'mechanical_prop': {},
            'thermo_prop': {},
            'structure': {},
            'meta_data': {}
        }

    def search_materials(self, elements: Optional[List[str]] = None, formula: Optional[str] = None,
                         mp_ids: Optional[List[str]] = None, max_results: int = 5,description=True) -> List[str]:
        """
        Search for materials by elements, chemical formula, or MP-IDs.

        Parameters
        ----------
        elements : list[str], optional
            List of chemical elements to search for (e.g., ["Mo", "S"]).
        formula : str, optional
            Chemical formula to search for (e.g., "MoS2").
        mp_ids : list[str], optional
            List of Materials Project IDs to search directly.
        max_results : int, default 5
            Maximum number of results to return.
        description : bool, default True
            If True, returns a list of dictionaries with basic material info.
            If False, returns only a list of material IDs.

        Returns
        -------
        list[str] or list[dict]
            List of material IDs or detailed descriptions for each result.

        Example
        -------
        >>> explorer = MaterialsProject_Explorer()
        >>> explorer.search_materials(formula="MoS2")
        [{'mp-id': 'mp-123', 'system': 'Mo-S', 'band_gap (eV)': 1.23, ...}]
        """
        material_ids = []

        if mp_ids:
            material_ids.extend(mp_ids)
            
            return list(set(material_ids))

        if formula:
            docs = self.mpr.materials.summary.search(formula=formula, num_chunks=1, chunk_size=max_results)
            material_ids.extend([d.material_id for d in docs])
            
            if description==True:
                rows = []
                for d in docs:
                    rows.append({
                        "mp-id": d.material_id.string,
                        "system": d.chemsys,
                        "band_gap (eV)": getattr(d, "band_gap", None),
                        "Ehull (eV/atom)": getattr(d, "energy_above_hull", None),
                        "magnetism": getattr(d, "magnetic_ordering", None),
                        "theoretical": getattr(d, "theoretical", None)
                    })
                
                return rows
                    

        if elements:
            docs = self.mpr.materials.summary.search(elements=elements, num_chunks=1, chunk_size=max_results)
            material_ids.extend([d.material_id for d in docs])
            
            if description==True:
                rows = []
                for d in docs:
                    rows.append({
                        "mp-id": d.material_id.string,
                        "system": d.chemsys,
                        "band_gap (eV)": getattr(d, "band_gap", None),
                        "Ehull (eV/atom)": getattr(d, "energy_above_hull", None),
                        "magnetism": getattr(d, "magnetic_ordering", None),
                        "theoretical": getattr(d, "theoretical", None)
                    })
                
                return rows
                



    def fetch_electronic_properties(self, mp_id: str) -> dict:

        """
        Fetch electronic properties for a given Materials Project ID.

        Parameters
        ----------
        mp_id : str
            Materials Project ID (e.g., "mp-123").

        Returns
        -------
        dict
            Electronic properties including band structure, density of states, bandgap, etc.

        Example
        -------
        >>> explorer.fetch_electronic_properties("mp-123")
        {'mp-id': 'mp-123', 'band_gap': 1.23, 'DOS': {...}, ...}
        """
        try:
            
            electronic_results={}
            
            electronic_results['mp-id']=mp_id
            
            
            
            docs = self.mpr.materials.electronic_structure.search(material_ids=[mp_id])
            doc=docs[0].dict()
            for key,value in doc.items():
                electronic_results[key]=value
                
        except Exception:
            electronic_results['electronic_structure']=None
                
        try:
            
            
            dos = self.mpr.get_dos_by_material_id(mp_id)

            dos_dict = dos.dict()
            for key,value in dos_dict.items():
                electronic_results[key]=value
                
        except Exception:
            electronic_results['DOS']=None
        
        try:

            bs = self.mpr.get_bandstructure_by_material_id(mp_id)
            bs_dict=bs.as_dict()
            for key,value in bs_dict.items():
                electronic_results[key]=value

        except Exception:
            electronic_results['bandstructure']=None
                

        self.results['electronic_prop']=electronic_results
            

        return electronic_results





    def fetch_mechanical_properties(self, mp_id: str) -> dict:
        """
        Fetch mechanical properties for a given Materials Project ID.

        Parameters
        ----------
        mp_id : str
            Materials Project ID.

        Returns
        -------
        dict or None
            Mechanical properties including elastic tensors, bulk/shear modulus, Poisson ratio.

        Example
        -------
        >>> explorer.fetch_mechanical_properties("mp-123")
        {'mp-id': 'mp-123', 'elastic_tensor': [...], ...}
        """
        try:
            
            docs = self.mpr.materials.elasticity.search(material_ids=[mp_id])
            doc = docs[0].dict()
            
            mechanical_prop=doc
            mechanical_prop['mp-id']=mp_id
            

            self.results['mechanical_prop']=doc


        except Exception as e:
            self.results['mechanical_prop']=None
            return None
        return mechanical_prop
    
    

    def fetch_thermodynamic_properties(self, mp_id: str) -> dict:
        """
        Fetch thermodynamic properties for a given Materials Project ID.

        Parameters
        ----------
        mp_id : str
            Materials Project ID.

        Returns
        -------
        dict or None
            Thermodynamic properties including formation energy, heat capacities, etc.

        Example
        -------
        >>> explorer.fetch_thermodynamic_properties("mp-123")
        {'mp-id': 'mp-123', 'formation_energy_per_atom': -2.34, ...}
        """
        try:
            docs = self.mpr.materials.thermo.search(material_ids=[mp_id])
            
            doc = docs[0].dict()
            
            thermo_prop=doc
            thermo_prop['mp-id']=mp_id
            
            self.results['thermo_prop']=thermo_prop

           
        except Exception as e:
            self.results['thermo_prop']=None
            return None

        return thermo_prop


    def fetch_structure(self, mp_id: str) -> dict:
        """
        Fetch structural properties for a given Materials Project ID.

        Parameters
        ----------
        mp_id : str
            Materials Project ID.

        Returns
        -------
        list of dict or None
            Structure data in pymatgen format.

        Example
        -------
        >>> explorer.fetch_structure("mp-123")
        [{'structure_pymatgen': Structure(...)}]
        """
        try:
            struct = self.mpr.get_structure_by_material_id(mp_id)
            

            self.results["structure"]=[{'mp-id':mp_id,'structure_pymatgen': struct}]
            
        except Exception as e:
            self.results["structure"]=None
            return None
        
        
    
        return [{'structure_pymatgen': struct}]
    
    
    def fetch_all_data(self,mp_id):
        """
        Fetch all available data (electronic, mechanical, thermodynamic, structure) for a material.

        Parameters
        ----------
        mp_id : str
            Materials Project ID.

        Returns
        -------
        dict
            Dictionary with keys: 'electronic_prop', 'mechanical_prop', 'thermo_prop', 'structure', 'meta_data'.

        Example
        -------
        >>> all_data = explorer.fetch_all_data("mp-123")
        >>> all_data['electronic_prop']['band_gap']
        1.23
        """
        electronic=self.fetch_electronic_properties(mp_id)
        mechanical=self.fetch_mechanical_properties(mp_id)
        thermo=self.fetch_thermodynamic_properties(mp_id)
        structure=self.fetch_structure(mp_id)
        
        return self.results

    def close(self):
        """
        Close the underlying MPRester session.

        Example
        -------
        >>> explorer.close()
        """
        self.mpr.session.close()
        
        
        
        
        
      













