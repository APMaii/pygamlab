

from jarvis.db.figshare import data

class Jarvis_Explorer:
    """
    Python interface for exploring the JARVIS-DFT dataset.

    This class allows searching materials by chemical formula, 
    and retrieving electronic, mechanical, thermodynamic, 
    and structural properties for entries.

    Attributes
    ----------
    all_entries : list[dict]
        List of entries matching the last search.
    results : dict
        Dictionary storing electronic, mechanical, thermodynamic, 
        structural, and metadata for the last fetched material.

    Examples
    --------
    >>> explorer = Jarvis_Explorer()
    >>> entries = explorer.search_materials("MoS2", dataset="dft_3d", max_results=3)
    >>> entries
    [{'JID': 'JVASP-123', 'formula': 'MoS2', 'spg': 'P63/mmc', ...}, ...]
    >>> explorer.fetch_electronic_properties("JVASP-123")
    {'JID': 'JVASP-123', 'Formula': 'MoS2', 'optb88vdw_bandgap': 1.23, ...}
    >>> all_data = explorer.fetch_all_data("JVASP-123")
    >>> all_data['electronic_prop']['optb88vdw_bandgap']
    1.23
    """
    def __init__(self):
        """
        Initialize the JARVIS Explorer instance.

        Attributes initialized:
        - all_entries: Stores search results
        - results: Dictionary to store properties fetched for a material
        """
        
        self.all_entries = None
        self.results = {
            'electronic_prop': {},
            'mechanical_prop': {},
            'thermo_prop': {},
            'structure': {},
            'meta_data': {}
        }
        

    def search_materials(self, formula=None,dataset="dft_3d", max_results=5,description=True):
        """
        Search the JARVIS dataset for materials matching a chemical formula.

        Parameters
        ----------
        formula : str, optional
            Chemical formula to search for (e.g., "MoS2"). If None, returns all entries.
        dataset : str, default "dft_3d"
            Name of the dataset to query.
        max_results : int, default 5
            Maximum number of entries to return.
        description : bool, default True
            If True, returns detailed description dictionaries for each entry.
            If False, returns only JIDs.

        Returns
        -------
        list[dict] or list[str] or None
            If description=True, returns a list of dicts containing basic material info:
                - JID
                - formula
                - space group
                - bandgap
                - formation energy
                - density
                - dimensionality
            If description=False, returns a list of JIDs.
            Returns None if no matches are found.

        Example
        -------
        >>> explorer = Jarvis_Explorer()
        >>> results = explorer.search_materials("MoS2", dataset="dft_3d", max_results=3)
        >>> results[0]['formula']
        'MoS2'
        """
        print(f"Loading JARVIS dataset: {dataset} ...")
        self.dataset = data(dataset)   # loads local dataset
    
        self.max_results = max_results

        if formula:
            results = [m for m in self.dataset if m.get("formula") == formula]
        else:
            results = self.dataset

        if not results:
            print("No results found.")
            return None

        self.all_entries = results[:max_results]
        
        if description==True:
        
            descriptions = []
            for entry in self.all_entries:
                jid = entry.get("jid")
                formula = entry.get("formula")
                spg = entry.get("spg_symbol")
                bandgap = entry.get("optb88vdw_bandgap")
                ehull = entry.get("formation_energy_peratom")
                dim = entry.get("dimensionality", "3D")
                density = entry.get("density")
                
                descriptions.append({'JID':jid , 'formula':formula,
                                     'spg':spg ,'bandgap':bandgap ,
                                     'ehull':ehull , 'dim':dim,
                                     'density':density})
            return descriptions
        
        else:
            jid_list=[]
            for entry in self.all_entries:
                jid = entry.get("jid")
                jid_list.append(jid)
            return jid_list
                

    def fetch_electronic_properties(self,jid):
        """
        Retrieve electronic properties for a given JARVIS ID.

        Parameters
        ----------
        jid : str
            JARVIS ID of the material.

        Returns
        -------
        dict or None
            Electronic properties including bandgaps, effective masses, 
            dielectric constants, magnetization, and superconducting Tc.
            Returns None if jid not found or search not performed.

        Example
        -------
        >>> explorer.fetch_electronic_properties("JVASP-123")
        {'JID': 'JVASP-123', 'Formula': 'MoS2', 'optb88vdw_bandgap': 1.23, ...}
        """
        if self.all_entries is None:
            print("First run search_materials().")
            return None
        
        # Find the entry with the matching JID
        entry = next((e for e in self.all_entries if e.get("jid") == jid), None)
    
        if entry is None:
            print(f"JID {jid} not found in current search results.")
            return None
        

        results={"JID": entry.get("jid"),
            "Formula": entry.get("formula"),
            "optb88vdw_bandgap": entry.get("optb88vdw_bandgap"),
            "mbj_bandgap": entry.get("mbj_bandgap"),
            "hse_gap": entry.get("hse_gap"),
            "avg_elec_mass": entry.get("avg_elec_mass"),
            "avg_hole_mass": entry.get("avg_hole_mass"),
            "epsx": entry.get("epsx"),
            "epsy": entry.get("epsy"),
            "epsz": entry.get("epsz"),
            "mepsx": entry.get("mepsx"),
            "mepsy": entry.get("mepsy"),
            "mepsz": entry.get("mepsz"),
            "slme": entry.get("slme"),
            "efg": entry.get("efg"),
            "max_efg": entry.get("max_efg"),
            "magmom_oszicar": entry.get("magmom_oszicar"),
            "magmom_outcar": entry.get("magmom_outcar"),
            "spillage": entry.get("spillage"),
            "Tc_supercon": entry.get("Tc_supercon")}
            
        self.results['electronic_prop']=results
        return results

    def fetch_mechanical_properties(self,jid):
        """
        Retrieve mechanical/elastic properties for a given JARVIS ID.

        Parameters
        ----------
        jid : str
            JARVIS ID of the material.

        Returns
        -------
        dict or None
            Mechanical properties including elastic tensor, bulk modulus, 
            shear modulus, Poisson ratio, and piezoelectric constants.
            Returns None if jid not found or search not performed.

        Example
        -------
        >>> explorer.fetch_mechanical_properties("JVASP-123")
        {'JID': 'JVASP-123', 'Formula': 'MoS2', 'elastic_tensor': [...], ...}
        """
        if self.all_entries is None:
            print("First run search_materials().")
            return None
        
        # Find the entry with the matching JID
        entry = next((e for e in self.all_entries if e.get("jid") == jid), None)
    
        if entry is None:
            print(f"JID {jid} not found in current search results.")
            return None
        

        results={ "JID": entry.get("jid"),
            "Formula": entry.get("formula"),
            "elastic_tensor": entry.get("elastic_tensor"),
            "bulk_modulus_kv": entry.get("bulk_modulus_kv"),
            "shear_modulus_gv": entry.get("shear_modulus_gv"),
            "poisson": entry.get("poisson"),
            "dfpt_piezo_max_eij": entry.get("dfpt_piezo_max_eij"),
            "dfpt_piezo_max_dij": entry.get("dfpt_piezo_max_dij")}
        
        self.results['mechanical_prop']=results
        return results

    def fetch_thermodynamic_properties(self,jid):
        """
        Retrieve thermodynamic properties for a given JARVIS ID.

        Parameters
        ----------
        jid : str
            JARVIS ID of the material.

        Returns
        -------
        dict or None
            Thermodynamic properties including Seebeck coefficients, 
            power factors, conductivity, heat capacities, Debye temperature,
            formation energy, and exfoliation energy.
            Returns None if jid not found or search not performed.

        Example
        -------
        >>> explorer.fetch_thermodynamic_properties("JVASP-123")
        {'JID': 'JVASP-123', 'Formula': 'MoS2', 'n-Seebeck': 200, ...}
        """
        if self.all_entries is None:
            print("First run search_materials().")
            return None
        
        
        # Find the entry with the matching JID
        entry = next((e for e in self.all_entries if e.get("jid") == jid), None)
    
        if entry is None:
            print(f"JID {jid} not found in current search results.")
            return None
        

        results={"JID": entry.get("jid"),
            "Formula": entry.get("formula"),
            "n-Seebeck": entry.get("n-Seebeck"),
            "p-Seebeck": entry.get("p-Seebeck"),
            "n-powerfact": entry.get("n-powerfact"),
            "p-powerfact": entry.get("p-powerfact"),
            "ncond": entry.get("ncond"),
            "pcond": entry.get("pcond"),
            "nkappa": entry.get("nkappa"),
            "pkappa": entry.get("pkappa"),
            "modes": entry.get("modes"),
            "max_ir_mode": entry.get("max_ir_mode"),
            "min_ir_mode": entry.get("min_ir_mode"),
            "formation_energy_peratom": entry.get("formation_energy_peratom"),
            "exfoliation_energy": entry.get("exfoliation_energy")}
            
        self.results['thermo_prop']=results
        return results

    def fetch_structure(self,jid):
        """
        Retrieve structural properties for a given JARVIS ID.

        Parameters
        ----------
        jid : str
            JARVIS ID of the material.

        Returns
        -------
        dict or None
            Structural properties including space group, atoms, dimensionality,
            density, crystal type, and links to raw data files.
            Returns None if jid not found or search not performed.

        Example
        -------
        >>> explorer.fetch_structure("JVASP-123")
        {'JID': 'JVASP-123', 'Formula': 'MoS2', 'spg_number': 194, ...}
        """
        if self.all_entries is None:
            print("First run search_materials().")
            return None
        
        
        # Find the entry with the matching JID
        entry = next((e for e in self.all_entries if e.get("jid") == jid), None)
    
        if entry is None:
            print(f"JID {jid} not found in current search results.")
            return None
        

        results={"JID": entry.get("jid"),
            "Formula": entry.get("formula"),
            "spg_number": entry.get("spg_number"),
            "spg_symbol": entry.get("spg_symbol"),
            "atoms": entry.get("atoms"),
            "dimensionality": entry.get("dimensionality"),
            "nat": entry.get("nat"),
            "density": entry.get("density"),
            "spg": entry.get("spg"),
            "crys": entry.get("crys"),
            "typ": entry.get("typ"),
            "xml_data_link": entry.get("xml_data_link"),
            "raw_files": entry.get("raw_files"),
            "reference": entry.get("reference"),
            "search": entry.get("search")}
            
        self.results['structure']=results
        return results
    
    
    
    def fetch_all_data(self,jid):
        
        """
        Fetch all available data (electronic, mechanical, thermodynamic, structural, metadata)
        for a given JARVIS ID.

        Parameters
        ----------
        jid : str
            JARVIS ID of the material.

        Returns
        -------
        dict
            Dictionary with keys:
            - electronic_prop
            - mechanical_prop
            - thermo_prop
            - structure
            - meta_data

        Example
        -------
        >>> all_data = explorer.fetch_all_data("JVASP-123")
        >>> all_data['electronic_prop']['optb88vdw_bandgap']
        1.23
        """
        electronic=self.fetch_electronic_properties(jid)
        mechanical=self.fetch_mechanical_properties(jid)
        thermo=self.fetch_thermodynamic_properties(jid)
        structure=self.fetch_structure(jid)
        
        
        entry = next((e for e in self.all_entries if e.get("jid") == jid), None)


        results={"JID": entry.get("jid"),
                "Formula": entry.get("formula"),
                "func": entry.get("func"),
                "kpoint_length_unit": entry.get("kpoint_length_unit"),
                "maxdiff_mesh": entry.get("maxdiff_mesh"),
                "maxdiff_bz": entry.get("maxdiff_bz"),
                "encut": entry.get("encut"),
                "optb88vdw_total_energy": entry.get("optb88vdw_total_energy"),
                "icsd": entry.get("icsd")}
            
        self.results['meta_data']=results
            
        return self.results
    
    
    
    
    
    
    

