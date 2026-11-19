from aflow import search, K




class Aflow_Explorer:
    """
    A unified interface to search and fetch various material properties from the AFLOW database.

    This class allows searching by species or chemical formula, and fetching electronic, mechanical,
    thermodynamic, and structural properties of materials using AFLOW keys.

    Attributes
    ----------
    all_entries : list
        Stores the results of the latest search query.
    results : dict
        Stores the latest fetched data for a material.

    Examples
    --------
    >>> explorer = Aflow_Explorer()
    >>> results = explorer.search_materials(formula="MoS2")
    >>> explorer.fetch_electronic_properties("aflow-123")
    {'formula': 'MoS2', 'bandgap': 1.23, ...}
    >>> all_data = explorer.fetch_all_data("aflow-123")
    >>> all_data['electronic_prop']['bandgap']
    1.23
    """
    def __init__(self):
        """
        Initialize the AFLOW Explorer.
        """
        
        self.all_entries = None
        self.results = {
            'electronic_prop': {},
            'mechanical_prop': {},
            'thermo_prop': {},
            'structure': {},
            'meta_data': {}
        }

    def search_materials(self, species=None, formula=None, max_results=5, batch_size=10,description=True):
        """
        Search for materials in AFLOW by species or chemical formula.

        Parameters
        ----------
        species : str or list[str], optional
            Element(s) to search for.
        formula : str, optional
            Chemical formula to search for.
        max_results : int, default 5
            Maximum number of results to return.
        batch_size : int, default 10
            Number of entries retrieved per batch from AFLOW.
        description : bool, default True
            If True, returns detailed information for each entry.
            If False, returns only AFLOW unique IDs (AUIDs).

        Returns
        -------
        list of dict
            List of results containing AUID and optionally detailed description.

        Example
        -------
        >>> explorer = Aflow_Explorer()
        >>> explorer.search_materials(formula="MoS2", max_results=3)
        [{'auid': 'aflow-123', 'prototype': 'MoS2', 'spacegroup_relax': 187, ...}]
        """
        # -------------------- Clean .select() --------------------
        self.batch_size=batch_size
        
        query = (
            search(batch_size=self.batch_size)
            .select(
                K.auid, K.compound, K.species, K.composition,

                # Electronic
                K.Egap, K.Egap_fit, K.Egap_type,
                K.delta_electronic_energy_convergence, K.delta_electronic_energy_threshold,
                K.ldau_TLUJ, K.dft_type,
                K.bader_atomic_volumes, K.bader_net_charges,
                K.spinD, K.spinF, K.spin_atom, K.spin_cell,
                K.scintillation_attenuation_length,

                # Mechanical
                K.ael_bulk_modulus_reuss, K.ael_bulk_modulus_voigt, K.ael_bulk_modulus_vrh,
                K.ael_shear_modulus_reuss, K.ael_shear_modulus_voigt, K.ael_shear_modulus_vrh,
                K.ael_poisson_ratio, K.ael_elastic_anisotropy,
                K.stress_tensor, K.forces,
                K.Pulay_stress, K.pressure, K.pressure_residual,

                # Thermodynamic
                K.agl_acoustic_debye, K.agl_debye,
                K.agl_gruneisen,
                K.agl_heat_capacity_Cp_300K, K.agl_heat_capacity_Cv_300K,
                K.agl_thermal_conductivity_300K, K.agl_thermal_expansion_300K,
                K.agl_bulk_modulus_isothermal_300K, K.agl_bulk_modulus_static_300K,
                K.eentropy_atom, K.eentropy_cell,
                K.enthalpy_atom, K.enthalpy_cell,
                K.enthalpy_formation_atom, K.enthalpy_formation_cell,
                K.energy_atom, K.energy_cell, K.energy_cutoff,
                K.entropic_temperature,

                # Structural
                K.Bravais_lattice_orig, K.Bravais_lattice_relax,
                K.Pearson_symbol_orig, K.Pearson_symbol_relax,
                K.lattice_system_orig, K.lattice_system_relax,
                K.lattice_variation_orig, K.lattice_variation_relax,
                K.spacegroup_orig, K.spacegroup_relax,
                K.sg, K.sg2,
                K.prototype, K.stoich, K.stoichiometry, K.geometry,
                K.natoms, K.nspecies, K.nbondxx, K.species_pp, K.species_pp_ZVAL, K.species_pp_version,
                K.positions_cartesian, K.positions_fractional,K.valence_cell_iupac, K.valence_cell_std,
                K.volume_atom, K.volume_cell,K.density
            )
        )
        
        self.batch_size = batch_size

        # -------------------- Filters --------------------
        if species:
            if isinstance(species, list):
                for el in species:
                    query = query.filter(K.species == el)
            else:
                query = query.filter(K.species == species)

        if formula:
            query = query.filter(K.compound == formula)

        all_entries = []
        try:
            for i, entry in enumerate(query):
                if i < max_results:
                    all_entries.append(entry)
                else:
                    break
            self.all_entries = all_entries
            
            #return all_entries
        
            if description==True:
                results=[]
                for entry in self.all_entries:
                    results.append({'auid':getattr(entry, 'auid', None),
                                    'prototype':getattr(entry, 'prototype', None),
                                    'spacegroup_relax':getattr(entry, 'spacegroup_relax', None),
                                    'dft_type':getattr(entry, 'dft_type', None),
                                    'spinD':getattr(entry, 'spinD', None),
                                    'spinF':getattr(entry, 'spinF', None),
                                    'enthalpy_formation_cell':getattr(entry, 'enthalpy_formation_cell', None),
                                    'natoms':getattr(entry, 'natoms', None)})
                return results

            else:
                results=[]
                for entry in self.all_entries:
                    results.append({'auid':getattr(entry, 'auid', None)})
                return results
                
            #    return
        except AssertionError:
            print("No results found for this query.")
            return

    # -------------------- Electronic properties --------------------
    def fetch_electronic_properties(self,auid):
        """
        Fetch electronic properties for a given AFLOW AUID.

        Parameters
        ----------
        auid : str
            AFLOW unique ID of the material.

        Returns
        -------
        dict
            Dictionary containing electronic properties such as bandgap, spin polarization, etc.

        Example
        -------
        >>> explorer.fetch_electronic_properties("aflow-123")
        {'formula': 'MoS2', 'bandgap': 1.23, ...}
        """
        if self.all_entries is None:
            print("First use search_materials function")
            return None
        
        #entry = next((e for e in self.all_entries if e.get("auid") == auid), None)
        
        for r in self.all_entries:
            if r.auid==auid:
                entry=r
                
                
                
        
        if entry is None:
            print(f"auid {auid} not found in current search results.")
            return None
        


        electronic_results={'formula': getattr(entry, 'compound', None),
            'bandgap': getattr(entry, 'Egap', None),
            'bandgap_fit': getattr(entry, 'Egap_fit', None),
            'bandgap_type': getattr(entry, 'Egap_type', None),
            'delta_elec_energy_convergence': getattr(entry, 'delta_electronic_energy_convergence', None),
            'delta_elec_energy_threshold': getattr(entry, 'delta_electronic_energy_threshold', None),
            'ldau_TLUJ': getattr(entry, 'ldau_TLUJ', None),
            'dft_type': getattr(entry, 'dft_type', None),
            'bader_atomic_volumes': getattr(entry, 'bader_atomic_volumes', None),
            'bader_net_charges': getattr(entry, 'bader_net_charges', None),
            'spinD': getattr(entry, 'spinD', None),
            'spinF': getattr(entry, 'spinF', None),
            'spin_atom': getattr(entry, 'spin_atom', None),
            'spin_cell': getattr(entry, 'spin_cell', None),
            'scintillation_attenuation_length': getattr(entry, 'scintillation_attenuation_length', None)}

        self.results['electronic_prop'] = electronic_results
        return electronic_results

    # -------------------- Mechanical properties --------------------
    def fetch_mechanical_properties(self,auid):
        """
        Fetch mechanical properties for a given AFLOW AUID.

        Parameters
        ----------
        auid : str
            AFLOW unique ID of the material.

        Returns
        -------
        dict
            Dictionary containing bulk modulus, shear modulus, elastic tensors, etc.
        """
        if self.all_entries is None:
            print("First use search_materials function")
            return None
        
        
        
        #entry = next((e for e in self.all_entries if e.get("auid") == auid), None)
        
        for r in self.all_entries:
            if r.auid==auid:
                entry=r
        
        if entry is None:
            print(f"auid {auid} not found in current search results.")
            return None
        

        mechanical_results={ 'formula': getattr(entry, 'compound', None),
            'bulk_modulus_reuss': getattr(entry, 'ael_bulk_modulus_reuss', None),
            'bulk_modulus_voigt': getattr(entry, 'ael_bulk_modulus_voigt', None),
            'bulk_modulus_vrh': getattr(entry, 'ael_bulk_modulus_vrh', None),
            'shear_modulus_reuss': getattr(entry, 'ael_shear_modulus_reuss', None),
            'shear_modulus_voigt': getattr(entry, 'ael_shear_modulus_voigt', None),
            'shear_modulus_vrh': getattr(entry, 'ael_shear_modulus_vrh', None),
            'poisson_ratio': getattr(entry, 'ael_poisson_ratio', None),
            'elastic_anisotropy': getattr(entry, 'ael_elastic_anisotropy', None),
            'stress_tensor': getattr(entry, 'stress_tensor', None),
            'forces': getattr(entry, 'forces', None),
            'Pulay_stress': getattr(entry, 'Pulay_stress', None),
            'pressure': getattr(entry, 'pressure', None),
            'pressure_residual': getattr(entry, 'pressure_residual', None)}

        self.results['mechanical_prop'] = mechanical_results
        return mechanical_results

    # -------------------- Thermodynamic properties --------------------
    def fetch_thermodynamic_properties(self,auid):
        """
        Fetch thermodynamic properties for a given AFLOW AUID.

        Parameters
        ----------
        auid : str
            AFLOW unique ID of the material.

        Returns
        -------
        dict
            Dictionary containing Debye temperature, heat capacity, enthalpy, and related properties.
        """
        if self.all_entries is None:
            print("First use search_materials function")
            return None
        
        
        #entry = next((e for e in self.all_entries if e.get("auid") == auid), None)
        
        for r in self.all_entries:
            if r.auid==auid:
                entry=r
        
        if entry is None:
            print(f"auid {auid} not found in current search results.")
            return None
        

        thermo_results={'formula': getattr(entry, 'compound', None),
            'acoustic_debye': getattr(entry, 'agl_acoustic_debye', None),
            'debye': getattr(entry, 'agl_debye', None),
            'gruneisen': getattr(entry, 'agl_gruneisen', None),
            'heat_capacity_Cp_300K': getattr(entry, 'agl_heat_capacity_Cp_300K', None),
            'heat_capacity_Cv_300K': getattr(entry, 'agl_heat_capacity_Cv_300K', None),
            'thermal_conductivity_300K': getattr(entry, 'agl_thermal_conductivity_300K', None),
            'thermal_expansion_300K': getattr(entry, 'agl_thermal_expansion_300K', None),
            'bulk_modulus_isothermal_300K': getattr(entry, 'agl_bulk_modulus_isothermal_300K', None),
            'bulk_modulus_static_300K': getattr(entry, 'agl_bulk_modulus_static_300K', None),
            'entropy_atom': getattr(entry, 'eentropy_atom', None),
            'entropy_cell': getattr(entry, 'eentropy_cell', None),
            'enthalpy_atom': getattr(entry, 'enthalpy_atom', None),
            'enthalpy_cell': getattr(entry, 'enthalpy_cell', None),
            'enthalpy_formation_atom': getattr(entry, 'enthalpy_formation_atom', None),
            'enthalpy_formation_cell': getattr(entry, 'enthalpy_formation_cell', None),
            'energy_atom': getattr(entry, 'energy_atom', None),
            'energy_cell': getattr(entry, 'energy_cell', None),
            'energy_cutoff': getattr(entry, 'energy_cutoff', None),
            'entropic_temperature': getattr(entry, 'entropic_temperature', None)}

        self.results['thermo_prop'] = thermo_results
        return thermo_results

    # -------------------- Structural properties --------------------
    def fetch_structure(self,auid):
        """
        Fetch structural properties for a given AFLOW AUID.

        Parameters
        ----------
        auid : str
            AFLOW unique ID of the material.

        Returns
        -------
        dict
            Dictionary containing lattice, space group, atomic positions, density, and valence info.
        """
        if self.all_entries is None:
            print("First use search_materials function")
            return None
        
       # entry = next((e for e in self.all_entries if e.get("auid") == auid), None)
        for r in self.all_entries:
            if r.auid==auid:
                entry=r
        
        if entry is None:
            print(f"auid {auid} not found in current search results.")
            return None
        

        structure_results={'formula': getattr(entry, 'compound', None),
            'Bravais_lattice_orig': getattr(entry, 'Bravais_lattice_orig', None),
            'Bravais_lattice_relax': getattr(entry, 'Bravais_lattice_relax', None),
            'Pearson_symbol_orig': getattr(entry, 'Pearson_symbol_orig', None),
            'Pearson_symbol_relax': getattr(entry, 'Pearson_symbol_relax', None),
            'lattice_system_orig': getattr(entry, 'lattice_system_orig', None),
            'lattice_system_relax': getattr(entry, 'lattice_system_relax', None),
            'lattice_variation_orig': getattr(entry, 'lattice_variation_orig', None),
            'lattice_variation_relax': getattr(entry, 'lattice_variation_relax', None),
            'spacegroup_orig': getattr(entry, 'spacegroup_orig', None),
            'spacegroup_relax': getattr(entry, 'spacegroup_relax', None),
            'sg': getattr(entry, 'sg', None),
            'sg2': getattr(entry, 'sg2', None),
            'prototype': getattr(entry, 'prototype', None),
            'stoich': getattr(entry, 'stoich', None),
            'stoichiometry': getattr(entry, 'stoichiometry', None),
            'geometry': getattr(entry, 'geometry', None),
            'natoms': getattr(entry, 'natoms', None),
            'nspecies': getattr(entry, 'nspecies', None),
            'nbondxx': getattr(entry, 'nbondxx', None),
            'composition': getattr(entry, 'composition', None),
            'compound': getattr(entry, 'compound', None),
            'species': getattr(entry, 'species', None),
            "species_pp": getattr(entry, "species_pp", None),
            "species_pp_ZVAL": getattr(entry, "species_pp_ZVAL", None),
            "species_pp_version": getattr(entry, "species_pp_version", None),
            "positions_cartesian": getattr(entry, "positions_cartesian", None),
            "positions_fractional": getattr(entry, "positions_fractional", None),

            # Valence info
            "valence_cell_iupac": getattr(entry, "valence_cell_iupac", None),
            "valence_cell_std": getattr(entry, "valence_cell_std", None),
            "volume_atom": getattr(entry, "volume_atom", None),
            "volume_cell": getattr(entry, "volume_cell", None),
            "density": getattr(entry, "density", None)}

        self.results['structure'] = structure_results
        return structure_results
    # -------------------- Aggregate all data --------------------
    def fetch_all_data(self,auid):
        """
        Fetch all available data (electronic, mechanical, thermodynamic, structure) for a material.

        Parameters
        ----------
        auid : str
            AFLOW unique ID of the material.

        Returns
        -------
        dict
            Dictionary with keys: 'electronic_prop', 'mechanical_prop', 'thermo_prop', 'structure'.

        Example
        -------
        >>> all_data = explorer.fetch_all_data("aflow-123")
        >>> all_data['electronic_prop']['bandgap']
        1.23
        """
        if self.all_entries is None:
            print('first use search_materials function')
            return None

        electro=self.fetch_electronic_properties(auid)
        mechanical=self.fetch_mechanical_properties(auid)
        thermo=self.fetch_thermodynamic_properties(auid)
        structure=self.fetch_structure(auid)
        
        return self.results
    
    

    







    
    
    
    
    
    