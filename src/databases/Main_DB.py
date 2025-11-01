from .COD import COD_Explorer
from .Material_projects import MaterialsProject_Explorer
from .Jarvis import Jarvis_Explorer
from .Aflow import Aflow_Explorer


class GAM_Explorer:
    """
    Unified interface to explore material properties from multiple databases:
    Materials Project, COD (Crystallography Open Database), Jarvis (NIST), and AFLOW.

    This class allows searching and fetching electronic, mechanical, thermodynamic,
    and structural properties from a selected backend or from all backends simultaneously.

    Parameters
    ----------
    backend : str, default='jarvis'
        Backend database to use. Options are:
        'material_project', 'cod', 'jarvis', 'aflow', 'all'.
    timeout : int, default=60
        Timeout in seconds for database queries (applies to COD).
    api_key : str, optional
        API key required for Materials Project access.
    dataset : str, default='dft_3d'
        Dataset to use when querying Jarvis.
    max_results : int, default=5 
        Maximum number of materials to return for search queries.
    batch_size : int, default=10
        Batch size for queries (applies to AFLOW).

    Attributes
    ----------
    backend : str
        Normalized backend name.
    _explorer : object
        Instance of the selected backend explorer.
    search_results : dict or list
        Stores results from the last search query.

    Examples
    --------
    # Example 1: Using Jarvis backend to search and fetch properties
    >>> gam = GAM_Explorer(backend='jarvis')
    >>> results = gam.search_materials(formula="MoS2")
    >>> electronic_props = gam.fetch_electronic_properties(gam_id="JARVIS-123")
    >>> mechanical_props = gam.fetch_mechanical_properties(gam_id="JARVIS-123")
    >>> thermo_props = gam.fetch_thermodynamic_properties(gam_id="JARVIS-123")
    >>> structure = gam.fetch_structure(gam_id="JARVIS-123")
    >>> all_data = gam.fetch_all_data(gam_id="JARVIS-123")

    # Example 2: Using Materials Project backend
    >>> gam_mp = GAM_Explorer(backend='material_project', api_key="MP_API_KEY")
    >>> results_mp = gam_mp.search_materials(elements=["Mo","S"])
    >>> electronic_mp = gam_mp.fetch_electronic_properties(gam_id="mp-123")

    # Example 3: Using AFLOW backend
    >>> gam_aflow = GAM_Explorer(backend='aflow')
    >>> results_aflow = gam_aflow.search_materials(species="MoS2")
    >>> mechanical_aflow = gam_aflow.fetch_mechanical_properties(gam_id="aflow-456")

    # Example 4: Using COD backend for structure data
    >>> gam_cod = GAM_Explorer(backend='cod')
    >>> results_cod = gam_cod.search_materials(formula="MoS2")
    >>> structure_cod = gam_cod.fetch_structure(gam_id="COD-789")

    # Example 5: Using all backends simultaneously
    >>> gam_all = GAM_Explorer(backend='all', api_key="MP_API_KEY")
    >>> results_all = gam_all.search_materials(formula="MoS2", elements=["Mo","S"], species="MoS2", mp_ids=["mp-123"])
    >>> all_props = gam_all.fetch_all_data(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456", code_id="COD-789")
    >>> all_electronic = gam_all.fetch_electronic_properties(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456")
    >>> all_mechanical = gam_all.fetch_mechanical_properties(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456")
    >>> all_thermo = gam_all.fetch_thermodynamic_properties(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456")
    >>> all_structures = gam_all.fetch_structure(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456", code_id="COD-789")
    """

    def __init__(self,backend='jarvis',timeout=60,api_key=None,dataset="dft_3d",
                 max_results=5,batch_size=10):
        
        self.backend=backend.lower().strip()
        self.timeout=timeout
        self.api_key=api_key
        self.dataset=dataset
        self.max_results=max_results
        self.batch_size=batch_size
        self.search_results=None
        


        if self.backend in ['material_project','material project' ,'materials_project','materials project' ,'mp','mpapi','mp-api','mp_api']:
            if api_key==None:
                raise ValueError("You must provide an API key for Materials Project.")
                
            self._backend='mp'
            self._explorer=MaterialsProject_Explorer(self.api_key)
                
                
        elif self.backend in ['cod','crystallography open database','crystallography','crystallography_open_database']:
            
            self._backend='cod'
            self._explorer=COD_Explorer(self.timeout)
            
            
        elif self.backend in ['nist','jarvis']:
            
            
            self._backend='jarvis'
            self._explorer=Jarvis_Explorer()
            

        elif self.backend in ['aflow','aflowlib','aflow ++']:
            
            self._backend='aflow'
            self._explorer=Aflow_Explorer()
            


        elif self.backend in ['all']:
            if api_key==None:
                raise ValueError("You must provide an API key for Materials Project.")
                         
            self._backend='all'
            self._explorer1=MaterialsProject_Explorer(self.api_key)
            self._explorer2=COD_Explorer(self.timeout)
            self._explorer3=Jarvis_Explorer()
            self._explorer4=Aflow_Explorer()
            
            
        else:
            raise ValueError("Options for Backend is [material_project,'crystallography open database','jarvis','aflow','all'")
            
            
        
        
    def search_materials(self,formula=None,elements=None,mp_ids=None,species=None,description=True):
        """
        Search for materials using the selected backend.

        Parameters
        ----------
        formula : str, optional
            Chemical formula to search for.
        elements : list, optional
            Elements to search for (only for Materials Project).
        mp_ids : list, optional
            Materials Project IDs (only for Materials Project backend).
        species : str or list, optional
            Species for AFLOW backend.
        description : bool, default=True
            If True, returns detailed descriptions. Otherwise, only IDs.

        Returns
        -------
        list or dict
            Search results. Dict is returned if backend='all'.

        Raises
        ------
        ValueError
            If backend-specific arguments are used incorrectly.
        """
        self.description=description
        if mp_ids!=None and self._backend!='mp':
            raise ValueError('only Material Project Backend has mp_ids argument')
            
        if species!=None and self._backend!='aflow':
            raise ValueError('only Aflow Backend has Species argument')
        
        if elements!=None and self._backend!='mp':
            raise ValueError('only Material Project Backend has Elements argument')
            
            
        
        if self._backend=='mp':
            self.search_results=self._explorer.search_materials(elements,formula,mp_ids,self.max_results,self.description)
            
        
            
        elif self._backend=='cod':
            self.search_results=self._explorer.search_materials(formula)

            
        elif self._backend=='jarvis':
            self.search_results=self._explorer.search_materials(formula,self.dataset,self.max_results,self.description)

        elif self._backend=='aflow':
            self.search_results=self._explorer.search_materials(species,formula,self.max_results,self.batch_size,self.description)

            
        elif self._backend=='all':
            self.search_results={}
            
            self.search_results['Material Project']=self._explorer1.search_materials(elements,formula,mp_ids,self.max_results,self.description)

            self.search_results['COD']=self._explorer2.search_materials(formula)
            
            self.search_results['Jarvis']=self._explorer3.search_materials(formula,self.dataset,self.max_results,self.description)
            self.search_results['Aflow']=self._explorer4.search_materials(species,formula,self.max_results,self.batch_size,self.description)
            

        return self.search_results

            
            
            

    def fetch_electronic_properties(self,gam_id=None,mp_id= None,juid=None , auid= None):
        """
        Fetch electronic properties for a material.

        Parameters
        ----------
        gam_id : str, optional
            Material ID for single backend (Jarvis, AFLOW, or Materials Project).
        mp_id : str, optional
            Materials Project ID (required when backend='all').
        juid : str, optional
            Jarvis ID (required when backend='all').
        auid : str, optional
            AFLOW ID (required when backend='all').

        Returns
        -------
        dict or None
            Electronic properties of the material. Returns a dictionary containing
            results from multiple backends if backend='all'. Returns None for COD backend.

        Raises
        ------
        ValueError
            If required IDs are not provided.
        NotImplementedError
            If COD backend is selected.

        Examples
        --------
        >>> gam = GAM_Explorer(backend='jarvis')
        >>> gam.fetch_electronic_properties(gam_id="JARVIS-123")
        {'band_gap': 1.5, 'eigenvalues': [...], ...}

        >>> gam = GAM_Explorer(backend='all', api_key="MP_API_KEY")
        >>> gam.fetch_electronic_properties(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456")
        {'Material Project': {...}, 'Jarvis': {...}, 'Aflow': {...}, 'COD': None}
        """
  
        if self._backend=='all':
            if mp_id==None or juid==None or auid==None:
                raise ValueError('You must insert your mp_id, juid and auid for fetching data')
                
            else:
                electronic_result={}
                electronic_result['Material Project']=self._explorer1.fetch_electronic_properties(mp_id)
                electronic_result['Jarvis']=self._explorer3.fetch_electronic_properties(juid)
                electronic_result['Aflow']=self._explorer4.fetch_electronic_properties(auid)
                electronic_result['COD']=None
            
        elif self._backend=='cod':
            electronic_result=None
            
            raise NotImplementedError(
            "COD backend does not support fetching electronic properties. "
            "Only structure data is available."
        )
            
        else:
            if gam_id==None:
                raise ValueError('You must insert your gam_id for fetching data')
                
            electronic_result=self._explorer.fetch_electronic_properties(gam_id)

        return electronic_result
        
        
        
    def fetch_mechanical_properties(self,gam_id=None,mp_id= None,juid=None , auid= None):
        """
        Fetch mechanical properties for a material.

        Parameters
        ----------
        gam_id : str, optional
            Material ID for single backend (Jarvis, AFLOW, or Materials Project).
        mp_id : str, optional
            Materials Project ID (required when backend='all').
        juid : str, optional
            Jarvis ID (required when backend='all').
        auid : str, optional
            AFLOW ID (required when backend='all').

        Returns
        -------
        dict or None
            Mechanical properties of the material. Returns a dictionary containing
            results from multiple backends if backend='all'. Returns None for COD backend.

        Raises
        ------
        ValueError
            If required IDs are not provided.
        NotImplementedError
            If COD backend is selected.

        Examples
        --------
        >>> gam = GAM_Explorer(backend='jarvis')
        >>> gam.fetch_mechanical_properties(gam_id="JARVIS-123")
        {'bulk_modulus_reuss': 120.5, 'bulk_modulus_voigt': 123.0, ...}

        >>> gam = GAM_Explorer(backend='all', api_key="MP_API_KEY")
        >>> gam.fetch_mechanical_properties(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456")
        {'Material Project': {...}, 'Jarvis': {...}, 'Aflow': {...}, 'COD': None}
        """
        
        if self._backend=='all':
            if mp_id==None or juid==None or auid==None:
                raise ValueError('You must insert your mp_id, juid and auid for fetching data')
                
            else:
                mechanical_result={}
                mechanical_result['Material Project']=self._explorer1.fetch_mechanical_properties(mp_id)
                mechanical_result['Jarvis']=self._explorer3.fetch_mechanical_properties(juid)
                mechanical_result['Aflow']=self._explorer4.fetch_mechanical_properties(auid)
                mechanical_result['COD']=None
            
        elif self._backend=='cod':
            mechanical_result=None
            raise NotImplementedError(
            "COD backend does not support fetching electronic properties. "
            "Only structure data is available."
        )
            
        else:
            if gam_id==None:
                raise ValueError('You must insert your gam_id for fetching data')
                
            mechanical_result=self._explorer.fetch_mechanical_properties(gam_id)

        return mechanical_result
        
        
    def fetch_thermodynamic_properties(self,gam_id=None,mp_id= None,juid=None , auid= None):
        """
        Fetch thermodynamic properties for a material.

        Parameters
        ----------
        gam_id : str, optional
            Material ID for single backend (Jarvis, AFLOW, or Materials Project).
        mp_id : str, optional
            Materials Project ID (required when backend='all').
        juid : str, optional
            Jarvis ID (required when backend='all').
        auid : str, optional
            AFLOW ID (required when backend='all').

        Returns
        -------
        dict or None
            Thermodynamic properties of the material. Returns a dictionary containing
            results from multiple backends if backend='all'. Returns None for COD backend.

        Raises
        ------
        ValueError
            If required IDs are not provided.
        NotImplementedError
            If COD backend is selected.

        Examples
        --------
        >>> gam = GAM_Explorer(backend='jarvis')
        >>> gam.fetch_thermodynamic_properties(gam_id="JARVIS-123")
        {'enthalpy_formation_cell': -5.4, 'entropy_cell': 12.3, ...}

        >>> gam = GAM_Explorer(backend='all', api_key="MP_API_KEY")
        >>> gam.fetch_thermodynamic_properties(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456")
        {'Material Project': {...}, 'Jarvis': {...}, 'Aflow': {...}, 'COD': None}
        """
        if self._backend=='all':
            if mp_id==None or juid==None or auid==None:
                raise ValueError('You must insert your mp_id, juid and auid for fetching data')
                
            else:
                thermo_result={}
                thermo_result['Material Project']=self._explorer1.fetch_thermodynamic_properties(mp_id)
                thermo_result['Jarvis']=self._explorer3.fetch_thermodynamic_properties(juid)
                thermo_result['Aflow']=self._explorer4.fetch_thermodynamic_properties(auid)
                thermo_result['COD']=None
            
        elif self._backend=='cod':
            thermo_result=None
            raise NotImplementedError(
            "COD backend does not support fetching electronic properties. "
            "Only structure data is available."
        )
            
        else:
            if gam_id==None:
                raise ValueError('You must insert your gam_id for fetching data')
                
            thermo_result=self._explorer.fetch_thermodynamic_properties(gam_id)

        return thermo_result
        
        
    def fetch_structure(self,gam_id=None,mp_id= None,juid=None , auid= None,code_id=None):
        """
        Fetch structural properties for a material.

        Parameters
        ----------
        gam_id : str, optional
            Material ID for single backend (Jarvis, AFLOW, or Materials Project).
        mp_id : str, optional
            Materials Project ID (required when backend='all').
        juid : str, optional
            Jarvis ID (required when backend='all').
        auid : str, optional
            AFLOW ID (required when backend='all').
        code_id : str, optional
            COD ID (required when backend='all').

        Returns
        -------
        dict or None
            Structural properties of the material. Returns a dictionary containing
            results from multiple backends if backend='all'.

        Raises
        ------
        ValueError
            If required IDs are not provided.

        Examples
        --------
        >>> gam = GAM_Explorer(backend='jarvis')
        >>> gam.fetch_structure(gam_id="JARVIS-123")
        {'spacegroup_relax': 194, 'Bravais_lattice_relax': 'hexagonal', ...}

        >>> gam = GAM_Explorer(backend='all', api_key="MP_API_KEY")
        >>> gam.fetch_structure(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456", code_id="COD-789")
        {'Material Project': {...}, 'COD': {...}, 'Jarvis': {...}, 'Aflow': {...}}
        """
        if self._backend=='all':
            if mp_id==None or juid==None or auid==None or code_id==None:
                raise ValueError('You must insert your mp_id, juid , auid and code_id for fetching data')
                
            else:
                structure_result={}
                structure_result['Material Project']=self._explorer1.fetch_structure(mp_id)
                structure_result['COD']=self._explorer2.fetch_structure(code_id)
                structure_result['Jarvis']=self._explorer3.fetch_structure(juid)
                structure_result['Aflow']=self._explorer4.fetch_structure(auid)

        else:
            if gam_id==None:
                raise ValueError('You must insert your gam_id for fetching data')
                
            structure_result=self._explorer.fetch_structure(gam_id)

        return structure_result
    
       
        
        
    def fetch_all_data(self,gam_id=None,mp_id= None,juid=None , auid= None,code_id=None):
        """
        Fetch all available properties (electronic, mechanical, thermodynamic, structure)
        for a material.

        Parameters
        ----------
        gam_id : str, optional
            Material ID for single backend (Jarvis, AFLOW, or Materials Project).
        mp_id : str, optional
            Materials Project ID (required when backend='all').
        juid : str, optional
            Jarvis ID (required when backend='all').
        auid : str, optional
            AFLOW ID (required when backend='all').
        code_id : str, optional
            COD ID (required when backend='all').

        Returns
        -------
        dict
            Dictionary of all properties. If backend='all', returns results
            from all backends; otherwise, returns results from the selected backend.

        Raises
        ------
        ValueError
            If required IDs are not provided.

        Examples
        --------
        >>> gam = GAM_Explorer(backend='jarvis')
        >>> gam.fetch_all_data(gam_id="JARVIS-123")
        {'electronic_prop': {...}, 'mechanical_prop': {...}, 'thermo_prop': {...}, 'structure': {...}}

        >>> gam = GAM_Explorer(backend='all', api_key="MP_API_KEY")
        >>> gam.fetch_all_data(mp_id="mp-123", juid="JARVIS-123", auid="aflow-456", code_id="COD-789")
        {'Material Project': {...}, 'COD': {...}, 'Jarvis': {...}, 'Aflow': {...}}
        """
        if self._backend=='all':
            if mp_id==None or juid==None or auid==None:
                raise ValueError('You must insert your mp_id, juid and auid for fetching data')
                
            else:
                all_result={}
                all_result['Material Project']=self._explorer1.fetch_all_data(mp_id)
                all_result['COD']=self._explorer2.fetch_all_data(code_id)
                all_result['Jarvis']=self._explorer3.fetch_all_data(juid)
                all_result['Aflow']=self._explorer4.fetch_all_data(auid)

        else:
            if gam_id==None:
                raise ValueError('You must insert your gam_id for fetching data')
                
            all_result=self._explorer.fetch_all_data(gam_id)

        return all_result
    
    
            
            
            
        
        
        
        
        
        
        
