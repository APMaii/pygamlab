
import math
import requests
from typing import List, Dict, Union
import warnings
import re
from collections import defaultdict

class COD_Explorer:
    """
    Python wrapper for the Crystallography Open Database (COD).

    Provides methods to search materials by chemical formula, 
    retrieve COD IDs, and fetch crystallographic information files (CIFs).

    Attributes
    ----------
    timeout : int
        Request timeout in seconds for HTTP queries.
    url : str
        Base URL of the COD website.
    api_url : str
        Endpoint URL for COD search queries.
    formula : str
        Last processed chemical formula in COD format.
    all_elements : list[str]
        List of all recognized chemical elements (1- and 2-letter symbols).

    Examples
    --------
    >>> explorer = COD_Explorer(timeout=30)
    >>> ids = explorer.search_materials("TiO2")
    >>> ids
    [900856, 900857, 901234]
    >>> structure = explorer.fetch_structure(ids[0])
    >>> print(structure['cif'][:200])  # Print first 200 chars of CIF
    """
    all_elements = [
        "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
        "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
        "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
        "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
        "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
    ]
    
    
    def __init__(self, timeout: int = 60):
        """
        Initialize a COD_Explorer instance.

        Parameters
        ----------
        timeout : int, optional
            Timeout in seconds for HTTP requests (default is 60).

        Example
        -------
        >>> explorer = COD_Explorer(timeout=30)
        """
        self.timeout = timeout
        self.url = "https://www.crystallography.net/"
        self.api_url = f"{self.url}/cod/result"
        self.formula=''
        
        
    def search_materials(self,formula):
        """
        Search COD for materials matching a chemical formula.

        Parameters
        ----------
        formula : str
            Chemical formula (e.g., 'TiO2', 'C6H6').

        Returns
        -------
        List[int]
            List of COD IDs corresponding to structures matching the formula.

        Example
        -------
        >>> explorer = COD_Explorer()
        >>> ids = explorer.search_materials("C6H6")
        >>> print(ids)
        [123456, 123457]
        """
        ids=self._get_cod_ids(formula)
        return ids
    

    def _get_cod_ids(self, formula: str) -> List[int]:
        """
        Internal method to query COD and return IDs for a given formula.

        Parameters
        ----------
        formula : str
            Chemical formula in standard notation.

        Returns
        -------
        List[int]
            List of COD IDs found.

        Raises
        ------
        requests.HTTPError
            If the COD server returns an error.

        Example
        -------
        >>> explorer = COD_Explorer()
        >>> explorer._get_cod_ids("H2O")
        [856789, 856790]
        """
        
        cod_formula=self._get_cod_format_formula(formula)
        #cod_formula = Composition(formula).hill_formula
        print('---------')
        print(cod_formula)
        params = {"formula": cod_formula, "format": "json"}
        response = requests.get(self.api_url, params=params, timeout=self.timeout)
        response.raise_for_status()
        return [int(entry["file"]) for entry in response.json()]
    
    
    def _corrected_formula(self, comp, ordered):
        """
        Convert element counts into a COD-compatible spaced formula.

        Parameters
        ----------
        comp : dict
            Mapping of element symbols to quantities.
        ordered : list[str]
            List of element symbols in the desired order (Hill system).

        Returns
        -------
        str
            Space-separated formula with integer quantities removed for 1s,
            e.g., 'C H6 O'.

        Example
        -------
        >>> comp = {'C': 6, 'H': 6, 'O': 1}
        >>> ordered = ['C', 'H', 'O']
        >>> explorer._corrected_formula(comp, ordered)
        'C H6 O'
        """
        parts = []
        for e in ordered:
            count = int(comp[e])
            if count == 0 or count == 1:
                parts.append(f"{e}")
            else:
                parts.append(f"{e}{count}")
        return " ".join(parts)

            
       
    def _get_cod_format_formula(self, formula: str) -> str:
        """
        Parse a chemical formula and return a COD-compatible format.

        The formula is parsed case-insensitively and supports 1- or 2-letter elements.
        The resulting formula is Hill-ordered and spaced, e.g., 'Ti O2', 'C6 H6 O'.

        Parameters
        ----------
        formula : str
            Input chemical formula (e.g., 'TiO2', 'C6H6').

        Returns
        -------
        str
            Formatted COD-compatible chemical formula.

        Raises
        ------
        ValueError
            If an unknown element is found in the formula.

        Example
        -------
        >>> explorer = COD_Explorer()
        >>> explorer._get_cod_format_formula("H2O")
        'H2 O'
        """
        formula = formula.strip().replace(" ", "").lower()
        comp = defaultdict(float)
        i = 0
        n = len(formula)

        while i < n:
            # Try two-letter element
            if i+1 < n and formula[i:i+2].capitalize() in self.all_elements:
                el = formula[i:i+2].capitalize()
                i += 2
                print(el)
            elif formula[i].capitalize() in self.all_elements:
                el = formula[i].capitalize()
                i += 1
                print(el)
            else:
                raise ValueError(f"Unknown element starting at '{formula[i:]}'")

            # Parse optional number
            num = ""
            while i < n and (formula[i].isdigit() or formula[i] == "."):
                num += formula[i]
                i += 1
            amt = float(num) if num else 1.0
            comp[el] += amt

        # Hill ordering
        els = list(comp.keys())
        if "C" in els:
            ordered = ["C"]
            if "H" in els:
                ordered.append("H")
            ordered += sorted([e for e in els if e not in {"C", "H"}])
        else:
            ordered = sorted(els)

        #return " ".join(f"{e}{self._formula_double_format(comp[e], ignore_ones=True)}" for e in ordered)
        
        final=self._corrected_formula(comp,ordered)
        self.formula=final
        return final
        
        #pre_final= " ".join(f"{e}{int(comp[e])}" for e in ordered)
        #final=pre_final.replace('1','').replace('0','')
        #return final

        

    def _get_cif_by_id(self, cod_id: int) -> str:
        """
        Retrieve the raw CIF text for a given COD ID.

        Parameters
        ----------
        cod_id : int
            COD ID of the desired crystal structure.

        Returns
        -------
        str
            Raw CIF text from the COD database.

        Raises
        ------
        requests.HTTPError
            If the COD server returns an error or the CIF is not found.

        Example
        -------
        >>> explorer = COD_Explorer()
        >>> cif_text = explorer._get_cif_by_id(856789)
        >>> print(cif_text[:200])
        """
        response = requests.get(f"{self.url}cod/{cod_id}.cif", timeout=self.timeout)
        response.raise_for_status()
        return response.text
    
    
    def fetch_structure(self,cod_id):
        """
        Fetch structure data for a given COD ID.

        Parameters
        ----------
        cod_id : int
            COD ID of the structure to fetch.

        Returns
        -------
        dict
            Dictionary containing:
            - 'formula': COD-formatted chemical formula
            - 'cod_id': COD ID of the structure
            - 'cif': Raw CIF content as a string

        Example
        -------
        >>> explorer = COD_Explorer()
        >>> ids = explorer.search_materials("TiO2")
        >>> structure = explorer.fetch_structure(ids[0])
        >>> print(structure['formula'])
        'Ti O2'
        """

        #structure_results=[]
        #ids=self._get_cod_ids(formula)
        
        #for specific_id in ids:
        structure_results={'formula' : self.formula,
                                  'cod_id' : cod_id,
                                  'cif': self._get_cif_by_id(cod_id)}
        
        return structure_results
    
    def fetch_all_data(self,cod_id):
        """
        Fetch all available data for a given COD ID.

        Parameters
        ----------
        cod_id : int
            COD ID of the structure to fetch.

        Returns
        -------
        dict
            Dictionary containing:
            - 'formula': COD-formatted chemical formula
            - 'cod_id': COD ID of the structure
            - 'cif': Raw CIF content as a string

        Notes
        -----
        Currently identical to `fetch_structure` but designed for future expansion
        to include more metadata fields from COD.

        Example
        -------
        >>> explorer = COD_Explorer()
        >>> data = explorer.fetch_all_data(856789)
        >>> print(data.keys())
        dict_keys(['formula', 'cod_id', 'cif'])
        """
        structure_results={'formula' : self.formula,
                                  'cod_id' : cod_id,
                                  'cif': self._get_cif_by_id(cod_id)}
        
        return structure_results
        
            
        
        
        

    
    
    









