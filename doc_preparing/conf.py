# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html






import os
import sys
sys.path.insert(0, os.path.abspath('..'))



# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PyGamLab'
copyright = '2025, Ali Pilehvar Meibody'
author = 'Ali Pilehvar Meibody'
release = '2.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


# Extensions
extensions = [
    "nbsphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",        # Google/NumPy style docstrings
    "sphinx.ext.viewcode",        # Links to source code
    "sphinx_autodoc_typehints"    # Type hint support
]

templates_path = ['_templates']
exclude_patterns = []

html_show_sourcelink = False
# Control TOC behavior
toc_object_entries = False
#toc_object_entries_show_parents = False


autodoc_mock_imports = ['mp_api', 'aflow', 'jarvis', 'emmet', 'mp_api.client']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Theme
html_theme = "furo"
html_static_path = ['_static']
html_theme_options = {
    "sidebar_hide_name": False,
    "light_logo": "logo_dark.png",  # Optional logo for light mode
    "dark_logo": "logo_dark.png",    # Optional logo for dark mode
    "navigation_with_keys": True,
    #"show_toc_level": 1,  # This limits the sidebar TOC depth
    #"top_of_page_button":None,
}

# Autodoc default options
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": False,      # Set True if you want _private members
    "show-inheritance": True,
}



#-----
# Disable execution of notebooks
nbsphinx_execute = "never"

# Optional: allow existing notebook outputs
nbsphinx_allow_errors = True
