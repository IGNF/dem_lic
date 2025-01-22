# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../src'))


project = 'demgen_lic'
copyright = '2025, Edmond SAINT DENIS'
author = 'Edmond SAINT DENIS'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # Pour supporter les docstrings au format NumPy/Google
    "sphinx.ext.githubpages"
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

html_theme_options = {
    "navbar_end": ["search-field", "navbar-icon-links"],
    "primary_sidebar_end": ["sidebar-ethical-ads"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/ESaint-Denis/dem_lic",
            "icon": "fab fa-github",
        },
    ],
}