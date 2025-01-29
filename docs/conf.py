# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../src'))


project = 'dem_lic'
copyright = '2025, Edmond SAINT DENIS'
author = 'Edmond SAINT DENIS'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "autodoc2",
    # "sphinx.ext.napoleon",
    # "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "myst_parser",
]

autodoc2_render_plugin = "myst"

autodoc2_packages = [
    "../src/dem_lic",  
]


# autodoc_default_options = {
#     "members": True,
#     "undoc-members": False,
#     "show-inheritance": True,
#     "inherited-members": True,
#     "special-members": "__init__",
#     "exclude-members": "__weakref__",
# }

# napoleon_use_param = True
# napoleon_google_docstring = True
# napoleon_numpy_docstring = True
# napoleon_include_init_with_doc = True
# napoleon_include_private_with_doc = False
# napoleon_include_special_with_doc = True




templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

# html_sidebars = {
#   "**": []
# }
html_theme_options = {
    # "show_nav_level": 0,  # Affiche les sous-niveaux dans la navigation
    # "navigation_depth": 1,  # Profondeur maximale des niveaux affichés
    # "show_toc_level": 0,  # Affiche les ancres internes des pages (titres)
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/ESaint-Denis/dem_lic",
            "icon": "fab fa-github",
        },
    ],
}