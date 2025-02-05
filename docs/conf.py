# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PCIpy'
copyright = '2025, Q. Baghi, J. Baker, E. Castelli'
author = 'Q. Baghi, J. Baker, E. Castelli'
release = '0.1'

# -- Path setup ---------------------------------------------------
# https://howto-sphinx.readthedocs.io/en/latest/api-modules.html
sys.path.insert(0, os.path.abspath('../pcipy'))


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# extensions = ['myst_parser',
              # "sphinx.ext.autosectionlabel",
              # "sphinx.ext.autodoc",
              # "sphinx.ext.viewcode",
              # "sphinx.ext.napoleon",
              # "sphinx.ext.mathjax",
# ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
