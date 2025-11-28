# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
# docs/source/conf.py
from __future__ import annotations

import os
import sys
from datetime import datetime

# Ensure package import works for autodoc
sys.path.insert(0, os.path.abspath("../.."))

project = "chemsmart"
author = "Zhang Lab, The Chinese University of Hong Kong"
copyright = f"{datetime.now():%Y}, {author}"
release = "0.1.11"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",  # Markdown support
    "sphinx.ext.autodoc",  # API from docstrings
    "sphinx.ext.napoleon",  # Google/NumPy docstrings
    "sphinx.ext.autosummary",  # summary tables + stub pages
    "sphinx.ext.intersphinx",  # cross-links to external docs
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",  # render type hints nicely
]

autosummary_generate = True

# MyST options (Markdown)
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3

# Cross-links to Python stdlib
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3",
        None,
    ),  # let Sphinx find objects.inv
}


# If heavy optional deps cause import errors on RTD, mock them here:
autodoc_mock_imports = [
    # "rdkit", "torch", "pyscf", ...
]

templates_path = ["_templates"]  # create or remove
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"  # matches your requirements.txt
html_static_path = ["_static"]  # create or remove
html_title = f"{project} {release}"
