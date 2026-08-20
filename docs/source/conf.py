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

from datetime import datetime

project = "chemsmart"
author = "Zhang Lab, The Chinese University of Hong Kong"
copyright = f"{datetime.now():%Y}, {author}"
try:
    from importlib import metadata as _metadata

    release = _metadata.version("chemsmart")
except Exception:  # pragma: no cover - docs built outside an install
    release = "unknown"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",  # Markdown support
    "sphinx.ext.intersphinx",  # cross-links to external docs
    "sphinx.ext.mathjax",
]

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


exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"  # matches your requirements.txt
html_static_path = ["_static"]  # create or remove
html_title = f"{project} {release}"
