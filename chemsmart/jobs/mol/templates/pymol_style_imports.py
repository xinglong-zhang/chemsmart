"""
Shared import helper for PyMOL style template scripts.
"""

import importlib.util
import os

SHARED_MODULE_NAME = "metallic_element_colors"


def load_metallic_element_colors():
    """Load metallic_element_colors.py from the template or job directory."""
    search_paths = []
    try:
        search_paths.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        pass
    search_paths.append(os.getcwd())

    for directory in search_paths:
        module_path = os.path.join(directory, f"{SHARED_MODULE_NAME}.py")
        if os.path.exists(module_path):
            spec = importlib.util.spec_from_file_location(
                SHARED_MODULE_NAME, module_path
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module

    raise ImportError(
        f"{SHARED_MODULE_NAME}.py not found in: {', '.join(search_paths)}"
    )
