import importlib_resources
import os

def get_resource_path(pkg, path):
    """Get the full path to a resource file within a package."""
    ref = importlib_resources.files(pkg) / path
    with importlib_resources.as_file(ref) as full_path:
        if os.path.exists(full_path):
            return str(full_path)
        else:
            raise FileNotFoundError(f"Resource file {path} not found in package {pkg}.")