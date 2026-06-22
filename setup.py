"""Setup script for GOAC.

Builds the Fortran extension modules (GOAC.f90 and ABCEwald.f90)
using numpy.f2py with the meson backend.
"""
import os
import sys
import subprocess
import shutil
import glob

import setuptools
from setuptools import setup
from setuptools.command.build_ext import build_ext


F2PY_EXTENSIONS = {
    "GOAC._goac_fortran": {
        "source": "GOAC/GOAC.f90",
        "dest_is_pkg": True,  # install inside GOAC/ package dir
    },
    "ABCEwald": {
        "source": "GOAC/ABCEwald.f90",
        "dest_is_pkg": False,  # install at top level
    },
}


class FortranBuildExt(build_ext):
    """Custom build command to compile Fortran sources via numpy.f2py.

    Overrides build_extension() so that setuptools never tries to compile
    the .f90 files itself — we delegate everything to numpy.f2py instead.
    """

    def run(self):
        # Ensure numpy is importable
        try:
            import numpy  # noqa: F401
        except ImportError:
            raise RuntimeError(
                "numpy is required to build the Fortran extensions. "
                "Install it first: pip install numpy"
            )
        # Do NOT call super().run() — we handle all extensions ourselves.
        self.build_all_fortran_extensions()

    def build_extension(self, ext):
        """Called once per Extension.  We skip setuptools' logic entirely."""
        pass  # handled by build_all_fortran_extensions

    def build_all_fortran_extensions(self):
        """Build every Fortran extension using f2py with the meson backend."""
        build_dir = os.path.abspath(self.build_temp)
        os.makedirs(build_dir, exist_ok=True)

        for ext_name, info in F2PY_EXTENSIONS.items():
            source = os.path.abspath(info["source"])
            module_name = ext_name.split(".")[-1]  # "_goac_fortran" or "ABCEwald"

            print("*" * 60, file=sys.stderr)
            print(f"Building {ext_name} from {os.path.basename(source)}...",
                  file=sys.stderr)
            print("*" * 60, file=sys.stderr)

            self._run_f2py(source=source, module_name=module_name,
                           work_dir=build_dir)

            so_file = self._find_so(module_name, build_dir)
            if not so_file:
                raise RuntimeError(
                    f"Could not find compiled {module_name} extension"
                )

            if info["dest_is_pkg"]:
                dest_dir = os.path.join(
                    os.path.abspath(self.build_lib), "GOAC"
                )
            else:
                dest_dir = os.path.abspath(self.build_lib)
            os.makedirs(dest_dir, exist_ok=True)
            dest = os.path.join(dest_dir, os.path.basename(so_file))
            shutil.copy2(so_file, dest)
            os.remove(so_file)
            print(f"  -> Installed to {dest}", file=sys.stderr)

    def _run_f2py(self, source, module_name, work_dir):
        """Run numpy.f2py to compile a Fortran file into a Python extension."""
        # f2py outputs the .so into the current working directory
        old_cwd = os.getcwd()
        os.chdir(work_dir)
        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "numpy.f2py",
                    "-c",
                    source,
                    "-m",
                    module_name,
                    "--backend",
                    "meson",
                    "--dep",
                    "openmp",
                ],
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                print(result.stdout, file=sys.stderr)
                print(result.stderr, file=sys.stderr)
                raise RuntimeError(
                    f"f2py failed for {os.path.basename(source)} "
                    f"(module: {module_name})"
                )
            # Print any warnings from stderr
            if result.stderr and "Warning" in result.stderr:
                print(result.stderr, file=sys.stderr)
        finally:
            os.chdir(old_cwd)

    def _find_so(self, module_name, search_dir):
        """Find a .so file matching the given module name."""
        pattern = os.path.join(search_dir, f"{module_name}*.so")
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
        # Also search one level deeper (meson puts output in subdirs sometimes)
        for root, dirs, files in os.walk(search_dir):
            for f in files:
                if f.startswith(module_name) and f.endswith(".so"):
                    return os.path.join(root, f)
        return None


# Define dummy extension modules so setuptools invokes build_ext.
# The actual Fortran compilation is handled by FortranBuildExt.
_dummy_ext_modules = [
    setuptools.Extension("GOAC._goac_fortran", sources=["GOAC/GOAC.f90"]),
    setuptools.Extension("ABCEwald", sources=["GOAC/ABCEwald.f90"]),
]

setup(
    cmdclass={"build_ext": FortranBuildExt},
    ext_modules=_dummy_ext_modules,
)
