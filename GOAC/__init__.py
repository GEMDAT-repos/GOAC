"""
GOAC: Global Optimization of Atomistic Configurations by Coulomb Energies.

This package provides solvers for optimizing Coulomb energies of crystal
structures with configurational disorder (partial site occupations).
"""

import sys
import os as _os

# ----- sys.path hack for implicit relative imports -----
# The existing Python files use imports like:
#   from IterationProblem import Iteration_Problem
# instead of the proper relative form:
#   from .IterationProblem import Iteration_Problem
#
# We add the package directory to sys.path so that these
# implicit-absolute imports resolve correctly at runtime.
_pkg_dir = _os.path.dirname(__file__)
if _pkg_dir not in sys.path:
    sys.path.insert(0, _pkg_dir)

# ----- Re-export Fortran functions from _goac_fortran -----
# GOAC.f90 compiles to a Python extension module named "_goac_fortran"
# (renamed to avoid clashing with this Python package directory).
#
# The existing solver code does:
#   from GOAC import greedy
#   from GOAC import random_samples, branch_n_bound, ...
#
# By importing everything from _goac_fortran into this __init__.py,
# those imports continue to work unchanged.

try:
    from GOAC._goac_fortran import (  # type: ignore[import-untyped]
        greedy,
        random_samples,
        branch_n_bound,
        local_minimizer,
        monte_carlo,
        remc,
        ga,
        rega,
        energy,
        occupied,
        solution_unique,
        valid_solutions,
    )
except ImportError as _exc:
    # Allow import even if native extensions are missing
    # (e.g. during documentation builds or on unsupported platforms).
    import warnings as _warnings
    _warnings.warn(
        f"Could not load compiled Fortran extension (_goac_fortran): {_exc}",
        ImportWarning,
    )


# ----- Package metadata -----
__author__ = "Konstantin Köster"
__copyright__ = "Copyright 2024, GOAC"
__credits__ = ["Konstantin Köster", "Tobias Binninger", "Payam Kaghazch"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = ""
__email__ = "p.kaghazchi@fz-juelich.de"
__status__ = "Development"
