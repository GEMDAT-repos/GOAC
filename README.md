# GOAC: Global Optimization of Atomistic Configurations by Coulomb Energies

> **Note:** This is a temporary fork with added `pyproject.toml` / `setup.py` to enable `pip install` directly from the repository. For the original GOAC code, documentation, and issue tracking, please refer to the upstream repository at **https://iffgit.fz-juelich.de/k.koester/goac**.

GOAC is a Python-based command line tool to approach optimization of Coulomb energies for crystall structures with configurational disorders, e.g., partial site occupations. GOAC's core consists of serveral Fortran based heurstics that allow for parallelization by OpenMP. Moreover, GOAC can be used to wrap Atomistic Configuration problems to existing software such as Gurobi. Thereby, GOAC is, to the best of our knowledege, at its publication (2024) among the most powerful software for optimizing atomistic configuration problems by Coulomb energies.

Details on how to use GOAC as command line tool or directly from own Python code can be found in the Documentation PDF file and scientific details are availabe in the corresponding publication under: https://doi.org/10.1038/s41524-025-01690-7

## Installation

### Install from prebuilt wheels (recommended)

Each [GitHub Release](https://github.com/GEMDAT-repos/GOAC/releases) ships prebuilt binary wheels for CPython 3.10–3.14 on Linux (x86_64), macOS (Apple Silicon) and Windows (64-bit), so no Fortran compiler is needed. Install GOAC from the release's assets, letting pip resolve the runtime dependencies (`numpy`, `pymatgen`, `gurobipy`) from PyPI:

```sh
pip install GOAC \
    --find-links https://github.com/GEMDAT-repos/GOAC/releases/expanded_assets/vX.Y.Z
```

Replace `vX.Y.Z` with the release tag you want. The `--find-links` URL is the release's asset listing, which pip resolves the GOAC wheel from while pulling the dependencies from PyPI as usual.

> **Note:** Do not add `--no-index`. That flag tells pip to ignore PyPI entirely, but the release page hosts only GOAC — not its dependencies — so the install would fail to resolve `numpy`/`pymatgen`/`gurobipy`. For a fully offline install you would need to additionally provide those dependencies (e.g. a second `--find-links` pointing at a local wheelhouse).
>
> If your platform has no matching wheel (e.g. an Intel Mac, for which no wheel is built), pip falls back to building the source distribution, which requires a Fortran compiler — see the manual build below.

### Build from source

Install a Python3 environment and meson on your system along with a Fortran compiler (e.g., gfortran) and a corresponding OpenMP library.

### Requirements:

The required external Python packages can be installed for example via pip by:
```sh
pip install -r requirements.txt
```

The requirements contain the licensed package Gurobi and the license should be checked before installation. GOAC can be used completely without Gurobi (do not using the Gurobi solver option) and relying on the internal solvers only. When using Gurobi for optimization (setting Gurobi solver), please follow the offical licensing options, including free academic licenses, under: https://www.gurobi.com/solutions/licensing/

### Compilation of Fortran Code

If all requirements are satisfied, clone the repository and run the following commands in your copy (in case you want to use another compiler adjust settings accodingly):
```sh
cd GOAC
FC="gfortran" python3 -m numpy.f2py -c GOAC.f90 -m GOAC --backend meson --dep openmp
FC="gfortran" python3 -m numpy.f2py -c ABCEwald.f90 -m ABCEwald --backend meson --dep openmp
cd ..
```

### Test of Installation

To test your installation, run the following command:
```sh
python3 GOAC/ -f LCO.cif -p "Li*:c=1.0" -p "Co*:c=3.5" -p "O*:c=-2.0" -s "random-mc" -n 4 -w 1
```
Depending on your system you should obtain an optimized cif "out-0.cif" and a file "out-summary.txt" within a few seconds. For more advanced usage of the GOAC code please have a look at the Documentation PDF file in this repository.


## Citing GOAC for your Research

If you find GOAC helpful for your research, please considere citing the packages GOAC is working with as well as citing the GOAC code directly by the following reference: https://doi.org/10.1038/s41524-025-01690-7


## License

GOAC is released under the MIT License. The terms of the license are as follows:

The MIT License (MIT) 

Copyright 2024 GOAC

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
