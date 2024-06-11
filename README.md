# GOAC: Global Optimization of Atomistic Arrangements by Coulomb Energies

GOAC is a Python-based command line tool to approach optimization of Coulomb energies for crystall structures with configurational disorders, e.g., partial site occupations. GOAC's core consists of serveral Fortran based heurstics that allow for parallelization by OpenMP. Moreover, GOAC can be used to wrap Atomistic Arrangement problems to existing software such as Gurobi. Details on how to use GOAC as command line tool or as Python package can be found in the Documentation PDF file and scientific details are availabe in the corresponding publication under: [DOI]

## Installation
Install a Python3.10 or higher Python-environment on your system along with a Fortran compiler (e.g., gfortran7.5.0) and a corresponding OpenMP library.

### Requirements:
Some Python package are required before installing GOAC.
```sh
pip install pymatgen
```
