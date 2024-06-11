# GOAC: Global Optimization of Atomistic Arrangements by Coulomb Energies

GOAC is a Python-based command line tool to approach optimization of Coulomb energies for crystall structures with configurational disorders, e.g., partial site occupations. GOAC's core consists of serveral Fortran based heurstics that allow for parallelization by OpenMP. Moreover, GOAC can be used to wrap Atomistic Arrangement problems to existing software such as Gurobi. Thereby, GOAC is, to the best of our knowledege, at its publication (2024) among the most powerful software for optimizing atomistic arrangement problems by Coulomb energies.

Details on how to use GOAC as command line tool or as Python package can be found in the Documentation PDF file and scientific details are availabe in the corresponding publication under: [DOI]

## Installation

Install a Python3.10, or higher, Python-environment on your system along with a Fortran compiler (e.g., gfortran7.5.0) and a corresponding OpenMP library.

### Requirements:

Next to standard Python packages such as time, warnings, re, copy, and optparse some more Python package are required before installing GOAC:
```sh
pip install joblib
pip install numpy
pip install pymatgen
pip install gurobipy
```

The functionalities of Gurobi required to run GOAC are freely available but corresponding licences should be checked. For using Gorbi for optimization, please follow the offical licensing options, including free academic licenses, under: https://www.gurobi.com/solutions/licensing/

### Compilation of Fortran Code

If all requirements are satisfied, clone the repository and run the following commands in your copy (in case you want to use another compiler adjust settings accodingly):
```sh
cd GOAC
python -m numpy.f2py -c --fcompiler=gfortran --f90flags='-fopenmp -Wall' -m GOAC.f90
cd ..
```

### Test of Installation

To test your installation, run the following command:
```sh
python GOAC/ -f LCO.cif -p "Li*=1.0" -p "Co*=3.5" -p "O*=-2.0" -s "random-mc" -n 4 -w 1 -c 4
```
Depending on your system you should obtain an optimized cif "out-0.cif" and a file "out-summary.txt" within a few seconds. For more advanced usage of the GOAC code please have a look at the Documentation PDF file in this repository.


## Citing GOAC for your Research

If you find GOAC helpful for your research, please considere citing the packages GOAC is working with as well as citing the GOAC code directly by the following reference:
[DOI]


## License

GOAC is released under the MIT License. The terms of the license are as follows:

The MIT License (MIT) Copyright (c) 2011-2012 MIT & LBNL

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
