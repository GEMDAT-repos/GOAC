__author__ = "Konstantin Köster"
__copyright__ = "Copyright 2024, GOAC"
__credits__ = ["Konstantin Köster", "Tobias Binninger", "Payam Kaghazch"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = ""
__email__ = "p.kaghazchi@fz-juelich.de"
__status__ = "Development"

from GOAC import greedy
import numpy as np
from Solver import Solver
from IterationProblem import Iteration_Problem
import time

class Greedy_Solver(Solver):
    def __init__(self, name:str, problem:Iteration_Problem, n=1, w=False):
        super().__init__(name=name, problem=problem, n=n, w=w)

    def initialize(self, options:str=None):
        self.opt = {}
        if options is not None:
            # set user specific opt-dict
            f = open(options, "r")
            for line in f.readlines():
                line = line.strip().split()
                if line[0].casefold() == "mc_kt":
                    if len(line) > 2:
                        self.opt[line[0].casefold()] = []
                        for i in range(1, len(line)):
                            self.opt[line[0].casefold()].append(float(line[i]))
                    else:
                        self.opt[line[0].casefold()] = float(line[1])
                else:
                    self.opt[line[0].casefold()] = float(line[1])
                    if abs(int(float(line[1])) - float(line[1])) < 10E-7:
                        self.opt[line[0].casefold()] = int(float(line[1]))
        else:
            self.opt["samples"] = 1
        if "tol" not in self.opt.keys():
            self.opt["tol"] = 1e-08

        self.species_nums = []
        self.species_occs = []
        self.shared_sites = []
        num_species = 0
        max_site_num = 0
        for it_site in self.problem.iterate_sites.values():
            for el in it_site.site.species.get_el_amt_dict().keys():
                num_species += 1
                if it_site.num > max_site_num:
                    max_site_num = it_site.num
        i = 0
        self.np_alpha = np.zeros((num_species, max_site_num))
        for it_site in self.problem.iterate_sites.values():
            start = i
            end = i+len(it_site.site.species.get_el_amt_dict().keys())
            for el, occ in it_site.site.species.get_el_amt_dict().items():
                current_shared = [False,]*num_species
                for l in range(start, end):
                    current_shared[l] = True
                self.shared_sites.append(current_shared)
                self.species_nums.append(it_site.num)
                self.species_occs.append(int(it_site.num*occ+0.5))
                for j in range(it_site.num):
                    self.np_alpha[i,j] = self.problem.alpha[i][j]
                i += 1

        self.shared_sites = np.array(self.shared_sites)
        self.species_occs = np.array(self.species_occs)
        self.species_nums = np.array(self.species_nums)


    def solve(self, out_name: str):
        #Solving the model by greedy fortran implementation
        start = time.time()
        best_energies, best_solutions = greedy(self.species_nums, self.species_occs, self.shared_sites,
                                               self.problem.const, self.np_alpha, self.problem.beta,
                                               self.opt["samples"], self.opt["tol"])
        out_energies = best_energies
        out_solutions = best_solutions

        print("Solver finished in " + str(round(time.time()-start,0)) + "s")

        data = sorted(zip(out_energies, out_solutions), key=lambda x: x[0])
        out_energies, out_solutions = list(zip(*data))

        # Getting n-best solutions, writing their cif files and a summary txt
        out = "Filename.cif \t\t E / eV\n"

        is_nan = False
        for sol in range(len(out_energies)):
            name = out_name + "-" + str(sol) + ".cif"
            if sol < self.write:
                if not self.problem.print_solution_cif_fotran(name=name, x=out_solutions[sol]):
                    is_nan = True
            else:
                name = "NaN"
            if not is_nan:
                if out_energies[sol] < 10E100:
                    out += name + " \t\t " + str(round(out_energies[sol], 5)) + "\n"
                else:
                    is_nan = True
            else:
                out += ("No more solutions to print. "
                        "To obtain more solutions, increase samples in the options file or "
                        "deactivate/change tol in the options file. tol -1 prints also solutions "
                        "that are same in energy.")
                break

        # Write summary tx with cif-file names and energies
        f = open(out_name + "-summary.txt", "w")
        f.write(out)
        f.close()
