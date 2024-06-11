__author__ = "Konstantin Köster"
__copyright__ = "Copyright 2024, GOAC"
__credits__ = ["Konstantin Köster", "Tobias Binninger", "Payam Kaghazch"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = ""
__email__ = "p.kaghazchi@fz-juelich.de"
__status__ = "Development"

from GOAC import greedy
from GOAC import branch_n_bound
from GOAC import local_minimizer
from GOAC import monte_carlo
from GOAC import remc
import numpy as np
from Solver import Solver
from IterationProblem import Iteration_Problem
import time

class Greedy_Solver(Solver):
    def __init__(self, name:str, problem:Iteration_Problem, cores:int, n=1, w=False):
        super().__init__(name=name, problem=problem, cores=cores, n=n, w=w)

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
            self.opt["stop_time"] = -1
            self.opt["stop_steps_no_improve"] = -1
            if self.name == "Greedy-MC":
                self.opt["mc_steps"] = 10000000
                self.opt["mc_kT"] = 0.5
                self.opt["mc_sim_an"] = 1.0
                self.opt["mc_sim_an_steps"] = 10000001
                self.opt["mc_write_steps"] = 100000
            elif self.name == "Greedy-SA":
                self.opt["mc_steps"] = 10000000
                self.opt["mc_kT"] = 5
                self.opt["mc_sim_an"] = 0.99
                self.opt["mc_sim_an_steps"] = 20000
                self.opt["mc_write_steps"] = 100000
            elif self.name == "Greedy-REMC":
                self.opt["mc_steps"] = 100000
                self.opt["re_repeat"] = 10
                self.opt["mc_kT"] = [0.1, 0.2, 0.5]

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
                                               self.n)
        if self.name == "Greedy-BB":
            lowest_energies, lowest_solutions = branch_n_bound(self.species_nums, self.species_occs, self.shared_sites,
                                                               best_solutions, self.problem.const, self.np_alpha,
                                                               self.problem.beta, self.opt["samples"])
        elif self.name == "Greedy-LM":
            lowest_energies, lowest_solutions = local_minimizer(self.species_nums, self.species_occs, self.shared_sites,
                                                                best_solutions, self.problem.const, self.np_alpha,
                                                                self.problem.beta, self.opt["samples"],
                                                                self.opt["stop_time"], self.opt["stop_steps_no_improve"])
        elif self.name == "Greedy-MC" or self.name == "Greedy-SA":
            lowest_energies, lowest_solutions, cs = monte_carlo(self.species_nums, self.species_occs, self.shared_sites,
                                                                best_solutions, self.problem.const, self.np_alpha,
                                                                self.problem.beta, self.opt["samples"],
                                                                self.opt["mc_steps"], self.opt["mc_kt"],
                                                                self.opt["mc_sim_an"], self.opt["mc_sim_an_steps"],
                                                                self.opt["mc_write_steps"], self.opt["stop_time"],
                                                                self.opt["stop_steps_no_improve"])
        elif self.name == "Greedy-REMC":
            lowest_energies, lowest_solutions = remc(self.species_nums, self.species_occs, self.shared_sites,
                                                     best_solutions, self.problem.const, self.np_alpha,
                                                     self.problem.beta, self.opt["samples"],
                                                     self.opt["mc_steps"], self.opt["re_repeat"],
                                                     self.opt["mc_kt"], self.opt["stop_time"],
                                                     self.opt["stop_steps_no_improve"])
        if self.name in ["Greedy-BB", "Greedy-LM", "Greedy-MC", "Greedy-SA", "Greedy-REMC"]:
            out_energies = np.zeros(max(self.n, self.write))
            shape = list(best_solutions.shape)
            shape[0] = max(self.n, self.write)
            out_solutions = np.zeros(tuple(shape))
            out_energies[:self.n] = lowest_energies
            out_solutions[:self.n, :] = lowest_solutions
            if self.opt["samples"] > 1:
                for j in range(self.opt["samples"]):
                    for i in range(self.n):
                        if lowest_energies[i, j] < np.max(out_energies):
                            out_solutions[np.argmax(out_energies)] = lowest_solutions[i, j]
                            out_energies[np.argmax(out_energies)] = lowest_energies[i, j]
        else:
            out_energies = best_energies
            out_solutions = best_solutions

        print("Solver finished in " + str(round(time.time()-start,0)) + "s")

        data = sorted(zip(out_energies, out_solutions), key=lambda x: x[0])
        out_energies, out_solutions = list(zip(*data))

        # Getting n-best solutions, writing their cif files and a summary txt
        out = "Filename.cif \t\t E / eV\n"

        for sol in range(len(out_energies)):
            name = out_name + "-" + str(sol) + ".cif"
            out += name + " \t\t " + str(round(out_energies[sol], 5)) + "\n"
            if sol < self.write:
                self.problem.print_solution_cif_fotran(name=name, x=out_solutions[sol])

        # Write summary tx with cif-file names and energies
        f = open(out_name + "-summary.txt", "w")
        f.write(out)
        f.close()
