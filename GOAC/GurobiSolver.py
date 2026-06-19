__author__ = "Konstantin Köster"
__copyright__ = "Copyright 2024, GOAC"
__credits__ = ["Konstantin Köster", "Tobias Binninger", "Payam Kaghazch"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = ""
__email__ = "p.kaghazchi@fz-juelich.de"
__status__ = "Development"

from Solver import Solver
from IterationProblem import Iteration_Problem
import gurobipy as gp
from gurobipy import GRB
import numpy as np


class Gurobi_Solver(Solver):
    def __init__(self, name:str, problem:Iteration_Problem, n=1, w=False):
        super().__init__(name=name, problem=problem, n=n, w=w)

    def initialize(self, options:str=None):
        opt = {}
        if options is not None:
            #set user specific opt-dict
            f = open(options, "r")
            for line in f.readlines():
                line = line.strip().split()
                opt[line[0]] = float(line[1])
                if abs(int(float(line[1]))-float(line[1])) < 10E-7:
                    opt[line[0]] = int(float(line[1]))
        elif self.name == "Gurobi-Heuristic":
            #set default-opt for heuristic
            opt['Presolve'] = 2
            opt['MIPFocus'] = 1
            opt['MIPGap'] = 0.0
            opt['MIPGapAbs'] = 0.0
            opt['Threads'] = 4
            opt['PoolSearchMode'] = 2
            opt['PoolSolutions'] = self.n
            opt['TimeLimit'] = 300
            opt['Heuristics'] = 0.99999
            opt['NoRelHeurTime'] = 300
        else:
            #set default-opt dict
            opt['Presolve'] = 2
            opt['MIPFocus'] = 2
            opt['MIPGap'] = 0.0
            opt['MIPGapAbs'] = 0.0
            opt['Threads'] = 4
            opt['PoolSearchMode'] = 2
            opt['PoolSolutions'] = self.n

        #Initialize Gurobi model
        self.model = gp.Model("Coulomb")
        #set options
        for key, value in opt.items():
            self.model.setParam(key, value)

        # Create gurobi variables of the problem
        # VARIABLES x_ij is occ [0,1] of site i in position j
        tuple_list = []
        i = -1
        for it_site in self.problem.iterate_sites.values():
            for el in it_site.site.species.get_el_amt_dict().keys():
                i += 1
                for j in range(it_site.num):
                    tuple_list.append((i,j))
        self.x = self.model.addVars(tuple_list, vtype=GRB.BINARY)

        #Add Constraints to gurobi model
        #Constraint that each iterative site has the correct integer occupation
        i = -1
        site_nums = []
        for it_site in self.problem.iterate_sites.values():
            for el, occ in it_site.site.species.get_el_amt_dict().items():
                i += 1
                site_nums.append(it_site.num)
                self.model.addConstr(self.x.sum(i,'*') == int(it_site.num*occ+0.5))

        #Constrain to avoid doubl occupancies in iterative site with multiple species
        i_counter = -1
        for it_site in self.problem.iterate_sites.values():
            i_counter += 1
            if len(it_site.site.species.get_el_amt_dict().keys()) > 1:
                shared_i = np.arange(i_counter,i_counter+len(it_site.site.species.get_el_amt_dict().keys()))
                i_counter += len(it_site.site.species.get_el_amt_dict().keys())-1
                for j in range(it_site.num):
                    if abs(np.sum(list(it_site.site.species.get_el_amt_dict().values()))-1) <= 10E-5:
                        self.model.addConstr(gp.quicksum(self.x[i,j] for i in shared_i) == 1)
                    else:
                        self.model.addConstr(gp.quicksum(self.x[i,j] for i in shared_i) <= 1)


        #set objective function for model
        objective = (self.problem.const
                     + gp.quicksum(self.x[i,j]*self.problem.alpha[i][j]
                                   for i in range(len(site_nums))
                                   for j in range(site_nums[i]))
                     + gp.quicksum(self.problem.beta[i,j,i,l]*self.x[i,j]*self.x[i,l]
                                   for i in range(len(site_nums))
                                   for j in range(site_nums[i])
                                   for l in range(j+1, site_nums[i]))
                     + gp.quicksum(self.problem.beta[i,j,k,l]*self.x[i,j]*self.x[k,l]
                                   for i in range(len(site_nums))
                                   for j in range(site_nums[i])
                                   for k in range(i+1, len(site_nums))
                                   for l in range(site_nums[k])))

        self.model.setObjective(objective, sense=GRB.MINIMIZE)

        #Add all updates to the model and print it
        self.model.update()
        self.model.printStats()
        self.model.write("model.mps")

    def solve(self, out_name:str):
        #Gurobi performs optimization
        self.model.optimize()

        #Getting n-best solutions, writing their cif files and a summary txt
        out = "Filename.cif \t\t E / eV\n"
        is_nan = False
        for sol in range(self.model.SolCount):
            self.model.setParam(GRB.Param.SolutionNumber, sol)
            x = self.model.Xn
            name = out_name + "-" + str(sol) + ".cif"
            if sol < self.write:
                if not self.problem.print_solution_cif_gurobi(name=name, x=x):
                    is_nan = True
            else:
                name = "NaN"
            if not is_nan:
                out += name + " \t\t " + str(round(self.model.PoolObjVal, 5)) + "\n"
            else:
                out += ("No more solutions to print. "
                        "To obtain more solutions, increase PoolSolutions in the "
                        "options file. Consider also that maybe not more solutions exist.")
                break

        #Write summary tx with cif-file names and energies
        f = open(out_name + "-summary.txt", "w")
        f.write(out)
        f.close()
