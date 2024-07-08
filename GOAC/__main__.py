__author__ = "Konstantin Köster"
__copyright__ = "Copyright 2024, GOAC"
__credits__ = ["Konstantin Köster", "Tobias Binninger", "Payam Kaghazch"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = ""
__email__ = "p.kaghazchi@fz-juelich.de"
__status__ = "Development"

import optparse
from IterationProblem import Iteration_Problem
from GurobiSolver import Gurobi_Solver
from GreedySolver import Greedy_Solver
from RandomSolver import Random_Solver

def input_parser():
    global cif_file
    global charges
    global fixed_sites
    global dry_run
    global cores
    global solver_name
    global solver_options
    global n_best_solutions
    global out_name
    global write
    global write_energy
    global write_part_struct
    global read_energy
    global supercell
    charges = {}
    fixed_sites = []
    p = optparse.OptionParser()
    p.add_option('--file', '-f', type="string", help="Read cif from FILENAME.", dest="cif_file")
    p.add_option('--property', '-p', action="append", type="string",
                 help="Add properties as 'Na*=1.0' for charges or 'Na1=fixed' to fix sites.", dest="properties")
    p.add_option('--dry', '-d', help="Set to perform dry run.", dest="dry", default=False, action='store_true')
    p.add_option('--solver_options', '-y', help="File with parameters for the set solver.",
                 dest="solver_options", default=None, type="string")
    p.add_option('--cores', '-c', help="Number of threads the code will use", dest="cores", default=1, type="int")
    p.add_option('--solver', '-s', help="Name of solver to solve given Coulomb problem", dest="solver",
                 default="", type="string")
    p.add_option('--number', '-n', help="The 'n' best soultions are returned", dest="n",
                 default=1, type="int")
    p.add_option('--out-name', '-o', help="The out-name is written in front of all output-files", dest="out_name",
                 default="out", type="string")
    p.add_option('--write', '-w', help="Set number of cif files that will be written", dest="write",
                 default=0, type="int")
    p.add_option('--write_energy', '-a', help="Write files with energy data.", dest="write_energy",
                 default=False, action='store_true')
    p.add_option('--write_part_struct', '-b', help="Write structures to calculate energy data",
                 dest="write_part_struct", default=False, action='store_true')
    p.add_option('--read_energy', '-r', help="Read energies from files.", dest="read_energy",
                 default=False, action='store_true')
    p.add_option('--supercell', '-e', help="Dimensions of supercell given as dxdxd.", dest="supercell",
                 default="1x1x1")
    options, arguments = p.parse_args()
    cif_file = options.cif_file
    dry_run = options.dry
    cores = options.cores
    solver_name = options.solver.casefold()
    solver_options = options.solver_options
    n_best_solutions = options.n
    out_name = options.out_name
    write = options.write
    write_energy = options.write_energy
    write_part_struct = options.write_part_struct
    read_energy = options.read_energy
    supercell = options.supercell.split("x")
    for i in range(len(supercell)):
        supercell[i] = float(supercell[i])

    #parse properties such as chrages and fixed sites
    if options.properties is None:
        return
    for property in options.properties:
        if "=" not in property:
            raise "Invalid property given. Properties must contain '='."
        property = property.split("=")
        if "fixed" == property[1]:
            fixed_sites.append(property[0])
        else:
            try:
                charges[property[0]] = float(property[1])
            except ValueError as verr:
                raise "Properties must be either 'fixed' for fixing sites or floats for charges."
    return


def main():

    line = ("##################################################################\n"
            "##################################################################\n"
            "#*#                                                            #*#\n"
            "#*#                   o-o    o-o     O      o-o                #*#\n" 
            "#*#                  o      o   o   / \    /                   #*#\n" 
            "#*#                  |  -o  |   |  o---o  O                    #*#\n" 
            "#*#                  o   |  o   o  |   |   \                   #*#\n" 
            "#*#                   o-o    o-o   o   o    o-o                #*#\n"
            "#*#                                                            #*#\n"
            "#*# Global Optimization of Atomistic Configurations by Coulomb #*#\n"
            "#*#                                                            #*#\n"
            "#*#  Konstantin Köster, Tobias Binninger, and Payam Kaghazchi  #*#\n"
            "#*#              DOI:                                          #*#\n"
            "#*#                                                            #*#\n"
            "##################################################################\n"
            "##################################################################\n")

    print(line)

    input_parser()
    problem = Iteration_Problem(cif_file=cif_file, fixed_sites=fixed_sites, charges=charges, supercell=supercell)
    print(problem)

    if dry_run:
        return

    if read_energy:
        problem.read_energy()
    else:
        problem.calc_coulomb_matrices(cores=cores, write_part_struct=write_part_struct)
    if write_energy:
        problem.write_energy()


    match solver_name:
        case "gurobi":
            solver = Gurobi_Solver(name="Gurobi", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "gurobi-heuristic":
            solver = Gurobi_Solver(name="Gurobi-Heuristic", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "greedy":
            solver = Greedy_Solver(name="Greedy", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "greedy-bb":
            solver = Greedy_Solver(name="Greedy-BB", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "greedy-lm":
            solver = Greedy_Solver(name="Greedy-LM", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "greedy-mc":
            solver = Greedy_Solver(name="Greedy-MC", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "greedy-sa":
            solver = Greedy_Solver(name="Greedy-SA", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "greedy-remc":
            solver = Greedy_Solver(name="Greedy-REMC", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random":
            solver = Random_Solver(name="Random", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-bb":
            solver = Random_Solver(name="Random-BB", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-lm":
            solver = Random_Solver(name="Random-LM", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-mc":
            solver = Random_Solver(name="Random-MC", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-sa":
            solver = Random_Solver(name="Random-SA", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-remc":
            solver = Random_Solver(name="Random-REMC", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-ga":
            solver = Random_Solver(name="Random-GA", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-rega":
            solver = Random_Solver(name="Random-REGA", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case "random-hybrid":
            solver = Random_Solver(name="Random-Hybrid", problem=problem, cores=cores, n=n_best_solutions, w=write)
        case _:
            raise "The given solver name is unknown. Chose from: [...]" # TODO fill solver


    solver.initialize(options=solver_options)
    solver.solve(out_name=out_name)


if __name__ == '__main__':
    main()
