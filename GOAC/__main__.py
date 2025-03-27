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
import numpy as np

def input_parser():
    global cif_file
    global charges
    global fixed_sites
    global dry_run
    global solver_name
    global solver_options
    global n_best_solutions
    global out_name
    global write
    global write_energy
    global write_part_struct
    global read_energy
    global supercell
    global scan_supercells
    charges = {}
    fixed_sites = []
    p = optparse.OptionParser()
    p.add_option('--file', '-f', type="string", help="Read cif from FILENAME.", dest="cif_file")
    p.add_option('--property', '-p', action="append", type="string",
                 help="Add properties as 'Na*=1.0' for charges or 'Na1=fixed' to fix sites.", dest="properties")
    p.add_option('--dry', '-d', help="Set to perform dry run.", dest="dry", default=False, action='store_true')
    p.add_option('--solver_options', '-y', help="File with parameters for the set solver.",
                 dest="solver_options", default=None, type="string")
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
    p.add_option('--scan_supercells', '-x', help="Scan all combinations of cell-sizes up to the specified supercell.",
                 dest="scan_supercells", default=False, action='store_true')
    p.add_option('--supercell', '-e', help="Dimensions of supercell given as dxdxd.", dest="supercell",
                 default="1x1x1")
    options, arguments = p.parse_args()
    cif_file = options.cif_file
    dry_run = options.dry
    solver_name = options.solver.casefold()
    solver_options = options.solver_options
    n_best_solutions = options.n
    out_name = options.out_name
    write = options.write
    write_energy = options.write_energy
    write_part_struct = options.write_part_struct
    read_energy = options.read_energy
    scan_supercells = options.scan_supercells
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
            "#*#       DOI: https://doi.org/10.48550/arXiv.2409.08808       #*#\n"
            "#*#                                                            #*#\n"
            "##################################################################\n"
            "##################################################################\n")

    print(line)

    global supercell
    input_parser()
    problems = []
    supercell_list = []
    if scan_supercells:
        for x in range(1, int(supercell[0])+1):
            for y in range(1, int(supercell[1])+1):
                for z in range(1, int(supercell[2])+1):
                    current_supercell = [x, y, z]
                    problem = Iteration_Problem(cif_file=cif_file, fixed_sites=fixed_sites, charges=charges,
                                                supercell=current_supercell)
                    problems.append(problem)
                    supercell_list.append(current_supercell)
    else:
        problem = Iteration_Problem(cif_file=cif_file, fixed_sites=fixed_sites, charges=charges, supercell=supercell)
        problems.append(problem)

    for i, problem in enumerate(problems):
        if scan_supercells:
            supercell_str = (str(supercell_list[i][0]) + "x" + str(supercell_list[i][1]) + "x" +
                             str(supercell_list[i][2]))
            current_out_name = out_name + "-" + supercell_str
        if not problem.occ_error:
            print(problem)
        else:
            if not scan_supercells:
                raise "It seems ocupations result in exactly xx.5 ions and mathematical rounding would cause over-occupation of the site. This case is not defined, please provide proper occupations."
            else:
                print("Skipping supercell " + supercell_str + " because occupations cannot be matched in cell shape.")
                continue

        if dry_run:
            return

        if write_part_struct:
            problem.write_partial_cifs()

        if read_energy:
            problem.read_energy()
        else:
            problem.calc_coulomb_matrices()

        if write_energy:
            problem.write_energy()

        if problem.total_combinations == 0:
            out = "Filename.cif \t\t E / eV\n"
            out += cif_file + " \t\t " + str(round(problem.const, 5)) + "\n"
            f = open(out_name + "-summary.txt", "w")
            f.write(out)
            f.close()
            print("Nothing to do here. The provided cif file has no iterative sites.")
            return

        match solver_name:
            case "gurobi":
                solver = Gurobi_Solver(name="Gurobi", problem=problem, n=n_best_solutions, w=write)
            case "gurobi-heuristic":
                solver = Gurobi_Solver(name="Gurobi-Heuristic", problem=problem, n=n_best_solutions, w=write)
            case "greedy":
                solver = Greedy_Solver(name="Greedy", problem=problem, n=n_best_solutions, w=write)
            case "random":
                solver = Random_Solver(name="Random", problem=problem, n=n_best_solutions, w=write)
            case "random-bb":
                solver = Random_Solver(name="Random-BB", problem=problem, n=n_best_solutions, w=write)
            case "random-lm":
                solver = Random_Solver(name="Random-LM", problem=problem, n=n_best_solutions, w=write)
            case "random-mc":
                solver = Random_Solver(name="Random-MC", problem=problem, n=n_best_solutions, w=write)
            case "random-sa":
                solver = Random_Solver(name="Random-SA", problem=problem, n=n_best_solutions, w=write)
            case "random-remc":
                solver = Random_Solver(name="Random-REMC", problem=problem, n=n_best_solutions, w=write)
            case "random-ga":
                solver = Random_Solver(name="Random-GA", problem=problem, n=n_best_solutions, w=write)
            case "random-rega":
                solver = Random_Solver(name="Random-REGA", problem=problem, n=n_best_solutions, w=write)
            case "random-hybrid":
                solver = Random_Solver(name="Random-Hybrid", problem=problem, n=n_best_solutions, w=write)
            case _:
                raise ("The given solver name is unknown. Chose from: Gurobi, Gurobi-Heuristic, Greedy, Random, Random-BB, Random-LM, Random-MC, Random-SA, Random-REMC, Random-GA, Random-REGA, Random-Hybrid")


        solver.initialize(options=solver_options)
        if scan_supercells:
            solver.solve(out_name=current_out_name)
        else:
            solver.solve(out_name=out_name)
    if scan_supercells:
        files = []
        energies = []
        for supercell in supercell_list:
            supercell_str = (str(supercell[0]) + "x" + str(supercell[1]) + "x" + str(supercell[2]))
            try:
                f = open(out_name + "-" + supercell_str + "-summary.txt", "r")
                for line in f.readlines():
                    line = line.strip().split()
                    if len(line) == 2:
                        files.append(line[0])
                        energies.append(round(float(line[1])/np.prod(supercell), 5))
                f.close()
            except IOError:
                print("Did not found file: " + out_name + "-" + supercell_str)
        if len(energies) > 0 and len(files) > 0:
            pairs = list(zip(energies, files))
            pairs.sort(key=lambda p: (p[0], np.prod(np.array(p[1].split("-")[1].split("x")).astype(np.int64))))

            out = "Filename.cif \t\t E / eV\n"
            #for e, f in sorted(zip(energies, files)):
            for e, f in pairs:
                out += f + " \t\t " + str(round(e, 5)) + "\n"
            f = open("Supercell-Scan-Report.txt", "w")
            f.write(out)
            f.close()
            print("Supercell Scan Report written.")

if __name__ == '__main__':
    main()
