from IterationProblem import Iteration_Problem

class Solver():
    def __init__(self, name:str, problem:Iteration_Problem, cores:int, n=1, w=False):
        self.name = name
        self.problem = problem
        self.cores = cores
        self.n = n
        self.write = w

        print("Creating solver of type: " + self.name)

    def initialize(self, options:str=None):
        """Function to map a coulomb problem to the specific solver
        and to set default or user-parameters for the solver"""
        return

    def solve(self, out_name:str):
        """Function that actually solves a coulomb problem and returns a structure+energy or a list of both"""
        return