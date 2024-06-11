import copy

from pymatgen.core import Structure
import warnings
import re
import numpy as np
import copy
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.cif import CifWriter
from joblib import Parallel, delayed
import time

class Iterate_Site():
    """Simple class to store  pymatgen site objects along with the number of appearances in the structure
    and coordinates of all site positions as list"""
    def __init__(self, site, num:int):
        self.site = site
        self.num = num
        self.coords = [site.coords,]
class Iteration_Problem():
    """Class to store a iteration probelm:
    1. Cif file is read
    2. iterative sites are identified and stored in iterate_sites
    3. Helper functions for assigning charges (get_oxidation_states),
    calculating ewald (calc_coulomb_matrices), and printing (__str__)"""

    def __init__(self, cif_file:str, fixed_sites:list, charges:dict, supercell:list):
        self.cif_file = cif_file
        self.fixed_sites = fixed_sites
        self.charges = charges

        #Read cif-file with pymatgen
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.struct = Structure.from_file(self.cif_file)
            self.struct = self.struct.make_supercell(supercell)

        #Identify iteratie sites
        self.iterate_sites = {}
        for site in self.struct.sites:
            #check if user gave site as fixed
            fixed = False
            for el, fra in site.species.get_el_amt_dict().items():
                #Create label for different elements as shared sites have only
                # label of last element+ number by default
                label_num = re.findall(r'\d+', site.label)[-1]
                label = el + label_num
                for l in self.fixed_sites:
                    #check for wildcard match
                    if "*" in l:
                        l = l.replace("*", "")
                        if l in label:
                            fixed = True
                    #check for exact match
                    else:
                        if l == label:
                            fixed = True
            #Identify sites that are not fixed and have combinatios to iterate over and store them
            if not fixed and (len(site.species.as_dict()) > 1 or
                              abs(list(site.species.as_dict().values())[0]-1.0) > 10e-5):
                #Raise error if occupations of one site are more than 1
                if sum(list(site.species.as_dict().values())) > 1.0:
                    raise "Site " + site.label + " has an occupation greater 1."
                if site.label not in self.iterate_sites.keys():
                    self.iterate_sites[site.label] = Iterate_Site(site, 1)
                else:
                    self.iterate_sites[site.label].num += 1
                    self.iterate_sites[site.label].coords.append(site.coords)


    def calc_coulomb_matrices(self, cores:int, write_part_struct:bool):
        """Function to calculate the different Coulomb contributions Const, Self, Alpha, Beta"""
        self.cores = cores

        def calc_const_energy():
            """Function to calculate constant energy of fixed sites"""
            #copy structure to make manipulations in copy
            struct = copy.deepcopy(self.struct)
            #delete all iterative sites
            delete_sites = []
            for i, site in enumerate(struct.sites):
                if site.label in self.iterate_sites.keys():
                    delete_sites.append(i)
            struct.remove_sites(delete_sites)

            #set oxidation states and set occupancies to 1.0 as partial occupations of fixed
            # sites are already accounted for by average oxidation states
            site_oxidations = self.get_oxidation_states(struct)
            for site in struct.sites:
                site.species = {list(site.species.get_el_amt_dict().keys())[0]: 1.0}
            struct.add_oxidation_state_by_site(site_oxidations)
            out = 0
            if len(struct.sites) > 1:
                #Calc ewald
                es = EwaldSummation(struct)
                if write_part_struct:
                    cifWriter = CifWriter(struct)
                    cifWriter.write_file("const.cif")
                out = es.total_energy
            self.fixed_struct = copy.deepcopy(struct)
            self.fixed_site_oxidations = copy.deepcopy(site_oxidations)
            return out

        def calc_self_and_alpha_energy():
            """Function to self energies of iterative sites along with their interactions to fixed sites"""
            # copy structure to make manipulations in copy
            struct = copy.deepcopy(self.struct)
            self_e = []
            alpha = []
            #Iterate over all iterative sites and all elements per site
            j = -1
            for l, it_site in self.iterate_sites.items():
                for el in it_site.site.species.get_el_amt_dict().keys():
                    j += 1
                    #create label for adding charge states
                    label_num = re.findall(r'\d+', l)[-1]
                    label = el + label_num
                    self_e.append([])
                    alpha.append([])
                    #Iterate over all site positions in parallel
                    def calc_s_a_parallel(i):
                        #create new structure with only one site
                        new_struct = Structure(lattice=struct.lattice, species=[{el:1.0,},], labels=[label,],
                                               coords=[it_site.coords[i],], coords_are_cartesian=True)
                        #add oxidation states
                        site_oxidations = self.get_oxidation_states(new_struct)
                        oxidation_state = site_oxidations[0]
                        new_struct.add_oxidation_state_by_site(site_oxidations)
                        #calc and append Ewald of structure
                        es = EwaldSummation(new_struct)
                        s = es.total_energy
                        if write_part_struct:
                            cifWriter = CifWriter(new_struct)
                            cifWriter.write_file("self-" + str(i) + ".cif")

                        #create structure with site and all fixed sites
                        a_struct = copy.deepcopy(self.fixed_struct)
                        a_struct.append(species={el:1.0,}, coords=it_site.coords[i], coords_are_cartesian=True)
                        #add oxidation states
                        site_oxidations = copy.deepcopy(self.fixed_site_oxidations)
                        site_oxidations.append(oxidation_state)
                        a_struct.add_oxidation_state_by_site(site_oxidations)
                        #calc and append Ewald of structure
                        es = EwaldSummation(a_struct)
                        a = es.total_energy-self.const
                        if write_part_struct:
                            cifWriter = CifWriter(a_struct)
                            cifWriter.write_file("alpha-" + str(j) + "-" + str(i) + ".cif")
                        return s, a

                    #Parallelization
                    #start = time.time()
                    results = Parallel(n_jobs=self.cores)(delayed(calc_s_a_parallel)(i) for i in range(it_site.num))
                    for i in range(len(results)):
                        self_e[-1].append(results[i][0])
                        alpha[-1].append(results[i][1])
                    #print(time.time()-start)

                    #Non-parallel code
                    #self_e.append([])
                    #alpha.append([])
                    #start = time.time()
                    #for i in range(it_site.num):
                    #    new_struct = Structure(lattice=struct.lattice, species=[{el:1.0,},], labels=[label,],
                    #                           coords=[it_site.coords[i],], coords_are_cartesian=True)
                    #    site_oxidations = self.get_oxidation_states(new_struct)
                    #    oxidation_state = site_oxidations[0]
                    #    new_struct.add_oxidation_state_by_site(site_oxidations)
                    #    es = EwaldSummation(new_struct)
                    #    self_e[-1].append(es.total_energy)

                    #    a_struct = copy.deepcopy(self.fixed_struct)
                    #    a_struct.append(species={el:1.0,}, coords=it_site.coords[i], coords_are_cartesian=True)
                    #    site_oxidations = copy.deepcopy(self.fixed_site_oxidations)
                    #    site_oxidations.append(oxidation_state)
                    #    a_struct.add_oxidation_state_by_site(site_oxidations)
                    #    es = EwaldSummation(a_struct)
                    #    alpha[-1].append(es.total_energy-self.const)

                    #print(time.time() - start)
                    #print(np.array(alpha[-1])-np.array(alpha[-2]))
            #cifWriter = CifWriter(new_struct)
            #cifWriter.write_file("test.cif")
            #cifWriter = CifWriter(a_struct)
            #cifWriter.write_file("test.cif")
            return self_e, alpha

        def calc_beta_energy():
            """Function to calculate energies between two iterative sites"""
            # copy structure to make manipulations in copy
            struct = copy.deepcopy(self.struct)
            #get dimension for beta matrix and creat beta matrix full of zeros
            max_j = -1
            num_num_sites = 0
            for it_site in self.iterate_sites.values():
                num_num_sites += len(it_site.site.species.get_el_amt_dict().keys())
                if max_j < it_site.num:
                    max_j = it_site.num
            beta = np.zeros((num_num_sites, max_j, num_num_sites, max_j))
            #beta2 = np.zeros((num_num_sites, max_j, num_num_sites, max_j))
            #iterate over all sites and elements per site
            i = -1
            for lab1, it_site1 in self.iterate_sites.items():
                for el1 in it_site1.site.species.get_el_amt_dict().keys():
                    #get one leable for charges
                    label_num = re.findall(r'\d+', lab1)[-1]
                    label1 = el1 + label_num
                    i += 1

                    #iterate over all positions of site in parallel
                    def calc_b_parallel(j):
                        beta_tmp = np.zeros((num_num_sites, max_j))
                        k = -1
                        #iterate over all sites and elemts per site
                        for lab2, it_site2 in self.iterate_sites.items():
                            for el2 in it_site2.site.species.get_el_amt_dict().keys():
                                #get other label for charges
                                label_num = re.findall(r'\d+', lab2)[-1]
                                label2 = el2 + label_num
                                k += 1
                                #iterate over all site positions
                                for l in range(it_site2.num):
                                    #check to only calculate one half of the beta matrix as interacions are symmetric
                                    if (k == i and l > j) or k > i:
                                        #create new strucuter with just two sites
                                        new_struct = Structure(lattice=struct.lattice,
                                                               species=[{el1:1.0,}, {el2:1.0,}],
                                                               labels=[label1, label2],
                                                               coords=[it_site1.coords[j], it_site2.coords[l]],
                                                               coords_are_cartesian=True)
                                        #add oxdation states
                                        site_oxidations = self.get_oxidation_states(new_struct)
                                        new_struct.add_oxidation_state_by_site(site_oxidations)
                                        #calc ewald for structure
                                        es = EwaldSummation(new_struct)
                                        beta_tmp[k,l] = es.total_energy - self.self_e[i][j] - self.self_e[k][l]
                                        if write_part_struct:
                                            cifWriter = CifWriter(new_struct)
                                            cifWriter.write_file("beta-" + str(i) + "-" + str(j) + "-" +
                                                                 str(k) + "-" + str(l) + ".cif")
                        return beta_tmp

                    #Parallelization
                    #start = time.time()
                    results = Parallel(n_jobs=self.cores)(delayed(calc_b_parallel)(j) for j in range(it_site1.num))
                    #merge parallel results
                    for j, result in enumerate(results):
                        beta[i,j,:,:] = result
                    #print(time.time() - start)

                    #Non-parallel code
                    #start = time.time()
                    #for j in range(it_site1.num):
                    #    k = -1
                    #    for lab2, it_site2 in self.iterate_sites.items():
                    #        for el2 in it_site2.site.species.get_el_amt_dict().keys():
                    #            label_num = re.findall(r'\d+', lab2)[-1]
                    #            label2 = el2 + label_num
                    #            k += 1
                    #            for l in range(it_site2.num):
                    #                if (k == i and l > j) or k > i:
                    #                    new_struct = Structure(lattice=struct.lattice,
                    #                                           species=[{el1:1.0,}, {el2:1.0,}],
                    #                                           labels=[label1, label2],
                    #                                           coords=[it_site1.coords[j], it_site2.coords[l]],
                    #                                           coords_are_cartesian=True)
                    #                    site_oxidations = self.get_oxidation_states(new_struct)
                    #                    new_struct.add_oxidation_state_by_site(site_oxidations)
                    #                    es = EwaldSummation(new_struct)
                    #                    beta2[i,j,k,l] = es.total_energy - self.self_e[i][j] - self.self_e[k][l]
                    #print(time.time() - start)
            #cifWriter = CifWriter(new_struct)
            #cifWriter.write_file("test.cif")
            #print(np.sum(beta-beta2))
            return beta


        start = time.time()
        print("Started Coulomb matrices calculation")
        const = calc_const_energy()
        #print(const)
        self.const = const
        self_e, alpha = calc_self_and_alpha_energy()
        #print(self_e)
        #print(alpha)
        self.self_e = self_e
        self.alpha = alpha
        beta = calc_beta_energy()
        # print(beta)
        self.beta = beta
        print("Finished Coulomb matrices calculations in: " + str(round(time.time()-start)) + "s")
        print()

    def get_oxidation_states(self, struct:Structure):
        """Function to obtain oxidation states per site.
        Required to set average oxidation states for fixed sites and to
        to allow for site-specific (instead of element-specific) charges by the user"""
        site_oxidations = []
        for site in struct.sites:
            oxidation_state = 0
            #Average oxidation state for all species in a site
            for el, fra in site.species.get_el_amt_dict().items():
                # Create label for different elements as shared sites have only label
                # of last element+number by default
                label_num = re.findall(r'\d+', site.label)[-1]
                label = el + label_num
                #Iterate over charges given by user
                for l, c in self.charges.items():
                    #Wildcarde charges
                    if "*" in l:
                        l = l.replace("*", "")
                        if l in label:
                            oxidation_state += c*fra
                    #Exact match charges
                    else:
                        if l == label:
                            oxidation_state += c*fra
            site_oxidations.append(oxidation_state)
        return site_oxidations

    def write_energy(self):
        f = open("const", "w")
        f.write(str(self.const))
        f.close()
        f = open("self_e", "w")
        g = open("alpha", "w")
        i = -1
        for l, it_site in self.iterate_sites.items():
            for el in it_site.site.species.get_el_amt_dict().keys():
                i += 1
                for j in range(it_site.num):
                    f.write(str(self.self_e[i][j]) + "\n")
                    g.write(str(self.alpha[i][j]) + "\n")
        f.close()
        g.close()
        f = open("beta", "w")
        i = -1
        for l, it_site in self.iterate_sites.items():
            for el in it_site.site.species.get_el_amt_dict().keys():
                i += 1
                for j in range(it_site.num):
                    k = -1
                    for lab2, it_site2 in self.iterate_sites.items():
                        for el2 in it_site2.site.species.get_el_amt_dict().keys():
                            k += 1
                            for l in range(it_site2.num):
                                if (k == i and l > j) or k > i:
                                    f.write(str(self.beta[i,j,k,l]) + "\n")
        f.close()

    def read_energy(self):
        print("Reading energies from files, any given oxidation states will be ignored!")
        f = open("const", "r")
        self.const = float(f.readline().strip())
        f.close()
        f = open("self_e", "r")
        g = open("alpha", "r")
        self.self_e = []
        self.alpha = []
        i = -1
        for l, it_site in self.iterate_sites.items():
            for el in it_site.site.species.get_el_amt_dict().keys():
                self.self_e.append([])
                self.alpha.append([])
                i += 1
                for j in range(it_site.num):
                    self.self_e[-1].append(float(f.readline().strip()))
                    self.alpha[-1].append(float(g.readline().strip()))
        f.close()
        g.close()
        f = open("beta", "r")
        max_j = -1
        num_num_sites = 0
        for it_site in self.iterate_sites.values():
            num_num_sites += len(it_site.site.species.get_el_amt_dict().keys())
            if max_j < it_site.num:
                max_j = it_site.num
        self.beta = np.zeros((num_num_sites, max_j, num_num_sites, max_j))
        i = -1
        for l, it_site in self.iterate_sites.items():
            for el in it_site.site.species.get_el_amt_dict().keys():
                i += 1
                for j in range(it_site.num):
                    k = -1
                    for lab2, it_site2 in self.iterate_sites.items():
                        for el2 in it_site2.site.species.get_el_amt_dict().keys():
                            k += 1
                            for l in range(it_site2.num):
                                if (k == i and l > j) or k > i:
                                    self.beta[i,j,k,l] = float(f.readline().strip())
        f.close()


    def __str__(self):
        out = "\n"
        all_combinations = []

        oxidation_states = {}
        for site in self.struct.sites:
            #Average oxidation state for all species in a site
            for el, fra in site.species.get_el_amt_dict().items():
                oxidation_state = 0
                # Create label for different elements as shared sites have only label
                # of last element+number by default
                label_num = re.findall(r'\d+', site.label)[-1]
                label = el + label_num
                #Iterate over charges given by user
                for l, c in self.charges.items():
                    #Wildcarde charges
                    if "*" in l:
                        l = l.replace("*", "")
                        if l in label:
                            oxidation_state += c
                    #Exact match charges
                    else:
                        if l == label:
                            oxidation_state += c
                    oxidation_states[label] = oxidation_state

        for label, oxidation_state in oxidation_states.items():
            line = "Using oxidation state of " + str(oxidation_state) + " for site " + label + "\n"
            out += line
        out += "\n"

        composition = {}
        for i, site in enumerate(self.struct.sites):
            for el, fra in site.species.get_el_amt_dict().items():
                if el not in composition.keys():
                    composition[el] = fra
                else:
                    composition[el] += fra
        line = ""
        total_charge = 0
        for el, fra in composition.items():
            line += el + str(int(fra+0.5)) + " "
            for l, c in self.charges.items():
                if el in l:
                    total_charge += int(fra+0.5)*c
                    break
        out += "The desired composition is: " + line + "\n"
        out += "The total charge of the supecell is: " + str(total_charge) + "\n\n"

        for label, it_site in self.iterate_sites.items():
            line = "Found iterative site " + label + " with:\n"
            out += line
            num_ions = []
            for el, occ in it_site.site.species.as_dict().items():

                line = ("Place " + str(int(occ*it_site.num+0.5)) + " ions of " + el + " in "
                        + str(it_site.num) + " positions.\n")
                out += line
                num_ions.append(int(occ*it_site.num+0.5))
            combinations = np.float128(1)
            for num_ion in num_ions:
                combinations *= np.float128(np.math.factorial(num_ion))
            if sum(num_ions) > it_site.num:
                raise("It seems ocupations result in exactly xx.5 ions and mathematical rounding would cause over-occupation of the site. This case is not defined, please provide proper occupations.")
            combinations *= np.float128(np.math.factorial(it_site.num-sum(num_ions)))
            combinations = np.log10(np.float128(np.math.factorial(it_site.num))/combinations)
            all_combinations.append(combinations)
            line = "Site has " + str(combinations) + " log10(combinations).\n\n"
            out += line
        out += "\n"
        total_combinations = np.float128(1)
        for i in range(len(all_combinations)):
            total_combinations *= np.power(10, all_combinations[i])
        total_combinations = np.log10(total_combinations)
        line = "Total log10(combinations) of the problem are: " + str(total_combinations) + ".\n"
        out += line

        return out

    def print_solution_cif_gurobi(self, name:str, x):
        """Function to write a cif file of a structure for a given gurobi solution x of the problem"""
        # copy structure to make manipulations in copy
        struct = copy.deepcopy(self.struct)
        # delete all iterative sites
        delete_sites = []
        for i, site in enumerate(struct.sites):
            if site.label in self.iterate_sites.keys():
                delete_sites.append(i)
        struct.remove_sites(delete_sites)

        #Iterate over all iterative sites and ions and positions
        i = -1
        var_counter = -1
        for it_site in self.iterate_sites.values():
            for el in it_site.site.species.as_dict().keys():
                i += 1
                for j in range(it_site.num):
                    var_counter += 1
                    #check if site is occupied in solution, var_counter needed as Grurobi output id 1D
                    if np.abs(x[var_counter]-1) < 10E-6:
                        struct.append(species={el: 1.0,}, coords=it_site.coords[j], coords_are_cartesian=True)

        #Write structure as cif
        cifWriter = CifWriter(struct)
        cifWriter.write_file(name)


    def print_solution_cif_fotran(self, name: str, x):
        """Function to write a cif file of a structure for a given Fortran solution x of the problem"""
        # copy structure to make manipulations in copy
        struct = copy.deepcopy(self.struct)
        # delete all iterative sites
        delete_sites = []
        for i, site in enumerate(struct.sites):
            if site.label in self.iterate_sites.keys():
                delete_sites.append(i)
        struct.remove_sites(delete_sites)

        # Iterate over all iterative sites and ions and positions
        i = -1
        for it_site in self.iterate_sites.values():
            for el in it_site.site.species.as_dict().keys():
                i += 1
                for j in range(it_site.num):
                    if x[i,j]:
                        struct.append(species={el: 1.0, }, coords=it_site.coords[j], coords_are_cartesian=True)

        # Write structure as cif
        cifWriter = CifWriter(struct)
        cifWriter.write_file(name)