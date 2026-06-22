__author__ = "Konstantin Köster"
__copyright__ = "Copyright 2024, GOAC"
__credits__ = ["Konstantin Köster", "Tobias Binninger", "Payam Kaghazch"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = ""
__email__ = "p.kaghazchi@fz-juelich.de"
__status__ = "Development"

import copy

from pymatgen.core import Structure
import warnings
import re
import numpy as np
import math
import copy
from pymatgen.io.cif import CifWriter
import time
import ABCEwald

def base_label(label):
    """Return the crystallographic site label without the supercell index suffix.

    Newer pymatgen (>= 2024.x) appends a unique "_<n>" suffix to each site label
    in Structure.make_supercell (e.g. "Li1" -> "Li1_29"), even for a [1,1,1] cell.
    GOAC groups symmetry-equivalent iterative sites by their shared label, so the
    suffix must be stripped to restore the original grouping behaviour."""
    return re.sub(r'_\d+$', '', label)

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
        self.occ_error = False

        #Read cif-file with pymatgen
        with warnings.catch_warnings(record=True) as recorded_warnings:
            warnings.simplefilter("ignore")
            self.struct = Structure.from_file(self.cif_file)
            self.struct = self.struct.make_supercell(supercell)
        if self.struct is None:
            for w in recorded_warnings:
                warnings.warn_explicit(message=w.message, category=w.category, filename=w.filename,
                                       lineno=w.lineno, source=w.source)
            raise ("Error while reading cif file. Please check the warnings above. If you believe your cif is correct,"
                   "consider to save your cif as POSCAR and re-save this POSCAR as cif by using e.g., VESTA."
                   "This might resolve this issue in case your cif is slightly corrupted.")

        #Identify iteratie sites
        self.iterate_sites = {}
        self.coords = []
        self.num_species_site = []
        self.const_sites = []
        q = []
        max_q = 0
        sigma = []
        out = "\n"
        all_combinations = []
        oxidation_states = {}
        composition = {}
        avg_el_charges = {}
        for i, site in enumerate(self.struct.sites):
            #check if user gave site as fixed
            self.coords.append(site.coords)
            fixed = False
            fras = []
            chgs = []
            sigmas = []
            for el, fra in site.species.get_el_amt_dict().items():
                if el not in composition.keys():
                    composition[el] = fra
                else:
                    composition[el] += fra

                #Create label for different elements as shared sites have only
                # label of last element+ number by default
                fras.append(fra)
                label_num = re.findall(r'\d+', base_label(site.label))[-1]
                label = el + label_num
                for l in self.fixed_sites:
                    #check for wildcard match
                    if "*" in l:
                        l = l.replace("*", "")
                        if l in label:
                            # Wildcard only for numbers in labels
                            # Else N* would match every N and Ni
                            if label.replace(l, "").isdigit():
                                fixed = True
                    #check for exact match
                    else:
                        if l == label:
                            fixed = True
                #iterater over charges given by user
                chg_set = False
                for l, c in self.charges.items():
                    #Wildcarde charges
                    if "*" in l:
                        l = l.replace("*", "")
                        if l in label:
                            # Wildcard only for numbers in labels
                            # Else N* would match every N and Ni
                            if label.replace(l, "").isdigit():
                                chgs.append(c[0])
                                sigmas.append(c[1])
                                chg_set = True
                    #Exact match charges
                    else:
                        if l == label:
                            chgs.append(c[0])
                            sigmas.append(c[1])
                            chg_set = True
                if not chg_set:
                    chgs.append(0)
                    sigmas.append(-1)
                oxidation_states[label] = [chgs[-1], sigmas[-1]]
                if el not in avg_el_charges.keys():
                    avg_el_charges[el] = [chgs[-1],]
                else:
                    avg_el_charges[el].append(chgs[-1])

            #Identify sites that are not fixed and have combinatios to iterate over and store them
            if not fixed and (len(site.species.as_dict()) > 1 or
                              abs(list(site.species.as_dict().values())[0]-1.0) > 10e-5):
                #Append vector for charges
                q.append(chgs)
                sigma.append(sigmas)
                if len(chgs) > max_q:
                    max_q = len(chgs)
                self.num_species_site.append(len(chgs))
                self.const_sites.append(False)
                #Raise error if occupations of one site are more than 1
                if sum(list(site.species.as_dict().values())) > 1.0:
                    raise "Site " + site.label + " has an occupation greater 1."
                it_label = base_label(site.label)
                if it_label not in self.iterate_sites.keys():
                    self.iterate_sites[it_label] = Iterate_Site(site, 1)
                else:
                    self.iterate_sites[it_label].num += 1
                    self.iterate_sites[it_label].coords.append(site.coords)
            else:
                #Append avg. charge to q
                q.append([np.sum(np.array(chgs)*np.array(fras)),])
                sigma.append([np.sum(np.array(sigmas)*np.array(fras)),])
                self.num_species_site.append(1)
                if max_q < 1:
                    max_q = 1
                self.const_sites.append(True)

        # Prepare q and sigma array for use with Fortran function
        self.q = np.zeros((len(self.coords), max_q))
        for i, q_site in enumerate(q):
            for j, q_species in enumerate(q_site):
                self.q[i,j] = q_species
        self.sigma = np.zeros((len(self.coords), max_q))
        for i, s_site in enumerate(sigma):
            for j, s_species in enumerate(s_site):
                self.sigma[i,j] = s_species

        self.coords = np.array(self.coords)
        self.num_species_site = np.array(self.num_species_site)
        self.const_sites = np.array(self.const_sites)

        self.species_site_map = np.zeros((len(self.coords), max_q, 2))
        i = -1
        for it_site in self.iterate_sites.values():
            for ii, el in enumerate(it_site.site.species.get_el_amt_dict().keys()):
                i += 1
                for j in range(it_site.num):
                    coord_id = -1
                    for c in range(len(self.coords)):
                        if np.linalg.norm(self.coords[c]-it_site.coords[j]) < 10E-07:
                            coord_id = c
                            break
                    if coord_id >= 0:
                        self.species_site_map[coord_id, ii, :] = [i+1,j+1]
        for label, oxidation_state in oxidation_states.items():
            line = "Using oxidation state of " + str(oxidation_state[0]) + " for site " + label + "\n"
            out += line
        out += "\n"

        if np.max(self.sigma) > 0:
            for label, oxidation_state in oxidation_states.items():
                if oxidation_state[1] < 0:
                    line = "Using point-charge for site " + label + "\n"
                else:
                    line = "Using smeared charge with width of " + str(oxidation_state[1]) + " Ang for site " + label + "\n"
                out += line
            out += "\n"

        line = ""
        total_charge = 0
        total_amounts = {}
        for el, fra in composition.items():
            line += el + str(int(fra + 0.5)) + " "
            total_amounts[el] = int(fra + 0.5)
            total_charge += np.sum(avg_el_charges[el])/len(avg_el_charges[el])*int(fra + 0.5)
        out += "The desired composition is: " + line + "\n"
        out += "The total charge of the supercell is: " + str(total_charge) + "\n\n"

        inf_site_combinations = False
        placed_amounts = {}
        for label, it_site in self.iterate_sites.items():
            line = "Found iterative site " + label + " with:\n"
            out += line
            num_ions = []
            for el, occ in it_site.site.species.get_el_amt_dict().items(): #as_dict().items():
                line = ("Place " + str(int(occ * it_site.num + 0.5)) + " ions of " + el + " in "
                        + str(it_site.num) + " positions.\n")
                out += line
                if not el in placed_amounts.keys():
                    placed_amounts[el] = int(occ * it_site.num + 0.5)
                else:
                    placed_amounts[el] += int(occ * it_site.num + 0.5)
                num_ions.append(int(occ * it_site.num + 0.5))

            if sum(num_ions) > it_site.num:
                self.occ_error = True
                print("---------------------------------------------------------------------------")
                print("WARNING: Occupancies of " + label + " could not be matched in your supercell.")
                print("---------------------------------------------------------------------------")
                print()
                break
            try:
                combinations = np.float128(1)
                for num_ion in num_ions:
                    combinations *= np.float128(math.factorial(num_ion))
                combinations *= np.float128(math.factorial(it_site.num - sum(num_ions)))
                combinations = np.log10(np.float128(math.factorial(it_site.num)) / combinations)
                all_combinations.append(combinations)
                line = "Site has " + str(combinations) + " log10(combinations).\n\n"
            except ValueError as e:
                if "Exceeds the limit" in str(e):
                    all_combinations.append(np.float128(1))
                    line = "Site has +INF log10(combinations).\n\n"
                    inf_site_combinations = True
                else:
                    raise e
            out += line
        out += "\n"
        #for el in placed_amounts.keys():
        #    if total_amounts[el] < placed_amounts[el]:
        #        self.occ_error =True
        #        print("---------------------------------------------------------------------------")
        #        print("WARNING: Occupancies of " + el + " could not be matched in your supercell.")
        #        print("---------------------------------------------------------------------------")
        #        print()
        #        break
        if not self.occ_error:
            if not inf_site_combinations:
                try:
                    total_combinations = np.float128(1)
                    for i in range(len(all_combinations)):
                        total_combinations *= np.power(10, all_combinations[i])
                    total_combinations = np.log10(total_combinations)
                    line = "Total log10(configurations) of the problem are: " + str(total_combinations) + ".\n"
                except ValueError as e:
                    if "Exceeds the limit" in str(e):
                        total_combinations = np.float128(1)
                        line = "Total log10(configurations) of the problem are: +INF.\n"
                    else:
                        raise e
            else:
                total_combinations = np.float128(1)
                line = "Total log10(configurations) of the problem are: +INF.\n"
            out += line
            self.total_combinations = total_combinations
            self.total_charge = total_charge
            self.out = out

    def calc_coulomb_matrices(self):
        """Function to calculate the different Coulomb contributions Const, Self, Alpha, Beta"""
        start = time.time()
        print("Started Coulomb matrices calculation")

        self.const, self.alpha, self.beta, self.a = ABCEwald.ewald.getabc(r=self.coords, q=self.q, sigma=self.sigma,
                                                                  lat=self.struct.lattice.as_dict()["matrix"],
                                                                  alpha=-1, cutoff_real=-1, cutoff_fourier=-1, acc=1E-10,
                                                                  species_site_map=self.species_site_map,
                                                                  num_species_site=self.num_species_site,
                                                                  num_species=np.max(self.species_site_map[:,:,0]),
                                                                  max_site_num=np.max(self.species_site_map[:,:,1]),
                                                                  const_sites=self.const_sites)

        print("Finished Coulomb matrices calculations in: " + str(round(time.time() - start, 1)) + "s")
        print()

        if np.abs(self.total_charge) > 1e-05:
            print("-------- WARNING - Charged Supercell --------")
            print("Calculating correction for charged supercell")
            ang2bohr = 0.52917721067121
            hartree2eV = 27.21138624598130
            self.const -= (ang2bohr*hartree2eV/2*np.pi/self.struct.volume/
                           self.a**2.0*self.total_charge**2.0)
            print()

    def write_energy(self):
        f = open("const", "w")
        f.write(str(self.const) + "\n")
        f.close()
        f = open("alpha", "w")
        i = -1
        for l, it_site in self.iterate_sites.items():
            for el in it_site.site.species.get_el_amt_dict().keys():
                i += 1
                for j in range(it_site.num):
                    f.write(str(self.alpha[i][j]) + "\n")
        f.close()
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
        f = open("alpha", "r")
        self.alpha = []
        i = -1
        for l, it_site in self.iterate_sites.items():
            for el in it_site.site.species.get_el_amt_dict().keys():
                self.alpha.append([])
                i += 1
                for j in range(it_site.num):
                    self.alpha[-1].append(float(f.readline().strip()))
        f.close()
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

    def write_partial_cifs(self):
        # copy structure to make manipulations in copy
        struct = copy.deepcopy(self.struct)
        # delete all iterative sites
        delete_sites = []
        for i, site in enumerate(struct.sites):
            if base_label(site.label) in self.iterate_sites.keys():
                delete_sites.append(i)
        struct.remove_sites(delete_sites)

        #write constant ions in cif file
        if len(struct.sites) > 1:
            with warnings.catch_warnings(record=True) as recorded_warnings:
                warnings.simplefilter("ignore")
                cifWriter = CifWriter(struct)
                cifWriter.write_file("const.cif")
        fixed_struct = copy.deepcopy(struct)

        i = -1
        for l, it_site in self.iterate_sites.items():
            for el in it_site.site.species.get_el_amt_dict().keys():
                i += 1
                # create label
                label_num = re.findall(r'\d+', l)[-1]
                label = el + label_num
                # Iterate over all site positions in parallel
                for j in range(it_site.num):
                    # create new structure with only one site
                    new_struct = Structure(lattice=struct.lattice, species=[{el: 1.0, }, ], labels=[label, ],
                                           coords=[it_site.coords[i], ], coords_are_cartesian=True)
                    with warnings.catch_warnings(record=True) as recorded_warnings:
                        warnings.simplefilter("ignore")
                        cifWriter = CifWriter(new_struct)
                        cifWriter.write_file("self-" + str(i) + "-" + str(j) + ".cif")

                    # create structure with site and all fixed sites
                    a_struct = copy.deepcopy(fixed_struct)
                    a_struct.append(species={el: 1.0, }, coords=it_site.coords[i], coords_are_cartesian=True)
                    with warnings.catch_warnings(record=True) as recorded_warnings:
                        warnings.simplefilter("ignore")
                        cifWriter = CifWriter(a_struct)
                        cifWriter.write_file("alpha-" + str(i) + "-" + str(j) + ".cif")

                    k=-1
                    # iterate over all sites and elemts per site
                    for lab2, it_site2 in self.iterate_sites.items():
                        for el2 in it_site2.site.species.get_el_amt_dict().keys():
                            # get other label for charges
                            label_num = re.findall(r'\d+', lab2)[-1]
                            label2 = el2 + label_num
                            k += 1
                            for l in range(it_site2.num):
                                # check to only calculate one half of the beta matrix as interacions are symmetric
                                if (k == i and l > j) or k > i:
                                    new_struct = Structure(lattice=struct.lattice,
                                                           species=[{el: 1.0, }, {el2: 1.0, }],
                                                           labels=[label, label2],
                                                           coords=[it_site.coords[j], it_site2.coords[l]],
                                                           coords_are_cartesian=True)
                                    with warnings.catch_warnings(record=True) as recorded_warnings:
                                        warnings.simplefilter("ignore")
                                        cifWriter = CifWriter(new_struct)
                                        cifWriter.write_file("beta-" + str(i) + "-" + str(j) + "-" +
                                                             str(k) + "-" + str(l) + ".cif")

    def __str__(self):
        return self.out

    def print_solution_cif_gurobi(self, name:str, x):
        """Function to write a cif file of a structure for a given gurobi solution x of the problem"""
        # copy structure to make manipulations in copy
        struct = copy.deepcopy(self.struct)
        # delete all iterative sites
        delete_sites = []
        for i, site in enumerate(struct.sites):
            if base_label(site.label) in self.iterate_sites.keys():
                delete_sites.append(i)
        struct.remove_sites(delete_sites)

        #Iterate over all iterative sites and ions and positions
        i = -1
        var_counter = -1
        added_species = False
        for it_site in self.iterate_sites.values():
            for el in it_site.site.species.as_dict().keys():
                i += 1
                for j in range(it_site.num):
                    var_counter += 1
                    #check if site is occupied in solution, var_counter needed as Grurobi output id 1D
                    if np.abs(x[var_counter]-1) < 10E-6:
                        added_species = True
                        struct.append(species={el: 1.0,}, coords=it_site.coords[j], coords_are_cartesian=True)

        #Write structure as cif
        if added_species:
            with warnings.catch_warnings(record=True) as recorded_warnings:
                warnings.simplefilter("ignore")
                cifWriter = CifWriter(struct)
                cifWriter.write_file(name)
        return added_species


    def print_solution_cif_fotran(self, name: str, x):
        """Function to write a cif file of a structure for a given Fortran solution x of the problem"""
        # copy structure to make manipulations in copy
        struct = copy.deepcopy(self.struct)
        # delete all iterative sites
        delete_sites = []
        for i, site in enumerate(struct.sites):
            if base_label(site.label) in self.iterate_sites.keys():
                delete_sites.append(i)
        struct.remove_sites(delete_sites)

        added_species = False
        # Iterate over all iterative sites and ions and positions
        i = -1
        for it_site in self.iterate_sites.values():
            for el in it_site.site.species.as_dict().keys():
                i += 1
                for j in range(it_site.num):
                    if x[i,j]:
                        struct.append(species={el: 1.0, }, coords=it_site.coords[j], coords_are_cartesian=True)
                        added_species = True
        # Write structure as cif
        if added_species:
            with warnings.catch_warnings(record=True) as recorded_warnings:
                warnings.simplefilter("ignore")
                cifWriter = CifWriter(struct)
                cifWriter.write_file(name)
        return added_species
