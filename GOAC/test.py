from CoulombToolbox import greedy
import numpy as np

num_it_sites = 4
num_species = 10
max_site_num = 24
species_nums = np.zeros(num_species)+4
species_occs = np.zeros(num_species)+1
shared_sites = np.array([[1,0,0,0,0,0,0,0,0,0], [2,0,0,0,0,0,0,0,0,0], [3,4,5,6,0,0,0,0,0,0], [7,8,9,10,0,0,0,0,0,0]])
const_e = -151.453
self_e = np.zeros((num_species, max_site_num))
a_e = np.zeros((num_species, max_site_num))
b_e = np.zeros((num_species, max_site_num, num_species, max_site_num))

shared_sites = []
for i in range(num_species):
    tmp = []
    for j in range(num_species):
        tmp.append(False)
    shared_sites.append(tmp)

shared_sites = np.array(shared_sites)

best_energy, best_solution = greedy(species_nums, species_occs, shared_sites, const_e, self_e, a_e, b_e)
print(best_solution)
