# cell parameters
fitness_threshold = 7.5
apoptosis_threshold =  5 #10 # threshold for the cell to go in apoptosis - should be decreased over time?
max_copy_plasmid_rate = 0#0.1
max_lose_plasmid_rate = 0#0.1

# protein dynamics parameters
alpha = 10
Kd = 1
n = 2
delta = 1

# fitness protein parameters
alpha_fitness = 2.5
Kd_fitness = 1
n_fitness = 2
delta_fitness = 0.2

# apoptosis protein parameters
alpha_apoptosis = 1
Kd_apoptosis = 1
n_apoptosis = 4
delta_apoptosis = 0.01

# find best neighbours
prob_find = 0.9 # probability of taking the best neighbour or randomly initializing a cell
