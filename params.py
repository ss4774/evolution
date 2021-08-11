# cell parameters
fitness_threshold = 7.5
apoptosis_threshold =  5 #10 # threshold for the cell to go in apoptosis - should be decreased over time?
max_copy_plasmid_rate = 0.1
max_lose_plasmid_rate = 0.1

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

# other parameters
neighbourhood = "moore"
prob_elitism = 0.6 # find best neighbours: probability of taking the best neighbour and copying its state after apoptosis
prob_cross = 0.3 # probability of crossover if find_best fails. If elitism and crossover fails, cell will be randomly initialised

