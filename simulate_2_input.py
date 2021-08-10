#import cell
#import models
import generate_population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from collections import defaultdict

from params import *

observables_global = ["in_1", "in_2", "eval"]
observables_local =  []#["fitness","apoptosis","y"]#, "fitness", "apoptosis"]

N = 5 # size of the lattice is N x N
N_inputs = 2 # number of inputs: in_1, in_2,...in_N_inputs
N_layers = 1 # number of internal layers
N_terms = 3 # number of terms per layer
iterations = 3 # number of learning iterations

t_end = 100
dt = 0.1
plot_resolution = 1

functions = {}


pop = generate_population.generate_cells(N=N, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms)
"""
for i in range(N):
    for j in range(N):
        print(pop[i,j].mod_degradation)
"""




df = pd.DataFrame(dtype=float)

for it in range(iterations):
    for t in np.arange(0, t_end+dt, dt):
        T = t + t_end*it

        if (0 <= t < 25):
            global_vars = {"in_1": 0,
                        "in_2": 0,
                        "eval": 0,
                        "learn": 1}
        elif (25 <= t < 50):
            global_vars = {"in_1": 0,
                        "in_2": 10,
                        "eval": 10,
                        "learn": 1}
        elif (50 <= t < 75):
            global_vars = {"in_1": 10,
                        "in_2": 0,
                        "eval": 10,
                        "learn": 1}
        else:
            global_vars = {"in_1": 10,
                        "in_2": 10,
                        "eval": 0,
                        "learn": 1}

        track = (t % plot_resolution) == 0

        if track:
            d = {'t':T}
            for obs in observables_global:
                d[obs] = global_vars[obs]
        
            functions[T] = defaultdict(int)

        
        for i in range(N):
            for j in range(N):
                if track:    
                    prefix= f'cell_{i},{j}_'
                    for obs in observables_local:
                        if obs in pop[i][j].state:
                            if obs == 'y':
                                ff = pop[i][j].decode_function()
                                if ff:
                                    d[prefix+ff] = pop[i][j].state[obs]
                            else:
                                d[prefix+obs] = pop[i][j].state[obs]                                                
                action = pop[i][j].update(dt, global_vars)
                if action:
                    if len(action) == 2:
                        plasmid_pair = action
                        generate_population.give_plasmid(pop, i, j, plasmid_pair)
                    elif action == "apoptosis":
                        # print("kill")
                        # find the best cell in the neighbourhood. If there isn't any, generate a new cell randomly                        
                        best_neigbhour = False
                        if np.random.random() < prob_find_best:
                            best_neigbhour = generate_population.find_best_neighbour(pop, i,j)
                            if best_neigbhour:
                                pop[i,j] = best_neigbhour.copy()
                                #print("apoptosis best")
                        if best_neigbhour == False:
                            if np.random.random() < prob_cross:
                                cross =  generate_population.cross_random_neighbours(pop, i,j)
                                if cross:
                                    pop[i,j] = cross

                        if (best_neigbhour == False) and (cross == False):
                            pop_ = generate_population.generate_cells(N=1, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms)
                            pop[i,j] = pop_[0,0]
                            #print("apoptosis rand")
                
                if track:
                    functions[T][pop[i,j].decode_function()] += 1

        if track:     
            df = df.append(d, ignore_index=True, sort=False)
            


df.to_csv('test.txt', index=False)

f = open('functions.txt', 'w')
for t in functions:
    NF = sorted(functions[t].items(), key=lambda x: x[1])
    f.write(f't={t}:')
    for n,func in NF:
        f.write(f'{n}:{func}; ')
    f.write("\n")
f.close()


df = df.set_index('t')
df.plot()
plt.show()


    


