#import cell
#import models
import generate_population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from collections import defaultdict

from parameter_values import *

observables_global = ["in_1", "in_2", "eval"]
observables_local =  []#["fitness","apoptosis","y"]#, "fitness", "apoptosis"]

N = 5 # size of the lattice is N x N
N_inputs = 2 # number of inputs: in_1, in_2,...in_N_inputs
N_layers = 0 # number of internal layers
N_terms = 2 # number of terms per layer
iterations = 3 # number of learning iterations

layered = False # if plasmid selection is performed by layers
ordered_crossover = True # if crossover takes the first half of plasmids from the first cell and second half from the second (ordered)

t_end = 100
dt = 0.1
plot_resolution = 1

functions = {}


pop = generate_population.generate_cells(N=N, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms, layered=layered)

"""
for i in range(N):
    for j in range(N):
        print(pop[i,j].mod_degradation)
"""


states = {"in_1":(0,0,1,1),
          "in_2":(0,1,0,1),
          "out":(0,1,1,0)}


df = pd.DataFrame(dtype=float)

for it in range(iterations):
    for t in np.arange(0, t_end+dt, dt):
        T = t + t_end*it

        t_end_div_4 = t_end/4

        if (0 <= t < t_end_div_4):
            global_vars = {"in_1": states["in_1"][0]*10,
                        "in_2": states["in_2"][0]*10,
                        "eval": states["out"][0]*10,
                        "learn": 1}
        elif (t_end_div_4 <= t < 2*t_end_div_4):
            global_vars = {"in_1": states["in_1"][1]*10,
                        "in_2": states["in_2"][1]*10,
                        "eval": states["out"][1]*10,
                        "learn": 1}
        elif (2*t_end_div_4 <= t < 3*t_end_div_4):
            global_vars = {"in_1": states["in_1"][2]*10,
                        "in_2": states["in_2"][2]*10,
                        "eval": states["out"][2]*10,
                        "learn": 1}
        else:
            global_vars = {"in_1": states["in_1"][3]*10,
                        "in_2": states["in_2"][3]*10,
                        "eval": states["out"][3]*10,
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
                        generate_population.apoptosis(pop, i, j, N_inputs, N_layers, N_terms, ordered_crossover=ordered_crossover)                        
                
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


    


