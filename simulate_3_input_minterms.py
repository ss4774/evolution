#import cell
#import models
import population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

observables_global = ["in_1", "in_2", "eval"]
observables_local =  ["fitness","apoptosis","y"]#, "fitness", "apoptosis"]

states = {"in_1":(0,0,0,0,10,10,10,10),
          "in_2":(0,0,10,10,0,0,10,10),
          "in_3":(0,10,0,10,0,10,0,10),
          "eval":(0,10,10,0,10,0,0,10)}

N = 10 # size of the lattice is N x N
N_inputs = 3 # number of inputs: in_1, in_2,...in_N_inputs
iterations = 3 # number of learning iterations

t_end = 200
dt = 0.1
plot_resolution = 1

functions = {}

p = population.population()
p.generate_cells_minterms(N, N_inputs)

df, functions = p.simulate(states, observables_local, observables_global, t_end, dt=dt, iterations = iterations, plot_resolution = plot_resolution)

df.to_csv('test.txt', index=False)

f = open('functions.txt', 'w')
for t in functions:
    NF = sorted(functions[t].items(), key=lambda x: x[1])
    f.write(f't={t};')
    for n,func in NF:
        f.write(f'{n}:{func}; ')
    f.write("\n")
f.close()


df = df.set_index('t')
df.plot()
plt.show()


    


