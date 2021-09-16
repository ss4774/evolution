#import cell
#import models
import population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

observables_global = ["in_1", "in_2", "eval"]
observables_local =  ["fitness","apoptosis","y"]#, "fitness", "apoptosis"]

states = {"in_1":(0,0,10,10),
          "in_2":(0,10,0,10),
          "eval":(0,10,10,0)}

N = 2 # size of the lattice is N x N
N_inputs = 2 # number of inputs: in_1, in_2,...in_N_inputs
N_layers = 1 # number of internal layers
N_terms = 3 # number of terms per layer
iterations = 3 # number of learning iterations

t_end = 100
dt = 0.1
plot_resolution = 1

functions = {}

"""
plasmids = [(["in_1"], "x11", "YES"),
            (["in_2"], "x11", "NOT"),
            (["in_1"], "x12", "NOT"),
            (["in_2"], "x12", "YES"),
            (["x11"], "y", "NOT"),
            (["x12"], "y", "NOT")]
"""
"""
plasmids = [(["in_1"], "not_in_1", "NOT"),
            (["in_2"], "not_in_2", "NOT"),
            (["in_1", "not_in_2"], "y", "AND"),
            (["not_in_1", "in_2"], "y", "AND")]
"""

plasmids = [(["in_1", "in_2"], "y", "AND01"),
            (["in_1", "in_2"], "y", "AND10")]


p = population.population()
p.generate_cells_function(N, plasmids)
p.generate_minterms=True

df, functions = p.simulate(states, observables_local, observables_global, t_end, dt=dt, iterations = iterations, plot_resolution = plot_resolution)

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


    


