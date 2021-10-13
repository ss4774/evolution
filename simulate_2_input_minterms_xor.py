#import cell
#import models
import population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import datetime

observables_global = ["in_1", "in_2", "eval"]
observables_local =  ["fitness","apoptosis","y"]#, "fitness", "apoptosis"]


func = '(~(in_1) & in_2) | (in_1 & ~(in_2))'
states = {"in_1":(0,0,10,10),
          "in_2":(0,10,0,10),
          "eval":(0,10,10,0)}

N = 2 # size of the lattice is N x N
N_inputs = 2 # number of inputs: in_1, in_2,...in_N_inputs
iterations = 3 # number of learning iterations

t_end = 100
dt = 0.5
plot_resolution = 1

functions = {}

p = population.population()
p.generate_cells_minterms(N, N_inputs)

#df, functions = p.simulate(states, observables_local, observables_global, t_end, dt=dt, iterations = iterations, plot_resolution = plot_resolution, track_states=False)
_, functions = p.simulate(states, observables_local, observables_global, t_end, dt=dt, iterations = iterations, plot_resolution = plot_resolution, track_states=False)

#df.to_csv('test.txt', index=False)

#f = open('examples\\functions2\\functions2.txt', 'w')

now = datetime.datetime.now()
yr,mth,day,h,m,s = now.year, now.month, now.day, now.hour, now.minute, now.second
f = open(f'examples\\functions_xor2\\functions2_{yr}_{mth}_{day}_{h}_{m}_{s}.txt', 'w')
if func:
    f.write("-1;")
    f.write(func)
    f.write("\n")
for t in functions:
    NF = sorted(functions[t].items(), key=lambda x: x[1])
    f.write(f't={t};')
    for n,func in NF:
        f.write(f'{n}:{func}; ')
    f.write("\n")
f.close()

print("Final functions: ", NF)
#df = df.set_index('t')
#df.plot()
#plt.show()


    


