#import cell
#import models
import population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import datetime

observables_global = None
observables_local =  None

N_inputs = 2 # number of inputs: in_1, in_2,...in_N_inputs
N_layers = 1 # number of internal layers
N_terms = 2 # number of terms per layer
iterations = 3 # number of learning iterations

repeats = 10 # how many repeats for a given scenario (function, lattice size)

t_end = 100
dt = 0.1
plot_resolution = 1


func_xor = '(~(in_1) & in_2) | (in_1 & ~(in_2))'
states_xor = {"in_1":(0,0,10,10),
          "in_2":(0,10,0,10),
          "eval":(0,10,10,0)}
label_xor = "xor"

func_and = '(in_1 & in_2)'
states_and = {"in_1":(0,0,10,10),
          "in_2":(0,10,0,10),
          "eval":(0,0,0,10)}
label_and = "and"

func_or = '(in_1 | in_2)'
states_or = {"in_1":(0,0,10,10),
          "in_2":(0,10,0,10),
          "eval":(0,10,10,10)}
label_or = "or"


func_all = (func_xor, func_or, func_and)
states_all = (states_xor, states_or, states_and)
labels = (label_xor, label_or, label_and)

#func_all = (func_or,)
#states_all = (states_or,)
#labels = (label_or,)

lattice_sizes = (2,5,7,10)
#lattice_sizes = (5,)

for func, states, label in zip(func_all, states_all, labels):
    for N in lattice_sizes:
        for _ in range(repeats): # 10 repetitions
  


            #N = 1 # size of the lattice is N x N

            p = population.population()
            p.generate_cells(N, N_inputs, N_layers, N_terms)

            _, functions = p.simulate(states, observables_local, observables_global, t_end, dt=dt, iterations = iterations, plot_resolution = plot_resolution, track_states=False)

            now = datetime.datetime.now()
            yr,mth,day,h,m,s = now.year, now.month, now.day, now.hour, now.minute, now.second
            f = open(f'tests\\functions_simple_{label}{N_inputs}\\{N*N}\\functions2_{yr}_{mth}_{day}_{h}_{m}_{s}.txt', 'w')
            if func:
                f.write("-1;")
                f.write(func)
                f.write("\n")
            for t in functions:
                NF = sorted(functions[t].items(), key=lambda x: x[1])
                f.write(f't={t};')
                for ff,nn in NF:
                    f.write(f'{ff}:{nn}; ')
                f.write("\n")
            f.close()

            print(f"Final functions for {label}: ", NF)

