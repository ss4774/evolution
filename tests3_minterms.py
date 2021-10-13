#import cell
#import models
import population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import datetime

observables_global = None
observables_local =  None

N_inputs = 3 # number of inputs: in_1, in_2,...in_N_inputs
iterations = 3 # number of learning iterations

repeats = 10 # how many repeats for a given scenario (function, lattice size)

t_end = 200
dt = 0.1
plot_resolution = 1


func_and = '(in_1 & in_2 & in_3)'
states_and = {"in_1":(0,0,0,0,10,10,10,10),
              "in_2":(0,0,10,10,0,0,10,10),
              "in_3":(0,10,0,10,0,10,0,10),
              "eval":(0,0,0,0,0,0,0,10)}
label_and = "and"

func_or = '(in_1 | in_2 | in_3)'
states_or = {"in_1":(0,0,0,0,10,10,10,10),
             "in_2":(0,0,10,10,0,0,10,10),
             "in_3":(0,10,0,10,0,10,0,10),
             "eval":(0,10,10,10,10,10,10,10)}
label_or = "or"

func_xor = '(~(in_1) & ~(in_2) & in_3) | (~(in_1) & in_2 & ~(in_3)) | (~(in_1) & ~(in_2) & in_3) | (in_1 & in_2 & in_3)'
states_xor = {"in_1":(0,0,0,0,10,10,10,10),
              "in_2":(0,0,10,10,0,0,10,10),
              "in_3":(0,10,0,10,0,10,0,10),
              "eval":(0,10,10,0,10,0,0,10)}
label_xor = "xor"

func_all = (func_or, func_and, func_xor)
states_all = (states_or, states_and, states_xor)
labels = (label_or, label_and, label_xor)

lattice_sizes = (2,5,7,10)

for func, states, label in zip(func_all, states_all, labels):
    for N in lattice_sizes:
        for _ in range(repeats): # 10 repetitions
  
            #N = 1 # size of the lattice is N x N

            p = population.population()
            p.generate_cells_minterms(N, N_inputs)

            _, functions = p.simulate(states, observables_local, observables_global, t_end, dt=dt, iterations = iterations, plot_resolution = plot_resolution, track_states=False)

            now = datetime.datetime.now()
            yr,mth,day,h,m,s = now.year, now.month, now.day, now.hour, now.minute, now.second
            f = open(f'tests\\functions_minterms_{label}{N_inputs}\\{N*N}\\functions2_{yr}_{mth}_{day}_{h}_{m}_{s}.txt', 'w')
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

