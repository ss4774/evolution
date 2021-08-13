#import cell
#import models
import generate_population

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

observables_global = ["eval", "in_1"]#, "in_2", "eval"]
observables_local = ["y", "fitness", "apoptosis"]



N = 1
N_inputs = 1
N_layers = 1
N_terms = 1

t_end = 50
dt = 0.1
plot_resolution = 1



pop = generate_population.generate_cells(N=N, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms)
"""
for i in range(N):
    for j in range(N):
        print(pop[i,j].mod_degradation)
"""




df = pd.DataFrame(dtype=float)

for t in np.arange(0, t_end+dt, dt):
    if t < 25:
        global_vars = {"in_1": 0,
                       "in_2": 0,
                       "eval": 0,
                       "learn": 1}
    elif t < 50:
        global_vars = {"in_1": 10,
                       "in_2": 0,
                       "eval": 10,
                       "learn": 1}
    elif t < 75:
        global_vars = {"in_1": 0,
                       "in_2": 10,
                       "eval": 10,
                       "learn": 1}
    else:
        global_vars = {"in_1": 10,
                       "in_2": 10,
                       "eval": 10,
                       "learn": 0}

    track = (t % plot_resolution) == 0

    if track:
        d = {'t':t}
        for obs in observables_global:
            d[obs] = global_vars[obs]
    
    
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
                                        
                    pop_ = generate_population.generate_cells(N=1, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms)
                    pop[i,j] = pop_[0,0]

    if track:
        df = df.append(d, ignore_index=True, sort=False)


df.to_csv('test.txt', index=False)

df = df.set_index('t')
df.plot()
plt.show()


    


