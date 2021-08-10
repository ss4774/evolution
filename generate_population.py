import cell
import models

from itertools import product
import numpy as np

from params import *

def give_plasmid(pop, i, j, plasmid_pair):
    plasmid_id, plasmid = plasmid_pair
    
    N = pop.shape[0]
    
    H = [(i-1) % N, (i+1) % N]
    V = [(j-1) % N, (j+1) % N]
    #candidates = []
    #for candidate in product(H,V):
    #    if pop[candidate].mode == 1:
    #        candidates.append(candidate)
    candidates = list(product(H,V))
    
    candidate = candidates[np.random.choice(range(len(candidates)), size=1, replace=False)[0]]

    #print(pop[candidate].decode_function())    
    #print(pop[candidate].plasmids)       
    if  pop[candidate].mode == 1:
        if 'params' in plasmid:
            params = plasmid['params'].copy
        else:
            params = None
        if 'mod_degradation' in plasmid:
            mod_degradation = plasmid['mod_degradation']
        else:
            mod_degradation=None

        pop[candidate].add_plasmid(plasmid_id, plasmid['inputs'].copy(), plasmid['output'], plasmid['f'], params, mod_degradation)        
    #print(pop[candidate].decode_function())    
    #print(pop[candidate].plasmids)    

def find_best_neighbour(pop, i,j, mode='apoptosis'):
    N = pop.shape[0]    
    if N == 1 or np.random.random() > prob_find:
        return False
        
    H = [(i-1) % N, (i+1) % N]
    V = [(j-1) % N, (j+1) % N]
    
    candidates = list(product(H,V))

    best_coords = [candidates[0]]
    best_val = pop[candidates[0]].state[mode]
    
    for candidate in candidates[1:]:
        val = pop[candidate].state[mode]
        if val == best_val:
            best_coords.append(candidate)        
        elif mode == 'apoptosis': # lower is better
            if val < best_val:
                best_coords = [candidate]
        elif mode == 'fitness': # higher is better
            if val > best_val:
                best_coords = [candidate]
    
    if len(best_coords) > 1:
        # if there are multiple equally good candidates just pick one randomly
        best_i, best_j = best_coords[np.random.choice(range(len(best_coords)), size=1, replace=False)[0]]
    else:
        best_i, best_j = best_coords[0]

    #print(pop[best_i, best_j].decode_function())

    return pop[best_i, best_j]


def generate_plasmid(inputs, output, f, params = None, mod_degradation=None):
    d = {}
    d['inputs'] = inputs
    d['output'] = output
    d['f'] = f
    if params:
        d['params'] = params
    if mod_degradation:
        d['mod_degradation'] = mod_degradation

    return d
        
def generate_plasmids(ins, outs, prefix):
    plasmids = {}
    plasmids_ids = []


    for i, (input,output) in enumerate(product(ins, outs)):
        label1 = f'{prefix}{i}_YES'
        plasmid1 = generate_plasmid([input], output, models.YES)
        plasmids[label1] = plasmid1
        plasmids_ids.append(label1)

        label2 = f'{prefix}{i}_NOT'
        plasmid2 = generate_plasmid([input], output, models.NOT)
        plasmids[label2] = plasmid2
        plasmids_ids.append(label2)
        
    return plasmids, plasmids_ids

def generate_program_plasmids(N_inputs, N_layers = 1, N_terms = 3):
    if type(N_terms) == int:
        N_terms = [N_terms]
        N_terms *= N_layers

    layers = {}
    plasmids = {}

    for i, n_terms in enumerate(N_terms):
        if i == 0:
            ins = [f'in_{j}' for j in range(1, N_inputs+1)]
        else:
            ins = [f'x_{i}_{j}' for j in range(1, n_terms+1)]

        outs = [f'x_{i+1}_{j}' for j in range(1, n_terms+1)]

        # input and middle layer operons
        Ps, P_ids = generate_plasmids(ins, outs, prefix = f'P{i}_')
        layers[f'layer{i}'] = P_ids
        plasmids.update(Ps)
        
    # output layer operons
    Ps, P_ids = generate_plasmids(outs, ["y"], prefix = f'P{i+1}_')
    layers['layer_out'] = P_ids
    plasmids.update(Ps)

    return plasmids, layers
   

def generate_cells(N, N_inputs, N_layers = 1, N_terms = 3, layered=False):
    plasmids, layers = generate_program_plasmids(N_inputs, N_layers, N_terms)  
    population = np.zeros((N,N), dtype=object)

    for i in range(N):      
        for j in range(N):

            C = cell.cell()
            
            """
            basic functions
            """
            # fitness
            C.add_basic_function("fitness_operon", ["eval", "y"], "fitness", models.EQU, params = (alpha_fitness, Kd_fitness, n_fitness), mod_degradation=delta_fitness)
            
            # apoptosis
            C.add_basic_function("apoptosis_operon", ["fitness"], "apoptosis", models.NOT, params=(alpha_apoptosis, Kd_apoptosis, n_apoptosis), mod_degradation = delta_apoptosis)

            """
            plasmids
            """
            # randomly choose plasmids from the available ones          
            cell_plasmids = []        
            
            if layered:
                for layer in layers:
                    #if layer == 'layer_out':
                    #    n = 1 # take only one operon from the output layer
                    #else:
                    n = np.random.randint(1, len(layers[layer])) # take a random number of plasmids from a layer
                    
                    cell_plasmids.extend(np.random.choice(layers[layer], size=n, replace=False))
            else:
               
                n = np.random.randint(1, len(plasmids)) # take a random number of plasmids from a layer
                cell_plasmids.extend(np.random.choice(list(plasmids.keys()), size=n, replace=False))
                

                
            for plasmid_id in cell_plasmids:
                plasmid = plasmids[plasmid_id]

                #if 'mod_degradation' in plasmid:
                #    C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'], mod_degradation = plasmid['mod_degradation'])
                #else:    
                C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])

            population[i,j] = C
            
    return population

        

if __name__ == "__main__":
    #print(generate_program_module(3, 2, 3, file_name="test.txt"))
    #generate_evaluation_module()
    generate_cells(2, 3, 2, 3)


    



