import cell
import models

from itertools import product
import numpy as np

from params import *

def get_neighbours(i,j,N):
    if neighbourhood == "moore":
        H = [(i-1) % N, (i+1) % N, i]
        V = [(j-1) % N, (j+1) % N, j]
        neighbours = set(product(H,V)) - {(i,j)}
        neighbours = list(neighbours)
    else:
        H = [(i-1) % N, (i+1) % N]
        V = [(j-1) % N, (j+1) % N]
        neighbours = list(product(H,V))
    
    return neighbours

# add a given plasmid (one) to a cell that is a neighbour of i,j in the population pop
def give_plasmid(pop, i, j, plasmid_pair):
    plasmid_id, plasmid = plasmid_pair
    
    N = pop.shape[0]
    if N <= 1:
        return
    #H = [(i-1) % N, (i+1) % N]
    #V = [(j-1) % N, (j+1) % N]
    #candidates = list(product(H,V))
    candidates = get_neighbours(i,j,N)
    
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

# add given plasmids to C (cell)
def give_plasmids_directed(C, plasmid_pairs): # C...cell to which the plasmids are copied; plasmid_pairs = [(plasmid_id1, plasmid1),(plasmid_id2, plasmid2),...]
    for plasmid_id, plasmid in plasmid_pairs:
        if 'params' in plasmid:
            params = plasmid['params'].copy
        else:
            params = None
        if 'mod_degradation' in plasmid:
            mod_degradation = plasmid['mod_degradation']
        else:
            mod_degradation=None

        C.add_plasmid(plasmid_id, plasmid['inputs'].copy(), plasmid['output'], plasmid['f'], params, mod_degradation)  


def apoptosis(pop, i, j, N_inputs, N_layers, N_terms):
    # find the best cell in the neighbourhood. If there isn't any, generate a new cell randomly                        
    elitism = False
    crossover = False
    r = np.random.random()

    if r < prob_elitism:
        elitism = find_best_neighbour(pop, i,j)
        if elitism:
            pop[i,j] = elitism.copy()
            
    elif r < prob_elitism + prob_cross:    
        crossover =  cross_random_neighbours(pop, i,j)
        if crossover:
            pop[i,j] = crossover

    if (elitism == False) and (crossover == False): # if elitism and crossover failes the cell is randomly initialised
        pop_ = generate_cells(N=1, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms)
        pop[i,j] = pop_[0,0]
        

    """
    best_neigbhour = find_best_neighbour(pop, i,j)
    if best_neigbhour:
        pop[i,j] = best_neigbhour.copy()
        #print("apoptosis best")
    else:
        pop_ = generate_cells(N=1, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms)
        pop[i,j] = pop_[0,0]
        #print("apoptosis rand")
    """
def cross_random_neighbours(pop, i,j):
    N = pop.shape[0]    
    if N == 1:
        return False

    candidates = get_neighbours(i,j,N)

    C1, C2 = np.random.choice(range(len(candidates)), size=2, replace=False)
    C1, C2 = pop[candidates[C1]], pop[candidates[C2]]
    
    #print(C1.plasmids.keys())
    #print(C2.plasmids.keys())


    C = C1.copy()   # copy the first cell

    # lose half of the existing plasmids of the first cell
    for i in range(len(C.plasmids)//2):
        C.lose_plasmid()

    # copy half of the plasmids of the second cell
    plasmids2 = np.array(list(C2.plasmids.items()))
    if len(plasmids2) > 1:
        plasmids2 = plasmids2[np.random.choice(range(len(plasmids2)), size=len(plasmids2)//2, replace=False)]

    
    give_plasmids_directed(C, plasmids2)

    #print(C.plasmids.keys())

    return C

def find_best_neighbour(pop, i,j, mode='apoptosis'):
    N = pop.shape[0]    
    if N == 1:
        return False
        
    candidates = get_neighbours(i,j,N)

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

def generate_program_plasmids_function(functions): # functions = [(in1, out1, func1), (in2, out2, func2)...]
    plasmids = {}
    plasmids_ids = []

    for input, output, function in functions:
        label = f'{output}={function}({input})'
        if function == "YES":
            plasmid = generate_plasmid([input], output, models.YES)
        elif function == "NOT":
            plasmid = generate_plasmid([input], output, models.NOT)
        else:
            print("Invalid option!")
            return

        plasmids[label] = plasmid
        plasmids_ids.append(label)


    return plasmids_ids, plasmids

def generate_cells_function(N, functions):
    _, plasmids = generate_program_plasmids_function(functions)

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
            # add all plasmids to a cell          
            for plasmid_id, plasmid in plasmids.items():
                C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])

            population[i,j] = C
            
    return population

        


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
    
    if N_layers == 0:
        ins = [f'in_{j}' for j in range(1, N_inputs+1)]
        i = 0
    else:
        ins = outs
    
    # output layer operons
    Ps, P_ids = generate_plasmids(ins, ["y"], prefix = f'P{i+1}_')
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


    



