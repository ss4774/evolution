import cell
import models

from itertools import product, combinations
import numpy as np

from parameter_values import *

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
    if pop[candidate].mode == 1:
        if max_plasmids and len(pop[candidate].plasmids) >= max_plasmids: # if this cell already has a maximal number of allowed plasmids
            return
       
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


def apoptosis(pop, i, j, N_inputs, N_layers, N_terms, ordered_crossover=False, amplify_inputs=False):
    # find the best cell in the neighbourhood. If there isn't any, generate a new cell randomly                        
    elitism = False
    crossover = False
    r = np.random.random()

    if r < prob_elitism:
        elitism = find_best_neighbour(pop, i,j)
        if elitism:
            pop[i,j] = elitism
            
    elif r < prob_elitism + prob_cross:    
        crossover =  cross_random_neighbours(pop, i,j, ordered=ordered_crossover)
        if crossover:
            pop[i,j] = crossover

    if (elitism == False) and (crossover == False): # if elitism and crossover failes the cell is randomly initialised
        pop_ = generate_cells(N=1, N_inputs=N_inputs, N_layers=N_layers, N_terms=N_terms, amplify_inputs=amplify_inputs)
        pop[i,j] = pop_[0,0]
        
def cross_random_neighbours(pop, i,j, ordered = False, duplicate=True):
    N = pop.shape[0]    
    if N == 1:
        return False

    candidates = get_neighbours(i,j,N)

    C1, C2 = np.random.choice(range(len(candidates)), size=2, replace=False)
    C1, C2 = pop[candidates[C1]], pop[candidates[C2]]

    C = C1.copy()   # copy the first cell

    if duplicate: # take all plasmids from both cells
        plasmids2 = list(C2.plasmids.items())
        give_plasmids_directed(C, plasmids2) 

    elif ordered: # take the first half from C1 and second half from C2
        plasmid_ids1 = sorted(list(C.plasmids.keys())) # order by alphabet
        plasmid_ids1 = plasmid_ids1[len(plasmid_ids1)//2:] # drop the second half of the plasmids
        for p_id in plasmid_ids1:
            C.lose_plasmid(p_id)
        
        plasmids2 = sorted(list(C2.plasmids.items())) # order by alphabet
        plasmids2 = plasmids2[len(plasmids2)//2:] # keep the second half of the plasmids       
        give_plasmids_directed(C, plasmids2)    
    
    else: # take some random half of plasmids from C1 and some random half of plasmids from C2
        # lose half of the existing plasmids of the first cell
        for i in range(len(C.plasmids)//2):
            C.lose_plasmid()
        # copy half of the plasmids of the second cell
        plasmids2 = np.array(list(C2.plasmids.items()))
        if len(plasmids2) > 1:
            plasmids2 = plasmids2[np.random.choice(range(len(plasmids2)), size=len(plasmids2)//2, replace=False)]  
        give_plasmids_directed(C, plasmids2)    

    C.state['apoptosis'] = 0
    return C

def find_best_neighbour(pop, i,j, criterium='apoptosis'):
    N = pop.shape[0]    
    if N == 1:
        return False
        
    candidates = get_neighbours(i,j,N)

    best_coords = [candidates[0]]
    best_val = pop[candidates[0]].state[criterium]
    
    for candidate in candidates[1:]:
        val = pop[candidate].state[criterium]
        if val == best_val:
            best_coords.append(candidate)        
        elif criterium == 'apoptosis': # lower is better
            if val < best_val:
                best_coords = [candidate]
        elif criterium == 'fitness': # higher is better
            if val > best_val:
                best_coords = [candidate]
    
    if len(best_coords) > 1:
        # if there are multiple equally good candidates just pick one randomly
        best_i,best_j = best_coords[np.random.choice(range(len(best_coords)), size=1, replace=False)[0]]
    else:
        best_i,best_j = best_coords[0]

    C = pop[best_i, best_j].copy()
    C.state['apoptosis'] = 0 

    return C


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

def generate_plasmids_input_amplification(inputs, prefix): # generate input layer plasmids to amplify the inputs using YES and NOT gates
    plasmids = {}
    plasmids_ids = []

    for i, input in enumerate(inputs):
        label1 = f'{prefix}{i}_YES'
        output = f'YES({input})'
        plasmid1 = generate_plasmid([input], output, models.YES)
        plasmids[label1] = plasmid1
        plasmids_ids.append(label1)

        label2 = f'{prefix}{i}_NOT'
        output = f'NOT({input})'
        plasmid2 = generate_plasmid([input], output, models.NOT)
        plasmids[label2] = plasmid2
        plasmids_ids.append(label2)

        
    return plasmids, plasmids_ids


def generate_plasmids_YES_NOT(ins, outs, prefix): # generate plasmids using only YES and NOT gates
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

def generate_plasmids_AND(ins, outs, prefix): # generate plasmids using only YES and NOT gates
    plasmids = {}
    plasmids_ids = []

    ins2 = combinations(ins,2)

    for i, (inputs,output) in enumerate(product(ins2, outs)):
        label = f'{prefix}{i}_AND'
        plasmid = generate_plasmid(list(inputs), output, models.AND)
        plasmids[label] = plasmid
        plasmids_ids.append(label)
        
    return plasmids, plasmids_ids



def generate_program_plasmids(N_inputs, N_layers = 1, N_terms = 3, gates = ["AND"], amplify_inputs = False): # generate program using only YES and NOT gates
    if type(N_terms) == int:
        N_terms = [N_terms]
        N_terms *= N_layers

    layers = {}
    plasmids = {}

    for i, n_terms in enumerate(N_terms):
        if i == 0:            
            ins = [f'in_{j}' for j in range(1, N_inputs+1)]
            if amplify_inputs:
                ins = [f"YES({x})" for x in ins] + [f"NOT({x})" for x in ins]
        else:
            ins = [f'x_{i}_{j}' for j in range(1, n_terms+1)]

        outs = [f'x_{i+1}_{j}' for j in range(1, n_terms+1)]

        # input and middle layer operons
        Ps, P_ids = generate_plasmids_YES_NOT(ins, outs, prefix = f'P{i}_')
        if "AND" in gates:
            Ps2, P_ids2 = generate_plasmids_AND(ins, outs, prefix = f'P{i}_')
            Ps.update(Ps2)
            P_ids.extend(P_ids2)
        layers[f'layer{i}'] = P_ids
        plasmids.update(Ps)
    
    if N_layers == 0:
        ins = [f'in_{j}' for j in range(1, N_inputs+1)]
        i = 0
    else:
        ins = outs
    
    # output layer operons
    Ps, P_ids = generate_plasmids_YES_NOT(ins, ["y"], prefix = f'P{i+1}_')
    if "AND" in gates:
        Ps2, P_ids2 = generate_plasmids_AND(ins, ["y"], prefix = f'P{i+1}_')
        Ps.update(Ps2)
        P_ids.extend(P_ids2)

    layers['layer_out'] = P_ids
    plasmids.update(Ps)

    return plasmids, layers
   

def generate_program_plasmids_function(functions): # functions = [(ins1, out1, func1), (ins2, out2, func2)...]
    plasmids = {}
    plasmids_ids = []

    for inputs, output, function in functions:
        label = f'{output}={function}({",".join(inputs)})'
        if function == "YES":
            plasmid = generate_plasmid(inputs, output, models.YES)
        elif function == "NOT":
            plasmid = generate_plasmid(inputs, output, models.NOT)
        elif function == "AND":
            plasmid = generate_plasmid(inputs, output, models.AND)   
        elif function == "OR":
            plasmid = generate_plasmid(inputs, output, models.OR)  
        elif function == "NOR":
            plasmid = generate_plasmid(inputs, output, models.NOR)
        else:
            print("Invalid option!")
            return

        plasmids[label] = plasmid
        plasmids_ids.append(label)


    return plasmids_ids, plasmids

def generate_cells_function(N, functions, N_inputs, amplify_inputs = False): # generate cells using a specific function - all cells have the same plasmids, the same function
    _, plasmids = generate_program_plasmids_function(functions)

    ins = [f'in_{j}' for j in range(1, N_inputs+1)]
    input_plasmids, _ = generate_plasmids_input_amplification(ins, "P(in)_")

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

            
            #input amplification
            if amplify_inputs:
                
                for plasmid_id, plasmid in input_plasmids.items():
                     C.add_basic_function(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])


            """
            plasmids
            """
            # add all plasmids to a cell          
            for plasmid_id, plasmid in plasmids.items():
                C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])

            population[i,j] = C
            
    return population



def generate_cells(N, N_inputs, N_layers = 1, N_terms = 3, layered=False, amplify_inputs = False):
    plasmids, layers = generate_program_plasmids(N_inputs, N_layers, N_terms, amplify_inputs = amplify_inputs)  
    population = np.zeros((N,N), dtype=object)

    ins = [f'in_{j}' for j in range(1, N_inputs+1)]
    input_plasmids, _ = generate_plasmids_input_amplification(ins, "P(in)_")

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


            #input amplification
            if amplify_inputs:                
                for plasmid_id, plasmid in input_plasmids.items():
                     C.add_basic_function(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])


            """
            plasmids
            """
            # randomly choose plasmids from the available ones          
            cell_plasmids = []     

             
            # if a specific number of plasmids should be taken from each layer
            if layered:
                for layer in layers:
                    #if layer == 'layer_out':
                    #    n = 1 # take only one operon from the output layer
                    #else:                   
                    n = np.random.randint(1, len(layers[layer])) # take a random number of plasmids from a layer
                    if max_plasmids_per_layer:
                        n = min(n, max_plasmids_per_layer)
                    cell_plasmids.extend(np.random.choice(layers[layer], size=n, replace=False))
            else: # if layers are ignored               
                n = np.random.randint(1, len(plasmids)) 
                if max_plasmids:
                    n = min(n, max_plasmids)
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


    



