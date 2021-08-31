import cell
import models

from itertools import product, combinations
from collections import defaultdict

import numpy as np

import parameter_values

import pandas as pd

class population_generator:
    def __init__(self,  neighbourhood = "moore", ordered_crossover=True):        
        self.neighbourhood = neighbourhood
        self.pop = None
        
        self.N = 0
        self.N_inputs = 2
        self.N_layers = 1
        self.N_terms = 2

        self.ordered_crossover=ordered_crossover
        self.generate_minterms = False
        
        self.possible_plasmids = None
        self.layers = None


    def get_neighbours(self, i,j):
        N = self.N
        neighbourhood = self.neighbourhood

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
    def give_plasmid(self, i, j, plasmid_pair):
        plasmid_id, plasmid = plasmid_pair
        
        N = self.N
        pop = self.pop
        

        if N <= 1:
            return
        #H = [(i-1) % N, (i+1) % N]
        #V = [(j-1) % N, (j+1) % N]
        #candidates = list(product(H,V))
        candidates = self.get_neighbours(i,j)
        
        candidate = candidates[np.random.choice(range(len(candidates)), size=1, replace=False)[0]]

        #print(pop[candidate].decode_function())    
        #print(pop[candidate].plasmids)       
        if pop[candidate].mode == 1:
            if parameter_values.max_plasmids and len(pop[candidate].plasmids) >= parameter_values.max_plasmids: # if this cell already has a maximal number of allowed plasmids
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
    def give_plasmids_directed(self, C, plasmid_pairs): # C...cell to which the plasmids are copied; plasmid_pairs = [(plasmid_id1, plasmid1),(plasmid_id2, plasmid2),...]
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


    def apoptosis(self, i, j):
        pop = self.pop
        
        # find the best cell in the neighbourhood. If there isn't any, generate a new cell randomly                        
        elitism = False
        crossover = False
        r = np.random.random()

        if r < parameter_values.prob_elitism:
            elitism = self.find_best_neighbour(i,j)
            if elitism:
                pop[i,j] = elitism
                
        elif r < parameter_values.prob_elitism + parameter_values.prob_cross:            
            crossover =  self.cross_random_neighbours(i,j)
            if crossover:
                pop[i,j] = crossover

        if (elitism == False) and (crossover == False): # if elitism and crossover failes the cell is randomly initialised
            if self.generate_minterms:
                C = self.generate_cell_minterms()            
            else:            
                C = self.generate_cell()
            
            pop[i,j] = C
            
    def cross_random_neighbours(self, i,j, duplicate=False):
        
        pop = self.pop
        N = self.N
        ordered = self.ordered_crossover

        if N == 1:
            return False


        candidates = self.get_neighbours(i,j)

        C1, C2 = np.random.choice(range(len(candidates)), size=2, replace=False)
        C1, C2 = pop[candidates[C1]], pop[candidates[C2]]

        C = C1.copy()   # copy the first cell

        if duplicate: # take all plasmids from both cells
            plasmids2 = list(C2.plasmids.items())
            self.give_plasmids_directed(C, plasmids2) 

        elif ordered: # take the first half from C1 and second half from C2
            plasmid_ids1 = sorted(list(C.plasmids.keys())) # order by alphabet
            plasmid_ids1 = plasmid_ids1[len(plasmid_ids1)//2:] # drop the second half of the plasmids
            for p_id in plasmid_ids1:
                C.lose_plasmid(p_id)
            
            plasmids2 = sorted(list(C2.plasmids.items())) # order by alphabet
            plasmids2 = plasmids2[len(plasmids2)//2:] # keep the second half of the plasmids       
            self.give_plasmids_directed(C, plasmids2)    
        
        else: # take some random half of plasmids from C1 and some random half of plasmids from C2
            # lose half of the existing plasmids of the first cell
            for i in range(len(C.plasmids)//2):
                C.lose_plasmid()
            # copy half of the plasmids of the second cell
            plasmids2 = np.array(list(C2.plasmids.items()))
            if len(plasmids2) > 1:
                plasmids2 = plasmids2[np.random.choice(range(len(plasmids2)), size=len(plasmids2)//2, replace=False)]  
            self.give_plasmids_directed(C, plasmids2)    

        C.state['apoptosis'] = 0
        return C

    def find_best_neighbour(self, i,j, criterium='apoptosis'):
        N = self.N
        pop = self.pop

        if N == 1:
            return False
            
        candidates = self.get_neighbours(i,j)

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


    def generate_plasmid(self, inputs, output, f, params = None, mod_degradation=None):
        d = {}
        d['inputs'] = inputs
        d['output'] = output
        d['f'] = f
        if params:
            d['params'] = params
        if mod_degradation:
            d['mod_degradation'] = mod_degradation

        return d
            
    def generate_plasmids_YES_NOT(self, ins, outs, prefix): # generate plasmids using only YES and NOT gates
        plasmids = {}
        plasmids_ids = []

        for i, (input,output) in enumerate(product(ins, outs)):
            label1 = f'{prefix}{i}_YES'
            plasmid1 = self.generate_plasmid([input], output, models.YES)
            plasmids[label1] = plasmid1
            plasmids_ids.append(label1)

            label2 = f'{prefix}{i}_NOT'
            plasmid2 = self.generate_plasmid([input], output, models.NOT)
            plasmids[label2] = plasmid2
            plasmids_ids.append(label2)
            
        return plasmids, plasmids_ids

    def generate_plasmids_AND(self, ins, outs, prefix): # generate plasmids using only YES and NOT gates
        plasmids = {}
        plasmids_ids = []

        ins2= combinations(ins,2)

        for i, (inputs,output) in enumerate(product(ins2, outs)):
            label = f'{prefix}{i}_AND'
            plasmid = self.generate_plasmid(list(inputs), output, models.AND)
            plasmids[label] = plasmid
            plasmids_ids.append(label)
            
        return plasmids, plasmids_ids

    def generate_cells_function(self, N, functions):
        self.N = N

        _, plasmids = self.generate_program_plasmids_function(functions)

        population = np.zeros((N,N), dtype=object)

        for i in range(N):      
            for j in range(N):

                C = cell.cell()
                
                """
                basic functions
                """
                # fitness
                C.add_basic_function("fitness_operon", ["eval", "y"], "fitness", models.EQU, params = (parameter_values.alpha_fitness, parameter_values.Kd_fitness, parameter_values.n_fitness), mod_degradation=parameter_values.delta_fitness)
                
                # apoptosis
                C.add_basic_function("apoptosis_operon", ["fitness"], "apoptosis", models.NOT, params=(parameter_values.alpha_apoptosis, parameter_values.Kd_apoptosis, parameter_values.n_apoptosis), mod_degradation = parameter_values.delta_apoptosis)

                """
                plasmids
                """
                # add all plasmids to a cell          
                for plasmid_id, plasmid in plasmids.items():
                    C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])

                population[i,j] = C
                
        self.pop = population

    def generate_program_plasmids(self, gates = ["AND"]): # generate program using only YES and NOT gates; 
        N_inputs = self.N_inputs
        N_layers = self.N_layers
        N_terms = self.N_terms

        if type(N_terms) == int:
            N_terms = [N_terms]
            N_terms *= N_layers

        layers = {}
        plasmids = {}

        for i, n_terms in enumerate(N_terms):
            if i == 0:
                ins = [f'in_{j}' for j in range(1, N_inputs+1)]
            else:
                if parameter_values.full_connectivity:
                    ins.extend([f'x_{i}_{j}' for j in range(1, n_terms+1)])
                else:
                    ins = [f'x_{i}_{j}' for j in range(1, n_terms+1)]

            outs = [f'x_{i+1}_{j}' for j in range(1, n_terms+1)]

            # input and middle layer operons
            Ps, P_ids = self.generate_plasmids_YES_NOT(ins, outs, prefix = f'P{i}_')
            if "AND" in gates:
                Ps2, P_ids2 = self.generate_plasmids_AND(ins, outs, prefix = f'P{i}_')
                Ps.update(Ps2)
                P_ids.extend(P_ids2)
            layers[f'layer{i}'] = P_ids
            plasmids.update(Ps)
        
        if N_layers == 0:
            ins = [f'in_{j}' for j in range(1, N_inputs+1)]
            i = 0
        else:
            if parameter_values.full_connectivity:
                ins += outs
            else:   
                ins = outs
        
        # output layer operons
        Ps, P_ids = self.generate_plasmids_YES_NOT(ins, ["y"], prefix = f'P{i+1}_')
        if "AND" in gates:
            Ps2, P_ids2 = self.generate_plasmids_AND(ins, ["y"], prefix = f'P{i+1}_')
            Ps.update(Ps2)
            P_ids.extend(P_ids2)

        layers['layer_out'] = P_ids
        plasmids.update(Ps)

        return plasmids, layers

    def generate_program_plasmids_minterms(self): # generate all minterms for a given number of inputs    
                
        N_inputs = self.N_inputs

        plasmids = {}
        plasmids_ids = []

        ins = [f'in_{j}' for j in range(1, N_inputs+1)]
        ins2= combinations(ins,2)

        # only for two input gates 
        # TODO: an arbitrary number of inputs
        for i, inputs in enumerate(ins2):
            inputs = list(inputs)
            
            label = f'P0_{i}_AND00'
            output = 'y'
            plasmid = self.generate_plasmid(inputs, output, models.AND00)
            plasmids[label] = plasmid
            plasmids_ids.append(label)

            label = f'P0_{i}_AND01'
            output = 'y'
            plasmid = self.generate_plasmid(inputs, output, models.AND01)
            plasmids[label] = plasmid
            plasmids_ids.append(label)

            label = f'P0_{i}_AND10'
            output = 'y'
            plasmid = self.generate_plasmid(inputs, output, models.AND10)
            plasmids[label] = plasmid
            plasmids_ids.append(label)

            label = f'P0_{i}_AND11'
            output = 'y'
            plasmid = self.generate_plasmid(inputs, output, models.AND11)
            plasmids[label] = plasmid
            plasmids_ids.append(label)

            
        return plasmids, plasmids_ids


    def generate_program_plasmids_function(self, functions): # functions = [(ins1, out1, func1), (ins2, out2, func2)...]
        plasmids = {}
        plasmids_ids = []

        for inputs, output, function in functions:
            label = f'{output}={function}({",".join(inputs)})'
            if function == "YES":
                plasmid = self.generate_plasmid(inputs, output, models.YES)
            elif function == "NOT":
                plasmid = self.generate_plasmid(inputs, output, models.NOT)
            elif function == "AND":
                plasmid = self.generate_plasmid(inputs, output, models.AND)   
            elif function == "AND00":
                plasmid = self.generate_plasmid(inputs, output, models.AND00)   
            elif function == "AND01":
                plasmid = self.generate_plasmid(inputs, output, models.AND01)               
            elif function == "AND10":
                plasmid = self.generate_plasmid(inputs, output, models.AND10)   
            elif function == "AND11":
                plasmid = self.generate_plasmid(inputs, output, models.AND11)   
            elif function == "OR":
                plasmid = self.generate_plasmid(inputs, output, models.OR)  
            elif function == "NOR":
                plasmid = self.generate_plasmid(inputs, output, models.NOR)
            else:
                print("Invalid option!")
                return

            plasmids[label] = plasmid
            plasmids_ids.append(label)


        return plasmids_ids, plasmids

    def generate_cell(self, layered=False):
        C = cell.cell()

        if self.possible_plasmids:
            plasmids, layers = self.possible_plasmids, self.layers             
        else:
            self.possible_plasmids, self.layers = self.generate_program_plasmids() 
            plasmids, layers = self.possible_plasmids, self.layers   

        """
        basic functions
        """
        # fitness
        C.add_basic_function("fitness_operon", ["eval", "y"], "fitness", models.EQU, params = (parameter_values.alpha_fitness, parameter_values.Kd_fitness, parameter_values.n_fitness), mod_degradation=parameter_values.delta_fitness)
        
        # apoptosis
        C.add_basic_function("apoptosis_operon", ["fitness"], "apoptosis", models.NOT, params=(parameter_values.alpha_apoptosis, parameter_values.Kd_apoptosis, parameter_values.n_apoptosis), mod_degradation = parameter_values.delta_apoptosis)

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
                if parameter_values.max_plasmids_per_layer:
                    n = min(n, parameter_values.max_plasmids_per_layer)
                cell_plasmids.extend(np.random.choice(layers[layer], size=n, replace=False))
        else: # if layers are ignored               
            n = np.random.randint(1, len(plasmids)) 
            if parameter_values.max_plasmids:
                n = min(n, parameter_values.max_plasmids)
            cell_plasmids.extend(np.random.choice(list(plasmids.keys()), size=n, replace=False))
    
        for plasmid_id in cell_plasmids:
            plasmid = plasmids[plasmid_id]

            #if 'mod_degradation' in plasmid:
            #    C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'], mod_degradation = plasmid['mod_degradation'])
            #else:    
            C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])

        return C

    def generate_cells(self, N, N_inputs, N_layers, N_terms):
        
        self.generate_minterms = False  
        
        self.N = N
        self.N_inputs = N_inputs
        self.N_layers = N_layers
        self.N_terms = N_terms
        
        self.possible_plasmids, self.layers = self.generate_program_plasmids() 

        population = np.zeros((N,N), dtype=object)

        for i in range(N):      
            for j in range(N):

                C = self.generate_cell()                
                population[i,j] = C
                
        self.pop = population

    def generate_cell_minterms(self):
        C = cell.cell()

        if self.possible_plasmids:
            plasmids = self.possible_plasmids
        else:
            plasmids, _ = self.generate_program_plasmids_minterms()
            self.possible_plasmids = plasmids
        
        """
        basic functions
        """
        # fitness
        C.add_basic_function("fitness_operon", ["eval", "y"], "fitness", models.EQU, params = (parameter_values.alpha_fitness, parameter_values.Kd_fitness, parameter_values.n_fitness), mod_degradation=parameter_values.delta_fitness)
        
        # apoptosis
        C.add_basic_function("apoptosis_operon", ["fitness"], "apoptosis", models.NOT, params=(parameter_values.alpha_apoptosis, parameter_values.Kd_apoptosis, parameter_values.n_apoptosis), mod_degradation = parameter_values.delta_apoptosis)

        """
        plasmids
        """
        # randomly choose plasmids from the available ones          
        cell_plasmids = []     
                                                                
        n = np.random.randint(1, len(plasmids)) 
        if parameter_values.max_plasmids:
            n = min(n, parameter_values.max_plasmids)
        cell_plasmids.extend(np.random.choice(list(plasmids.keys()), size=n, replace=False))
    
        for plasmid_id in cell_plasmids:
            plasmid = plasmids[plasmid_id]

            C.add_plasmid(plasmid_id, plasmid['inputs'], plasmid['output'], plasmid['f'])

        return C        

    def generate_cells_minterms(self, N, N_inputs):
        self.generate_minterms = True # for latter calls of generation (e.g., when apoptosis occurs)
        self.N = N
        self.N_inputs = N_inputs
        
        population = np.zeros((N,N), dtype=object)

        plasmids, _ = self.generate_program_plasmids_minterms()
        self.possible_plasmids = plasmids

        for i in range(N):      
            for j in range(N):

                C = self.generate_cell_minterms()
                    
                population[i,j] = C
                
        self.pop = population

    def simulate(self, states, observables_local, observables_global, t_end, dt=0.1, iterations = 1, plot_resolution = 1):

        N = self.N
        pop = self.pop

        df = pd.DataFrame(dtype=float)
        functions = {}

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
                                self.give_plasmid(i, j, plasmid_pair)
                            elif action == "apoptosis":
                                self.apoptosis(i, j)                        
                        
                        if track:
                            functions[T][pop[i,j].decode_function()] += 1

                if track:     
                    df = df.append(d, ignore_index=True, sort=False)        

        return df, functions

if __name__ == "__main__":
    #generate_cells(2, 3, 2, 3)
    pass


        



