import numpy as np
import models

import copy

from params import *

from sympy import symbols
from sympy.logic.boolalg import to_cnf

class cell:
    def __init__(self):
        self.plasmids = {}
        self.basic_functions = {}
        self.state = {}
        self.mod_degradation = {} # proteins with a modified degradation rates

        self.mode = 1 # 0 ... not learning; 1 ... learning; 2 ... optimized

    def copy(self):
        c = cell()

        c.plasmids = copy.deepcopy(self.plasmids)
        c.basic_functions = copy.deepcopy(self.basic_functions)
        c.state = copy.deepcopy(self.state)
        c.mod_degradation = copy.deepcopy(self.mod_degradation)

        c.mode = self.mode

        return c



    """
    changing functionalities of a cell
    """    
    # plasmids and functions - operons
    # dictionary: 
    # - key = label
    # - value = dictionary with keys
    # -- inputs
    # -- output
    # -- f
    def add_operon(self, label, inputs, output, function, add_to, params=None, mod_degradation=None):
        d = {}
        d['inputs'] = inputs
        d['output'] = output
        d['f'] = function        
        
        if params:
            d['params'] = params

        if mod_degradation:
            self.mod_degradation[output] = mod_degradation
            d['mod_degradation'] = mod_degradation

        add_to[label] = d

        for x in inputs + [output]:
            if x not in self.state:
                self.state[x] = 0

    def add_basic_function(self, label, inputs, output, function, params=None, mod_degradation=None):
        self.add_operon(label, inputs, output, function, self.basic_functions, params=params, mod_degradation=mod_degradation)

    def add_plasmid(self, label, inputs, output, function, params=None, mod_degradation=None):
        self.add_operon(label, inputs, output, function, self.plasmids, params=params, mod_degradation=mod_degradation)

    def lose_plasmid(self):
        if self.plasmids:
            plasmid_ids = list(self.plasmids.keys())
            plasmid_id = np.random.choice(plasmid_ids, size=1, replace=False)[0]
            plasmid = self.plasmids[plasmid_id]
            del self.plasmids[plasmid_id]
                            
            if 'mod_degradation' in plasmid: # if output of this plasmid has a modified degradation we will potentially delete it from mod_degradation. We need to check if it has modified degradation also due to other plasmids
                output = plasmid['output']
                for p in self.plasmids.values():
                    if (output == p['output']) and 'mod_degradation' in p:
                        break
                else:
                    del self.mod_degradation[output]

    def copy_plasmid(self, rate, dt):      
        if self.plasmids:
            if np.random.random() < rate*dt:
                plasmid_ids = list(self.plasmids.keys())
                plasmid_id = np.random.choice(plasmid_ids, size=1, replace=False)[0]
                return (plasmid_id, self.plasmids[plasmid_id])
        

    """
    state update functions
    """
    def update_operons(self, d_state, operons):
        for _, operon in operons.items():
            f = operon['f']
            input_labels = operon['inputs']
            inputs = []
            for i in input_labels:
                inputs.append(self.state[i])          
            
            if 'params' in operon:
                d_state[operon['output']] += f(*inputs, operon['params'])    
            else:
                d_state[operon['output']] += f(*inputs)

    def update_plasmids(self, d_state):
        self.update_operons(d_state, self.plasmids)

    def update_basic_functions(self, d_state):
        self.update_operons(d_state, self.basic_functions)

    def update(self, dt=0.1, global_vars={}):
        d_state = {}
        self.state.update(global_vars)
        
        # go through the state and apply degradation to all species
        for x, conc in self.state.items():
            if x in self.mod_degradation:
                d_state[x] = models.degrade(conc, self.mod_degradation[x])             
            else:
                d_state[x] = models.degrade(conc)   

        # go trough all the plasmids and update the state
        self.update_plasmids(d_state)

        # go through all the basic_functions and update the state
        self.update_basic_functions(d_state)       

        # apply the changes
        for x, dx in d_state.items():
            self.state[x] += dx * dt
                 
        # switching between cell modes (0 ... not learning; 1 ... learning; 2 ... optimized)
        if self.mode > 0:
            if ('learn' in self.state) and (self.state['learn'] <= 0):                
                self.mode = 0              
            elif ('fitness' in self.state) and (self.state['fitness'] > fitness_threshold):                
                self.mode = 2
            elif ('fitness' in self.state) and (self.state['fitness'] < fitness_threshold):
                self.mode = 1        

        # learning functions
        if self.mode == 1:
            # delete plasmids        
            #lose_plasmid_rate = max_lose_plasmid_rate * np.exp(-(apoptosis_threshold-self.state['apoptosis']))    # lose plasmid rate is proportional with apoptosis signal            
            lose_plasmid_rate = max_lose_plasmid_rate * np.exp(-self.state['fitness'])    # lose plasmid rate is proportional with apoptosis signal            
            if np.random.random() < lose_plasmid_rate*dt:
                self.lose_plasmid()

            # apoptosis - if apoptosis should be triggered, the update function will return a string "apoptosis"    
            if ('apoptosis' in self.state) and (self.state['apoptosis'] > apoptosis_threshold):
                return 'apoptosis'

        # copy plasmid - if a plasmid will be copied, the update function will return a plasmid information to copy to another cell
        if self.mode > 0:
            #copy_plasmid_rate = max_copy_plasmid_rate * np.exp(-self.state['apoptosis']) # copy plasmid rate is inversely proportional with apoptosis signal
            copy_plasmid_rate = max_copy_plasmid_rate * np.exp(max(0, -(fitness_threshold - self.state['fitness']))) # copy plasmid rate is inversely proportional with apoptosis signal
            p = self.copy_plasmid(copy_plasmid_rate, dt)
            if p:
                return p



    """
    other functions
    """
    def decode_function(self):
        functions = {}            
        plasmids_ids = sorted(list(self.plasmids.keys()))

        inputs = []

        for p in plasmids_ids:
            plasmid = self.plasmids[p]
            output = plasmid['output']
            input = plasmid['inputs'][0]            
            inputs.append(input)


            if input in functions:
                input = functions[input]
            
            if 'NOT' in p:
                func = f'NOT({input})'
            else:
                func = input

            if output in functions:
                functions[output] += f' OR {func}'
            else:
                functions[output] = func
        if 'y' in functions:
            syms = symbols(",".join(inputs))
            if type(syms) != 'tuple':
                syms = (syms,)
            
            for input, sym in zip(inputs, syms):
                locals()[input] = sym

            f = functions['y']
            f = f.replace("OR", "|").replace("NOT", "~")
            f = to_cnf(f, True)
            f = str(f)

            return f"y={f}"

if __name__ == "__main__":
    pass

