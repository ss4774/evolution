import population
import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from multiprocessing import Pool, cpu_count


observables_global = None
observables_local =  None

iterations = 3 # number of learning iterations
repeats = 1 # how many repeats for a given scenario (function, lattice size)
t_end = 200
dt = 0.1
plot_resolution = 1
lattice_sizes = (2,5,7,10) # size of lattice is NxN


def prepare_data(N_inputs):
	if N_inputs == 2:
		func_and = '(in_1 & in_2)'
		states_and = {"in_1":(0,0,10,10),
				"in_2":(0,10,0,10),
				"eval":(0,0,0,10)}

		func_or = '(in_1 | in_2)'
		states_or = {"in_1":(0,0,10,10),
				"in_2":(0,10,0,10),
				"eval":(0,10,10,10)}

		func_xor = '(~(in_1) & in_2) | (in_1 & ~(in_2))'
		states_xor = {"in_1":(0,0,10,10),
				"in_2":(0,10,0,10),
				"eval":(0,10,10,0)}
	elif N_inputs == 3:
		func_and = '(in_1 & in_2 & in_3)'
		states_and = {"in_1":(0,0,0,0,10,10,10,10),
					"in_2":(0,0,10,10,0,0,10,10),
					"in_3":(0,10,0,10,0,10,0,10),
					"eval":(0,0,0,0,0,0,0,10)}

		func_or = '(in_1 | in_2 | in_3)'
		states_or = {"in_1":(0,0,0,0,10,10,10,10),
					"in_2":(0,0,10,10,0,0,10,10),
					"in_3":(0,10,0,10,0,10,0,10),
					"eval":(0,10,10,10,10,10,10,10)}

		func_xor = '(~(in_1) & ~(in_2) & in_3) | (~(in_1) & in_2 & ~(in_3)) | (~(in_1) & ~(in_2) & in_3) | (in_1 & in_2 & in_3)'
		states_xor = {"in_1":(0,0,0,0,10,10,10,10),
					"in_2":(0,0,10,10,0,0,10,10),
					"in_3":(0,10,0,10,0,10,0,10),
					"eval":(0,10,10,0,10,0,0,10)}

	label_and = "and"
	label_or = "or"
	label_xor = "xor"

	func_all = (func_or, func_and, func_xor)
	states_all = (states_or, states_and, states_xor)
	labels = (label_or, label_and, label_xor)

	return func_all, states_all, labels


def work(mode, func, states, label, N, N_inputs):
	p = population.population()
	if mode == "minterms":
		p.generate_cells_minterms(N, N_inputs)
	elif mode == "simple":
		N_layers = 1 # number of internal layers
		N_terms = 2 # number of terms per layer
		p.generate_cells(N, N_inputs, N_layers, N_terms)
	else:
		raise Exception("Incorrect mode: " + mode)

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


def work_consecutively(mode, func_all, states_all, labels, N_inputs):
	for func, states, label in zip(func_all, states_all, labels):
		for N in lattice_sizes:
			for _ in range(repeats):
				work(mode, func, states, label, N, N_inputs)


def prepare_parallel_args(mode, func_all, states_all, labels, N_inputs):
	args = []
	for func, states, label in zip(func_all, states_all, labels):
		for N in lattice_sizes:
			for _ in range(repeats):
				args.append((mode, func, states, label, N, N_inputs))
	return args


def work_parallel(mode, func_all, states_all, labels, N_inputs):
	pool = Pool(cpu_count()) # Set number of processes to the maximum possible value
	args = prepare_parallel_args(mode, func_all, states_all, labels, N_inputs)
	pool.starmap(work, args)


def convergence2(mode, parallel=False):
	func_all, states_all, labels = prepare_data(2) # Prepare data for specific number of intputs
	if parallel:
		return work_parallel(mode, func_all, states_all, labels, 2)
	return work_consecutively(mode, func_all, states_all, labels, 2)


def convergence3(mode, parallel=False):
	func_all, states_all, labels = prepare_data(3) # Prepare data for specific number of intputs
	if parallel:
		return work_parallel(mode, func_all, states_all, labels, 3)
	return work_consecutively(mode, func_all, states_all, labels, 3)


def convergence2_minterms():
	convergence2("minterms")


def convergence2_simple():
	convergence2("simple", True)


def convergence3_minterms():
	convergence3("minterms")


if __name__ == "__main__":
	convergence2_simple()
