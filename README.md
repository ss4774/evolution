# Synthetic evolution of bacterial populations

This repository is complementing the paper *Programmable evolution of bacterial populations*.

## Main Files

* [`cell.py`](cell.py): a Python class implementing a synthetic cell.
* [`population.py`](population.py): functionalities to generate a synthetic population, distribute the functionalities among the cells in the population, and run the simulations.
* [`models.py`](models.py): Hill functions and protein degradation models composing the models of logic functions.
* [`convergence2_simple.py`](convergence2_simple.py): analysis of convergence for 2-input logic functions using YES and NOT operons.
* [`convergence2_minterms.py`](convergence2_minterms.py): analysis of convergence for 2-input logic functions using operons encoding minterms.
* [`convergence3_minterms.py`](convergence3_minterms.py): analysis of convergence for 3-input logic functions using operons encoding minterms.
* [`plot_convergence.ipynb`](plot_convergence.ipynb): plot the graphs visualising the convergence of solutions.
* [`parameter_values.py`](parameter_values.py): paramemeter values used in simulations.
* [`simulate_2_input.py`](simulate_2_input.py): an example of a simulation of the evolution process using YES and NOT operons.
* [`simulate_2_input_minterms.py`](simulate_2_input_minterms.py): an example of a simulation of the evolution process using operons encoding minterms.
* [`simulate_2_input_fixed_and.py`](simulate_2_input_fixed_and.py): an example of a simulation in which AND function is preset as well as optimised and thus sustained throughout the simulation.
* [`simulate_2_input_fixed_or.py`](simulate_2_misimulate_2_input_fixed_or.py): an example of a simulation in which OR function is preset as well as optimised and thus sustained throughout the simulation.
* [`simulate_2_input_fixed_xor.py`](simulate_2_input_fixed_xor.py): an example of a simulation in which XOR function is preset as well as optimised and thus sustained throughout the simulation.



## Dependencies
The code can be used in a combination with Python 3 programming environment with an additional installation of the following libraries:
* numpy
* matplotlib
* pandas
* sympy

[//]: # (## How to cite this work)
[//]: # (Please cite this work as:)

[//]: # (TODO)

[//]: # (The paper is available at TODO)

## Contact
Please direct your questions and comments to [miha.moskon@fri.uni-lj.si](mailto:miha.moskon@fri.uni-lj.si)

[//]: # (## References)

[//]: # (TODO)