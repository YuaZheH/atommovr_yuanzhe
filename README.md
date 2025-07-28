# atommovr: a simulation framework for atom rearrangement algorithms.

Authors: Nikhil Harle*, Bo-Yu Chen*, Bob Bao, Hannes Bernien.

This is a Python-based framework for comparing algorithms for rearranging neutral atoms in tweezer arrays or optical lattices. Developed in the Bernien Lab at UChicago.

Contents:
- Tutorial: `tutorial_clean.ipynb`.
- Implementations of individual algorithms: `atommover.algorithms/`
- Stochastic loading: `atommover.utils.core.py`
- AOD commands and implementing moves: `atommover.utils.move_utils.py`

Installation instructions:
- To make the virtual environment, please navigate to this folder and run the following command: `conda env create -f environment.yml`.

