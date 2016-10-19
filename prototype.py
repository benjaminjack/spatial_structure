#! /usr/bin/env python3

import numpy as np

# Number patch rows
patch_rows = 10
# Number patch columns
patch_columns = 10
# Probability of placing a cell in a patch upon initialization
cell_prob = 0.1

# Make a matrix of cells, where 1 is one cell, and 0 is no cells
cells = np.random.choice([0, 1],
                         size = [patch_rows, patch_columns],
                         p = [1 - cell_prob, cell_prob])


init_phage = np.random.choice([0, 1],
                              size = [patch_rows, patch_columns],
                              p = [0.1, 0.9])

infecting_phage = np.zeros([patch_rows, patch_columns])

for row, col in np.ndindex(cells.shape):
    if cells[row, col] != 0 and init_phage[row, col] != 0:
        infecting_phage[row, col] = 1

free_phage = np.zeros([patch_rows, patch_columns])

print(infecting_phage)
