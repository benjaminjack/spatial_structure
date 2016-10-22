#! /usr/bin/env python3


import numpy as np

class Patch:
    """
    A single patch in the simulation.
    """
    def __init__(self, row, col, length, cells, phages):
        self.cells = cells
        self.phages = phages
        self.infected_cells = 0
        self.replicated_cells = 0
        self.row = row
        self.col = col
        self.length = length
        self.time = 0

    def execute_until(self, time):
        # Calculate propensities for all reactions
        prop = np.array([0.3 * self.cells * self.phages])
        # Sort propensities
        prop.sort()
        # Find total propensity of the system
        prop_sum = prop.sum()
        # Calculate time until next reaction will occur
        dt = -(1/prop_sum)*np.log(np.random.rand())
        # First check to see if any infected cells should burst
        # If not, we continue with gillespie algorithm
        # Randomly select next reaction to occur
        next_reaction = np.random.multinomial(1, prop/prop_sum).nonzero()[0][0]
        # Update reactant counts based on which reaction is occuring
        print(next_reaction)
        if next_reaction == 0:
            self.infect_cell()
        self.time += dt

    def replicate_cells(self):
        """
        Reaction 0
        """
        pass

    def infect_cell(self):
        """
        Reaction 1
        """
        self.cells -= 1
        self.phages -= 1
        self.infected_cells += 1

    def burst_cell(self):
        """
        Reaction 2
        """
        pass

    def __str__(self):
        return "Phage: " + str(self.phages) + "\nCells: " + str(self.cells) + " " + str(self.time)
