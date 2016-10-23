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
        self.goo = 5
        self.row = row
        self.col = col
        self.length = length
        self.time = 0

    def execute_until(self, time):
        while (self.time < time):

            # Calculate propensities for all reactions
            prop = np.array([0.3 * self.cells * self.phages,
                             0.1 * self.infected_cells,
                             0.05 * self.cells,
                             0.05 * self.goo * self.phages])

            # Find total propensity of the system
            prop_sum = prop.sum()

            # End the simulation if propensities are 0 and no cells are infected
            if prop_sum == 0:
                break

            # Calculate time until next reaction will occur
            dt = -(1/prop_sum)*np.log(np.random.rand())
            # Update current time in this patch
            self.time += dt
            # If we've gone past the time for the next diffusion event, break
            # out of simulation and do not execute next reaction
            if self.time > time:
                break

            # Randomly select next reaction to occur
            next_reaction = np.random.multinomial(1,
                                                  prop/prop_sum).nonzero()[0][0]
            # Execute next reaction
            if next_reaction == 0:
                self.infect_cell()
            elif next_reaction == 1:
                self.burst_cell()
            elif next_reaction == 2:
                self.replicate_cell()
            elif next_reaction == 3:
                self.kill_phage()

    def replicate_cell(self):
        """
        Reaction 0
        """
        self.replicated_cells += 1

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
        self.phages += 40
        self.infected_cells -= 1

    def kill_phage(self):
        self.phages -= 1

    def __str__(self):
        return "Phage: " + str(self.phages) + "\nCells: " + str(self.cells) + " " + str(self.time)

class Plate:

    def __init__(self, rows, cols, length):
        self.rows = rows
        self.cols = cols
        self.length = length
        self.phage_diffusion = 0.1
        self.phage_iter = 1
        self.cell_diffusion = 0.05
        self.cell_iter = 1
        self.time = 0
        self.matrix = np.empty([rows, cols], dtype = object)
        for i in range(rows):
            for j in range(cols):
                cells = np.random.randint(0,2)
                phages = np.random.randint(0,10)
                self.matrix[i][j] = Patch(i, j, length, cells, phages)

    def iterate(self):
        dt_phage = (1/4)*((self.length**2)/self.phage_diffusion)*self.phage_iter
        dt_cell = (1/4)*((self.length**2)/self.cell_diffusion)*self.cell_iter
        if dt_phage < dt_cell:
            self.time += dt_phage
            for i in range(self.rows):
                for j in range(self.cols):
                    self.matrix[i][j].execute_until(self.time)
                    self.matrix[i][j].time = self.time
            self.diffuse_phages()
            self.phage_iter += 1
        else:
            self.time += dt_cell
            for i in range(self.rows):
                for j in range(self.cols):
                    self.matrix[i][j].execute_until(self.time)
                    self.matrix[i][j].time = self.time
            self.diffuse_cells()
            self.cell_iter += 1

    def diffuse_phages(self):
        for i in range(self.rows):
            for j in range(self.cols):
                for p in range(self.matrix[i][j].phages):
                    direction = np.random.randint(5)
                    if direction == 0:
                        self.move_phage(i, j, i + 1, j)
                    elif direction == 1:
                        self.move_phage(i, j, i - 1, j)
                    elif direction == 2:
                        self.move_phage(i, j, i, j + 1)
                    elif direction == 3:
                        self.move_phage(i, j, i, j - 1)


    def move_phage(self, i, j, i2, j2):
        if i2 < self.rows and j2 < self.cols and i2 >= 0 and j2 >= 0:
            self.matrix[i2][j2].phages += 1
            self.matrix[i][j].phages -= 1

    def move_cell(self, i, j, i2, j2):
        if i2 < self.rows and j2 < self.cols and i2 >= 0 and j2 >= 0:
            if self.matrix[i2][j2].cells < 3:
                self.matrix[i2][j2].cells += 1
            self.matrix[i][j].replicated_cells -= 1

    def diffuse_cells(self):
        for i in range(self.rows):
            for j in range(self.cols):
                for p in range(self.matrix[i][j].replicated_cells):
                    direction = np.random.randint(5)
                    if direction == 0:
                        self.move_cell(i, j, i + 1, j)
                    elif direction == 1:
                        self.move_cell(i, j, i - 1, j)
                    elif direction == 2:
                        self.move_cell(i, j, i, j + 1)
                    elif direction == 3:
                        self.move_cell(i, j, i, j - 1)
                    elif direction == 4:
                        self.move_cell(i, j, i, j)

    def __str__(self):
        phages = ""
        cells = ""
        for i in range(self.rows):
            for j in range(self.cols):
                cells += str(self.matrix[i][j].cells) + "\t"
                phages += str(self.matrix[i][j].phages) + "\t"
            cells += "\n"
            phages += "\n"
        out = "Phages: \n" + phages
        out += "\nCells: \n" + cells + "\n"
        return out
