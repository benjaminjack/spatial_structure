#! /usr/bin/env python3

import numpy as np

class Patch:
    """
    A single patch in the simulation.

    Each patch in the simulation is treated as it's own set of reactions
    executed using the gillespie algorithm. Reactions are executed until the
    next diffusion event, determined by the Plate class, occurs.

    Attributes:
        row: Plate row where this patch is located.
        col: Plate column where this patch is located.
        length: Length and width of this patch (in same units rate constants).
        cells: Number of uninfected cells in patch.
        infected_cells: Number of infected cells in patch which are waiting to
            lyse.
        replicated_cells: Number of cells which have just replicated and are
            waiting to diffuse.
        phages: Number of phages in patch.
        goo: Number of goo particles in patch.
        burst_size: Number of phage progeny that in infected cell releases when
            it lyses.
        k_replicate: Replication rate constant.
        k_infect: Infection rate constant.
        k_lysis: Lysis rate constant.
        k_goo: Rate constant of phage infecting goo particles.
        time: Current time within the patch. This time may not always match the
            Plate time because each patch is an independent set of reactions
            being executed by the gillespie algorithm.

    """
    def __init__(self, row, col, length, cells, phages, goo, burst_size,
                 k_replicate, k_infect, k_lysis, k_goo):
        """
        Inits Patch with row and column positions, dimensions, number of cells,
        number of phage, number of goo particles, and associated rate constants.
        """
        self.cells = cells
        self.phages = phages
        self.infected_cells = 0
        self.replicated_cells = 0
        self.k_replicate = k_replicate
        self.k_lysis = k_lysis
        self.k_infect = k_infect
        self.k_goo = k_goo
        self.goo = goo
        self.burst_size = burst_size
        self.row = row
        self.col = col
        self.length = length
        self.time = 0

    def execute_until(self, time):
        """
        Execute reactions following gillespie algorithm until a given time.

        Args:
            time: Time at which reactions should stop executing.
        """
        while (self.time < time):
            # Calculate propensities for all reactions
            prop = np.array([self.k_infect * self.cells * self.phages,
                             self.k_lysis * self.infected_cells,
                             self.k_replicate * self.cells,
                             self.k_goo * self.goo * self.phages])

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
        Increase newly-replicated cell count by 1.
        """
        self.replicated_cells += 1

    def infect_cell(self):
        """
        Infect a cell. Increase infected cells by 1, decrease cells and phages
        by 1 each.
        """
        self.cells -= 1
        self.phages -= 1
        self.infected_cells += 1

    def burst_cell(self):
        """
        Lyse an infected cell and add new phage particles.
        """
        self.phages += self.burst_size
        self.infected_cells -= 1

    def kill_phage(self):
        """
        Kill a phage particle if it falls into a goo particle.
        """
        self.phages -= 1

    def __str__(self):
        """
        Get the current number of phages, cells, and gillespie simulation time.
        For debugging purposes.
        """
        return "Phage: " + str(self.phages) + "\nCells: " + \
            str(self.cells) + "\nTime:" + str(self.time)

class Plate:
    """
    A plate that contains patches in a grid. The plate controls the diffusion
    of particles between patches.

    Attributes:
        rows: Rows in the grid.
        cols: Columns in the grid.
        length: Length and width of each patch within grid.
        phage_diffusion: Diffusion constant for phage particles.
        cell_diffusion: Diffusion constant for newly-replicated cells.
        max_cell_density: Maximum number of cells per patch.
    """
    def __init__(self, rows, cols, length, phage_diffusion, cell_diffusion,
                 goo, max_cell_density, burst_size, k_infect, k_lysis, k_goo,
                 k_replicate):
        """
        Inits Plate class by constructing a grid of patches.
        """
        self.rows = rows
        self.cols = cols
        self.length = length
        self.phage_diffusion = phage_diffusion
        self.phage_iter = 1 # Phage diffusion events
        self.cell_diffusion = cell_diffusion
        self.cell_iter = 1 # Cell diffusion events
        self.time = 0
        self.max_cell_density = max_cell_density
        # Grid of patches
        self.matrix = np.empty([rows, cols], dtype = object)
        for i in range(rows):
            for j in range(cols):
                # Randomly distribute phage and cell particles
                cells = np.random.randint(0,2)
                phages = np.random.randint(0,10)
                self.matrix[i][j] = Patch(i, j, length, cells, phages, goo,
                                          burst_size, k_replicate, k_infect,
                                          k_lysis, k_goo)

    def iterate(self):
        """
        Execute a diffusion event.
        """
        # Determine which should diffuse first: phages or cells
        dt_phage = (1/4)*((self.length**2)/self.phage_diffusion)*self.phage_iter
        dt_cell = (1/4)*((self.length**2)/self.cell_diffusion)*self.cell_iter
        if dt_phage < dt_cell:
            # Phages diffuse first...
            # Update system time
            self.time += dt_phage
            # Iterate over all patches and execute reactions until the time to
            # the next diffusion event.
            for i in range(self.rows):
                for j in range(self.cols):
                    self.matrix[i][j].execute_until(self.time)
                    self.matrix[i][j].time = self.time
            # Diffuse phages
            self.diffuse_phages()
            self.phage_iter += 1
        else:
            # Cells diffuse first
            self.time += dt_cell
            for i in range(self.rows):
                for j in range(self.cols):
                    self.matrix[i][j].execute_until(self.time)
                    self.matrix[i][j].time = self.time
            self.diffuse_cells()
            self.cell_iter += 1

    def diffuse_phages(self):
        """
        For each phage particle, move the particle up, down, left, right, or
        keep it in the same position, according to a uniform distribution.
        """
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
        """
        Move a phage particle, while checking for boundaries of plate.

        Args:
            i: initial row
            j: initial col
            i2: new row
            j2: new col
        """
        if i2 < self.rows and j2 < self.cols and i2 >= 0 and j2 >= 0:
            # Only move phage particle if it is within bounds of plate.
            self.matrix[i2][j2].phages += 1
            self.matrix[i][j].phages -= 1

    def move_cell(self, i, j, i2, j2):
        """
        Move a cell, while checking that the patch is not at the maximum cell
        density and that the new patch is within plate boundaries. Newly-
        replicated cells are only moved once and then converted to normal cells.
        """
        if i2 < self.rows and j2 < self.cols and i2 >= 0 and j2 >= 0:
            if self.matrix[i2][j2].cells < self.max_cell_density:
                self.matrix[i2][j2].cells += 1
            # If a patch is at max cell density, and the replicated cell is at
            # the edge of the plate, or can't move anywhere, it is effectively
            # destroyed.
            self.matrix[i][j].replicated_cells -= 1

    def diffuse_cells(self):
        """
        Diffuse each newly-replicated cell up, down, left, or right, according
        to a uniform distribution.
        """
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
        """
        Generate matrix of cell and phage counts for debugging and
        visualization purposes.
        """
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

def main():
    """
    Define a plate and iterate through simulation until a certain time point.
    """
    my_plate = Plate(10, # rows
                     10, # cols
                     1, # patch length/width
                     0.1, # cell diffusion
                     0.1, # phage diffusion
                     5, # goo particles
                     3, # max cell density
                     40, # burst size
                     0.3, # infection rate constant
                     0.01, # lysis rate constant
                     0.05, # phage-goo interaction rate constant
                     0.02 # replication rate constant
                     )
    while (my_plate.time < 50):
        my_plate.iterate()

    print(my_plate)

if __name__ == "__main__":
    main()
