#! /usr/bin/env python3

import numpy as np

class Cell:
    """
    A class that defines a Cell that can replicate and become infected by a
    Phage.
    """
    def __init__(self, row, col, replication_time = 10):
        self.replication_time = replication_time
        self.infected = False
        self.attached_phage = 0
        self.age = 1
        self.row = row
        self.col = col

    def replicate(self):
        # Randomly number between 0 - 9
        move = self.to_coords(np.random.randint(9))
        return (self.row + move[0], self.col + move[1])

    def to_coords(self, num):
        """
        Convert a number between 1 and 9 to a set of relative coordinates.
        """
        return {
            0: (0, 0),
            1: (0, 1),
            2: (0, -1),
            3: (1, 0),
            4: (-1, 0),
            5: (-1, 1),
            6: (1, 1),
            7: (1, -1),
            8: (-1, -1)
        }[num]

class Plate:
    """
    A class that defines a Plate object containing a matrix of patches.
    """
    def __init__(self, rows, cols, cell_fill, phage_fill, max_cell_density):
        # Empty for now because cells replicate at a fixed interval
        self.rows = rows
        self.cols = cols
        self.phage_fill = phage_fill
        self.cell_fill = cell_fill
        self.max_cell_density = max_cell_density
        self.indices = np.ndindex((rows, cols))

        # I'll worry about nutrients later
        # self.nutrient_matrix = np.zeros([rows, cols])
        # Populate cells
        self.populate_cells()
        self.populate_phages()

    def populate_cells(self):
        # A matrix to track the number of cells in each patch
        self.cell_matrix = np.random.choice(
            [0, 1],
            size = [self.rows, self.cols],
            p = [1 - self.cell_fill, self.cell_fill]
        )
        # List to store cell objects
        self.cell_list = []
        # Populate cell list with cell objects
        for row, col in self.indices:
            if self.cell_matrix[row, col] == 1:
                new_cell = Cell(row, col)
                self.cell_list.append(new_cell)

    def populate_phages(self):
        # A matrix to track the number of phages in each patch
        self.phage_matrix = np.random.choice(
            [0, 1],
            size = [self.rows, self.cols],
            p = [1 - self.phage_fill, self.phage_fill]
        )

    def iterate(self):
        # First iterate over Cells
        new_cells = []
        for cell in self.cell_list:
            if cell.infected == True:
                if cell.lysis_time == 0:
                    cell.infected == False
                    self.phage_matrix[cell.row, cell.col] += 40
                else:
                    cell.lysis_time -= 1
            elif cell.age < cell.replication_time:
                cell.age += 1
            else:
                (new_row, new_col) = cell.replicate()
                if self.is_occupied(new_row, new_col) == False:
                    new_cell = Cell(new_row, new_col)
                    new_cells.append(new_cell)
                    self.cell_matrix[new_row, new_col] += 1
                    cell.age = 1
        self.cell_list.extend(new_cells)

        # Now iterate over phage matrix and diffuse phage
        new_phage_matrix = np.zeros([self.rows, self.cols], dtype='int')
        for row, col in np.ndindex(self.phage_matrix.shape):
            for i in range(self.phage_matrix[row, col]):
                move = self.to_coords(np.random.randint(9))
                new_row = row + move[0]
                new_col = col + move[1]
                if new_row < self.rows and new_row >=0 and new_col < self.cols and new_col >= 0:
                    new_phage_matrix[new_row, new_col] += 1
                else:
                    new_phage_matrix[row, col] += 1
        self.phage_matrix = new_phage_matrix

        for row, col in np.ndindex(self.phage_matrix.shape):
            cell_count = self.cell_matrix[row, col]
            phage_count = self.phage_matrix[row, col]
            if cell_count > 0 and phage_count > 0:
                self.infect_cell(row, col)



    def infect_cell(self, row, col):
        for cell in self.cell_list:
            if cell.row == row and cell.col == col:
                if cell.infected == False:
                    cell.infected = True
                    cell.lysis_time = 10
                    break
        self.cell_matrix[row, col] -= 1

    def is_occupied(self, row, col):
        """
        Returns True if a given set of coordinates extends outside of the plate
        or if a given patch has reached its maximum cell-carrying capacity.
        """
        if row > self.rows - 1 or row < 0 or col > self.cols - 1 or col < 0:
            return True
        if self.cell_matrix[row, col] >= self.max_cell_density:
            return True
        return False

    def to_coords(self, num):
        """
        Convert a number between 1 and 9 to a set of relative coordinates.
        """
        return {
            0: (0, 0),
            1: (0, 1),
            2: (0, -1),
            3: (1, 0),
            4: (-1, 0),
            5: (-1, 1),
            6: (1, 1),
            7: (1, -1),
            8: (-1, -1)
        }[num]

    def __str__(self):
        out = "Cell matrix: \n" + str(self.cell_matrix)
        out += "\n\nPhage matrix: \n" +str(self.phage_matrix)
        return out

my_plate = Plate(20, 20, 0.1, 0.1, 3)

print(my_plate)

for i in range(100):
    my_plate.iterate()

print(my_plate)
