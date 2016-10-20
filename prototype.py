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

class Phage:
    """
    A class that defines a Phage particle that attaches to and infects Cell
    objects.
    """
    def __init__(self, row, col, lysis_time):
        self.row = row
        self.col = col
        self.lysis_time = lysis_time

    def diffuse(self):
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
        # List to store phage objects
        self.phage_list = []
        # Populate phage list with Phage objects
        for row, col in self.indices:
            if self.phage_matrix[row, col] == 1:
                new_phage = Phage(row, col)
                self.cell_list.append(new_phage)

    def iterate(self):
        # First iterate over Cells
        new_cells = []
        for cell in self.cell_list:
            if cell.age < cell.replication_time:
                cell.age += 1
            else:
                (new_row, new_col) = cell.replicate()
                if self.is_occupied(new_row, new_col) == False:
                    new_cell = Cell(new_row, new_col)
                    new_cells.append(new_cell)
                    self.cell_matrix[new_row, new_col] += 1
                    cell.age = 1
        self.cell_list.extend(new_cells)
        # Now iterate over phages
        new_phage = []
        for phage in self.phage_list:
            # Check if phage is in same patch as a cell
            if cell_matrix[phage.row, phage.column] > 0:
                


    def is_occupied(self, row, col):
        """
        Returns True if a given set of coordinates extends outside of the plate
        or if a given patch has reached its maximum cell-carrying capacity.
        """
        if row > self.rows - 1 or col > self.cols - 1:
            return True
        if self.cell_matrix[row, col] >= self.max_cell_density:
            return True
        return False

    def __str__(self):
        return "Cell matrix: \n" + str(self.cell_matrix)


my_plate = Plate(10, 10, 0.1, 3)

for i in range(1000):
    my_plate.iterate()

print(my_plate)
