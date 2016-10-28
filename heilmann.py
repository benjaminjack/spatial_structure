import numpy as np

delta_t = 1
k_phage_jumping = 10 # lambda
burst_size = 80 # Beta
replication_time = 300 # T
k_degradation = 10**-2 # delta
k_infection = 10**-3 # alpha

rows = 10
cols = 10

phage = np.random.multinomial(size=(rows, cols), dtype='int')
cells = np.zeros((rows, cols), dtype='int')
lysis = np.ones((rows, cols), dtype='int')*-1
replication = np.ones((rows, cols), dtype='int')*-1

p_decay = 1 - np.exp(-k_degradation * delta_t)

time = 0

def iterate():
    for i in range(rows):
        for j in range(cols):
            if infected_cells[i][j] == time:
                # Lyse cell
                infected_cells[i][j] == -1
                phage[i][j] += burst_size
