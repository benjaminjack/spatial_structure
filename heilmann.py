import numpy as np

delta_t = 1
k_diffuse = 10 # lambda
burst_size = 80 # Beta
replication_time = 300 # T
lysis_time = 10
k_degradation = 10**-2 # delta
k_infection = 10**-3 # alpha

rows = 10
cols = 10

phage = np.random.choice([0,100], p = [0.9, 0.1], size=(rows, cols))
cells = np.random.choice([0,1], p = [0.9, 0.1], size=(rows, cols))
lysis = np.ones((rows, cols), dtype='int')*-1
replication = np.ones((rows, cols), dtype='int')*-1 + 101*cells

p_decay = 1 - np.exp(-k_degradation * delta_t)
p_diffuse = 1 - np.exp(-k_diffuse * delta_t)

time = 0

def iterate(time):
    for i in range(rows):
        for j in range(cols):
            if lysis[i][j] == time:
                # Lyse cell
                lysis[i][j] == -1
                phage[i][j] += burst_size
                cells[i][j] == 0
            if replication[i][j] == time:
                # replicate cell randomly eventually
                if (i+1) < rows and cells[i+1][j] == 0:
                    cells[i+1][j] += 1
                    replication[i+1][j] = time + replication_time
                replication[i][j] = time + replication_time
            # Calculate probability of infection
            if phage[i][j] > 0 and cells[i][j] > 0 and lysis[i][j] == -1:
                p_infect = 1 - np.exp(-phage[i][j]*k_infection*(p_decay/k_degradation))
                infect = np.random.choice([0,1], p = [1 - p_infect, p_infect])
                if infect == 1:
                    lysis[i][j] = infect * lysis_time + time
            dead_phage = np.random.binomial(phage[i][j], p_decay)
            # print(dead_phage)
            phage[i][j] -= dead_phage
            to_diffuse = np.random.binomial(phage[i][j], p_decay)
            directions = np.random.multinomial(to_diffuse, [1/4]*4, size = 1)
            if i + 1 < rows and i - 1 >= 0:
                phage[i+1][j] += directions[0][0]
                phage[i-1][j] += directions[0][1]
                phage[i][j] -= (directions[0][0] + directions[0][1])
            if j + 1 < cols and (j - 1) >= 0:
                phage[i][j+1] += directions[0][2]
                phage[i][j-1] += directions[0][3]
                phage[i][j] -= (directions[0][2] + directions[0][3])


while time < 1000:
    iterate(time)
    time += delta_t
