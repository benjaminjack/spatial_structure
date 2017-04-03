import numpy as np

delta_t = 1  # time step in minutes
k_diffuse = 10  # lambda
burst_size = 80  # Beta
replication_time = 300  # T, in time steps
lysis_time = 10  # in time steps
decay_time = 10  # dead cell decay time in times steps, set to 0 for no decay
# k_degradation = 10**-2  # delta
k_eps = 10**-2  # rate that phage are lost to EPS
k_infection = 10**-3  # alpha

rows = 10
cols = 10

# Grid of phage
phage = np.random.choice([0, 100], p=[0.9, 0.1], size=(rows, cols))
# Grid of live cells
cells = np.random.choice([0, 1], p=[0.9, 0.1], size=(rows, cols))
# Grid of EPS
eps = np.random.choice([0, 1], p=[0.9, 0.1], size=(rows, cols))
# Grid of lysis timers (can also tell us where infected cells are)
lysis = np.ones((rows, cols), dtype='int')*-1
# Grid of replication timers
replication = np.ones((rows, cols), dtype='int')*-1 + 101*cells
# Grid of debris/dead cell timers
debris = np.ones((rows, cols), dtype='int')*-1

# Probability of phage particle decaying during given time step
# p_decay = 1 - np.exp(-k_degradation * delta_t)
# Probability of phage particle diffusing to neighboring patch during given
# time step
p_diffuse = 1 - np.exp(-k_diffuse * delta_t)

time = 0


def iterate(time):
    random_rows = list(range(0, rows))
    np.random.shuffle(random_rows)
    random_cols = list(range(0, cols))
    np.random.shuffle(random_cols)

    for i in random_rows:
        for j in random_cols:
            if debris[i][j] == time and decay_time > 0:
                # Reset timers for decaying dead cells
                debris[i][j] = -1
            if lysis[i][j] == time:
                # Check if it's time for cell to lyse
                lysis[i][j] = -1
                phage[i][j] += burst_size
                cells[i][j] = 0
                debris[i][j] = time + decay_time
                replication[i][j] = -1
            if replication[i][j] == time:
                # Replicate cells and diffuse randomly
                possible_moves = []
                # Construct a list of possible moves
                if (i+1) < rows and cells[i+1][j] == 0:
                    possible_moves.append((i+1, j))
                if (i-1) >= 0 and cells[i-1][j] == 0:
                    possible_moves.append((i-1, j))
                if (j+1) < cols and cells[i][j+1] == 0:
                    possible_moves.append((i, j+1))
                if (j-1) >= 0 and cells[i][j-1] == 0:
                    possible_moves.append((i, j-1))
                if len(possible_moves) > 0:
                    # Randomly select from list of possible moves
                    move_index = np.random.choice(len(possible_moves))
                    move = possible_moves[move_index]
                    cells[move[0]][move[1]] += 1
                    # Reset replication timer for daughter cell
                    replication[move[0]][move[1]] = time + replication_time
                # Reset replication counter for original cell
                replication[i][j] = time + replication_time
            # Calculate probability of infection and phage death
            if phage[i][j] > 0:
                # p_infect = 1 - np.exp(-phage[i][j] * k_infection *
                #                      (p_decay/k_degradation))
                p_infect = 1 - np.exp(-cells[i][j] * k_infection)
                p_eps = 1 - np.exp(-eps[i][j] * k_eps)
                p_debris = 1 - np.exp(-(debris[i][j] > 0) * k_infection)
                infect = np.random.multinomial(
                    phage[i][j],
                    [p_infect,
                     p_eps + p_debris,
                     1 - p_infect - p_eps - p_debris]
                )
                if infect[0] > 0:
                    phage[i][j] -= infect[0]  # Lose phage
                    if lysis[i][j] == -1:
                        # Only primary infections reset lysis timer
                        lysis[i][j] = lysis_time + time
                        # Infected cells can't replicate
                        replication[i][j] = -1
                elif infect[1] > 0:
                    # Lose phage to EPS or debris
                    phage[i][j] -= infect[1]
                if cells[i][j] > 0:
                    print(p_infect, p_eps, p_debris, infect)

            # Randomly diffuse phage
            # First calculate how many phage will diffuse
            to_diffuse = np.random.binomial(phage[i][j], p_diffuse)
            # Next calculate which direction they will diffuse in
            directions = np.random.multinomial(to_diffuse, [1/4]*4)
            if i + 1 < rows:
                phage[i+1][j] += directions[0]
                phage[i][j] -= directions[0]
            if i - 1 >= 0:
                phage[i-1][j] += directions[1]
                phage[i][j] -= directions[1]
            if j + 1 < cols:
                phage[i][j+1] += directions[2]
                phage[i][j] -= directions[2]
            if (j - 1) >= 0:
                phage[i][j-1] += directions[3]
                phage[i][j] -= directions[3]


# Iterate simulation for 1000 minutes (delta_t = 1 min)
while time < 1000:
    iterate(time)
    time += delta_t
