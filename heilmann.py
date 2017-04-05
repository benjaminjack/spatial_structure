import numpy as np

delta_t = 1  # time step in minutes
k_diffuse = 2.5  # lambda
burst_size = 2  # Beta
# replication_time = 100  # T, in time steps
k_replicate = 0.01
lysis_time = 20  # in time steps
decay_time = 0  # dead cell decay time in times steps, set to 0 for no decay
# k_degradation = 10**-2  # delta
k_eps = 0.35/100  # rate that phage are lost to EPS
k_infection = 0.25/100  # alpha

rows = 20
cols = 20

p_eps = 0.35
p_infect = 0.25
p_debris = 0.25
p_replicate = 0.01
p_diffuse = 0.025

# Grid of phage
phage = np.random.choice([0, burst_size], p=[0.7, 0.3], size=(rows, cols))
print(phage)
# Grid of live cells
cells = np.random.choice([0, 1], p=[0.7, 0.3], size=(rows, cols))
# Grid of EPS
eps = np.random.choice([0, 1], p=[0.9, 0.1], size=(rows, cols))
# Grid of lysis timers (can also tell us where infected cells are)
lysis = np.ones((rows, cols), dtype='int')*-1
# Grid of replication timers
# replication = np.ones((rows, cols), dtype='int')*-1 + 101*cells
# Grid of debris/dead cell timers
debris = np.ones((rows, cols), dtype='int')*-1

# Probability of phage particle decaying during given time step
# p_decay = 1 - np.exp(-k_degradation * delta_t)
# Probability of phage particle diffusing to neighboring patch during given
# time step
# p_diffuse = 1 - np.exp(-k_diffuse * delta_t)

time = 0
burn_in = 7000

output = {
    'burst_total': 0,
    'replication_total': 0,
    'eps_total': np.sum(eps),
    'lost_to_debris_or_eps': 0,
    'lost_to_primary_infection': 0,
    'lost_to_secondary_infection': 0,
    'alphab': 0
}


def iterate(time):
    random_rows = list(range(0, rows))
    np.random.shuffle(random_rows)
    random_cols = list(range(0, cols))
    np.random.shuffle(random_cols)

    for i in random_rows:
        for j in random_cols:
            if cells[i][j] == 0 and lysis[i][j] > 0:
                print(cells, lysis)
                raise RuntimeError("Error")
            if debris[i][j] == time and decay_time > 0:
                # Reset timers for decaying dead cells
                debris[i][j] = -1
            if lysis[i][j] == time:
                # Check if it's time for cell to lyse
                lysis[i][j] = -1
                phage[i][j] += burst_size
                cells[i][j] = 0
                if decay_time > 0:
                    debris[i][j] = time + decay_time
                # replication[i][j] = -1
                if time > burn_in:
                    output['burst_total'] += 1
            if cells[i][j] > 0 and lysis[i][j] == -1:
                # p_replicate = 1 - np.exp(-cells[i][j] * k_replicate)
                replicate = np.random.binomial(cells[i][j], p_replicate)
                if replicate == 1:
                    # Replicate cells and diffuse randomly
                    possible_moves = []
                    # Construct a list of possible moves
                    if (i+1) < rows:
                        down = i + 1
                    else:
                        down = 0
                    if (j+1) < cols:
                        right = j + 1
                    else:
                        right = 0

                    if cells[down][j] == 0:
                        possible_moves.append((down, j))
                    if cells[i-1][j] == 0:
                        possible_moves.append((i-1, j))
                    if cells[i][right] == 0:
                        possible_moves.append((i, right))
                    if cells[i][j-1] == 0:
                        possible_moves.append((i, j-1))
                    if cells[i-1][j-1] == 0:
                        possible_moves.append((i-1, j-1))
                    if cells[down][right] == 0:
                        possible_moves.append((down, right))
                    if cells[down][i-1] == 0:
                        possible_moves.append((down, i-1))
                    if cells[j-1][right] == 0:
                        possible_moves.append((j-1, right))
                    if len(possible_moves) > 0:
                        # Randomly select from list of possible moves
                        move_index = np.random.choice(len(possible_moves))
                        move = possible_moves[move_index]
                        cells[move[0]][move[1]] += 1
                        # Reset replication timer for daughter cell
                        # replication[move[0]][move[1]] = time + replication_time
                    # Reset replication counter for original cell
                    # replication[i][j] = time + replication_time
                    if time > burn_in:
                        output['replication_total'] += 1
            # Calculate probability of infection and phage death
            if phage[i][j] > 0:
                # p_infect = 1 - np.exp(-cells[i][j] * k_infection)
                # p_eps = 1 - np.exp(-eps[i][j] * k_eps)
                # p_debris = 1 - np.exp(-(debris[i][j] > 0) * k_infection)
                p_infect_adj = cells[i][j] * p_infect
                p_eps_adj = eps[i][j] * p_eps
                p_debris_adj = (debris[i][j] > 0) * p_debris
                infect = np.random.multinomial(
                    phage[i][j],
                    [p_infect_adj,
                     p_eps_adj + p_debris_adj,
                     1 - p_infect_adj - p_eps_adj - p_debris_adj]
                )
                if infect[0] > 0:
                    phage[i][j] -= infect[0]  # Lose phage
                    # if time > burn_in:
                    #     output['lost_to_secondary_infection'] += (infect[0] - 1)
                    if lysis[i][j] == -1:
                        # Only primary infections reset lysis timer
                        lysis[i][j] = lysis_time + time
                        # Infected cells can't replicate
                        # replication[i][j] = -1
                        if time > burn_in:
                            output['lost_to_primary_infection'] += infect[0]
                    else:
                        if time > burn_in:
                            output['lost_to_secondary_infection'] += infect[0]
                elif infect[1] > 0:
                    # Lose phage to EPS or debris
                    phage[i][j] -= infect[1]
                    if time > burn_in:
                        output['lost_to_debris_or_eps'] += infect[1]

            # Randomly diffuse phage
            # First calculate how many phage will diffuse
            to_diffuse = np.random.binomial(phage[i][j], p_diffuse)
            # Next calculate which direction they will diffuse in
            directions = np.random.multinomial(to_diffuse, [1/4]*4)

            if (i+1) < rows:
                down = i + 1
            else:
                down = 0
            if (j+1) < cols:
                right = j + 1
            else:
                right = 0

            phage[down][j] += directions[0]
            phage[i][j] -= directions[0]

            phage[i-1][j] += directions[1]
            phage[i][j] -= directions[1]

            phage[i][right] += directions[2]
            phage[i][j] -= directions[2]

            phage[i][j-1] += directions[3]
            phage[i][j] -= directions[3]


# Iterate simulation for 1000 minutes (delta_t = 1 min)
while time < 10000:
    if time % 1000 == 0:
        print(time)
    iterate(time)
    time += delta_t
    if time > burn_in:
        output['alphab'] += (burst_size*(np.sum(cells)-np.sum(lysis > 0))*p_infect)/(p_infect*(np.sum(cells) + np.sum(debris > 0)) + p_eps*output["eps_total"])

output['phage_final'] = np.sum(phage)
output['total_cells_final'] = np.sum(cells)
output['alphab_avg'] = output['alphab']/3000



print(output)
