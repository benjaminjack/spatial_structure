import numpy as np
import sys

def iterate(time,
            delta_t,
            p_diffuse,
            burst_size,
            p_replicate,
            lysis_time,
            decay_time,
            p_eps,
            p_infect,
            p_debris,
            random_eps,
            rows,
            cols,
            fill_phage,
            fill_cells,
            fill_eps,
            burn_in,
            phage,
            cells,
            lysis,
            debris,
            random_rows,
            random_cols,
            eps,
            output
            ):
    # Shuffle rows and columns so they're not traversed in the same order on
    # each iteration
    np.random.shuffle(random_rows)
    np.random.shuffle(random_cols)

    temp_phage_with_eps = 0
    temp_infection_with_eps = 0
    temp_infection_total = 0

    for i in random_rows:
        for j in random_cols:
            if cells[i][j] == 0 and lysis[i][j] > 0:
                # Sanity check
                print(cells, lysis)
                raise RuntimeError("Cell counts and lysis timers are out of "
                                   " sync.")
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
                # Record lysis event
                if time > burn_in:
                    output['burst_total'] += 1
            # Replicate cells (only those not infected)
            if cells[i][j] > 0 and lysis[i][j] == -1:
                replicate = np.random.binomial(cells[i][j], p_replicate)
                if replicate == 1:
                    # Replicate cells and diffuse randomly
                    possible_moves = []
                    # Check to see if cell is at edge of grid, and have the
                    # coordinates wrap around such that phage and cells move
                    # on a torus.
                    if (i+1) < rows:
                        down = i + 1
                    else:
                        down = 0
                    if (j+1) < cols:
                        right = j + 1
                    else:
                        right = 0
                    # Daughter cell can move into one of 8 surrounding cells
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
                    # Record replication event
                    if time > burn_in:
                        output['replication_total'] += 1
            # Process infections
            if phage[i][j] > 0:
                # Adjust probabilites based on presence or absence of cells,
                # eps, or debris (dead cells)
                p_infect_adj = cells[i][j] * p_infect
                p_eps_adj = eps[i][j] * p_eps
                p_debris_adj = (debris[i][j] > 0) * p_debris
                # Randomly distribute phage
                infect = np.random.multinomial(
                    phage[i][j],
                    [p_infect_adj,
                     p_eps_adj,
                     p_debris_adj,
                     1 - p_infect_adj - p_eps_adj - p_debris_adj]
                )
                if infect[0] > 0:
                    # Infection will occur
                    phage[i][j] -= infect[0]  # Lose phage
                    if lysis[i][j] == -1:
                        # Only primary infections reset lysis timer
                        lysis[i][j] = lysis_time + time
                        # Record primary infection event
                        if time > burn_in:
                            output['lost_to_primary_infection'] += infect[0]
                            if eps[i][j] > 0:
                                temp_infection_with_eps += 1
                            temp_infection_total += 1
                    else:
                        # This cell is already infected
                        if time > burn_in:
                            output['lost_to_secondary_infection'] += infect[0]
                elif infect[1] > 0:
                    # Lose phage to EPS
                    phage[i][j] -= infect[1]
                    if time > burn_in:
                        output['lost_to_eps'] += infect[1]
                elif infect[2] > 0:
                    # Lose phage to debris
                    phage[i][j] -= infect[2]
                    if time > burn_in:
                        output['lost_to_debris'] += infect[2]

            # Randomly diffuse phage
            # First calculate how many phage will diffuse
            to_diffuse = np.random.binomial(phage[i][j], p_diffuse)
            # Next calculate which direction they will diffuse in
            directions = np.random.multinomial(to_diffuse, [1/4]*4)
            # Permute boundaries as in cell replication
            if (i+1) < rows:
                down = i + 1
            else:
                down = 0
            if (j+1) < cols:
                right = j + 1
            else:
                right = 0
            # Phage only move orthogonallly
            phage[down][j] += directions[0]
            phage[i][j] -= directions[0]

            phage[i-1][j] += directions[1]
            phage[i][j] -= directions[1]

            phage[i][right] += directions[2]
            phage[i][j] -= directions[2]

            phage[i][j-1] += directions[3]
            phage[i][j] -= directions[3]

            # Record some values
            if eps[i][j] > 0 and phage[i][j] > 0:
                temp_phage_with_eps += phage[i][j]


    # Calculate and record alpha-b values and other info
    if time > burn_in:
        healthy_cells = np.sum(cells) - np.sum(lysis > 0)
        total_cells = np.sum(cells) + np.sum(debris > 0)
        output['alphab'] += (burst_size * p_infect * healthy_cells) / \
            ((p_infect*total_cells) + p_eps * output["eps_total"])

        temp = (np.where(cells > 0, cells, 1) + np.where(eps > 0, eps, 1) + np.where(lysis < 0, lysis, 1))
        output['cells_with_eps'] += np.sum(temp == 3)/healthy_cells
        try:
            output['phage_with_eps'] += temp_phage_with_eps/np.sum(phage)
            output['infection_with_eps'] += temp_infection_with_eps/temp_infection_total
        except:
            pass


def run_simulation(eps,
                   burst,
                   decay=0,
                   random_eps=True,
                   sim_time=10000,
                   burn_in=7000):

    params = dict(
        delta_t = 1,  # time step in minutes
        p_diffuse = 0.025,  # probability of phage diffusing to an orthogonal patch
        burst_size = burst,  # Beta
        p_replicate = 0.01,  # probability of cell dividing
        lysis_time = 20,  # in time steps
        decay_time = decay,  # dead cell decay time in times steps, set to 0 for no decay
        p_eps = 0.35,  # probability that phage are lost to EPS
        p_infect = 0.25,  # probability of infection
        p_debris = 0.25,  # probability of infecting dead cells

        random_eps = random_eps,

        rows = 20,
        cols = 20,

        fill_phage = 0.3,  # proportion of phage at initialization
        fill_cells = 0.3,  # proportion of cells at initialization
        fill_eps = eps,  # proportion of eps at initialization

        burn_in = burn_in,  # how many iterations should we wait before recording data?
        # Randomly place phage in grid
    )

    patches = params['rows'] * params['cols']
    phage_count_init = int(params['fill_phage'] * patches)
    phage_empties_init = patches - phage_count_init
    cell_count_init = int(params['fill_cells'] * patches)
    cell_empties_init = patches - cell_count_init
    eps_count_init = int(params['fill_eps'] * patches)
    eps_empties_init = patches - eps_count_init

    params['phage'] = np.random.choice(
                        [params['burst_size']]*phage_count_init + \
                        [0]*phage_empties_init,
                        replace=False, size=(params['rows'], params['cols'])
                        )
    # Randomly place live cells
    params['cells'] = np.random.choice(
                        [1]*cell_count_init + [0]*cell_empties_init,
                        replace=False, size=(params['rows'], params['cols']))
    # Grid of lysis timers (can also tell us where infected cells are)
    params['lysis'] = np.ones((params['rows'], params['cols']), dtype='int')*-1
    # Grid of debris/dead cell timers
    params['debris'] = np.ones((params['rows'], params['cols']), dtype='int')*-1
    params['random_rows'] = list(range(0, params['rows']))
    params['random_cols'] = list(range(0, params['cols']))
    # Grid of EPS
    if not params['random_eps']:
        # if eps is deterministic, patches from left to right, top to bottom
        params['eps'] = np.zeros((params['rows'], params['cols']))
        for i in range(0, params['rows']):
            for j in range(0, params['cols']):
                if np.sum(params['eps']) < eps_count_init:
                    params['eps'][i][j] = 1
    else:
        params['eps'] = np.random.choice([1]*eps_count_init + \
                                         [0]*eps_empties_init,
                                         replace=False,
                                         size=(params['rows'], params['cols']))

    params['output'] = {
        'burst_total': 0,
        'replication_total': 0,
        'eps_total': np.sum(params['eps']),
        'lost_to_eps': 0,
        'lost_to_debris': 0,
        'lost_to_primary_infection': 0,
        'lost_to_secondary_infection': 0,
        'alphab': 0,
        'cells_with_eps': 0,
        'phage_with_eps': 0,
        'infection_with_eps': 0
    }

    time = 0

    while time < sim_time:
        iterate(time, **params)
        time += params['delta_t']

    try:
        params['output']['phage_final'] = np.sum(params['phage'])
        params['output']['total_cells_final'] = np.sum(params['cells'])
        params['output']['alphab_avg'] = params['output']['alphab'] / \
                                         (sim_time - burn_in)
        params['output']['cells_with_eps'] = params['output']['cells_with_eps'] / (sim_time - burn_in)
        params['output']['phage_with_eps'] = params['output']['phage_with_eps'] / (sim_time - burn_in)
        params['output']['infection_with_eps'] = params['output']['infection_with_eps'] / (sim_time - burn_in)
        params['output']['p2c'] = params['output']['lost_to_primary_infection'] / \
            params['output']['burst_total']
        params['output']['p2ic'] = params['output']['lost_to_secondary_infection']/\
            params['output']['burst_total']
        params['output']['p2eps'] = params['output']['lost_to_eps']/\
            params['output']['burst_total']
        params['output']['p2debris'] = params['output']['lost_to_debris']/\
            params['output']['burst_total']
    except:
        params['output']['phage_final'] = np.sum(params['phage'])
        params['output']['total_cells_final'] = np.sum(params['cells'])
        params['output']['alphab_avg'] = params['output']['alphab'] / \
                                         (sim_time - burn_in)
        params['output']['p2c'] = 0
        params['output']['p2ic'] = 0
        params['output']['p2eps'] = 0
        params['output']['p2debris'] = 0
    return params['output']


def main():
    burst_sizes = [2, 6, 10, 20, 40, 60]
    eps_list = [0.1, 0.3, 0.6, 0.9]
    print("EPS\tburst\talpha-b\tp->c\tp->ic\tp->eps\tC:E/C\tP:E/P\tI:E/I")
    # output = run_simulation(0.1, 2)
    for burst in burst_sizes:
        for eps in eps_list:
            for i in range(3):
                output = run_simulation(eps, burst)
                print("{:}\t{:}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}".format(eps, burst, output['alphab_avg'], output['p2c'], output['p2ic'], output['p2eps'], output['cells_with_eps'], output['phage_with_eps'], output['infection_with_eps']))
                sys.stdout.flush()


if __name__ == "__main__":
    main()
