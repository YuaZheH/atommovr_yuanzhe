import numpy as np
import copy
import time
from collections import deque
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from atommover.utils.move_utils import Move
from atommover.utils.AtomArray import AtomArray
from atommover.algorithms.Algorithm_class import Algorithm
def generate_cost_matrix(current_positions, target_positions):
    num_atoms = len(current_positions)
    num_targets = len(target_positions)
    cost_matrix = np.zeros((num_atoms, num_targets))

    for i, current in enumerate(current_positions):
        for j, target in enumerate(target_positions):
            cost_matrix[i, j] = np.sqrt((current[0] - target[0])**2 + (current[1] - target[1])**2)
    return cost_matrix
def move_atom_and_show_grid(grid, start, end):
    #Initialize current position
    current_pos = start
    path = []
    
    while current_pos != end:
        path, current_pos = bfs_move_atom(grid, current_pos, end, path)
        
    path = flatten_tuple(path)[::-1]
    grid, path = generate_decomposed_move_set(grid, path)

    return path
def bfs_move_atom(grid, start, end, prev_path):
    queue = deque([(start[0], start[1], [(start[0], start[1])])]) #Use the queue to record current position and path
    visited = set() #Record the visited positions
    visited.add((start[0], start[1]))

    #Start finding the path
    while queue:
        current_row, current_col, path = queue.popleft() #Update current position

        #If we arrive end point, return the path
        if (current_row, current_col) == end:
            if len(prev_path) > 0:
                prev_path = prev_path, path
            else:
                prev_path = path

            return prev_path, end

        #Explore the next step (based on current position and end point)
        len_path = len(path) - 1
        dr = 1 if end[0] > path[len_path][0] else (-1 if end[0] < path[len_path][0] else 0) 
        dc = 1 if end[1] > path[len_path][1] else (-1 if end[1] < path[len_path][1] else 0)
        new_row, new_col = current_row + dr, current_col + dc

        #Check if there is an obstacle there (If no, start from this new point to find next step)
        if (new_row, new_col) not in visited and grid[new_row][new_col] == 0:
            visited.add((new_row, new_col))
            queue.append((new_row, new_col, path + [(new_row, new_col)])) 

    #If there is an obstacle on the path, we decompose the path: start->obstacle->target
    #Define the obstacle position
    obstacle = (path[len_path][0] + dr, path[len_path][1] + dc)

    # Update the move in path until obstacle
    if len(prev_path) > 0:
        prev_path = prev_path, path + [obstacle]
    else:
        prev_path = path + [obstacle]

    #[Path between start and obstacle] + [Path between obstacle and end]
    return prev_path, obstacle
def flatten_tuple(nested_tuple):
    # This function will flatten a nested tuple of lists into a single tuple of lists
    result = []
    
    def recursive_flatten(element):
        if isinstance(element, tuple):
            # If the element is a tuple, apply recursion to each item
            for item in element:
                recursive_flatten(item)
        elif isinstance(element, list):
            # If the element is a list, append it to the result
            result.append(tuple(element))

    # Start the recursion with the entire nested tuple
    recursive_flatten(nested_tuple)
    
    # Convert the list of tuples into a single tuple
    return tuple(result)

    left_eject, bot_eject, top_eject, right_eject = [], [], [], []
    len_x = len(matrix[0])
    len_y = len(matrix)
    for x in range(len_x):
        for y in range(len_y):
            if matrix[x][y] == 1 and target_config[x][y] == 0:
                if x >= y and x < len_y - y:
                    left_eject.append((x, y))
                elif x >= y and x >= len_y - y:
                    bot_eject.append((x, y))
                elif x < y and x <= len_y - y:
                    top_eject.append((x, y))
                elif x <= y and x >= len_y - y:
                    right_eject.append((x, y))
    return left_eject, bot_eject, top_eject, right_eject
def generate_decomposed_move_set(grid, path):
    decomposed_move_set = []

    #Iterate all path segments (((a1,b1), (a2,b2), (a3, b3), (a4,b4)), ((c1,d1), (c2,d2)))
    try:
        for segmant in path:
            from_row, from_col = segmant[0]
            for coordinate in segmant:
                to_row, to_col = coordinate
                # To exclude the frist move
                if from_row != to_row or from_col != to_col:
                    decomposed_move_set.append([Move(from_row, from_col, to_row, to_col)])
                    grid[from_row][from_col] = 0
                    grid[to_row][to_col] = 1
                    from_row, from_col = to_row, to_col
            # decomposed_move_set.append(segmant_moves)
    except IndexError:
        return grid, []
        
    return grid, decomposed_move_set
def define_current_and_target(matrix, target_config):
    current_positions = [(x, y) for x in range(len(matrix[0])) for y in range(len(matrix)) if matrix[x][y] == 1 if target_config[x][y] == 0] #NKH this should in theory not change anything...
    target_positions = [(x, y) for x in range(len(matrix[0])) for y in range(len(matrix)) if target_config[x][y] == 1 if matrix[x][y] == 0] #same here
    return current_positions, target_positions
def Hungarian_algorithm_original(atom_arrays: np.ndarray, target_config: np.ndarray, do_ejection: bool = False, final_size: list = []):
    move_set = []
    matrix = copy.deepcopy(atom_arrays)

    if len(final_size) == 0:
        final_size = [0, len(matrix[0])-1, 0, len(matrix)-1]

    time_1 = time.time()
    current_positions, target_positions = define_current_and_target(matrix, target_config)    #Define target positions for the center square in a matrix.
    time_2 = time.time()
    print(f"Defining current and target positions took {time_2 - time_1} seconds.")   
    cost_matrix = generate_cost_matrix(current_positions, target_positions)    # Generate the cost matrix using the current atom positions and the target positions
    time_3 = time.time()
    print(f"Generating cost matrix took {time_3 - time_2} seconds.")
    row_ind, col_ind = linear_sum_assignment(cost_matrix)    # row_ind and col_ind are arrays of indices indicating the optimal assignment
    time_4 = time.time()
    print(f"Solving linear sum assignment took {time_4 - time_3} seconds.")                     
    paired_indices = sorted(zip(row_ind, col_ind), key=lambda x: x[1])    # Pair up row_ind and col_ind and sort by col_ind

    if paired_indices:
        # Unzip the sorted pairs if paired_indices is not empty
        sorted_row_ind, sorted_col_ind = zip(*paired_indices)
    else:
        # Assign default values if paired_indices is empty
        sorted_row_ind, sorted_col_ind = [], []
    prepared_assignments = [(current_positions[i], target_positions[j]) for i, j in zip(sorted_row_ind, sorted_col_ind)]

    for start, target in prepared_assignments:
        Hungarian_move = []
        Hungarian_move = move_atom_and_show_grid(matrix, start, target)
        move_set.extend(Hungarian_move)

    #Optional ejection argument
    if do_ejection:
        eject_moves, eject_config = ejection(matrix, target_config, final_size)
        move_set.extend(eject_moves)
    else:
        eject_config = copy.deepcopy(matrix)

    success_flag = Algorithm.get_success_flag(eject_config.reshape(np.shape(target_config)), target_config, do_ejection=do_ejection, n_species = 1)
    time_5 = time.time()    
    print(f"returning the results took {time_5 - time_4} seconds.")     
    print(f"Total time for Hungarian algorithm: {time_5 - time_1} seconds.")          
    return eject_config, move_set, success_flag

import time
from atommover.utils import PhysicalParams

 

# Start timing
start_time = time.time()

# specifying parameters
load_prob = 0.5 # float; the probability that an individual site will be loaded
array_length = 30 # int; number of rows (or cols) of the square array

params = PhysicalParams(loading_prob = load_prob)

# Time array initialization
init_start = time.time()
atom_array_1 = AtomArray([array_length, array_length], 1, params)
init_time = time.time() - init_start
print(f"Array initialization time: {init_time:.4f} seconds")

# Time tweezer loading
load_start = time.time()
atom_array_1.load_tweezers()
load_time = time.time() - load_start
print(f"Tweezer loading time: {load_time:.4f} seconds")

# Time target generation
target_start = time.time()
atom_array_1.generate_target()
target_time = time.time() - target_start
print(f"Target generation time: {target_time:.4f} seconds")

# Time visualization
viz_start = time.time()
# atom_array_1.image()
# atom_array_1.plot_target_config()
viz_time = time.time() - viz_start
print(f"Visualization time: {viz_time:.4f} seconds")


# Time move calculation
moves_start = time.time()
move_sequence = Hungarian_algorithm_original(atom_array_1.matrix[:,:,0], atom_array_1.target)[1]
moves_time = time.time() - moves_start
print(f"Move sequence calculation time: {moves_time:.4f} seconds")

# Total time
total_time = time.time() - start_time
print(f"\nTotal execution time: {total_time:.4f} seconds")
print(f"Number of moves: {len(move_sequence) if move_sequence else 0}")
print(move_sequence)
print(move_sequence[1][0].from_row)