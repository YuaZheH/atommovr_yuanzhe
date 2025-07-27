## Blue lab algorithm
import copy

from atommover.utils.core import *
from atommover.utils.move_utils import *
from atommover.algorithms.Algorithm_class import Algorithm
from atommover.algorithms.source.ejection import ejection
from atommover.algorithms.source.Hungarian_works import flatten_tuple, bfs_move_atom


def bluelab(init_config: np.ndarray, target_config: np.ndarray, do_ejection : bool = False, final_size: list = []):
    """
    Implements the balance and compress algorithm used in the Bernien blue lab.
    """
    bluelab_success_flag = False
    matrix = copy.deepcopy(init_config)
    move_list = []

    # Counts moves components in the blue lab algorithm
    balance_moves_term = 0
    compress_moves_term = 0
    ejection_moves_term = 0

    if len(final_size) == 0:
        final_size = [0, len(matrix) - 1, 0, len(matrix[0]) - 1]

    # 1. Balance
    balance_config, move_list = balance(init_config, target_config, final_size[0], final_size[1], move_list)
    balance_moves_term += len(move_list)

    # 2. Compress
    compress_config, move_list = one_d_compress(balance_config, target_config, move_list)
    compress_moves_term += len(move_list) - balance_moves_term
    final_config = copy.deepcopy(compress_config)

    # 3 Eject to certain geoemetry
    if do_ejection:
        eject_moves, final_config = ejection(compress_config, target_config, [0, len(matrix) - 1, 0, len(matrix[0]) - 1])
        move_list.extend(eject_moves)
        ejection_moves_term = len(eject_moves)
        # 3.1 Check if the configuration is the same as the target configuration
        # if np.array_equal(final_config, target_config):
        #     bluelab_success_flag = True
    else:
        ejection_moves_term = 0
        
        # 3.2 Check if the configuration (inside range of target) the same as the target configuration
        # effective_config = np.multiply(final_config, target_config)
        # if np.array_equal(effective_config, target_config):
        #     bluelab_success_flag = True
    bluelab_success_flag = Algorithm.get_success_flag(final_config,target_config, do_ejection=do_ejection, n_species=1)    

    return final_config, move_list, bluelab_success_flag#, [balance_moves_term, compress_moves_term, ejection_moves_term]

def one_d_compress(matrix, target_config, move_list):
    
    start_col_index = np.min(np.where(target_config == 1)[1])
    end_col_index = np.max(np.where(target_config == 1)[1])

    for row_ind in range(len(matrix)):
        source_atom_pos = []
        target_atom_pos = []
        corresponding_col_index = start_col_index
        max_distance = 0

        for col_ind in range(len(matrix[0])):
            if matrix[row_ind][col_ind] == 1:
                source_atom_pos.append([row_ind, col_ind])
                if corresponding_col_index <= end_col_index:
                    target_atom_pos.append([row_ind, corresponding_col_index])
                else:
                    target_atom_pos.append([row_ind, col_ind])

                if abs(corresponding_col_index - col_ind) > max_distance:
                    max_distance = abs(corresponding_col_index - col_ind)

                if corresponding_col_index <= end_col_index:
                    corresponding_col_index += 1
                else:
                    break

        
        i = 0
        while i <= max_distance:
            moves_in_scan = []
            moved_index = []
            temp_filled_ind = []
            temp_empty_ind = []

            for index in range(len(source_atom_pos)):
                if source_atom_pos[index] != target_atom_pos[index]:
                    move_direction = 1 if source_atom_pos[index][1] < target_atom_pos[index][1] else -1
                    next_pos = source_atom_pos[index][1] + move_direction
                    # Check if the target position is empty before moving
                    if (matrix[row_ind][next_pos] == 0 or next_pos in temp_empty_ind) and next_pos not in temp_filled_ind:
                        move = Move(row_ind, source_atom_pos[index][1], row_ind, next_pos)
                        moved_index.append(index)
                        temp_filled_ind.append(next_pos)
                        temp_empty_ind.append(source_atom_pos[index][1])
                        source_atom_pos[index][1] = next_pos  # Update source position towards target
                        moves_in_scan.append(move)
            
            for index in range(len(source_atom_pos))[::-1]:
                if source_atom_pos[index] != target_atom_pos[index] and index not in moved_index:
                    move_direction = 1 if source_atom_pos[index][1] < target_atom_pos[index][1] else -1
                    next_pos = source_atom_pos[index][1] + move_direction
                    # Check if the target position is empty before moving
                    if (matrix[row_ind][next_pos] == 0 or next_pos in temp_empty_ind) and next_pos not in temp_filled_ind:
                        move = Move(row_ind, source_atom_pos[index][1], row_ind, next_pos)
                        temp_filled_ind.append(next_pos)
                        temp_empty_ind.append(source_atom_pos[index][1])
                        source_atom_pos[index][1] = next_pos  # Update source position towards target
                        moves_in_scan.append(move)
                    elif matrix[row_ind][next_pos] == 1:
                        source_atom_pos[index][1] = next_pos
        
                        
            if len(moves_in_scan) > 0:
                move_list.append(moves_in_scan)
                matrix, _ = move_atoms(matrix, moves_in_scan)
            i += 1

    return matrix, move_list

"""
# Save version of individual blue lab algo
def one_d_compress(matrix, move_list):

    # Identify total atom numbers and required number in a row
    for row_ind in range(len(matrix[0])):
        source_atom_pos = []
        target_atom_pos = []
        corresponding_col_index = 0 ##If you want to contract atoms to middle range, you can add a row_min here
        max_distance = 0
        
        
        for col_ind in range(len(matrix)):
            if matrix[row_ind][col_ind] == 1:
                source_atom_pos.append([row_ind,col_ind])
                target_atom_pos.append([row_ind,corresponding_col_index])
                if abs(corresponding_col_index - col_ind) > max_distance:
                    max_distance = abs(corresponding_col_index - col_ind)
                corresponding_col_index += 1

        i = 0
        while i in range(max_distance):
            moves_in_scan = []
            for index in range(len(source_atom_pos)):
                if source_atom_pos[index] != target_atom_pos[index]:
                    move = Move(row_ind, source_atom_pos[index][1],  row_ind,source_atom_pos[index][1]-1)
                    source_atom_pos[index][1] -= 1
                    moves_in_scan.append(move)
            if len(moves_in_scan) > 0:
                move_list.append(moves_in_scan)
                matrix, _ = move_atoms(matrix, moves_in_scan)
            i+=1

    return matrix, move_list
"""

def balance(init_config, target_config, row_min, row_max, move_list):
    matrix = copy.deepcopy(init_config)
    balance_config = copy.deepcopy(matrix)
    col_min = 0
    col_max = len(matrix[0]) - 1
    # Calculate the number of rows in the submatrix
    row_nums = row_max - row_min + 1
    
    # If there is only one row, return
    if row_nums == 1:
        return matrix, move_list
    
    # Calculate the middle row index
    middle_row = row_min + (row_nums // 2) - 1
    
    # Calculate the total number of 1s in the submatrix S[i:j+1]
    #N_total = sum(sum(row) for row in matrix[row_min:row_max + 1])

    # Calculate the required number of 1s in the submatrix S[i:m+1]
    n_req = int(np.sum(target_config[row_min:middle_row + 1, col_min:col_max + 1]))
    n_req_c = int(np.sum(target_config[middle_row + 1:row_max + 1, col_min:col_max + 1]))
    #print(f"n_req: {n_req}")

    # Calculate the actual number of 1s in the submatrix S[i:m+1]
    N_S_i_m = int(np.sum(matrix[row_min:middle_row + 1, col_min:col_max + 1]))
    N_S_i_m_c = int(np.sum(matrix[middle_row + 1:row_max + 1, col_min:col_max + 1]))
    
    if N_S_i_m > n_req and N_S_i_m_c < n_req_c:
        # Shift excess 1s down from S[i:m+1] to S[m+1:j+1]
        excess_atoms = N_S_i_m - n_req
        balance_config, move_list = down_move(matrix, target_config, excess_atoms, row_min, middle_row, row_max, col_min, col_max, move_list)
        #return matrix, move_list

    elif N_S_i_m < n_req and N_S_i_m_c > n_req_c:
        # Shift required 1s up from S[m+1:j+1] to S[i:m+1]
        needed_atoms = n_req - N_S_i_m
        balance_config, move_list = up_move(matrix, target_config, needed_atoms, row_min, middle_row, row_max, col_min, col_max, move_list)
        #return matrix, move_list

    # Recursively balance the submatrices
    balance_config, move_list = balance(balance_config, target_config, row_min, middle_row, move_list)
    balance_config, move_list = balance(balance_config, target_config, middle_row + 1, row_max, move_list)

    return balance_config, move_list

def down_move(matrix, target_config, excess_atoms, row_min, middle_row, row_max, col_min, col_max, move_list):
    #Initialize the move bound
    source_row = middle_row
    target_row = middle_row + 1
    normalize_row = 0
    balance_row_count = target_row - source_row
    stuff = 0
    brute_force_start = []
    brute_force_end = []

    #while excess_atoms >  0:
    while row_max >= source_row >= row_min and excess_atoms > 0:
        moves_in_scan = []

        for col in range(col_min,col_max + 1):
            if matrix[source_row + normalize_row, col] == 1 and matrix[source_row + normalize_row + 1, col] == 0 and matrix[target_row, col] == 0 and excess_atoms > 0:
                move = Move(source_row + normalize_row, col, source_row + normalize_row + 1, col)
                #print(f'{move.from_row}, {move.from_col} -> {move.to_row}, {move.to_col}')
                moves_in_scan.append(move)
                
                if balance_row_count == 1:
                    excess_atoms -= 1

        if len(moves_in_scan) > 0:
            move_list.append(moves_in_scan)
            matrix, _ = move_atoms(matrix, moves_in_scan)
    
        #TODO: Solve the stuff case later
        if sum(matrix[target_row,col_min:col_max+1]) == len(matrix[target_row,col_min:col_max+1]) and target_row + stuff < len(matrix) and excess_atoms > 0:
            for shift in range(stuff, -1, -1):
                moves_in_scan = []
                for col in range(col_min,col_max+1):
                    if matrix[target_row + shift, col] == 1 and matrix[target_row + shift + 1, col] == 0:
                        move = Move(target_row + shift, col, target_row + shift + 1, col)
                        #print(f'{move.from_row}, {move.from_col} -> {move.to_row}, {move.to_col}')
                        moves_in_scan.append(move)

                if len(moves_in_scan) > 0:
                    move_list.append(moves_in_scan)
                    matrix, _ = move_atoms(matrix, moves_in_scan)
                
            stuff += 1
            source_row = middle_row
            balance_row_count = target_row - source_row
            normalize_row = 0
        else:
            balance_row_count -= 1
            normalize_row += 1

        if balance_row_count == 0:
            source_row -= 1
            balance_row_count = target_row - source_row
            normalize_row = 0

            ## Slide up move to achieve balance
            # moves_in_scan = []
            # for col in range(col_min, col_max + 1):
            #     if matrix[source_row, col] == 1 and matrix[source_row + 1, col - 1] == 0 and matrix[target_row, col] == 0 and excess_atoms > 0:
    if excess_atoms > 0:
        #Find unbalance atoms
        for col_ind in range(col_min, col_max + 1):
            for row_ind in range(row_min, row_max + 1):
                if matrix[row_ind][col_ind] == 1 and target_config[row_ind][col_ind] == 0:
                    brute_force_start.append((row_ind,col_ind))
                if matrix[row_ind][col_ind] == 0 and target_config[row_ind][col_ind] == 1:
                    brute_force_end.append((row_ind,col_ind))
        
        for ind in range(len(brute_force_start)):
            if ind < len(brute_force_end):
                path = flatten_tuple(bfs_move_atom(matrix, brute_force_start[ind], brute_force_end[ind]))
                path = path[::-1]

            #Iterate all path segments (((a1,b1), (a2,b2)), ((c1,c2), (d1,d2)))
                for item in path:
                    current_pos = item[0] 
                    for next_pos in item:
                        if np.array_equal(current_pos, next_pos):
                            pass
                        else:
                            if matrix[current_pos[0]][current_pos[1]] == 1 and matrix[next_pos[0]][next_pos[1]] == 0:
                                matrix[current_pos[0]][current_pos[1]]=0
                                matrix[next_pos[0]][next_pos[1]]=1
                                move_list.append([Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1])])
                            current_pos = next_pos
    

    return matrix, move_list


def up_move(matrix, target_config, needed_atoms, row_min, middle_row, row_max, col_min, col_max, move_list):
    # Initialize the move bound
    source_row = middle_row + 1
    target_row = middle_row 
    normalize_row = 0
    balance_row_count = source_row - target_row
    stuff = 0
    brute_force_start = []
    brute_force_end = []

    # while needed_atoms > 0:
    while row_max >= source_row >= row_min and needed_atoms > 0:
        moves_in_scan = []

        for col in range(col_min,col_max+1):
            if matrix[source_row + normalize_row, col] == 1 and matrix[source_row + normalize_row - 1, col] == 0 and matrix[target_row, col] == 0 and needed_atoms > 0:
                move = Move(source_row + normalize_row, col, source_row + normalize_row - 1, col)
                #print(f'{move.from_row}, {move.from_col} -> {move.to_row}, {move.to_col}')
                moves_in_scan.append(move)
                
                if balance_row_count == 1:
                    needed_atoms -= 1

        if len(moves_in_scan) > 0:
            move_list.append(moves_in_scan)
            matrix, _ = move_atoms(matrix, moves_in_scan)
    
        if sum(matrix[target_row,col_min:col_max+1]) == len(matrix[target_row,col_min:col_max+1]) and target_row > 0 and needed_atoms > 0:
            for shift in range(stuff, -1, -1):
                moves_in_scan = []
                for col in range(col_min,col_max+1):
                    if matrix[target_row - shift, col] == 1 and matrix[target_row -shift - 1, col] == 0:
                        move = Move(target_row - shift, col, target_row - shift - 1, col)
                        #print(f'{move.from_row}, {move.from_col} -> {move.to_row}, {move.to_col}')
                        moves_in_scan.append(move)

                if len(moves_in_scan) > 0:
                    move_list.append(moves_in_scan)
                    matrix, _ = move_atoms(matrix, moves_in_scan)

            stuff += 1
            source_row = middle_row + 1
            balance_row_count = source_row - target_row
            normalize_row = 0
        else:
            balance_row_count -= 1
            normalize_row -= 1

        if balance_row_count == 0:
            source_row += 1
            balance_row_count = source_row - target_row
            normalize_row = 0
            
        #if source_row > row_max and needed_atoms > 0:
            ## Slide up move to achieve balance
            #print(f"Up move, Lattice1:{row_min}, {middle_row}, {col_min}, {col_max}")
            #print(f"Up move, Lattice2:{middle_row+1}, {row_max}, {col_min}, {col_max}")
    if needed_atoms > 0:
        #Find unbalance atoms
        for col_ind in range(col_min, col_max + 1):
            for row_ind in range(row_min, row_max + 1):
                if matrix[row_ind][col_ind] == 1 and target_config[row_ind][col_ind] == 0:
                    brute_force_start.append((row_ind,col_ind))
                if matrix[row_ind][col_ind] == 0 and target_config[row_ind][col_ind] == 1:
                    brute_force_end.append((row_ind,col_ind))
        for ind in range(len(brute_force_start)):
            if ind < len(brute_force_end):
                path = flatten_tuple(bfs_move_atom(matrix, brute_force_start[ind], brute_force_end[ind]))
                path = path[::-1]

            #Iterate all path segments (((a1,b1), (a2,b2)), ((c1,c2), (d1,d2)))
            for item in path:
                current_pos = item[0] 
                for next_pos in item:
                    if np.array_equal(current_pos, next_pos):
                        pass
                    else:
                        if matrix[current_pos[0]][current_pos[1]] == 1 and matrix[next_pos[0]][next_pos[1]] == 0:
                            matrix[current_pos[0]][current_pos[1]]=0
                            matrix[next_pos[0]][next_pos[1]]=1
                            move_list.append([Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1])])
                        current_pos = next_pos
                
    return matrix, move_list