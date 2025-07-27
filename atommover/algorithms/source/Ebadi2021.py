# This code implements the parallel rearrangement protocol described in [Nature 595, pp. 272-232](https://www.nature.com/articles/s41586-021-03582-4/figures/8).

from atommover.utils.core import *
from atommover.utils.move_utils import *
from atommover.algorithms.source.ejection import ejection
from atommover.algorithms.Algorithm_class import Algorithm

import copy


def parallel_row_col_rearrangement(init_config: np.ndarray, target_config: np.ndarray, max_rearrangement_cycles: int = 50, max_presorting_cycles: int = 5, do_ejection: bool = False) -> tuple[np.ndarray, list, bool]:
    """
        Generates a list of moves that implements the parallel rearrangement protocol described in [Ebadi et al. 2021](https://www.nature.com/articles/s41586-021-03582-4/figures/8).
    """
    move_list = []

    # 1. identify deficient columns and perform horizontal sorting (presorting)
    presorting_count = 0
    horizontal_config = copy.deepcopy(init_config)

    # Counts moves components in the parallel rearrangment
    presorting_moves_term = 0
    ejection_moves_term = 0
    colsort_moves_term = 0

    # Presorting
    while presorting_count < max_presorting_cycles:
        # counting atoms in each column
        column_counts = count_atoms_in_columns(horizontal_config)

        # finding columns with insufficient atoms
        deficient_cols = suff_col_list(target_config,column_counts)

        # calculating moves for presort
        horizontal_config, presorting_moves = horizontal_sorting(horizontal_config, target_config, deficient_cols, column_counts)
        presorting_moves_term += len(presorting_moves)
        move_list.extend(presorting_moves)

        # checking to see if the moves satisfy the column requirements
        horizontal_flag = horizontal_check(horizontal_config, target_config)
        if not horizontal_flag:
            presorting_count += 1
        else:
            break
    
    # If presorting fail, stop running the algorithm
    if not horizontal_flag:
        print('presorting fail, not enough atoms in all the columns.')
        return horizontal_config, move_list, horizontal_flag #, [0, 0, 0]
    else:
        pass
    
    # 2. Ejection of excess atoms
    if do_ejection:
        eject_config, eject_moves = eject_excess_atoms(horizontal_config, target_config, eject_all = True)
        ejection_moves_term += len(eject_moves)
        move_list.extend(eject_moves)
    else:
        # eject_config = copy.deepcopy(horizontal_config)
        eject_config, eject_moves = eject_excess_atoms(horizontal_config, target_config, eject_all = False)
        ejection_moves_term += len(eject_moves)
        move_list.extend(eject_moves)

    # 3. Parallel sorting within columns
    final_config, colsort_moves, _ = parallel_vertical_sort(eject_config, target_config, max_cycles = max_rearrangement_cycles)
    colsort_moves_term += len(colsort_moves)
    move_list.extend(colsort_moves)

    success_flag = Algorithm.get_success_flag(final_config, target_config, do_ejection=do_ejection, n_species = 1)
    
    # Statistics moves components (presorting, ejection, colsort) in the parallel rearrangment
    # all_moves_component = [presorting_moves_term, ejection_moves_term, colsort_moves_term]
    return final_config, move_list, success_flag #, all_moves_component


def suff_col_list(target_config, column_counts) -> list:
    sufficient_cols = []
    deficient_cols = []
    for col, count in enumerate(column_counts, start=0):
        if count >= np.sum(target_config[:,col]):
            sufficient_cols.append(col)
        else:
            deficient_cols.append(col)
    return deficient_cols


def horizontal_sorting(init_matrix, target_config, deficient_cols, column_counts):
    matrix = copy.deepcopy(init_matrix)
    horiz_moves_list = [] # a chronological list of moves or moves in parallel
    for col in deficient_cols:
        #left boundary condition
        if col == 0:
            surplus_left = 0
        else:
            surplus_left = column_counts[col - 1] - np.sum(target_config[:,col-1]) #NKH

        #right boundary condition
        if col == len(matrix[0]) - 1:
            surplus_right = 0
        else:
            surplus_right = column_counts[col + 1] - np.sum(target_config[:,col+1]) #NKH

        #find source column
        if surplus_left > surplus_right and surplus_left > 0:
            source_col = -1 
        elif surplus_left < surplus_right and surplus_right > 0:
            source_col = 1 
        elif surplus_left == surplus_right and surplus_left > 0:
            source_col = -1
        else:
            source_col = 0
            # print(f"No neighbor of column {col} has surplus atoms.")

        #find available rows
        empty_rows = [row for row in range(len(matrix)) if matrix[row][col] == 0]
        available_rows = [row for row in empty_rows if matrix[row][col + source_col] == 1]

        #perform horizontal move
        if source_col != 0:
            if np.sum(target_config[:,col+source_col]) == 0:
                unlimited_presorting = True
            else:
                unlimited_presorting = False
            matrix, moves_list = move_col(matrix, col, empty_rows, available_rows, np.sum(target_config[:,col]) - column_counts[col], surplus_left, surplus_right, source_col, unlimited_presorting) #NKH
            if len(moves_list) > 0:
                matrix, _ = move_atoms(matrix, moves_list)
                horiz_moves_list.append(moves_list)
    return matrix, horiz_moves_list

def move_col(matrix, col, empty_rows, available_rows, needed_atoms, surplus_left, surplus_right, col_shift, unlimited_presorting = False):
    """This code implements a set of horizontal moves to one column from a column on its left (`col_shift = 1`) or right (`col_shift = -1`)"""
    moves = []
    for row in available_rows:
        # move one atom from neighbor column
        move = Move(row, col+col_shift, row, col)
        # move = ((row,col+col_shift),(row,col))
        moves.append(move)
        
        needed_atoms -= 1
        if col_shift == -1:
            surplus_left -= 1
        else:
            surplus_right -= 1

        # stop if the number of atoms is sufficient
        if needed_atoms == 0 and unlimited_presorting == False:
            break

    return matrix, moves

def horizontal_check(matrix, target_config):
    column_counts = count_atoms_in_columns(matrix) # Count the number of atoms in each column
    all_columns_sufficient = all(count >= np.sum(target_config[:,col]) for col, count in enumerate(column_counts)) #NKH # Check if every column has at least 5 atoms

    return all_columns_sufficient

def eject_excess_atoms(init_matrix, target_config, eject_all = False):
    """ This code ejects excess atoms by scanning downward in
        the array (mimicking the action of a vertical AOD),
        and moving the lowest atom from each column containing
        excess atoms.
        
        This process is repeated until all excess
        atoms are ejected."""
    matrix = np.array(copy.deepcopy(init_matrix))
    move_list = [] # list of all the moves in the ejection process

    # Find the leftmost atom in the target configuration
    if not eject_all:
        left_eject_index = np.where(np.any(target_config == 1, axis=1))[0][0]
        right_eject_index = np.where(np.any(target_config == 1, axis=1))[0][-1]
    else:
        left_eject_index = 0
        right_eject_index = len(init_matrix)

    while np.sum(matrix[:,left_eject_index:right_eject_index+1]) > np.sum(target_config[:,left_eject_index:right_eject_index+1]):
        # Iterate through each column and count the number of atoms
        lowest_atom_dict = {}
        for col in range(len(matrix[0])):
            if np.sum(matrix[:,col]) > np.sum(target_config[:,col]):
                row = find_lowest_atom_in_col(matrix[:,col])
                lowest_atom_dict[col] = row
        cols_to_move = []
        moves_in_scan = [] # list of all the moves in one downward scan
        row = min(lowest_atom_dict.values())
        while row <= len(matrix)-1:
            cols_to_move.extend([key for key, value in lowest_atom_dict.items() if value == row])
            for col in cols_to_move:
                if col >= left_eject_index and col <= right_eject_index:
                    move = Move(row, col, row+1, col)
                    # move = ((row,col),(row+1, col))
                    moves_in_scan.append(move)
            row += 1
            if len(moves_in_scan) > 0:
                move_list.append(moves_in_scan)
                matrix, _ = move_atoms(matrix, moves_in_scan)
                moves_in_scan = []
    return matrix, move_list


def parallel_vertical_sort(init_matrix, target_config, max_cycles = 20):
    """ Implements parallel sorting within columns, as described in [Ebadi et al. 2021](https://www.nature.com/articles/s41586-021-03582-4#Sec8)."""
    matrix = np.array(copy.deepcopy(init_matrix))

    rearrangement_cycle = 0
    success_flag = False
    move_list = []

    # Find the leftmost atom in the target configuration
    left_eject_index = np.where(np.any(target_config == 1, axis=1))[0][0]
    right_eject_index = np.where(np.any(target_config == 1, axis=1))[0][-1]

    # Execute Parallel Sorting until Defect-free or until max number of cycles
    while success_flag == False:
        if rearrangement_cycle == max_cycles:
            return matrix, move_list, success_flag

        # Moving atoms upward, iterating row by row
        for i in range(1,len(matrix))[::-1]:
            up_moves = []
            row = matrix[i]
            # check if the atom should be moved up, and if so note down the move
            for j, spot in enumerate(row):
                num_target_spots = np.sum(target_config[:i,j])
                num_atoms = np.sum(matrix[:i,j])
                if spot == 0 or num_target_spots <= num_atoms or matrix[i-1,j] == 1:
                    pass
                else:
                    # up_moves.append(((i,j),(i-1,j)))
                    if j >= left_eject_index and j <= right_eject_index:
                        up_moves.append(Move(i,j,i-1,j))
                    else:
                        pass
            if len(up_moves) > 0:
                matrix, _ = move_atoms(matrix,up_moves)
                move_list.append(up_moves)
        
        # Moving atoms downward, iterating row by row
        for _i in range(len(matrix)-1):
            down_moves = []
            row = matrix[_i]
            # check if the atom should be moved down, and if so note down the move
            for _j, spot in enumerate(row):
                num_target_spots = np.sum(target_config[_i+1:,_j])
                num_atoms = np.sum(matrix[_i+1:,_j])
                if spot == 0 or num_target_spots <= num_atoms or matrix[_i+1,_j] == 1:
                    pass
                else:
                    # down_moves.append(((_i,_j),(_i+1,_j)))
                    if _j >= left_eject_index and _j <= right_eject_index:
                        down_moves.append(Move(_i,_j,_i+1,_j))
                    else:
                        pass
            if len(down_moves) > 0:
                matrix, _ = move_atoms(matrix,down_moves)
                move_list.append(down_moves)
        rearrangement_cycle += 1

        if np.array_equal(matrix[left_eject_index:right_eject_index, :], target_config[left_eject_index:right_eject_index, :]):
            success_flag = True
        
    return matrix, move_list, success_flag

def find_leftmost_atom_in_array(array):
    for row in array:
        if 1 in row:
            return row.index(1)
    return -1  # If no "1" is found

def find_rightmost_atom_in_array(array):
    for row in array:
        if 1 in row[::-1]:
            return len(row) - 1 - row[::-1].index(1)
    return -1  # If no "1" is found