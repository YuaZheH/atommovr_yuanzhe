from atommover.utils import *
from atommover.algorithms.source.parallel_Hungarian import generate_AOD_cmds

def generalized_balance_dual_species(rbcs_array: AtomArray):
    """
    Implement and extend the generalized balance algorithm in the Bernien lab for dual species arrays
    """
    generalized_balance_success_flag = False
    arrays = copy.deepcopy(rbcs_array)
    move_list = []
    
    final_size = [0, len(arrays.matrix[0])-1, 0, len(arrays.matrix)-1]

    row_min = final_size[0]
    row_max = final_size[1]
    col_min = final_size[2]
    col_max = final_size[3]

    balance_config = copy.deepcopy(arrays)
    balance_config, move_list = row_balance(balance_config, balance_config.target, row_min, row_max, col_min, col_max, move_list, 0)
    balance_moves_term = len(move_list)
    ejection_moves_term = 0

    # if do_ejection:
    #     eject_moves, eject_config = ejection(balance_config, target_config, final_size)
    #     move_list.extend(eject_moves)
    #     ejection_moves_term = len(eject_moves)
        
    # if np.array_equal(eject_config, target_config):
    #     generalized_balance_success_flag = True
    
    return balance_config, move_list, generalized_balance_success_flag, [balance_moves_term, ejection_moves_term]

def row_balance(atom_arrays: AtomArray, target_config: np.ndarray, row_min, row_max, col_min, col_max, move_list, recursive_flag):
    """
    Balance the number of atoms of top/down sublattice
    """
    # 1. Top_Bottom Lattice Balance
    row_nums = row_max - row_min + 1
    matrix = atom_arrays.matrix

    # if row_nums == 1 and recursive_flag == 0:
    #     return col_balance(matrix, target_config, row_min, row_max, col_min, col_max, move_list, 1)

    # 2. Left_Right Lattice Balance
    middle_row = row_min + (row_nums // 2) - 1

    # 3-1. Calculate the required number of atoms in the submatrix S[i:m+1]
    n_req_Rb = int(np.sum(target_config[row_min:middle_row + 1, col_min:col_max + 1, 0]))
    n_req_Cs = int(np.sum(target_config[row_min:middle_row + 1, col_min:col_max + 1, 1]))
    n_req_c_Rb = int(np.sum(target_config[middle_row + 1:row_max + 1, col_min:col_max + 1, 0]))
    n_req_c_Cs = int(np.sum(target_config[middle_row + 1:row_max + 1, col_min:col_max + 1, 1]))

    # 3-2. Calculate the actual number of atoms in the submatrix S[i:m+1]
    N_S_i_m_Rb = int(np.sum(matrix[row_min:middle_row + 1, col_min:col_max + 1, 0]))
    N_S_i_m_Cs = int(np.sum(matrix[row_min:middle_row + 1, col_min:col_max + 1, 1]))
    N_S_i_m_c_Rb = int(np.sum(matrix[middle_row + 1:row_max + 1, col_min:col_max + 1, 0]))
    N_S_i_m_c_Cs = int(np.sum(matrix[middle_row + 1:row_max + 1, col_min:col_max + 1, 1]))

    # 3-3. atom_difference_Rb is positive when upper half part has more Rb atoms than required
    if N_S_i_m_Rb > n_req_Rb and N_S_i_m_c_Rb < n_req_c_Rb:
        atom_difference_Rb = min(N_S_i_m_Rb - n_req_Rb, n_req_c_Rb - N_S_i_m_c_Rb)
    elif N_S_i_m_Rb < n_req_Rb and N_S_i_m_c_Rb > n_req_c_Rb:
        atom_difference_Rb = -min(n_req_Rb - N_S_i_m_Rb, N_S_i_m_c_Rb - n_req_c_Rb)
    # Don't move atoms if one side does not have enough atoms or the other side does not need atoms
    else:
        atom_difference_Rb = 0
    
    # 3-4. atom_difference_Cs is positive when upper half part has more Cs atoms than required
    if N_S_i_m_Cs > n_req_Cs and N_S_i_m_c_Cs < n_req_c_Cs:
        atom_difference_Cs = min(N_S_i_m_Cs - n_req_Cs, n_req_c_Cs - N_S_i_m_c_Cs)
    elif N_S_i_m_Cs < n_req_Cs and N_S_i_m_c_Cs > n_req_c_Cs:
        atom_difference_Cs = -min(n_req_Cs - N_S_i_m_Cs, N_S_i_m_c_Cs - n_req_c_Cs)
    else:
        atom_difference_Cs = 0

    # 4. Implement the moves of Rb and Cs atoms together
    if atom_difference_Rb != 0 or atom_difference_Cs != 0:
        move_list, balance_config = vertical_move(atom_arrays, target_config, atom_difference_Rb, atom_difference_Cs, row_min, middle_row, row_max, col_min, col_max, move_list)
    else:
        balance_config = atom_arrays
    
    balance_config, move_list = col_balance(balance_config, target_config, row_min, middle_row, col_min, col_max, move_list, recursive_flag)
    balance_config, move_list = col_balance(balance_config, target_config, middle_row + 1, row_max, col_min, col_max, move_list, recursive_flag)
        
    return balance_config, move_list

def col_balance(atom_arrays: AtomArray, target_config: np.ndarray, row_min, row_max, col_min, col_max, move_list, recursive_flag):
    """
    Balance the number of atoms of left/right sublattice
    """
    # 1. Left_Right Lattice Balance
    col_nums = col_max - col_min + 1
    matrix = atom_arrays.matrix

    # 2. Top_Bottom Lattice Balance
    middle_col = col_min + (col_nums // 2) - 1

    # Calculate the required number of 1s in the submatrix S[i:m+1]
    n_req_Rb = int(np.sum(target_config[row_min:row_max + 1, col_min:middle_col + 1, 0]))
    n_req_Cs = int(np.sum(target_config[row_min:row_max + 1, col_min:middle_col + 1, 1]))
    n_req_c_Rb = int(np.sum(target_config[row_min:row_max + 1, middle_col + 1:col_max + 1, 0]))
    n_req_c_Cs = int(np.sum(target_config[row_min:row_max + 1, middle_col + 1:col_max + 1, 1]))
    
    # Calculate the actual number of 1s in the submatrix S[i:m+1]
    N_S_i_m_Rb = int(np.sum(matrix[row_min:row_max + 1, col_min:middle_col + 1, 0]))
    N_S_i_m_Cs = int(np.sum(matrix[row_min:row_max + 1, col_min:middle_col + 1, 1]))
    N_S_i_m_c_Rb = int(np.sum(matrix[row_min:row_max + 1, middle_col + 1:col_max + 1, 0]))
    N_S_i_m_c_Cs = int(np.sum(matrix[row_min:row_max + 1, middle_col + 1:col_max + 1, 1]))
    
    # atom_difference_Rb is positive when left half part has more Rb atoms than required
    if N_S_i_m_Rb > n_req_Rb and N_S_i_m_c_Rb < n_req_c_Rb:
        atom_difference_Rb = min(N_S_i_m_Rb - n_req_Rb, n_req_c_Rb - N_S_i_m_c_Rb)
    elif N_S_i_m_Rb < n_req_Rb and N_S_i_m_c_Rb > n_req_c_Rb:
        atom_difference_Rb = -min(n_req_Rb - N_S_i_m_Rb, N_S_i_m_c_Rb - n_req_c_Rb)
    else:
        atom_difference_Rb = 0
    
    # atom_difference_Cs is positive when left half part has more Cs atoms than required
    if N_S_i_m_Cs > n_req_Cs and N_S_i_m_c_Cs < n_req_c_Cs:
        atom_difference_Cs = min(N_S_i_m_Cs - n_req_Cs, n_req_c_Cs - N_S_i_m_c_Cs)
    elif N_S_i_m_Cs < n_req_Cs and N_S_i_m_c_Cs > n_req_c_Cs:
        atom_difference_Cs = -min(n_req_Cs - N_S_i_m_Cs, N_S_i_m_c_Cs - n_req_c_Cs)
    else:
        atom_difference_Cs = 0
    
    # Implement the moves of Rb and Cs atoms together. NB: We may move left/right Rb atoms while move right/left Cs atoms
    if atom_difference_Rb != 0 or atom_difference_Cs != 0:
        move_list, balance_config = horizontal_move(atom_arrays, target_config, atom_difference_Rb, atom_difference_Cs, row_min, row_max, col_min, middle_col, col_max, move_list)
    else:
        balance_config = atom_arrays
    
    return balance_config, move_list

def vertical_move(atom_arrays_copy, target_config, Rb_atoms, Cs_atoms, row_min, middle_row, row_max, col_min, col_max, move_list):
    ## Determine the move direction for Rb and Cs atoms
    atom_arrays = copy.deepcopy(atom_arrays_copy)
    if Rb_atoms > 0:
        # Implement down_move function in single-species case
        move_direction_Rb = 1
        source_row_Rb = middle_row
        target_row_Rb = middle_row + 1
    # Otherwise, move atoms to up (decrease row index)
    else:
        # Implement up_move function in single-species case
        move_direction_Rb = -1
        source_row_Rb = middle_row + 1
        target_row_Rb = middle_row

    normalize_row_Rb = 0
    balance_row_count_Rb = (target_row_Rb - source_row_Rb)*move_direction_Rb
    stuff_Rb = 0

    if Cs_atoms > 0:
        move_direction_Cs = 1
        source_row_Cs = middle_row
        target_row_Cs = middle_row + 1
    else:
        move_direction_Cs = -1
        source_row_Cs = middle_row + 1
        target_row_Cs = middle_row

    normalize_row_Cs = 0
    balance_row_count_Cs = (target_row_Cs - source_row_Cs)*move_direction_Cs
    stuff_Cs = 0

    # Implement Rb move function independently
    start_row = source_row_Rb+normalize_row_Rb
    end_row = source_row_Rb+normalize_row_Rb+move_direction_Rb

    # Move Rb and Cs atoms
    # Check the boundary conditions
    while abs(Rb_atoms) > 0 and row_max >= start_row >= row_min and row_max >= end_row >= row_min:
        moves_in_scan_Rb = []
        for col in range(col_min, col_max+1):
            # Check if the source column has Rb atoms and the target column is empty
            if atom_arrays.matrix[source_row_Rb + normalize_row_Rb, col, 0] == 1 and np.sum(atom_arrays.matrix[source_row_Rb + normalize_row_Rb + move_direction_Rb, col, :]) == 0 and np.sum(atom_arrays.matrix[target_row_Rb, col, :]) == 0 and abs(Rb_atoms) > 0:
                move = Move(source_row_Rb + normalize_row_Rb, col, source_row_Rb + normalize_row_Rb + move_direction_Rb, col)
                moves_in_scan_Rb.append(move)

                if balance_row_count_Rb == 1:
                    Rb_atoms -= move_direction_Rb

        if len(moves_in_scan_Rb) > 0:
            move_list.append(moves_in_scan_Rb)
            _ = atom_arrays.move_atoms(moves_in_scan_Rb)
        
        # Check the boundary conditions
        if 0 <= target_row_Rb + move_direction_Rb <= len(atom_arrays.matrix)-1:
            # Check if the target column is blocked
            if np.sum(atom_arrays.matrix[target_row_Rb,col_min:col_max+1,:]) == len(atom_arrays.matrix[target_row_Rb,col_min:col_max,0]) and 0<=target_row_Rb+stuff_Rb <= len(atom_arrays.matrix)-1:
                for shift in range(stuff_Rb, -1, -1):
                    moves_in_scan_Rb = []
                    for col in range(col_min, col_max+1):
                        if np.sum(atom_arrays.matrix[target_row_Rb+move_direction_Rb*shift,col,:]) == 1 and np.sum(atom_arrays.matrix[target_row_Rb+move_direction_Rb*(shift+1),col,:]) == 0:
                            move = Move(target_row_Rb+move_direction_Rb*shift, col, target_row_Rb+move_direction_Rb*(shift+1), col)
                            moves_in_scan_Rb.append(move)

                if len(moves_in_scan_Rb) > 0:
                    move_list.append(moves_in_scan_Rb)
                    _ = atom_arrays.move_atoms(moves_in_scan_Rb)
                stuff_Rb += 1
                source_row_Rb = middle_row
                balance_row_count_Rb = (target_row_Rb - source_row_Rb)*move_direction_Rb
                normalize_row_Rb = 0
            else:
                balance_row_count_Rb -= 1
                normalize_row_Rb += move_direction_Rb
        else:
            break

        if balance_row_count_Rb == 0:
            source_row_Rb -= move_direction_Rb
            balance_row_count_Rb = (target_row_Rb - source_row_Rb)*move_direction_Rb
            normalize_row_Rb = 0
    
    while abs(Cs_atoms) > 0 and row_max >= start_row >= row_min and row_max >= end_row >= row_min:
        moves_in_scan_Cs = []
        for col in range(col_min, col_max+1):
            # Check if the source column has Cs atoms and the target column is empty
            if atom_arrays.matrix[source_row_Cs + normalize_row_Cs, col, 1] == 1 and np.sum(atom_arrays.matrix[source_row_Cs + normalize_row_Cs + move_direction_Cs, col, :]) == 0 and np.sum(atom_arrays.matrix[target_row_Cs, col, :]) == 0 and abs(Cs_atoms) > 0:
                move = Move(source_row_Cs + normalize_row_Cs, col, source_row_Cs + normalize_row_Cs + move_direction_Cs, col)
                moves_in_scan_Cs.append(move)

                if balance_row_count_Cs == 1:
                    Cs_atoms -= move_direction_Cs
        
        if len(moves_in_scan_Cs) > 0:
            move_list.append(moves_in_scan_Cs)
            _ = atom_arrays.move_atoms(moves_in_scan_Cs)

        if 0 <= target_row_Cs + move_direction_Cs <= len(atom_arrays.matrix)-1:
            if np.sum(atom_arrays.matrix[target_row_Cs,col_min:col_max+1,:]) == len(atom_arrays.matrix[target_row_Cs,col_min:col_max,0]) and 0<=target_row_Cs+stuff_Cs <= len(atom_arrays.matrix)-1:
                for shift in range(stuff_Cs, -1, -1):
                    moves_in_scan_Cs = []
                    for col in range(col_min, col_max+1):
                        if np.sum(atom_arrays.matrix[target_row_Cs+move_direction_Cs*shift,col,:]) == 1 and np.sum(atom_arrays.matrix[target_row_Cs+move_direction_Cs*(shift+1),col,:]) == 0:
                            move = Move(target_row_Cs+move_direction_Cs*shift, col, target_row_Cs+move_direction_Cs*(shift+1), col)
                            moves_in_scan_Cs.append(move)

                if len(moves_in_scan_Cs) > 0:
                    move_list.append(moves_in_scan_Cs)
                    _ = atom_arrays.move_atoms(moves_in_scan_Cs)
                stuff_Cs += 1
                source_row_Cs = middle_row
                balance_row_count_Cs = (target_row_Cs - source_row_Cs)*move_direction_Cs
                normalize_row_Cs = 0
            else:
                balance_row_count_Cs -= 1
                normalize_row_Cs += move_direction_Cs
        else:
            break

        if balance_row_count_Cs == 0:
            source_row_Cs -= move_direction_Cs
            balance_row_count_Cs = (target_row_Cs - source_row_Cs)*move_direction_Cs
            normalize_row_Cs = 0

    return move_list, atom_arrays

def horizontal_move(atom_arrays_copy, target_config, Rb_atoms, Cs_atoms, row_min, row_max, col_min, middle_col, col_max, move_list):
    ## Determine the move direction for Rb and Cs atoms
    # If more atoms in the left site, move atoms to right (increase column index)
    atom_arrays = copy.deepcopy(atom_arrays_copy)

    if Rb_atoms > 0:
        # Implement right_move function in single-species case
        move_direction_Rb = 1
        source_col_Rb = middle_col
        target_col_Rb = middle_col + 1
    # Otherwise, move atoms to left (decrease column index)
    else:
        # Implement left_move function in single-species case
        move_direction_Rb = -1
        source_col_Rb = middle_col + 1
        target_col_Rb = middle_col
    normalize_col_Rb = 0
    balance_col_count_Rb = (target_col_Rb - source_col_Rb)*move_direction_Rb
    stuff_Rb = 0

    if Cs_atoms > 0:
        move_direction_Cs = 1
        source_col_Cs = middle_col
        target_col_Cs = middle_col + 1
    else:
        move_direction_Cs = -1
        source_col_Cs = middle_col + 1
        target_col_Cs = middle_col

    normalize_col_Cs = 0
    balance_col_count_Cs = (target_col_Cs - source_col_Cs)*move_direction_Cs
    stuff_Cs = 0

    # Implement Rb move function independently
    start_col = source_col_Rb+normalize_col_Rb
    end_col = source_col_Rb+normalize_col_Rb+move_direction_Rb

    # Move Rb and Cs atoms
    # Check the boundary conditions
    while abs(Rb_atoms) > 0 and col_max >= start_col >= col_min and col_max >= end_col >= col_min:
        moves_in_scan_Rb = []
        for row in range(row_min, row_max+1):
            # Check if the source column has Rb atoms and the target column is empty
            if atom_arrays.matrix[row, source_col_Rb + normalize_col_Rb, 0] == 1 and np.sum(atom_arrays.matrix[row, source_col_Rb + normalize_col_Rb + move_direction_Rb,:]) == 0 and np.sum(atom_arrays.matrix[row, target_col_Rb, :]) == 0 and abs(Rb_atoms) > 0:
                move = Move(row, source_col_Rb + normalize_col_Rb, row, source_col_Rb + normalize_col_Rb + move_direction_Rb)
                moves_in_scan_Rb.append(move)

                if balance_col_count_Rb == 1:
                    Rb_atoms -= move_direction_Rb

        if len(moves_in_scan_Rb) > 0:
            move_list.append(moves_in_scan_Rb)
            _ = atom_arrays.move_atoms(moves_in_scan_Rb)
        
        # Check the boundary conditions
        if 0 <= target_col_Rb + move_direction_Rb <= len(atom_arrays.matrix)-1:
            # Check if the target column is blocked
            if np.sum(atom_arrays.matrix[row_min:row_max+1,target_col_Rb,:]) == len(atom_arrays.matrix[row_min:row_max+1,target_col_Rb,0]) and 0<=target_col_Rb+stuff_Rb <= len(atom_arrays.matrix)-1:
                for shift in range(stuff_Rb, -1, -1):
                    moves_in_scan_Rb = []
                    for row in range(row_min, row_max+1):
                        if np.sum(atom_arrays.matrix[row, target_col_Rb+move_direction_Rb*shift,:]) == 1 and np.sum(atom_arrays.matrix[row, target_col_Rb+move_direction_Rb*(shift+1),:]) == 0:
                            move = Move(row, target_col_Rb+move_direction_Rb*shift, row, target_col_Rb+move_direction_Rb*(shift+1))
                            moves_in_scan_Rb.append(move)

                if len(moves_in_scan_Rb) > 0:
                    move_list.append(moves_in_scan_Rb)
                    _ = atom_arrays.move_atoms(moves_in_scan_Rb)
                stuff_Rb += 1
                source_col_Rb = middle_col
                balance_col_count_Rb = (target_col_Rb - source_col_Rb)*move_direction_Rb
                normalize_col_Rb = 0
            else:
                balance_col_count_Rb -= 1
                normalize_col_Rb += move_direction_Rb
        else:
            break

        if balance_col_count_Rb == 0:
            source_col_Rb -= move_direction_Rb
            balance_col_count_Rb = (target_col_Rb - source_col_Rb)*move_direction_Rb
            normalize_col_Rb = 0
    
    while abs(Cs_atoms) > 0 and col_max >= start_col >= col_min and col_max >= end_col >= col_min:
        moves_in_scan_Cs = []
        for row in range(row_min, row_max+1):
            # Check if the source column has Cs atoms and the target column is empty
            if atom_arrays.matrix[row, source_col_Cs + normalize_col_Cs, 1] == 1 and np.sum(atom_arrays.matrix[row, source_col_Cs + normalize_col_Cs + move_direction_Cs,:]) == 0 and np.sum(atom_arrays.matrix[row, target_col_Cs, :]) == 0 and abs(Cs_atoms) > 0:
                move = Move(row, source_col_Cs + normalize_col_Cs, row, source_col_Cs + normalize_col_Cs + move_direction_Cs)
                moves_in_scan_Cs.append(move)

                if balance_col_count_Cs == 1:
                    Cs_atoms -= move_direction_Cs
        
        if len(moves_in_scan_Cs) > 0:
            move_list.append(moves_in_scan_Cs)
            _ = atom_arrays.move_atoms(moves_in_scan_Cs)
        
        if 0 <= target_col_Cs + move_direction_Cs <= len(atom_arrays.matrix)-1:
            if np.sum(atom_arrays.matrix[row_min:row_max+1,target_col_Cs,:]) == len(atom_arrays.matrix[row_min:row_max+1,target_col_Cs,0]) and 0<=target_col_Cs+stuff_Cs <= len(atom_arrays.matrix)-1:
                for shift in range(stuff_Cs, -1, -1):
                    moves_in_scan_Cs = []
                    for row in range(row_min, row_max+1):
                        if np.sum(atom_arrays.matrix[row, target_col_Cs+move_direction_Cs*shift,:]) == 1 and np.sum(atom_arrays.matrix[row, target_col_Cs+move_direction_Cs*(shift+1),:]) == 0:
                            move = Move(row, target_col_Cs+move_direction_Cs*shift, row, target_col_Cs+move_direction_Cs*(shift+1))
                            moves_in_scan_Cs.append(move)

                if len(moves_in_scan_Cs) > 0:
                    move_list.append(moves_in_scan_Cs)
                    _ = atom_arrays.move_atoms(moves_in_scan_Cs)
                stuff_Cs += 1
                source_col_Cs = middle_col
                balance_col_count_Cs = (target_col_Cs - source_col_Cs)*move_direction_Cs
                normalize_col_Cs = 0
            else:
                balance_col_count_Cs -= 1
                normalize_col_Cs += move_direction_Cs
        else:
            break
    
        if balance_col_count_Cs == 0:
            source_col_Cs -= move_direction_Cs
            balance_col_count_Cs = (target_col_Cs - source_col_Cs)*move_direction_Cs
            normalize_col_Cs = 0
    return move_list, atom_arrays