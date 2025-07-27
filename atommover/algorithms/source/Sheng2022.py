import copy
import time
import numpy as np
from sys import maxsize
from collections import deque

from atommover.utils.core import *
from atommover.utils.move_utils import Move
from atommover.utils.AtomArray import AtomArray


def first_dual_species_rearrange(rbcs_arrays: AtomArray):
    """
    Generates a list of moves that implements the dual-species rearrangement protocol described in [Sheng et al. 2022](https://journals.aps.org/prl/article/10.1103/PhysRevLett.128.083202/figures/2/large).
    """
    move_list = []

    arrays = copy.deepcopy(rbcs_arrays)
    # print("Initial arrays:")
    # arrays.image()
    # Initialize the information for atom arrays to group unfinished traps and identify M-, R- type atoms
    array_info = prepare_array_info(arrays)

    # 1. Fill innermost empty site (type-M)
    #innermost_M_config, innermost_M_moves, array_info = fill_innermost_empty_sites(arrays, array_info, 'M')
    fill_empty_config, fill_empty_moves, array_info = fill_innermost_empty_sites(arrays, array_info, 'M')
    move_list.extend(fill_empty_moves)

    # 2. Move out the rest of outermost M type atoms and fill all remaining sites
    if len(array_info['Type_M_Rb']+array_info['Type_M_Cs']) > 0:
        remove_M_config, remove_M_moves, array_info = remove_M_atoms(fill_empty_config, array_info)
        move_list.extend(remove_M_moves)
    else:
        remove_M_config = copy.deepcopy(fill_empty_config)

    # # 3. Fill the remaining empty sites (including the first round empty and M-type caused vacancies)
    # final_config, remaining_empty_moves = fill_remaining_empty_sites(remove_M_config, array_info)
    # move_list.extend(remaining_empty_moves)
    array_info = prepare_array_info(remove_M_config)
    final_config, final_moves, array_info = fill_innermost_empty_sites(remove_M_config, array_info, 'R')
    move_list.extend(final_moves)

    # putting moves in right format
    final_move_list = []
    for move in move_list:
        final_move_list.append([move])

    return final_config, final_move_list, array_info

### Helper functions for fill_innermost_empty_sites ###
def prepare_array_info(arrays: AtomArray) -> dict:
    """
    Unified function used to prepare array info
    """
    array_info = group_unfinished_atoms(arrays, find_outermost_atoms(arrays, classified_unfinished_atoms(arrays)))
    return array_info


#Input: matrix configuration, target; Output: list of groups of unfinished traps
def classified_unfinished_atoms(arrays: AtomArray) -> list:
    """
    Function used to identify traps in dual-species atom arrays.
    """
    array_info = {
        'Rb_fin': np.transpose(np.where(np.multiply(arrays.matrix[:,:,0], arrays.target_Rb) == 1)).tolist(),
        'Cs_fin': np.transpose(np.where(np.multiply(arrays.matrix[:,:,1], arrays.target_Cs) == 1)).tolist(),
        'Rb_wr': np.transpose(np.where(np.multiply(arrays.matrix[:,:,0], arrays.target_Cs) == 1)).tolist(),
        'Cs_wr': np.transpose(np.where(np.multiply(arrays.matrix[:,:,1], arrays.target_Rb) == 1)).tolist(),
        'Rb_emp': np.transpose(np.where((np.sum(arrays.matrix, axis=2) == 0) & (arrays.target_Rb == 1))).tolist(),
        'Cs_emp': np.transpose(np.where((np.sum(arrays.matrix, axis=2) == 0) & (arrays.target_Cs == 1))).tolist(),
        'Type_M_Rb': np.transpose(np.where(np.multiply(arrays.matrix[:,:,0], arrays.target_Cs) == 1)).tolist(),
        'Type_M_Cs': np.transpose(np.where(np.multiply(arrays.matrix[:,:,1], arrays.target_Rb) == 1)).tolist(),
        'Type_R_Rb': [],
        'Type_R_Cs': [],
        'left': min(np.where(np.any(arrays.target_Cs == 1, axis=1))[0][0], np.where(np.any(arrays.target_Rb == 1, axis=1))[0][0]),
        'right': max(np.where(np.any(arrays.target_Cs == 1, axis=1))[0][-1], np.where(np.any(arrays.target_Rb == 1, axis=1))[0][-1]),
        'top': min(np.where(np.any(arrays.target_Cs == 1, axis=0))[0][0], np.where(np.any(arrays.target_Rb == 1, axis=0))[0][0]),
        'bot': max(np.where(np.any(arrays.target_Cs == 1, axis=0))[0][-1], np.where(np.any(arrays.target_Rb == 1, axis=0))[0][-1]),
        'boundary_coordinate': [],
        'outermost_unfinished_sites': [],
        'unfinished_groups': [],
        'fail_flag': False
    }

    for i in range(len(arrays.matrix)):
        for j in range(len(arrays.matrix[0])):
            if [i,j] not in (array_info['Rb_fin'] + array_info['Rb_wr']) and arrays.matrix[i][j][0] == 1:
                array_info['Type_R_Rb'].append([i,j])
            if [i,j] not in (array_info['Cs_fin']+array_info['Cs_wr']) and arrays.matrix[i][j][1] == 1:
                array_info['Type_R_Cs'].append([i,j])

    return array_info

## Find the outermost unfinished atoms ##
def find_outermost_atoms(arrays: AtomArray, array_info: dict):
    # Find the boundary of the target configuration
    for i in range(array_info['left'], array_info['right']+1):
        array_info['boundary_coordinate'].append([i, array_info['top']])
        array_info['boundary_coordinate'].append([i, array_info['bot']])

    # Exclude the corner atoms we have counted
    for j in range(array_info['top']+1, array_info['bot']):
        array_info['boundary_coordinate'].append([array_info['left'], j])
        array_info['boundary_coordinate'].append([array_info['right'], j])

    # Find the outermost unfinished atoms
    for boundary_point in array_info['boundary_coordinate']:
        if boundary_point not in (array_info['Rb_fin'] + array_info['Cs_fin']):
            array_info['outermost_unfinished_sites'].append(boundary_point)
    
    return array_info

def group_unfinished_atoms(arrays: AtomArray, array_info: dict):
    visited_sites = np.zeros(arrays.matrix.shape[:2])
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    first_round_flag = True

    # Mark finished sites as visited
    for site in (array_info['Rb_fin']+array_info['Cs_fin']):
        visited_sites[site[0]][site[1]] = 1

    limit_cycle = 0

    while np.array_equal(visited_sites, (arrays.target_Cs + arrays.target_Rb)) == False:
        group_iteration = 0

        if limit_cycle == 10:
            break
        limit_cycle += 1
        # Iterate through all the outermost unfinished sites
        for group_start_pos in array_info['outermost_unfinished_sites']:

            # Initialize the queue and visited sites
            if first_round_flag:
                queue = deque([(group_start_pos[0], group_start_pos[1])])
                visited_sites[group_start_pos[0]][group_start_pos[1]] = 1
                group = [(group_start_pos[0], group_start_pos[1])]
            else:
                queue = deque(array_info['unfinished_groups'][group_iteration])

            # Start searching for the group
            while queue:
                found_atoms = False
                i, j = queue.popleft()
                for di, dj in directions:
                    ni, nj = i + di, j + dj
                    # Check if the site is within the target configuration and has not been visited
                    if array_info['left'] <= ni <= array_info['right'] and array_info['top'] <= nj <= array_info['bot'] and visited_sites[ni][nj] == 0 and [ni, nj] not in array_info['boundary_coordinate'] and [ni, nj] not in (array_info['Rb_fin'] + array_info['Cs_fin']):
                        queue.append((ni, nj))
                        visited_sites[ni][nj] = 1
                        if first_round_flag:
                            group.append((ni, nj))
                        else:
                            array_info['unfinished_groups'][group_iteration].append((ni, nj))
                        found_atoms = True
                        break
                if found_atoms:
                    break
            group_iteration += 1
            if first_round_flag:
                array_info['unfinished_groups'].append(group)
        first_round_flag = False
    return array_info
### End of helper functions ###

# Fill the innermost empty sites of dual-species arrays
def fill_innermost_empty_sites(arrays: AtomArray, array_info: dict, source_types: str = 'M', max_limit_cycle=20):
    """ 
    Function used to fill the innermost empty sites of dual-species atom arrays.

    NB: this differs from the Sheng et al PRL 2022 implementation in the following way:
    - On the second paragraph of the right column of pg 3, they specify that when they create new empty sites that used to be type-M atoms,
      they immediately iterate through the newly created empty sites and fill them with other type-M atoms.
      We do not do this, as we want to fill the empty sites that are closest to the middle first, to avoid creating vacancies.
    """
    #Initialize variables
    innermost_arrays = copy.deepcopy(arrays)
    innermost_moves = []
    limit_cycle = 0

    # For the process filling innermost empty sites (type-M), target is to use all available type-M atoms
    if source_types == 'M':
        control_sequence_len = len(array_info['Type_M_Rb'] + array_info['Type_M_Cs'])

    # For the process filling innermost empty sites (type-R), target is to use all available type-R atoms
    if source_types == 'R':
        control_sequence_len = len(array_info['Type_R_Rb'] + array_info['Type_R_Cs'])

    # Iterate through the unfinished groups and find their innermost empty site
    while control_sequence_len > 0:
        # Build a list to store all innermost empty sites' coordinates
        innermost_empty_traps = []
        array_info['fail_flag'] = False
        array_info['unfinished_groups'] = []
        array_info = group_unfinished_atoms(innermost_arrays, array_info)

        # Limit the cycle to avoid infinite loop
        if limit_cycle == max_limit_cycle:
            if source_types == 'M':
                source_types = 'R'
                limit_cycle = 0
            else:
                break
        limit_cycle += 1

        # Iterate through the unfinished groups and find their innermost empty site
        for group in array_info['unfinished_groups']:
            # Iterate the sites through the unfinished traps in the group
            for unfinished_trap in group[::-1]:
                if list(unfinished_trap) in array_info['Rb_emp']:
                    innermost_empty_traps.append([unfinished_trap, 'Rb'])
                    break
                elif list(unfinished_trap) in array_info['Cs_emp']:
                    innermost_empty_traps.append([unfinished_trap, 'Cs'])
                    break
        
        # Sort the innermost empty trapas based on the distance to the center of the arrays
        innermost_empty_traps.sort(key = lambda atom: (atom[0][0] - (array_info['top']+array_info['bot'])/2)**2 + (atom[0][1] - (array_info['left']+array_info['right'])/2)**2)

        # Fill each innermost empty site with the nearest type-M or type_R atoms
        for [empty_site, atom_type] in innermost_empty_traps:
            moves_in_scan = []
            # Update the source atoms list based on the source type
            if source_types == 'M':
                source_atoms = array_info['Type_M_{}'.format(atom_type)]
            if source_types == 'R':
                source_atoms = array_info['Type_R_{}'.format(atom_type)]

            # If there is no type-M atom, continue to the next empty site
            if len(source_atoms) == 0:
                break
            else:
                nearest_atom = min(source_atoms, key=lambda atom: (atom[0] - empty_site[0])**2 + (atom[1] - empty_site[1])**2)
                source_atoms.sort(key = lambda atom: (atom[0] - empty_site[0])**2 + (atom[1] - empty_site[1])**2)
                nearest_atoms = copy.deepcopy(source_atoms)

            path_found_flag = False

            while not path_found_flag:
                try:
                    nearest_atom = nearest_atoms[0]
                    path, path_found_flag = bfs_find_path(innermost_arrays, nearest_atom, empty_site)
                except IndexError:
                    break
                # Transform path to move objects
                if path_found_flag:
                    path = flatten_tuple(path)
                    path = path[::-1]
                    for item in path:
                        current_pos = item[0] 
                        for next_pos in item:
                            if np.array_equal(current_pos, next_pos):
                                pass
                            else:
                                moves_in_scan.append(Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1]))
                                current_pos = next_pos
                else:
                    nearest_atoms.remove(nearest_atom)

            # Implement the moves and update the arrays
            if len(moves_in_scan) > 0:
                # Apply the moves to the arrays
                innermost_moves.extend(moves_in_scan)
                _ = innermost_arrays.move_atoms(moves_in_scan)

                if source_types == 'M':
                    if atom_type == 'Cs':
                        array_info['Type_M_Cs'].remove(nearest_atom)
                        array_info['Cs_emp'].remove(list(empty_site))
                        array_info['Rb_emp'].append(nearest_atom)
                        array_info['Cs_fin'].append(list(empty_site))
                        
                    elif atom_type == 'Rb':
                        array_info['Type_M_Rb'].remove(nearest_atom)
                        array_info['Rb_emp'].remove(list(empty_site))
                        array_info['Cs_emp'].append(nearest_atom)
                        array_info['Rb_fin'].append(list(empty_site))
                    
                # Because the source atom is R-type, we do not need to append it into the Cs_empty_coordinate
                if source_types == 'R':
                    if atom_type == 'Cs':
                        array_info['Type_R_Cs'].remove(nearest_atom)
                        array_info['Cs_emp'].remove(list(empty_site))
                        array_info['Cs_fin'].append(list(empty_site))
                    elif atom_type == 'Rb':
                        array_info['Type_R_Rb'].remove(nearest_atom)
                        array_info['Rb_emp'].remove(list(empty_site))
                        array_info['Rb_fin'].append(list(empty_site))

            #  # If we use M-type atom, find a new M-type atom to fill new empty site
            # if source_types == 'M':
            #     innermost_arrays, moves_in_scan, array_info = refill_source_atoms(innermost_arrays, nearest_atom, array_info, 'Rb')
    
        # Update the control sequence length
        if source_types == 'M':
            control_sequence_len = len(array_info['Type_M_Rb']+array_info['Type_M_Cs'])
            if control_sequence_len == 0:
                source_types = 'R'
                limit_cycle = 0
                
        if source_types == 'R':
            control_sequence_len = len(array_info['Rb_emp'] + array_info['Cs_emp'])
    
    array_info['non_isolated_remaining_M'] = array_info['Type_M_Rb'] + array_info['Type_M_Cs']

    return innermost_arrays, innermost_moves, array_info

def fill_innermost_empty_sites_old(arrays: AtomArray, array_info: dict, source_types: str, max_limit_cycle=10):
    """ 
    Function used to fill the innermost empty sites of dual-species atom arrays.

    NB: this differs from the Sheng et al PRL 2022 implementation in the following way:
    - On the second paragraph of the right column of pg 3, they specify that when they create new empty sites that used to be type-M atoms,
      they immediately iterate through the newly created empty sites and fill them with other type-M atoms.
      We do not do this, as we want to fill the empty sites that are closest to the middle first, to avoid creating vacancies.
    """
    #Initialize variables
    innermost_arrays = copy.deepcopy(arrays)
    innermost_moves = []
    limit_cycle = 0

    # For the process filling innermost empty sites (type-M), target is to use all available type-M atoms
    if source_types == 'M':
        control_sequence_len = len(array_info['Type_M_Rb'] + array_info['Type_M_Cs'])

    # For the process filling innermost empty sites (type-R), target is to use all available type-R atoms
    if source_types == 'R':
        control_sequence_len = len(array_info['Type_R_Rb'] + array_info['Type_R_Cs'])

    # Iterate through the unfinished groups and find their innermost empty site
    while control_sequence_len > 0:
        # Build a list to store all innermost empty sites' coordinates
        innermost_empty_traps_Rb = []
        innermost_empty_traps_Cs = []

        # Limit the cycle to avoid infinite loop
        if limit_cycle == max_limit_cycle:
            break
        limit_cycle += 1

        # Iterate through the unfinished groups and find their innermost empty site
        for group in array_info['unfinished_groups']:
            # Iterate the sites through the unfinished traps in the group
            for unfinished_trap in group[::-1]:
                if list(unfinished_trap) in array_info['Rb_emp']:
                    innermost_empty_traps_Rb.append(unfinished_trap)
                    break
                elif list(unfinished_trap) in array_info['Cs_emp']:
                    innermost_empty_traps_Cs.append(unfinished_trap)
                    break
        
        innermost_empty_traps_Rb.sort(key = lambda atom: (atom[0] - len(arrays.matrix)/2)**2 + (atom[1] - len(arrays.matrix[0])/2)**2)
        innermost_empty_traps_Cs.sort(key = lambda atom: (atom[0] - len(arrays.matrix)/2)**2 + (atom[1] - len(arrays.matrix[0])/2)**2)

        # Fill each innermost empty site with the nearest type-M or type_R atoms
        for empty_site in innermost_empty_traps_Rb:
            moves_in_scan = []
            # Update the source atoms list based on the source type
            if source_types == 'M':
                source_atoms = array_info['Type_M_Rb']
            if source_types == 'R':
                source_atoms = array_info['Type_R_Rb']

    
            # If there is no type-M atom, continue to the next empty site
            if len(source_atoms) == 0:
                break
            else:
                nearest_atom = min(source_atoms, key=lambda atom: abs(atom[0] - empty_site[0]) + abs(atom[1] - empty_site[1]))
                source_atoms.sort(key = lambda atom: (atom[0] - empty_site[0])**2 + (atom[1] - empty_site[1])**2)
                nearest_atoms = copy.deepcopy(source_atoms)

            path_found_flag = False
            while not path_found_flag:
                try:
                    nearest_atom = nearest_atoms[0]
                    path, path_found_flag = bfs_find_path(innermost_arrays, nearest_atom, empty_site)
                except IndexError:
                    break
                # Transform path to move objects
                if path_found_flag:
                    path = flatten_tuple(path)
                    path = path[::-1]
                    for item in path:
                        current_pos = item[0] 
                        for next_pos in item:
                            if np.array_equal(current_pos, next_pos):
                                pass
                            else:
                                moves_in_scan.append(Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1]))
                                current_pos = next_pos

                else:
                    nearest_atoms.remove(nearest_atom)

            # # If there is no route to the empty site, find the second nearest atom
            # if len(moves_in_scan) == 0 and path_found_flag == False:
            #     # Find the second nearest atom
            #     temp_source_atom = nearest_atom
            #     source_atoms.remove(nearest_atom)
            #     if len(source_atoms) > 0:
            #         nearest_atom = min(source_atoms, key=lambda atom: abs(atom[0] - empty_site[0]) + abs(atom[1] - empty_site[1]))
            #         path, path_found_flag = bfs_find_path(innermost_arrays, nearest_atom, empty_site)

            #         if path_found_flag:
            #             # Transform path to move objects
            #             path = flatten_tuple(path)
            #             path = path[::-1]
            #             for item in path:
            #                 current_pos = item[0] 
            #                 for next_pos in item:
            #                     if np.array_equal(current_pos, next_pos):
            #                         pass
            #                     else:
            #                         moves_in_scan.append(Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1]))
            #                         current_pos = next_pos
            #     source_atoms.append(temp_source_atom)


            # Implement the moves and update the arrays
            if len(moves_in_scan) > 0 and path_found_flag == True:
                # Apply the moves to the arrays
                innermost_moves.extend(moves_in_scan)
                _ = innermost_arrays.move_atoms(moves_in_scan)
                moves_in_scan = []

                if source_types == 'M':
                    array_info['Type_M_Rb'].remove(nearest_atom)
                    array_info['Rb_emp'].remove(list(empty_site))
                    array_info['Cs_emp'].append(nearest_atom)
                    array_info['Rb_fin'].append(list(empty_site))
                    
                # Because the source atom is R-type, we do not need to append it into the Cs_empty_coordinate
                if source_types == 'R':
                    array_info['Type_R_Rb'].remove(nearest_atom)
                    array_info['Rb_emp'].remove(list(empty_site))
                    array_info['Rb_fin'].append(list(empty_site))

            # If we use M-type atom, find a new M-type atom to fill new empty site
            if source_types == 'M':
                innermost_arrays, moves_in_scan, array_info = refill_source_atoms(innermost_arrays, nearest_atom, array_info, 'Cs')
            
            if len(moves_in_scan) > 0:
                innermost_moves.extend(moves_in_scan)
                moves_in_scan = []

        for empty_site in innermost_empty_traps_Cs:
            moves_in_scan = []
            # Update the source atoms list based on the source type
            if source_types == 'M':
                source_atoms = array_info['Type_M_Cs']
            if source_types == 'R':
                source_atoms = array_info['Type_R_Cs']

            # If there is no type-M atom, continue to the next empty site
            if len(source_atoms) == 0:
                break
            else:
                nearest_atom = min(source_atoms, key=lambda atom: abs(atom[0] - empty_site[0]) + abs(atom[1] - empty_site[1]))
                source_atoms.sort(key = lambda atom: (atom[0] - empty_site[0])**2 + (atom[1] - empty_site[1])**2)
                nearest_atoms = copy.deepcopy(source_atoms)

            path_found_flag = False
            while not path_found_flag:
                try:
                    nearest_atom = nearest_atoms[0]
                    path, path_found_flag = bfs_find_path(innermost_arrays, nearest_atom, empty_site)
                except IndexError:
                    break
                # Transform path to move objects
                if path_found_flag:
                    path = flatten_tuple(path)
                    path = path[::-1]
                    for item in path:
                        current_pos = item[0] 
                        for next_pos in item:
                            if np.array_equal(current_pos, next_pos):
                                pass
                            else:
                                moves_in_scan.append(Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1]))
                                current_pos = next_pos

                else:
                    nearest_atoms.remove(nearest_atom)

            # Implement the moves and update the arrays
            if len(moves_in_scan) > 0:
                # Apply the moves to the arrays
                innermost_moves.extend(moves_in_scan)
                _ = innermost_arrays.move_atoms(moves_in_scan)
                moves_in_scan = []

                if source_types == 'M':
                    array_info['Type_M_Cs'].remove(nearest_atom)
                    array_info['Cs_emp'].remove(list(empty_site))
                    array_info['Rb_emp'].append(nearest_atom)
                    array_info['Cs_fin'].append(list(empty_site))
                    
                # Because the source atom is R-type, we do not need to append it into the Cs_empty_coordinate
                if source_types == 'R':
                    array_info['Type_R_Cs'].remove(nearest_atom)
                    array_info['Cs_emp'].remove(list(empty_site))
                    array_info['Cs_fin'].append(list(empty_site))

             # If we use M-type atom, find a new M-type atom to fill new empty site
            if source_types == 'M':
                innermost_arrays, moves_in_scan, array_info = refill_source_atoms(innermost_arrays, nearest_atom, array_info, 'Rb')
            
            if len(moves_in_scan) > 0:
                innermost_moves.extend(moves_in_scan)
                moves_in_scan = []

        # Update the control sequence length
        if source_types == 'M':
            control_sequence_len = len(array_info['Type_M_Rb']+array_info['Type_M_Cs'])
        if source_types == 'R':
            control_sequence_len = len(array_info['Rb_emp'] + array_info['Cs_emp'])
    
    return innermost_arrays, innermost_moves, array_info


def refill_source_atoms(arrays: AtomArray, refill_site: list, array_info: dict, needed_species: str): # Needed type is the type of the empty site
    refill_moves = []
    empty_site = refill_site

    if needed_species == 'Rb':
        source_atoms = array_info['Type_M_Rb']
    if needed_species == 'Cs':
        source_atoms = array_info['Type_M_Cs']

    # If there is no type-M atom, continue to the next empty site
    if len(source_atoms) == 0:
        return arrays, [], array_info
    else:
        nearest_atom = min(source_atoms, key=lambda atom: abs(atom[0] - empty_site[0]) + abs(atom[1] - empty_site[1]))
    
    path, path_found_flag = bfs_find_path(arrays, nearest_atom, refill_site)

    if path_found_flag:
        # Transform path to move objects
        path = flatten_tuple(path)
        path = path[::-1]
        for item in path:
            current_pos = item[0] 
            for next_pos in item:
                if np.array_equal(current_pos, next_pos):
                    pass
                else:
                    refill_moves.append(Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1]))
                    current_pos = next_pos
    
    # Implement the moves and update the arrays
    if len(refill_moves) > 0:
        # Apply the moves to the arrays
        _ = arrays.move_atoms(refill_moves)

        if needed_species == 'Rb':
            array_info['Type_M_Rb'].remove(nearest_atom)
            array_info['Rb_emp'].remove(list(empty_site))
            array_info['Cs_emp'].append(nearest_atom)
            array_info['Rb_fin'].append(list(empty_site))
        elif needed_species == 'Cs':
            array_info['Type_M_Cs'].remove(nearest_atom)
            array_info['Cs_emp'].remove(list(empty_site))
            array_info['Rb_emp'].append(nearest_atom)
            array_info['Cs_fin'].append(list(empty_site))
    
    return arrays, refill_moves, array_info

def remove_M_atoms(arrays: AtomArray, array_info: dict):
    """
    Remove remaining M-type atoms from the arrays.
    """
    remove_M_moves = []
    outer_M = []
    isolated_M_inside = []

    # 1. Identify remaining M-type atoms
    for M_type_atom in (array_info['Type_M_Rb']+array_info['Type_M_Cs']):
        if M_type_atom in array_info['Type_M_Rb']:
            outer_M.append([M_type_atom, 'Rb'])
        elif M_type_atom in array_info['Type_M_Cs']:
            outer_M.append([M_type_atom, 'Cs'])

    # 2. Move out the outermost M-type atoms (from out to in)
    if len(outer_M) > 0:
        outer_M.sort(key = lambda atom: (atom[0][0] - (array_info['top']+array_info['bot'])/2)**2 + (atom[0][1] - (array_info['left']+array_info['right'])/2)**2)
        outer_M = outer_M[::-1]

    # Find nearest available target empty site to put obstacles (empty; not in the path)
    emp_targets = []
    for i in range(len(arrays.matrix)):
        for j in range(len(arrays.matrix[0])):
            if not ((array_info['top'] <= i <= array_info['bot']) and (array_info['left'] <= j <= array_info['right'])):
                if arrays.matrix[i][j][0] == 0 and arrays.matrix[i][j][1] == 0 and (i,j):
                    emp_targets.append([i, j])

    for [atom, atom_type] in outer_M:
        moves_in_scan = []
        # Find the nearest empty site to the atom
        emp_targets.sort(key = lambda emp: (emp[0] - atom[0])**2 + (emp[1] - atom[1])**2)

        for emp_target in emp_targets:
            path, path_found_flag = bfs_find_path(arrays, atom, emp_target)

            if path_found_flag and ((emp_target[0] - (array_info['top']+array_info['bot'])/2)**2 + (emp_target[1] - (array_info['top']+array_info['bot'])/2)**2 - (atom[0] - (array_info['top']+array_info['bot'])/2)**2 - (atom[1] - (array_info['left']+array_info['right'])/2)**2) > 3:
                path = flatten_tuple(path)
                path = path[::-1]
                for item in path:
                    current_pos = item[0] 
                    for next_pos in item:
                        if np.array_equal(current_pos, next_pos):
                            pass
                        else:
                            moves_in_scan.append(Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1]))
                            current_pos = next_pos
                empty_site = emp_target
                break

        # Implement the moves and update the arrays
        if len(moves_in_scan) > 0 and path_found_flag:
            # Apply the moves to the arrays
            remove_M_moves.extend(moves_in_scan)
            _ = arrays.move_atoms(moves_in_scan)

            if atom_type == 'Cs':
                array_info['Type_M_Cs'].remove(atom)
                array_info['Rb_emp'].append(list(empty_site))
                array_info['Type_R_Cs'].append(list(empty_site))
                
            elif atom_type == 'Rb':
                array_info['Type_M_Rb'].remove(atom)
                array_info['Cs_emp'].append(atom)
                array_info['Type_R_Rb'].append(list(empty_site))

            outer_M.remove([atom, atom_type])
                
    return arrays, remove_M_moves, array_info

def bfs_find_path(arrays: AtomArray, start: list, end: list):
    # Initialize the queue
    start = tuple(start)
    end = tuple(end)
    queue = deque([(start, [(start[0], start[1])])])
    visited = set()
    path_found_flag = False

    # Mark the start node as visited
    visited.add(start)
    directions = [(0,1), (0,-1), (1,0), (-1,0)]

    # Storing the minimum time to reach at destination
    time = 0
    minTime = maxsize

    while queue:
        # Get the current node position for queue
        current_node, path = queue.popleft()

        if time >= minTime:
            return [], path_found_flag

        if current_node == end:
            minTime = min(time, minTime)
            queue.append((current_node, path))
            path_found_flag = True
            return path, path_found_flag
                
        # Visit the neighbors of the current node
        for direction in directions:
            dr, dc = direction[0], direction[1]
            
            # Check if the coordinates are within the matrix
            if current_node[0]+dr < 0 or current_node[0]+dr >= len(arrays.matrix) or current_node[1]+dc < 0 or current_node[1]+dc >= len(arrays.matrix[0]):
                continue
            else:
                # Check if the node is visited and there is obstacle
                if arrays.matrix[current_node[0]+dr][current_node[1]+dc][0] == 0 and arrays.matrix[current_node[0]+dr][current_node[1]+dc][1]==0 and (current_node[0]+dr, current_node[1]+dc) not in visited:
                    new_node = (current_node[0]+dr, current_node[1]+dc)
                    visited.add((current_node[0]+dr, current_node[1]+dc))
                    queue.append((new_node, path + [new_node]))
        time += 1
    
    # No feasible path found
    return [], path_found_flag

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

# Find a path between R-type atoms and outermost unfinished sites; Output: path, obstacle coordinates
def min_obstacle_path(arrays: AtomArray, start: list, end: list, array_info: dict):
    start = tuple(start)
    end = tuple(end)
    min_obstacle = maxsize
    queue = deque([(start, [(start[0], start[1])], [])])
    visited = set()
    path_found_flag = False

    # Mark the start node as visited
    visited.add(start)
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    while queue:
        current_node, path, obstacle_list = queue.popleft()
        obstacle_count = len(obstacle_list)

        if current_node == end:
            if obstacle_count < min_obstacle:
                min_obstacle = obstacle_count
                best_path = path
                best_obstacle_list = obstacle_list
                path_found_flag = True
                visited.remove(end)
            continue

        for direction in directions:
            dr, dc = direction
            new_node = (current_node[0] + dr, current_node[1] + dc)

            if (new_node[0] < 0 or new_node[0] >= len(arrays.matrix) or
                    new_node[1] < 0 or new_node[1] >= len(arrays.matrix[0])):        
                continue
            else:
                if array_info['left'] <= new_node[0] <= array_info['right'] and array_info['top'] <= new_node[1] <= array_info['bot']:
                    if new_node[0] == end[0] or new_node[1] == end[1] and np.sum(arrays.matrix[new_node[0],new_node[1],:]) == 0:
                        pass
                    else:
                        continue
            
            if new_node not in visited:
                visited.add(new_node)
                new_obstacle_list = obstacle_list.copy()
                if arrays.matrix[new_node[0]][new_node[1]][0] == 1 or arrays.matrix[new_node[0]][new_node[1]][1] == 1:
                    new_obstacle_list.append(new_node)
                queue.append((new_node, path + [new_node], new_obstacle_list))

    if path_found_flag:
        return best_path, path_found_flag, best_obstacle_list
    else:
        return [], path_found_flag, []

# Put the move list and path here, return a new move list according to path
def transform_path_to_move_object(path: list, moves: list) -> list:
    path = flatten_tuple(path)
    path = path[::-1]

    # Transform path into move object
    for item in path:
        current_pos = item[0]
        for next_pos in item:
            if np.array_equal(current_pos, next_pos):
                pass
            else:
                moves.append(Move(current_pos[0], current_pos[1], next_pos[0], next_pos[1]))
                current_pos = next_pos
    return moves