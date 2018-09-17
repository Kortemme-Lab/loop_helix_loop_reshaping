import os
import json

import numpy as np

import pyrosetta
from pyrosetta import rosetta

from loop_helix_loop_reshaping import simple_pose_moves
from loop_helix_loop_reshaping import pose_analysis

def insert_loop_helix_loop_unit(pose, insertion_points, insertion_id, lhl_unit):
    '''Insert a loop helix loop unit to a given insertion point.
    The lhl_unit is dictionary with phis, psis, omegas and cutpoint.
    '''
    # Run the terminal version of the function if the insertion point is on the terminal
    
    if insertion_points[insertion_id]['start'] == 0 or insertion_points[insertion_id]['stop'] == pose.size() + 1:
        insert_terminal_loop_helix_loop_unit(pose, insertion_points, insertion_id, lhl_unit)
        return

    # Remove the existing residues within the insertion point

    simple_pose_moves.remove_insertion_residues(pose, insertion_points, [insertion_id])

    # Add the residues

    n_pre_cut_residues = lhl_unit['cutpoint'] - 1
    n_post_cut_residues = len(lhl_unit['phis']) - 2 - n_pre_cut_residues

    simple_pose_moves.insert_alas(pose, insertion_points[insertion_id]['stop'], n_post_cut_residues, insert_after=False, reset_fold_tree=True)
    simple_pose_moves.insert_alas(pose, insertion_points[insertion_id]['start'], n_pre_cut_residues, insert_after=True, reset_fold_tree=True)

    # Apply the torsions

    simple_pose_moves.apply_linker(pose, lhl_unit, insertion_points[insertion_id]['start'])
    simple_pose_moves.update_insertion_points(insertion_points, insertion_id, len(lhl_unit['phis']) - 2 )

def insert_terminal_loop_helix_loop_unit(pose, insertion_points, insertion_id, lhl_unit):
    '''Insert a loop helix loop unit to a terminal of a protein and
    trim the overhanging loop.
    '''
    # Remove the existing residues within the insertion point

    simple_pose_moves.remove_insertion_residues(pose, insertion_points, [insertion_id])

    # Insert residues and apply the torsions 

    lhl_length = len(lhl_unit['phis']) - 2

    if insertion_points[insertion_id]['start'] == 0: # N-term
        simple_pose_moves.insert_alas(pose, 1, lhl_length + 1, insert_after=False, reset_fold_tree=True, fold_tree_root=pose.size())
        simple_pose_moves.apply_linker(pose, lhl_unit, 1)
        
        # Remove the one extra residue which was added for torsion application
        simple_pose_moves.delete_region(pose, 1, 1)

        # Remove the overhanging loop
        dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()

        overhanging_loop_length = 0
        for i in range(lhl_length):
            if dssp_str[i] != 'H':
                overhanging_loop_length += 1
            else:
                break
       
        if overhanging_loop_length > 0:
            simple_pose_moves.delete_region(pose, 1, overhanging_loop_length)
        simple_pose_moves.update_insertion_points(insertion_points, insertion_id, lhl_length - overhanging_loop_length)

    elif insertion_points[insertion_id]['stop'] == pose.size() + 1: # C-term
        simple_pose_moves.insert_alas(pose, pose.size(), lhl_length + 1, insert_after=True, reset_fold_tree=True, fold_tree_root=1) 
        simple_pose_moves.apply_linker(pose, lhl_unit, insertion_points[insertion_id]['start'])

        # Remove the one extra residue which was added for torsion application
        simple_pose_moves.delete_region(pose, pose.size(), pose.size())

        # Remove the overhanging loop
        dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()

        overhanging_loop_length = 0
        for i in range(lhl_length):
            if dssp_str[-1 - i] != 'H':
                overhanging_loop_length += 1
            else:
                break

        if overhanging_loop_length > 0:
            simple_pose_moves.delete_region(pose, pose.size() - overhanging_loop_length + 1, pose.size())
        simple_pose_moves.update_insertion_points(insertion_points, insertion_id, lhl_length - overhanging_loop_length)

def check_clashes_between_lhl_units(pose, insertion_points):
    '''Check the clashes between LHL units.
    Return True if there are no clashes between inserted LHL units.
    '''
    # Mutate all inserted residues to VAL

    inserted_residues = [i for ip in insertion_points for i in range(ip['start'] + 1, ip['stop'])]
    simple_pose_moves.mutate_residues(pose, inserted_residues, ['VAL'] * len(inserted_residues))

    # Check clashes between LHL units

    for i in range(len(insertion_points)):
        residues_i = [k for k in range(insertion_points[i]['start'] + 1, insertion_points[i]['stop'])]
        
        for j in range(i + 1, len(insertion_points)):
            residues_j = [k for k in range(insertion_points[j]['start'] + 1, insertion_points[j]['stop'])]

            if not pose_analysis.check_clashes_between_groups(pose, residues_i, residues_j):
                return False

    return True

def screen_compatible_loop_helix_loop_units(output_dir, pose, insertion_points, lhl_units, num_jobs, job_id, symmetric_lists=None, max_num_to_screen=None):
    '''Get LHL units that are compatible with each other.
    Dump pdb files of the compatible structures and a json
    file for the insertion points.
    Args:
        symmetric_lists : A list of list. Insertion_points_ids within an inner list will be symmetric.
    '''
    def index_to_combination(num_units_for_each_sym_list, index):
        '''Convert an index of a combination to the combination.'''
        combination = []

        for n in num_units_for_each_sym_list:
            combination.append(index % n)
            index = index // n 

        return combination

    # Update the symmetric_lists

    if symmetric_lists is None:
        symmetric_lists = [[i] for i in range(len(insertion_points))]

    # Get the number of LHL units for each symmetric list

    num_units_for_each_sym_list = [len(lhl_units[i]) for i in range(len(symmetric_lists))]

    num_all_combinations = 1
    for n in num_units_for_each_sym_list:
        num_all_combinations *= n

    print('There are {0} combinations of lhl units.'.format(num_all_combinations))

    # Update the max number of combinations to screen

    if max_num_to_screen:
        max_num_to_screen = min(max_num_to_screen, num_all_combinations)
    else:
        max_num_to_screen = num_all_combinations

    # Screen the LHL units

    for i in range(max_num_to_screen):
        if i % num_jobs != job_id:
            continue
        
        c = index_to_combination(num_units_for_each_sym_list, i)

        # Apply the LHL unit

        for sym_list_id in range(len(c)):
            
            for insertion_id in symmetric_lists[sym_list_id]:
            
                insert_loop_helix_loop_unit(pose, insertion_points, insertion_id, lhl_units[sym_list_id][c[sym_list_id]])
       
        # Check clashes

        if not check_clashes_between_lhl_units(pose, insertion_points):
            continue
        
        # Dump the outputs

        pose.dump_pdb(os.path.join(output_dir, 'model_{0}.pdb.gz'.format(i)))
        
        with open(os.path.join(output_dir, 'insertion_points_{0}.json'.format(i)), 'w') as f:
            json.dump(insertion_points, f)
