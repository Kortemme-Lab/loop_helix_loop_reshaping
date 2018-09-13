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

def screen_compatible_loop_helix_loop_units(output_dir, pose, insertion_points, lhl_units, num_jobs, job_id):
    '''Get LHL units that are compatible with each other.
    Dump pdb files of the compatible structures and a json
    file for the insertion points.
    '''
    
    def get_all_combinations(num_units_at_each_position, current_id):
        '''Get all combinations recursively.
        Return a list of combinations.
        '''
        if current_id == len(num_units_at_each_position) - 1:
            return [[i] for i in range(num_units_at_each_position[current_id])]
        else:
            low_level_combinations = get_all_combinations(num_units_at_each_position, current_id + 1)
            combinations = []
            
            for i in range(num_units_at_each_position[current_id]):
                for c in low_level_combinations:
                    combinations.append([i] + c)

            return combinations

    # Get the indices for all combinations of LHL units

    num_units_at_each_position = [len(x) for x in lhl_units]
    combinations = get_all_combinations(num_units_at_each_position, 0)

    print('There are {0} combinations of lhl units.'.format(len(combinations)))

    # Screen the LHL units

    for i, c in enumerate(combinations):

        # Apply the LHL unit

        for insertion_id in range(len(c)):
            insert_loop_helix_loop_unit(pose, insertion_points, insertion_id, lhl_units[insertion_id][c[insertion_id]])
       
        # Check clashes

        if not check_clashes_between_lhl_units(pose, insertion_points):
            continue
        
        # Dump the outputs

        pose.dump_pdb(os.path.join(output_dir, 'model_{0}.pdb.gz'.format(i)))
        
        with open(os.path.join(output_dir, 'insertion_points_{0}.json'.format(i)), 'w') as f:
            json.dump(insertion_points, f)
