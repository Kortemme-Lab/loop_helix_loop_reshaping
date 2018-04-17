#!/usr/bin/env python3
'''Generate the pos files required for fuzzball matching.
Usage:
    ./generate_pos_files_for_matching.py path_to_pdbs
'''
import sys
import os

import pyrosetta
from pyrosetta import rosetta


def update_fixed_region_positions(fixed_region_positions, reshaped_region_stop, removed_region_stop):
    '''Update the fixed region positions due to the change of length.'''
    diff = reshaped_region_stop - removed_region_stop
    
    return [i if i < removed_region_stop else i + diff for i in fixed_region_positions]

def find_reshaped_residues_that_point_to_fixed_positions(pose, fixed_region_positions, reshaped_region_start, reshaped_region_stop):
    '''Find reshaped residues that point to the fixed region
    positions.'''
    reshaped_positions = set()
    
    for i in range(reshaped_region_start, reshaped_region_stop + 1):
        if pose.residue(i).name3() == 'GLY':continue
        
        v1 = pose.residue(i).xyz('CB') - pose.residue(i).xyz('CA')
        for j in fixed_region_positions:
            v2 = pose.residue(j).nbr_atom_xyz() - pose.residue(i).xyz('CA')
   
            if v1.dot(v2) > 0:
                reshaped_positions.add(i)
                break

    return list(reshaped_positions)

def print_pymol_selection_for_residues(pose, residues):
    '''Print the pymol selection command for the residues.'''
    res_commands = ['(c. {0} and res {1})'.format(pose.pdb_info().chain(i), pose.pdb_info().number(i))
                    for i in residues]

    print('sele ' + ' or '.join(res_commands))

def generage_pose_file_for_1abe(input_pdb_path):
    pdb_name = os.path.basename(input_pdb_path)
    pdb_name = pdb_name[:-4] if pdb_name.endswith('.pdb') else pdb_name[:-7]
    pdb_name_splited = pdb_name.split('_')

    reshaped_region_start = 141  
    reshaped_region_stop = int(pdb_name_splited[-1]) - 1
    removed_region_stop = 167

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb_path)

    fixed_region_positions = [109, 110, 113, 117, 121, 124, 129, 131, 134, 136, 138, 140, 168, 170, 172, 199, 201, 203]

    new_fixed_region_positions = update_fixed_region_positions(fixed_region_positions,
            reshaped_region_stop, removed_region_stop)

    reshaped_positions = find_reshaped_residues_that_point_to_fixed_positions(pose, new_fixed_region_positions,
            reshaped_region_start, reshaped_region_stop)

    pos_residues = [str(i) for i in new_fixed_region_positions + reshaped_positions]

    output_pos_file_path = os.path.join(os.path.dirname(input_pdb_path), pdb_name + '.pos')
    with open(output_pos_file_path, 'w') as f:
        f.write(' '.join(pos_residues))

if __name__ == '__main__':
    path_to_pdbs = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1
   
    pyrosetta.init()
    
    target_files = [f for f in os.listdir(path_to_pdbs) if f.endswith('.pdb') or f.endswith('.pdb.gz')]
            
    for i, f in enumerate(target_files):
        if i % num_jobs == job_id:
            generage_pose_file_for_1abe(os.path.join(path_to_pdbs, f))
