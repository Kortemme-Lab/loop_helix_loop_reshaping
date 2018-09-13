#!/usr/bin/env python3

import os
import sys
import json

import pyrosetta
from pyrosetta import rosetta

import loop_helix_loop_reshaping as LHLR


def screen_compatible_loop_helix_loop_units(data_path, lhl_units_path, input_pdb, input_insertion_points_file, num_jobs, job_id):
    
    # Load insertion points
    
    with open(input_insertion_points_file, 'r') as f:
        insertion_points = json.load(f)
    
    # Load the lhl_units
    
    lhl_units = []

    for i in range(len(insertion_points)):
        lhl_units.append([])

        for lf in os.listdir(lhl_units_path):
            if lf.startswith('selected_lhl_units_{0}'.format(i)):
                with open(os.path.join(lhl_units_path, lf), 'r') as f:
                    lhl_units[i] += json.load(f)
   
    # Load the pose
    
    pose = rosetta.core.import_pose.pose_from_file(input_pdb)

    LHLR.screen_compatible_loop_helix_loop_units.screen_compatible_loop_helix_loop_units(
            data_path, pose, insertion_points, lhl_units, num_jobs, job_id)
    


if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1
   
    pyrosetta.init(options='-mute all')

    lhl_units_path = 'data/test_screen_single_insertion_loop_helix_loop_units'
    input_pdb = 'test_inputs/2lv8_cleaned.pdb'
    input_insertion_points_file = 'test_inputs/2lv8_insertion_points.json'

    screen_compatible_loop_helix_loop_units(data_path, lhl_units_path, 
            input_pdb, input_insertion_points_file, num_jobs, job_id)

