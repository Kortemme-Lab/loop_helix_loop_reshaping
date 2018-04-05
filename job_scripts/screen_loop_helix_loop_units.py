#!/usr/bin/env python3

import os
import sys
import json

import pyrosetta
from pyrosetta import rosetta

import loop_helix_loop_reshaping as LHLR

def screen(data_path, input_path, input_pdb, lhl_start, lhl_stop, num_jobs, job_id):
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)
  
    front_linker_dbs = []
    back_linker_dbs = []

    for f in os.listdir(input_path):
        if f.endswith('_front.json'):
            with open(os.path.join(input_path, f), 'r') as f:
                front_linker_dbs.append(json.load(f))

        elif f.endswith('_back.json'):
            with open(os.path.join(input_path, f), 'r') as f:
                back_linker_dbs.append(json.load(f))

    LHLR.build_loop_helix_loop_unit.screen_all_loop_helix_loop_units(data_path, pose, lhl_start, lhl_stop, 
            front_linker_dbs, back_linker_dbs, num_jobs, job_id)


if __name__ == '__main__':
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1
   
    pyrosetta.init(options='-extra_res_fa test_inputs/CRO.params -mute all')

    #screen(data_path, 'data/test_linker_selection', 'test_inputs/sfgfp_modified_from_5dph_cleaned.pdb',
    #        124, 138, num_jobs, job_id)

    screen(data_path, 'data/test_loop_selection_1abe', 'test_inputs/1abe_cleaned.pdb',
            141, 167, num_jobs, job_id)
