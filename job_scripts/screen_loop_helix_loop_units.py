#!/usr/bin/env python2

import os
import sys
import json

import pyrosetta
from pyrosetta import rosetta

import loop_helix_loop_reshaping as LHLR


if __name__ == '__main__':

    pyrosetta.init(options='-extra_res_fa test_inputs/CRO.params -mute all')

    pose = rosetta.core.pose.Pose()
    input_pdb ='test_inputs/sfgfp_modified_from_5dph_cleaned.pdb'
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)
   
    with open('data/test_linker_selection/selected_linkers_124_138_front.json', 'r') as f:
        front_linker_dbs = [json.load(f)]

    with open('data/test_linker_selection/selected_linkers_124_138_back.json', 'r') as f:
        back_linker_dbs = [json.load(f)]

    LHLR.build_loop_helix_loop_unit.screen_all_loop_helix_loop_units('debug', pose, 124, 138, front_linker_dbs, back_linker_dbs, 'debug/test.json')
