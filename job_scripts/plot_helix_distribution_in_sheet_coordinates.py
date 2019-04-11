#!/usr/bin/env python3
'''Plot the distribution of helices in the beta sheet coordinates.
'''

import numpy as np
import matplotlib.pyplot as plt

import pyrosetta
from pyrosetta import rosetta

def xyz_to_np_array(xyz):
    '''Convert an xyz vector to a numpy array.'''
    return np.array([xyz.x, xyz.y, xyz.z])

def get_sheet_ca_positions(pose, sheet_residues):
    '''Get CA positions of a beta sheet'''
    sheet_ca_positions = []

    for strand in sheet_residues:
        strand_ca_positions = []

        for seqpos in strand:
            if not (seqpos is None):
                xyz = pose.residue(seqpos).xyz('CA')
                strand_ca_positions.append(xyz_to_np_array(xyz))
            else:
                strand_ca_positions.append(None)

        sheet_ca_positions.append(strand_ca_positions)

    return sheet_ca_positions

def get_sheet_residue_local_frames(pose, sheet_residues, strand_directions, z_ref_residue):
    '''Get local frames for each beta sheet residues.
    The x direction is determined by the strand direction.
    The y direction is determined by neighbor strands.
    The z direction is determined by positions of CA and CB atoms.
    '''
    sheet_res_frames = []
    
    for i, strand in enumerate(sheet_residues):
        strand_res_frames = []

        for j, seqpos in enumerate(strand):
            if not (seqpos is None):
                # Get the Z direction
                
                ca_xyz = pose.residue(seqpos).xyz('CA')
            
                if pose.residue(seqpos).name3() == 'GLY':
                    cb_xyz = pose.residue(seqpos).xyz('2HA')
                else: 
                    cb_xyz = pose.residue(seqpos).xyz('CB')
           
                z_axis_no_sign = xyz_to_np_array((cb_xyz - ca_xyz).normalized())
                z_sign = 1 if (0 == (j - z_ref_residue) %2) else -1 
                z_axis = z_sign * z_axis_no_sign

                # Get the Y direction

                y_axis_no_sign = xyz_to_np_array((pose.residue(seqpos).xyz('C') 
                           - pose.residue(seqpos).xyz('N')).normalized())
                
                y_sign = 1 if strand_directions[i] else -1
                y_axis_before_correction = y_sign * y_axis_no_sign
                y_axis_before_correction = y_axis_before_correction - np.dot(y_axis_before_correction, z_axis) * z_axis
        
                y_axis = y_axis_before_correction / np.linalg.norm(y_axis_before_correction) 

                # Get the X direction

                x_axis = np.cross(y_axis, z_axis)

                strand_res_frames.append([x_axis, y_axis, z_axis])

            else:
                strand_res_frames.append(None)

        sheet_res_frames.append(strand_res_frames)
    
    return sheet_res_frames

def project_point_to_sheet(sheet_ca_positions, sheet_res_frames, point):
    '''Project a point onto the sheet coordinates.
    Use a gaussian function to weight the contributations
    of different local frames.
    '''
    ref_vertical_length = 3.3
    ref_horizontal_length = 4.6

    ca_p_linear = [ca for strand in sheet_ca_positions for ca in strand if not (ca is None)]
    
    res_frames_linear = [np.array(f) for strand in sheet_res_frames for f in strand if not (f is None)]
    
    sheet_coords_linear = [np.array([i * ref_horizontal_length, j * ref_vertical_length])
            for i in range(len(sheet_ca_positions)) for j in range(len(sheet_ca_positions[0]))
            if not (sheet_ca_positions[i][j] is None)]

    # Assign weights to sheet residues

    dists = [np.linalg.norm(point - ca) for ca in ca_p_linear]
    dist_scale = 5 ** 2
    weights = [np.exp(-d * d / dist_scale) for d in dists]
    total_weight = sum(weights)
    weights = [w / total_weight for w in weights]

    # Get the average positions

    mean_position = sum(weights[i] * ca_p_linear[i] for i in range(len(ca_p_linear)))
    mean_sheet_coord = sum(weights[i] * sheet_coords_linear[i] for i in range(len(sheet_coords_linear)))
    mean_frame = sum(weights[i] * res_frames_linear[i] for i in range(len(res_frames_linear)))

    # Get the projected sheet coordinates

    diff = point - mean_position
    
    x = mean_sheet_coord[0] + np.dot(diff, mean_frame[0])
    y = mean_sheet_coord[1] + np.dot(diff, mean_frame[1])

    return np.array([x, y])

def plot_the_underlying_sheet(sheet_ca_positions, sheet_res_frames):
    '''Plot the underlying sheet.'''
    X = []
    Y = []

    for i in range(len(sheet_ca_positions)):
        for j in range(len(sheet_ca_positions[i])):
            
            if sheet_ca_positions[i][j] is None:
                continue

            sheet_coord = project_point_to_sheet(sheet_ca_positions, sheet_res_frames, sheet_ca_positions[i][j])
            
            X.append(sheet_coord[0])
            Y.append(sheet_coord[1])

    plt.scatter(X, Y)
    plt.show()

def plot_test(sheet_ca_positions, sheet_res_frames, points):
    X = []
    Y = []

    for p in points:
        sheet_coord = project_point_to_sheet(sheet_ca_positions, sheet_res_frames, p)
        X.append(sheet_coord[0])
        Y.append(sheet_coord[1])

    plt.scatter(X, Y)
    #plt.show()


if __name__ == '__main__':

    pyrosetta.init()

    pose = rosetta.core.import_pose.pose_from_file('/home/xingjie/Softwares/scripts/loop_helix_loop_reshaping/data/test_screen_compatible_loop_helix_loop_units/model_104.pdb.gz')

    # Define residues in a beta sheet. The residues should be aligned.
    # Residues that cannot be aligned should be set to None.

    sheet_residues = [
            [None,None,None,27,28,29,30,None],
            [None,1,2,3,4,5,6,7],
            [57,58,59,60,61,62,63,64],
            [None,None,None,86,87,88,89,90],
            ]
    strand_directions = [True, True, True, True]
    
    # The residue on a strand that determines the Z-direction
    # This value must be set correctly such that as the indices
    # increase, the projected x, y coordinates also increase.
    z_ref_residue = 3 

    sheet_ca_positions = get_sheet_ca_positions(pose, sheet_residues)

    print(sheet_ca_positions)

    sheet_res_frames = get_sheet_residue_local_frames(pose, sheet_residues, strand_directions, z_ref_residue)

    print(sheet_res_frames)

    ca_points = [xyz_to_np_array(pose.residue(i).xyz('CA')) for i in range(1, pose.size() + 1)]
    #ca_points = [xyz_to_np_array(pose.residue(i).xyz('CA')) for i in range(69, 82)]
    #ca_points = [xyz_to_np_array(pose.residue(i).xyz('CA')) for i in [69, 73, 77, 81]]
    #ca_points = [np.array([i, i, i]) for i in range(10)]

    plot_test(sheet_ca_positions, sheet_res_frames, ca_points)
    plot_the_underlying_sheet(sheet_ca_positions, sheet_res_frames)
