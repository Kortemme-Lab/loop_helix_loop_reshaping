import os
import numpy as np

import pyrosetta
from pyrosetta import rosetta

from loop_helix_loop_reshaping import simple_pose_moves
from loop_helix_loop_reshaping import constraint
from loop_helix_loop_reshaping import pose_analysis
from loop_helix_loop_reshaping import filters


def apply_ideal_helix(pose, start, helix_length):
    '''Apply an ideal helix.'''
    for i in range(start, start + helix_length):
        pose.set_phi(i, -57)
        pose.set_psi(i, -47)
        pose.set_omega(i, 180)

def prepare_pose_for_lhl_screen(pose, lhl_start, lhl_stop, front_linker_length, back_linker_length):
    '''Trim the pose for lhl screen.
    A 10 residue helix will be appended to each of the linker.
    Set up the fold tree properly.
    '''
    # Remove the LHL region

    simple_pose_moves.delete_region(pose, lhl_start + 1, lhl_stop - 1)

    # Add the linker

    simple_pose_moves.insert_alas(pose, lhl_start + 1, back_linker_length + 10, insert_after=False, reset_fold_tree=True)
    simple_pose_moves.insert_alas(pose, lhl_start, front_linker_length + 10, insert_after=True, reset_fold_tree=True)
    rosetta.core.pose.correctly_add_cutpoint_variants(pose)

    # Set the appended residues to be helices

    apply_ideal_helix(pose, lhl_start + front_linker_length + 1, 20)

def residue_bb_rmsd_with_cutoff(pose, res1, res2, ca_cutoff=3):
    '''Calculate the backbone RMSD between two residues.
    If the CA distance is above a cutoff. Return infinity.
    '''
    if pose.residue(res1).xyz('CA').distance(pose.residue(res2).xyz('CA')) > ca_cutoff:
        return float('inf')

    bb_atms = ['N', 'CA', 'C']
    
    P1 = [pose.residue(res1).xyz(a) for a in bb_atms]
    P2 = [pose.residue(res2).xyz(a) for a in bb_atms]

    return np.sqrt(sum(P1[i].distance_squared(P2[i]) for i in range(len(P1))) / len(P1))

def residue_direction_vector(pose, res):
    '''Return the normalized vector pointing from N to C.'''
    return (pose.residue(res).xyz('C') - pose.residue(res).xyz('N')).normalized()

def close_helix_by_minimization(pose, movable_region_start, movable_region_end, helix_start, helix_end):
    '''Close a gap inside a helix by minimization.
    Return true if the gap could be closed.
    '''
    # Make a clone of poly ALA pose for minimization

    simple_pose_moves.mutate_pose_to_single_AA(pose, 'ALA')
    rosetta.core.pose.correctly_add_cutpoint_variants(pose)

    # Set hydrogen bond constraints for the helix

    pose.constraint_set().clear()
    helix_hbs = [(i + 4, i) for i in range(helix_start, helix_end - 3)]
    constraint.add_constraints_to_pose(pose, constraint.get_bb_hbond_constraint(helix_hbs))

    # Set score function

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn.set_weight(rosetta.core.scoring.base_pair_constraint, 1) #H-bond constraint

    # Set movemap

    mm = rosetta.core.kinematics.MoveMap()
    
    for i in range(movable_region_start, movable_region_end + 1):
        mm.set_bb(i, True)
   
    # Set the minimization mover

    min_opts = rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True )
    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.movemap(mm)
    min_mover.min_options(min_opts)

    # Close the chain

    for chainbreak_weight in [0.5, 1, 5, 10]:
        sfxn.set_weight(rosetta.core.scoring.chainbreak, chainbreak_weight)
        min_mover.score_function(sfxn)
        min_mover.apply(pose)

    chainbreak_energy = pose.energies().total_energies()[rosetta.core.scoring.chainbreak] 
    if chainbreak_energy > 0.2:
        return False

    # Minimize without constraints
    
    sfxn.set_weight(rosetta.core.scoring.base_pair_constraint, 0)
    min_mover.score_function(sfxn)
    min_mover.apply(pose)
    
    return True

def trim_helix_and_connect(original_pose, movable_region_start, movable_region_end, helix_start, helix_end, trim_start, trim_end, num_res_clashes_tolerance=0):
    '''Trim a part of a helix and connect the rest of the helix.
    Return true if a good connection could be made.
    '''
    # Trim the helix
    
    simple_pose_moves.delete_region(original_pose, trim_start, trim_end) 
    rosetta.core.pose.correctly_add_cutpoint_variants(original_pose)

    movable_region_end -= trim_end - trim_start + 1
    helix_end -= trim_end - trim_start + 1
  
    # Close the gap

    pose = original_pose.clone()
    
    if not close_helix_by_minimization(pose, movable_region_start, movable_region_end, helix_start, helix_end):
        return False
    
    # Check the secondary structure

    dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()
    for ss in dssp_str[helix_start - 1: helix_end]:
        if ss != 'H': return False

    # Check helix hbond scores

    helix_hb_scores = filters.get_helix_hbond_scores(pose, helix_start, helix_end)
    if len(helix_hb_scores) == 0 or max(helix_hb_scores) > -0.7:
        return False 

    # Check contact degrees of the helix

    helix_residues = list(range(helix_start, helix_end + 1))
    movable_residues = list(range(movable_region_start, movable_region_end + 1))
    not_movable_residues = list(range(1, movable_region_start)) + list(range(movable_region_end + 1, pose.size() + 1))
    
    contact_degrees = pose_analysis.contact_degrees_of_residues(pose, helix_residues, not_movable_residues)
    if np.median(contact_degrees) < 1:
        return False

    # Check clashes

    simple_pose_moves.mutate_residues(pose, helix_residues + not_movable_residues, ['VAL'] * len(helix_residues + not_movable_residues))

    num_res_clashes = pose_analysis.count_res_clashes_between_groups(pose, movable_residues, not_movable_residues, ignore_atom_beyond_cb=False)

    if num_res_clashes > num_res_clashes_tolerance:
        return False

    # Check buried unsatisfied hbonds

    buriedunsat_increase = pose_analysis.num_buried_unsatisfied_hbonds_change_upon_add_residues(
            pose, movable_region_start, movable_region_end)

    if buriedunsat_increase > 4:
        return False

    # Apply the torsions to the original_pose

    for i in range(movable_region_start, movable_region_end + 1):
        original_pose.set_phi(i, pose.phi(i))
        original_pose.set_psi(i, pose.psi(i))
        original_pose.set_omega(i, pose.omega(i))

    return True

def test_linker_pairs(pose, front_linker, back_linker, front_linker_start, back_linker_start, num_res_clashes_tolerance=0):
    '''Test if a pair of linkers could be used for form a loop-helix-loop unit.
    Return None if the pair of linkers cannot support a helix. Otherwise
    return (new_pose, reshaped_region_start, reshaped_region_stop).
    '''
    front_linker_length = len(front_linker['phis']) - 2
    back_linker_length = len(back_linker['phis']) - 2
   
    # Apply the linkers

    simple_pose_moves.apply_linker(pose, front_linker, front_linker_start)
    simple_pose_moves.apply_linker(pose, back_linker, back_linker_start)
   
    # Test if the gap is bridged

    front_helix = list(range(front_linker_start + 1 + front_linker_length, front_linker_start + 1 + front_linker_length + 10))
    back_helix = list(range(back_linker_start - 9, back_linker_start + 1))

    # Find residue pair RMSDs. The last residue of the back helix is ignored
    # because it would be problematic if it is trimed.

    res_pair_rmsds = [(r1, r2, residue_bb_rmsd_with_cutoff(pose, r1, r2)) for r1 in front_helix for r2 in back_helix[:-1]]

    # Find residues with the lowest direction difference

    res_pair_direction_dots = [(r[0], r[1], residue_direction_vector(pose, r[0]).dot(residue_direction_vector(pose, r[1]))) 
            for r in res_pair_rmsds if r[2] < 3]
    if len(res_pair_direction_dots) == 0: return None

    max_res_pair_direction_dot = max(res_pair_direction_dots, key=lambda x : x[2])
    if max_res_pair_direction_dot[2] < 0.5: return None
    
    # Trim the helix and bridge the gap
    
    new_pose = pose.clone()
    if not trim_helix_and_connect(new_pose, front_linker_start, back_helix[-1] + back_linker_length + 1, front_helix[0], 
            back_helix[-1], max_res_pair_direction_dot[0] + 1, max_res_pair_direction_dot[1], num_res_clashes_tolerance=num_res_clashes_tolerance):
        return None
   
    reshaped_region_start = front_linker_start
    reshaped_region_stop = back_linker_start + 1 + back_linker_length - (max_res_pair_direction_dot[1] - max_res_pair_direction_dot[0]) 
   
    # The local id of the cutpoint residue. A cutpoint will be made after this residue
    # The id of the start of LHL region (the residue before the front loop) is 1.
    cutpoint_residue = max_res_pair_direction_dot[0] - front_linker_start + 1

    return new_pose, reshaped_region_start, reshaped_region_stop, cutpoint_residue

def screen_loop_helix_loop_units_for_fixed_linker_length(output_dir, original_pose, lhl_start, lhl_stop, front_db, back_db, num_jobs, job_id, max_num_success=None,
        num_res_clashes_tolerance=0):
    '''Screen loop helix loop units for fixed length linkers.
    Return the torsions of the selected LHL units.
    '''
    # Prepare the pose

    pose = original_pose.clone()

    front_linker_length = len(front_db[0]['phis']) - 2
    back_linker_length = len(back_db[0]['phis']) - 2
    prepare_pose_for_lhl_screen(pose, lhl_start, lhl_stop, front_linker_length, back_linker_length) 

    front_linker_start = lhl_start
    back_linker_start = lhl_start + front_linker_length + 20

    # Define all tasks

    tasks = [(i, j) for i in range(len(front_db)) for j in range(len(back_db))]

    # Randomly shuffle the tasks of not all the tasks are screens

    if max_num_success:
        np.random.shuffle(tasks)

    # Run the tasks

    selected_lhl_units = []

    num_success = 0

    for i, task in enumerate(tasks):
        if job_id == i % num_jobs:
            if (i // num_jobs) % 100 == 0:
                print('Built {0} models after screening {1}/{2} pairs of linkers.'.format(num_success, i // num_jobs, len(tasks) // num_jobs + 1)) 
            
            test_result = test_linker_pairs(pose, front_db[task[0]], back_db[task[1]], 
                    front_linker_start, back_linker_start, num_res_clashes_tolerance=num_res_clashes_tolerance)
            
            if test_result is None: continue
                
            new_pose, reshaped_region_start, reshaped_region_stop, cutpoint_residue = test_result
           
            lhl_unit = {
                    'phis':[new_pose.phi(i) for i in range(reshaped_region_start, reshaped_region_stop + 1)],
                    'psis':[new_pose.psi(i) for i in range(reshaped_region_start, reshaped_region_stop + 1)],
                    'omegas':[new_pose.omega(i) for i in range(reshaped_region_start, reshaped_region_stop + 1)],
                    'cutpoint':cutpoint_residue    
                    } 
            selected_lhl_units.append(lhl_unit)

            num_success += 1

            #if num_success > 100:exit() ###DEBUG

            if max_num_success:
                if num_success >= max_num_success:
                    break

    return selected_lhl_units
            
def screen_all_loop_helix_loop_units(output_dir, pose, lhl_start, lhl_stop, front_linker_dbs, back_linker_dbs, num_jobs=1, job_id=0, max_num_success_each_db_pair=None,
        num_res_clashes_tolerance=0):
    '''Screen all loop helix loop units and record all possible designs.
    Return the torsions of the selected LHL units.
    '''
    # Do some input checking

    assert(1 < lhl_start < lhl_stop < pose.size())

    selected_lhl_units = []

    # Try all the combinations

    for front_db in front_linker_dbs:
        for back_db in back_linker_dbs:
        
            selected_lhl_units += screen_loop_helix_loop_units_for_fixed_linker_length(output_dir, pose, lhl_start, lhl_stop, front_db, back_db, num_jobs, job_id, 
                    max_num_success=max_num_success_each_db_pair, num_res_clashes_tolerance=num_res_clashes_tolerance)

    return selected_lhl_units
