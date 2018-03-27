import pyrosetta
from pyrosetta import rosetta

import simple_pose_moves
import pose_analysis


def prepare_linker_selection(pose, linker_length, lhl_start, lhl_stop, front_linker):
    '''Prepare the pose for the linker selection. Mutate the pose int ALA,
    GLY and PRO. Remove the residues in the loop-helix-loop region and add 
    residues for the linker. Setup the fold tree properly.

    NOTE: the inserted residues will be linker_length + 1

    Args:
        front_linker: True if the linker is in the front of the loop-helix-loop
            unit
    '''
    # Remove the LHL region

    simple_pose_moves.delete_region(pose, lhl_start, lhl_stop)

    # Add the linker

    if front_linker:
        simple_pose_moves.insert_alas(pose, lhl_start - 1, linker_length + 1, insert_after=True, reset_fold_tree=True)
    else:
        simple_pose_moves.insert_alas(pose, lhl_start, linker_length + 1, insert_after=False, reset_fold_tree=True)

    # Mutate the pose

    for i in range(1, pose.size() + 1):
        if (not pose.residue(i).is_protein()) or pose.residue(i).name3() in ['GLY', 'PRO']:
            continue

        simple_pose_moves.mutate_residues(pose, [i], ['ALA'])

    rosetta.core.pose.correctly_add_cutpoint_variants(pose)

def apply_linker(pose, linker, start):
    '''Apply a linker to a pose.'''
    
    for i in range(len(linker['phis'])):
        seqpos = start + i
        pose.set_phi(seqpos, linker['phis'][i]) 
        pose.set_psi(seqpos, linker['psis'][i]) 
        pose.set_omega(seqpos, linker['omegas'][i])

def select_non_clashing_linkers(pose, candidate_linkers, linker_start):
    '''Select linkers that do not clash with the scaffold.

    NOTE: The torsions from linker_start - 1 to linker_start + linker_length
    will be reset.
    '''
    linker_residues = list(range(linker_start, linker_start + len(candidate_linkers[0]['phis']) - 2))
    pose_residues = list(range(1, pose.size() + 1))

    selected_linkers = []

    for i, linker in enumerate(candidate_linkers):
        apply_linker(pose, linker, linker_start - 1)

        if pose_analysis.check_clashes_between_groups(pose, linker_residues, pose_residues):
            selected_linkers.append(linker)

        if 0 == (i + 1) % 100 or i + 1 == len(candidate_linkers):
            print 'Tested {0} candidates. Selected {1}/{2} linkers'.format(i + 1, len(selected_linkers), len(candidate_linkers))

    return selected_linkers

