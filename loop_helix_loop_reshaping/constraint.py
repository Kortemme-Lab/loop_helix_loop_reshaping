import numpy as np

import pyrosetta
from pyrosetta import rosetta

def get_bb_hbond_constraint(hbonds):
    '''Get a list of backbone hydrogen bond constraints. An Hbond is defined as
    (donor_res, acceptor_res)
    '''
    c_list = []
    
    hf_pair = rosetta.core.scoring.func.HarmonicFunc(2, 0.5)
    hf_angle = rosetta.core.scoring.func.HarmonicFunc(np.pi, 0.3)

    for donor_res, acceptor_res in set(hbonds):
        atomN = rosetta.core.id.NamedAtomID("N", donor_res)
        atomH = rosetta.core.id.NamedAtomID("H", donor_res)
        atomO = rosetta.core.id.NamedAtomID("O", acceptor_res)

        # Borrow the base_pair_constraint scoring term to constrain the Hbonds.
        # This can cause problem when one want to use the true base_pair_constraint.
        # Should add a dedicated term for hbond constraint when moving the function
        # to C++.

        c_list.append(rosetta.core.scoring.constraints.NamedAtomPairConstraint(
		atomH, atomO, hf_pair, rosetta.core.scoring.base_pair_constraint))
        c_list.append(rosetta.core.scoring.constraints.NamedAngleConstraint(
		atomN, atomH, atomO, hf_angle, rosetta.core.scoring.base_pair_constraint))
    
    return c_list

def add_constraints_to_pose(pose, constraint_list):
    '''Add a list of constraints to pose.'''
    cset = pose.constraint_set().clone() if pose.constraint_set() else rosetta.core.scoring.constraints.ConstraintSet()
    
    for c in constraint_list:
        cset.add_constraint(c)
    
    pose.constraint_set(cset)

def get_angle_constraint(pose, residues, angle, sd):
    '''Get an angle constraint between three residues.
    The angles are in radians. 
    '''
    hf = rosetta.core.scoring.func.HarmonicFunc(angle,  sd)
    a_ids = [rosetta.core.id.AtomID(pose.residue(i).atom_index("CA"), i) for i in residues]

    return rosetta.core.scoring.constraints.AngleConstraint(
            a_ids[0], a_ids[1], a_ids[2], hf, rosetta.core.scoring.angle_constraint)


