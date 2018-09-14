import pyrosetta
from pyrosetta import rosetta

from loop_helix_loop_reshaping import simple_pose_moves


def check_clashes_between_groups(pose, residues1, residues2, ignore_atom_beyond_cb=False):
    '''Return true if residues in two groups do not clash with each other.'''
    scale_factor = 0.36
    cb_id = 5

    for i in residues1:
        res_i = pose.residue(i)
        if res_i.is_virtual_residue(): continue
        
        for j in residues2:
            if i == j or i == j + 1 or i == j - 1:
                continue

            res_j = pose.residue(j)
            if res_j.is_virtual_residue(): continue

            # Rough check

            if res_i.xyz(1).distance(res_j.xyz(1)) > 10:
                continue

            # Detailed check
            
            atm_stop_i = min(cb_id, res_i.nheavyatoms()) if ignore_atom_beyond_cb else res_i.nheavyatoms()
            atm_stop_j = min(cb_id, res_j.nheavyatoms()) if ignore_atom_beyond_cb else res_j.nheavyatoms()

            for ai in range(1, atm_stop_i + 1):
                for aj in range(1, atm_stop_j + 1):

                    vi = res_i.xyz(ai)
                    vj = res_j.xyz(aj)
                    ri = res_i.atom_type(ai).lj_radius()
                    rj = res_j.atom_type(aj).lj_radius()

                    if (vi - vj).length_squared() < scale_factor * ((ri + rj) ** 2):
                        return False

    return True

def count_res_clashes_between_groups(pose, residues1, residues2, ignore_atom_beyond_cb=False):
    '''Return the number of residue clashes between two groups of residues'''
    scale_factor = 0.36
    cb_id = 5

    num_clashes = 0

    for i in residues1:
        res_i = pose.residue(i)
        if res_i.is_virtual_residue(): continue
        
        for j in residues2:
            if i == j or i == j + 1 or i == j - 1:
                continue

            res_j = pose.residue(j)
            if res_j.is_virtual_residue(): continue

            # Rough check

            if res_i.xyz(1).distance(res_j.xyz(1)) > 10:
                continue

            # Detailed check
            
            atm_stop_i = min(cb_id, res_i.nheavyatoms()) if ignore_atom_beyond_cb else res_i.nheavyatoms()
            atm_stop_j = min(cb_id, res_j.nheavyatoms()) if ignore_atom_beyond_cb else res_j.nheavyatoms()

            res_clash = False

            for ai in range(1, atm_stop_i + 1):
                for aj in range(1, atm_stop_j + 1):

                    vi = res_i.xyz(ai)
                    vj = res_j.xyz(aj)
                    ri = res_i.atom_type(ai).lj_radius()
                    rj = res_j.atom_type(aj).lj_radius()

                    if (vi - vj).length_squared() < scale_factor * ((ri + rj) ** 2):
                        res_clash = True

            if res_clash:
                num_clashes += 1

    return num_clashes

def contact_degrees_of_residues(pose, target_residues, residues_to_calc_against, cutoff_distance=10):
    '''Calculate the contact degrees of target residues against
    another group of residues. For one residue, the contact degree
    is the number of residues in residues_to_calc_against to which
    the CA-CA distance is within a given cutoff_distance.

    Return:
        A list of contact degrees for all target residues
    '''
    contact_degrees = []

    for r1 in target_residues:
        contact_degree = 0

        for r2 in residues_to_calc_against:
            if not pose.residue(r2).is_protein():
                continue
            
            if pose.residue(r1).xyz('CA').distance(pose.residue(r2).xyz('CA')) < cutoff_distance:
                contact_degree += 1

        contact_degrees.append(contact_degree)

    return contact_degrees

def num_buried_unsatisfied_hbonds(pose):
    '''Return the number of buried unsatisfied hbonds for a pose.'''
    buhf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<BuriedUnsatHbonds name="buriedunsat" cutoff="5" use_legacy_options="true" jump_number="0" />')
    return buhf.report_sm(pose)

def num_buried_unsatisfied_hbonds_change_upon_add_residues(pose, start, stop):
    '''Return the number change of buried unsatisfied hhonds
    upon adding a range of residues.
    '''
    pose_without_additions = pose.clone()
    simple_pose_moves.delete_region(pose_without_additions, start, stop)

    n_unsat_without_addition = num_buried_unsatisfied_hbonds(pose_without_additions)
    n_unsat = num_buried_unsatisfied_hbonds(pose)

    #print('num unsat: {0} - {1} = {2}'.format(n_unsat, n_unsat_without_addition, n_unsat - n_unsat_without_addition)) ###DEBUG

    return n_unsat - n_unsat_without_addition


