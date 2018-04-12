import pyrosetta
from pyrosetta import rosetta


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
            if pose.residue(r1).xyz('CA').distance(pose.residue(r2).xyz('CA')) < cutoff_distance:
                contact_degree += 1

        contact_degrees.append(contact_degree)

    return contact_degrees

