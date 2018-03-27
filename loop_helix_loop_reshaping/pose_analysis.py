import pyrosetta
from pyrosetta import rosetta


def check_clashes_between_groups(pose, residues1, residues2):
    '''Return true if residues in two groups do not clash with each other.'''
    scale_factor = 0.36

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

            for ai in range(1, res_i.nheavyatoms() + 1):
                for aj in range(1, res_j.nheavyatoms() + 1):

                    vi = res_i.xyz(ai)
                    vj = res_j.xyz(aj)
                    ri = res_i.atom_type(ai).lj_radius()
                    rj = res_j.atom_type(aj).lj_radius()

                    if (vi - vj).length_squared() < scale_factor * ((ri + rj) ** 2):
                        return False

    return True
