import pyrosetta
from pyrosetta import rosetta


def find_bb_hbonds(pose):
    '''Find backbone HBonds of the pose. The N at first residue and C at last residue
    are ignored.
    Return:
        A map from pairs of (donor_residue, acceptor_residue) to HBonds.
    '''
    bb_hbonds = {}
    hbset = rosetta.core.scoring.hbonds.HBondSet(pose)

    for i in range(1, hbset.nhbonds() + 1):
        if hbset.hbond(i).don_hatm_is_protein_backbone() and hbset.hbond(i).acc_atm_is_protein_backbone():
            hb = (int(hbset.hbond(i).don_res()), int(hbset.hbond(i).acc_res()))
            
            if hb[0] > 1 and hb[1] < pose.size(): # Ignore the N-term N and C-term C
                bb_hbonds[hb] = hbset.hbond(i)
    
    return bb_hbonds

def get_helix_hbond_scores(pose, helix_start, helix_stop):
    '''Return a list of helix hbond scores.'''
    # Get the backbone hbonds  

    bb_hbonds_map = find_bb_hbonds(pose)
    hb_acceptors = list(range(helix_start, helix_stop - 3))
    hb_energies = []

    for acc in hb_acceptors: 
        key = (acc + 4, acc) 
        if not key in bb_hbonds_map.keys():
            hb_energies.append(1.0)
        else:
            hb_energies.append(bb_hbonds_map[key].energy())
    
    return hb_energies
