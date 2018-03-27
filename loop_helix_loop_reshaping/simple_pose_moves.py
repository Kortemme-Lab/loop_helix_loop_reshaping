import pyrosetta
from pyrosetta import rosetta

def delete_region(pose, start, stop):
    '''Delete a region of the pose.'''
    if stop > pose.size() or start < 1 or start > stop:
        return
    pose.delete_residue_range_slow(start, stop)
    pose.conformation().detect_disulfides()

def insert_alas(pose, position, length, insert_after=True, reset_fold_tree=True):
    '''Insert a poly-ALA peptide before or after a given position.,
    Set the fold tree to have a cutpoint before or after inserted residues.
    '''
    assert(1 < position < pose.size())
    
    # Set the fold tree with a single cutpoint

    if reset_fold_tree:
        cutpoint = position if insert_after else position - 1
        ft = rosetta.core.kinematics.FoldTree()
        ft.add_edge(1, pose.size(), 1)
        ft.add_edge(1, cutpoint, -1)
        ft.add_edge(pose.size(), cutpoint + 1, -1)
        pose.fold_tree(ft)

    # Append the residues

    residue_type_set = pose.residue_type_set_for_pose()
    new_rsd = rosetta.core.conformation.ResidueFactory.create_residue( residue_type_set.name_map("ALA") )
   
    for i in range(length):
        if insert_after:
            pose.conformation().safely_append_polymer_residue_after_seqpos(new_rsd, position + i, True)
            pose.set_omega(position + i, 180)
        else:
            pose.conformation().safely_prepend_polymer_residue_before_seqpos(new_rsd, position, True)
            pose.set_omega(position, 180)

    if insert_after:
        rosetta.core.conformation.idealize_position(position + length, pose.conformation())
        rosetta.core.conformation.idealize_position(position + length + 1, pose.conformation())
    else:
        rosetta.core.conformation.idealize_position(position - 1, pose.conformation())
        rosetta.core.conformation.idealize_position(position, pose.conformation())

def mutate_residues(pose, res_list, aa_list):
    '''Mutate a list of residues. The list of AAs could
    either be 1 letter code or 3 letter code.
    '''
    aa_name_map = {'A':'ALA', 'P':'PRO', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET',
                   'F':'PHE', 'Y':'TYR', 'W':'TRP', 'S':'SER', 'T':'THR', 'C':'CYS',
                   'K':'LYS', 'R':'ARG', 'H':'HIS', 'D':'ASP', 'E':'GLU', 'N':'ASN',
                   'Q':'GLN', 'G':'GLY'}


    mutater = rosetta.protocols.simple_moves.MutateResidue()
    for i in range(len(res_list)):
        name = aa_list[i] if len(aa_list[i]) == 3 else aa_name_map[aa_list[i]]
        mutater.set_res_name(name)
        mutater.set_target(res_list[i])
        mutater.apply(pose)

def mutate_pose_to_single_AA(pose, aa_name, keep_aa=['GLY', 'PRO']):
    '''Mutate the residues in a pose to a single type of AA.
    Keep the residues in the keep_aa list.
    '''
    for i in range(1, pose.size() + 1):
        if (not pose.residue(i).is_protein()) or pose.residue(i).name3() in keep_aa:
            continue

        mutate_residues(pose, [i], [aa_name])

def apply_linker(pose, linker, start):
    '''Apply a linker to a pose.'''
    
    for i in range(len(linker['phis'])):
        seqpos = start + i
        pose.set_phi(seqpos, linker['phis'][i]) 
        pose.set_psi(seqpos, linker['psis'][i]) 
        pose.set_omega(seqpos, linker['omegas'][i])

