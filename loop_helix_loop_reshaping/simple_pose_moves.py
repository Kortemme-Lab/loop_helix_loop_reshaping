import pyrosetta
from pyrosetta import rosetta

def remove_cutpoint_variants(pose):
    '''Remove all cutpoint variants from the pose.'''
    for i in range(1, pose.size() + 1):
        rosetta.core.pose.remove_variant_type_from_pose_residue(pose,
                rosetta.core.chemical.CUTPOINT_LOWER, i)
        rosetta.core.pose.remove_variant_type_from_pose_residue(pose,
                rosetta.core.chemical.CUTPOINT_UPPER, i)

def delete_region(pose, start, stop):
    '''Delete a region of the pose.'''
    if stop > pose.size() or start < 1 or start > stop:
        return
    pose.delete_residue_range_slow(start, stop)
    pose.conformation().detect_disulfides()

def insert_alas(pose, position, length, insert_after=True, reset_fold_tree=True, fold_tree_root=1):
    '''Insert a poly-ALA peptide before or after a given position.,
    Set the fold tree to have a cutpoint before or after inserted residues.
    '''
    assert(1 <= position <= pose.size())

    # Set the fold tree with a single cutpoint

    def sub_fold_tree_add_edges_no_jump(ft, root, start, stop):
        '''Add edges to a sub-fold-tree that does not have
        and jumps.'''
        if start < root:
            ft.add_edge(root, start, -1)
        if stop > root:
            ft.add_edge(root, stop, -1)

    if reset_fold_tree:
        cutpoint = position if insert_after else position - 1
        ft = rosetta.core.kinematics.FoldTree()
        
        if fold_tree_root <= cutpoint and cutpoint < pose.size():
            sub_root = pose.size()
            ft.add_edge(fold_tree_root, sub_root, 1)
            sub_fold_tree_add_edges_no_jump(ft, sub_root, cutpoint + 1, pose.size())
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root, 1, cutpoint)
        
        elif fold_tree_root > cutpoint and cutpoint > 0:
            sub_root = 1
            ft.add_edge(fold_tree_root, sub_root, 1)
            sub_fold_tree_add_edges_no_jump(ft, sub_root, 1, cutpoint)
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root, cutpoint + 1, pose.size())

        else:
            sub_fold_tree_add_edges_no_jump(ft, fold_tree_root,  1, pose.size())
        
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
        
        if position + length + 1 <= pose.size():
            rosetta.core.conformation.idealize_position(position + length + 1, pose.conformation())
    else:
        if position - 1 > 0:
            rosetta.core.conformation.idealize_position(position - 1, pose.conformation())
        rosetta.core.conformation.idealize_position(position, pose.conformation())

def mutate_residues(pose, res_list, aa_list, protein_only=True):
    '''Mutate a list of residues. The list of AAs could
    either be 1 letter code or 3 letter code.
    '''
    aa_name_map = {'A':'ALA', 'P':'PRO', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET',
                   'F':'PHE', 'Y':'TYR', 'W':'TRP', 'S':'SER', 'T':'THR', 'C':'CYS',
                   'K':'LYS', 'R':'ARG', 'H':'HIS', 'D':'ASP', 'E':'GLU', 'N':'ASN',
                   'Q':'GLN', 'G':'GLY'}


    mutater = rosetta.protocols.simple_moves.MutateResidue()
    for i in range(len(res_list)):
        if protein_only and (not pose.residue(res_list[i]).is_protein()):
            continue

        name = aa_list[i] if len(aa_list[i]) == 3 else aa_name_map[aa_list[i]]
        mutater.set_res_name(name)
        mutater.set_target(res_list[i])
        mutater.apply(pose)

def mutate_pose_to_single_AA(pose, aa_name, keep_aa=['GLY', 'PRO']):
    '''Mutate the residues in a pose to a single type of AA.
    Keep the residues in the keep_aa list.
    '''
    keep_aa += [aa_name]

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

def update_insertion_points(insertion_points, insertion_id, new_length):
    '''Update the insertion point.'''
    old_length = insertion_points[insertion_id]['stop'] - insertion_points[insertion_id]['start'] - 1
    shift = new_length - old_length

    for i in range(len(insertion_points)):
        for k in ['start', 'stop']:
            if insertion_points[i][k] > insertion_points[insertion_id]['start']:
                insertion_points[i][k] += shift

def remove_insertion_residues(pose, insertion_points, insertion_ids=None):
    '''Remove the residues within insertion points and
    update the inerstion points.'''
    
    if insertion_ids is None:
        insertion_ids = range(len(insertion_points))
    
    for i in insertion_ids: 
        delete_region(pose, insertion_points[i]['start'] + 1, insertion_points[i]['stop'] - 1)
        update_insertion_points(insertion_points, i, 0)


