
def sasa_interface(item, receptor_selection, ligand_selection):
    
    from numpy import hstack
    from simtk.unit import nanometer
    from mdtraj import shrake_rupley
    
    receptor_list_atoms = item.topology.select(receptor_selection)
    ligand_list_atoms = item.topology.select(ligand_selection)
    complex_list_atoms = hstack((receptor_list_atoms, ligand_list_atoms))
    
    receptor = item.atom_slice(receptor_list_atoms)
    ligand = item.atom_slice(ligand_list_atoms)
    complex = item.atom_slice(complex_list_atoms)
    
    sasa_atoms_receptor = shrake_rupley(receptor, mode = 'atom')
    sasa_atoms_ligand = shrake_rupley(ligand, mode = 'atom')
    sasa_atoms_complex = shrake_rupley(complex, mode = 'atom')

    sasa_total_receptor = sasa_atoms_receptor.sum()
    sasa_total_ligand = sasa_atoms_ligand.sum()
    sasa_total_complex = sasa_atoms_complex.sum()

    sasa_total = sasa_total_receptor + sasa_total_ligand - sasa_total_complex
    sasa_total = sasa_total * nanometer**2
    
    return sasa_total


def sasa_angstroms(sasa_original):

    from simtk.unit import angstroms
    
    sasa_ang = sasa_original.in_units_of(angstroms**2)
    
    return sasa_ang


def sasa_is_drogable(sasa_original):

    from simtk.unit import angstroms

    if sasa_original > 800 * angstroms**2:
        return True

    else:
        return False


def interface_is_drogable(item, receptor_selection, ligand_selection):
    
    sasa_auxiliar = sasa_interface(item, receptor_selection, ligand_selection)
    return sasa_is_drogable(sasa_auxiliar)







