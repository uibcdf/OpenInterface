from simtk.unit import angstroms as _angstroms

def residues_interface(item, receptor_selection, ligand_selection, threshold = 15 * _angstroms):
    """residues_interface(item, receptor_selection, ligand_selection, threshold = 15 * _angstroms)

    List of residues defining the interface.

    This function defines the interface as the list of aminoacids of the receptor at a distance below a given threshold,
    and as the list of aminoacids of the ligand ...

    See reference Pepito Journal of whatever 2018.

    Parameters
    ----------

    item: mdtraj.Trajectory
        Molecular model as a mdtraj.Trajectory object containing the protein complex.

    receptor_selection: str
        Selection string of the receptor in item in mdtraj syntaxis (see: web de mdtraj con ejemplos de seleccion).

    ligand_selection: str
        Selection string of the ligand in item in mdtraj syntaxis (see: web de mdtraj con ejemplos de seleccion).

    threshold: quantity, default=15*angstroms
        Distance threshold defining the aminoacids belonging to the interface.

    Returns
    -------

    List:
        List of residue index of receptor interface

    List:
        List of residue index of ligand interface

    Examples
    --------
    .. highlight:: python

    >>> import pynterpred
    >>> import mdtraj as md

    >>> TPI = md.load('tpi.pdb')
    >>> list_interf_rec, list_interf_lig = pynterpred.b.residues_interface(TPI, "chainid 0", "chainid 1")

    >>> print("Residues in receptor:", list_interf_rec)
    >>> print("Residues in ligand:", list_interf_lig)

    See Also
    --------

    Notes
    -----

    Todo
    ----

    Warning
    -------

    """


    from simtk.unit import angstrom, nanometers
    import numpy as np
    
    list_atoms_CA_receptor = item.topology.select(receptor_selection+' and name CA')
    list_atoms_CA_ligand = item.topology.select(ligand_selection+' and name CA')
    num_CAs_receptor = len(list_atoms_CA_receptor)
    num_CAs_ligand = len(list_atoms_CA_ligand)
        
    # Para el receptor            
    CAs_elegidos_rec = []
    
    for CA_i_receptor in list_atoms_CA_receptor:
    
        si_es = False
        for CA_i_ligando in list_atoms_CA_ligand:
            
            pos_CA_i_receptor = item.xyz[0,CA_i_receptor,:]
            pos_CA_i_ligando = item.xyz[0,CA_i_ligando,:]
            
            distancia = np.sqrt( (pos_CA_i_receptor[0]-pos_CA_i_ligando[0])**2 +
                                 (pos_CA_i_receptor[1]-pos_CA_i_ligando[1])**2 +
                                 (pos_CA_i_receptor[2]-pos_CA_i_ligando[2])**2
                               ) * nanometers
        
            if distancia <= threshold:
                si_es = True
            
        if si_es == True:
            CAs_elegidos_rec.append(CA_i_receptor)
        
    lista_indices_residuos_elegidos_rec = []

    for ii in CAs_elegidos_rec:
        jj = item.topology.atom(ii).residue.index
        lista_indices_residuos_elegidos_rec.append(jj)
    
    # Para el ligando
    
    CAs_elegidos_lig = []
    
    for CA_i_ligand in list_atoms_CA_ligand:
    
        si_es = False
        for CA_i_receptor in list_atoms_CA_receptor:
        
            pos_CA_i_ligand = item.xyz[0,CA_i_ligand,:]
            pos_CA_i_receptor = item.xyz[0,CA_i_receptor,:]
        
            distancia = np.sqrt( (pos_CA_i_ligand[0]-pos_CA_i_receptor[0])**2 +
                                 (pos_CA_i_ligand[1]-pos_CA_i_receptor[1])**2 +
                                 (pos_CA_i_ligand[2]-pos_CA_i_receptor[2])**2
                                ) * nanometers
        
            if distancia <= threshold:
                si_es = True
                
        if si_es == True:
            CAs_elegidos_lig.append(CA_i_ligand)
    
    lista_indices_residuos_elegidos_lig = []

    for ii in CAs_elegidos_lig:
        jj = item.topology.atom(ii).residue.index
        lista_indices_residuos_elegidos_lig.append(jj) 

    return lista_indices_residuos_elegidos_rec, lista_indices_residuos_elegidos_lig






