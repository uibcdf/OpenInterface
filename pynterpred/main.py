from simtk.unit import angstroms


def interface (item, receptor=None, ligand=None, threshold_1st_shell=6*angstroms, threshold_2nd_shell=9*angstroms):

    """interface(item, receptor=None, ligand=None, threshold_1st_shell=6*angstroms, threshold_2nd_shell=9*angstroms)

    Caracterización de interfaces.

    XXXXXXXXXXXXXXXXXXXXXXXXXXXX
    xxxxxxxxxxxxxxxxxxxxxxxx
    xxxxxxxxxxxxxxxxxxxxxxxxx
    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    xxxxxxxxxxxxxxxxxxxxxxxxxxxx

    Parameters
    ----------

    item: molecular model
        Molecular model as a pdb file, pdb id, mdtraj object... sobre el que se calcula la interface

    receptor: str, default=None
       definicion del receptor.

    
    ligand: str, default=None
       definicion del ligando.

    threshold_1st_shell: quantity, default=6*angstroms
        Distancia limite que define la primera capa


    threshold_2nd_shell: quantity, default=9*angstroms
        Distancia limite que define la segunda capa
    
    
    Returns
    -------

    List of lists
        Lista de los residuos que componen las capas.

    Examples
    --------


    See Also
    --------

    Notes
    -----

    Todo
    ----

    Warning
    -------

    """


    return item


