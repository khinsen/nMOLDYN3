"""This modules implements the functions and procedures that are related to chemistry.

Functions:
    * detect_hb_atoms     : detects the putative H-Bond acceptor, donor and hydrogens in a MMTK chemical object
    * belongToAnAmine     : determines whether an atom is part of an amine group.
    * belongToAHydroxy    : determines whether an atom is part of a hydroxy group.
    * belongToAMethyl     : determines whether an atom is part of a mathyl group.
    * belongToAThiol      : determines whether an atom is part of a thiol group.
    * hierarchizeUniverse : hierarchizes the molecules of the universe.
"""

# The python distribution modules
import os

from Scientific import N

# The MMTK modules
from MMTK import Atom, AtomCluster, Molecule
from MMTK.NucleicAcids import NucleotideChain
from MMTK.ParticleProperties import ParticleScalar
from MMTK.Proteins import PeptideChain, Protein

# The nMOLDYN modules
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error

Symbols = ["h","he",\
"li","be","b","c","n","o","f","ne",\
"na","mg","al","si","p","s","cl","ar",\
"k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga","ge","as","se","br","kr",\
"rb","sr","y","zr","nb","mo","tc","ru","rh","pd","ag","cd","in","sn","sb","te","i","xe",\
"cs","ba","la",\
"ce","pr","nd","pm","sm","eu","gd","tb","dy","ho","er","tm","yb","lu",\
"hf","ta","w","re","os","ir","pt","au","hg","tl","pb","bi","po","at","rn",\
"fr","ra","ac",\
"th","pa","u","np","pu","am","cm","bk","cf","es","fm","md","no","lr",\
"rf","db","sg","bh","hs","mt","ds","rg","uub"]

M_to_Symbol = {
 "1.0"  : "h", 
 "4.0"  : "he", 
 "6.9"  : "li", 
 "9.0"  : "be", 
 "10.8" : "b",
 "12.0" : "c",
 "14.0" : "n",
 "16.0" : "o",
 "19.0" : "f",
 "20.2" : "ne",
 "23.0" : "na",
 "24.3" : "mg",
 "27.0" : "al",
 "28.1" : "si",
 "31.0" : "p",
 "32.1" : "s",
 "35.5" : "cl",
 "39.1" : "k",
 "39.9" : "ar",
 "40.1" : "ca",
 "45.0" : "sc",
 "47.9" : "ti",
 "50.9" : "v",
 "52.0" : "cr",
 "54.9" : "mn",
 "55.8" : "fe",
 "58.7" : "ni",
 "58.9" : "co",
 "63.5" : "cu",
 "65.4" : "zn",
 "69.7" : "ga",
 "72.6" : "ge",
 "74.9" : "as",
 "79.0" : "se",
 "79.9" : "br",
 "83.8" : "kr",
 "85.5" : "rb",
 "87.6" : "sr",
 "88.9" : "y",
 "91.2" : "zr",
 "92.9" : "nb",
 "95.9" : "mo",
 "98.0" : "tc",
"101.1" : "ru",
"102.9" : "rh",
"106.4" : "pd",
"107.9" : "ag",
"112.4" : "cd",
"114.8" : "in",
"118.7" : "sn",
"121.8" : "sb",
"126.9" : "i",
"127.6" : "te",
"131.3" : "xe",
"132.9" : "cs",
"137.3" : "ba",
"138.9" : "la",
"140.1" : "ce",
"140.9" : "pr",
"144.2" : "nd",
"145.0" : "pm",
"150.4" : "sm",
"152.0" : "eu",
"157.3" : "gd",
"158.9" : "tb",
"162.5" : "dy",
"164.9" : "ho",
"167.3" : "er",
"168.9" : "tm",
"173.0" : "yb",
"175.0" : "lu",
"178.5" : "hf",
"180.9" : "ta",
"183.8" : "w",
"186.2" : "re",
"190.2" : "os",
"192.2" : "ir",
"195.1" : "pt",
"197.0" : "au",
"200.6" : "hg",
"204.4" : "tl",
"207.2" : "pb",
"209.0" : "bi",
"209.0" : "po",
"210.0" : "at",
"222.0" : "rn",
"223.0" : "fr",
"226.0" : "ra",
"227.0" : "ac",
"231.0" : "pa",
"232.0" : "th",
"237.0" : "np",
"238.0" : "u",
"243.0" : "am",
"244.0" : "pu",
"247.0" : "cm",
"247.0" : "bk",
"251.0" : "cf",
"252.0" : "es",
"257.0" : "fm",
"258.0" : "md",
"259.0" : "no",
"261.0" : "rf",
"262.0" : "lr",
"262.0" : "db",
"264.0" : "bh",
"266.0" : "sg",
"268.0" : "mt",
"277.0" : "hs"
}

# Dictionnary associating tuples of residue names (values) to their corresponding chemical family (key). 
residusChemFamily = {'acidic'      : ('Asp','Glu'),
                     'aliphatic'   : ('Ile','Leu','Val'),
                     'aromatic'    : ('His','Phe','Trp','Tyr'),
                     'basic'       : ('Arg','His','Lys'),
                     'charged'     : ('Arg','Asp','Glu','His','Lys'),
                     'hydrophobic' : ('Ala','Cys','Cyx','Gly','His','Ile','Leu','Lys','Met','Phe','Thr','Trp','Tyr','Val'),
                     'polar'       : ('Arg','Asn','Asp','Cys','Gln','Glu','His','Lys','Ser','Thr','Trp','Tyr'),
                     'small'       : ('Ala','Asn','Asp','Cys','Cyx','Gly','Pro','Ser','Thr','Val')}

def detect_hb_atoms_object(universe, heavy_hydrogen_cutoff = 0.11, objects = None):
    """Crude Detection of the putative H-Bond acceptor, donor and hydrogens in a MMTK chemical object.
    For object for which the connectivity is available (Protein, PeptideChain, NucleotideChain, 
    Molecule), the algorithm uses it otherwise (AtomCluster) it uses distance criteria.
    
    @param universe: the universe in which to search for H-bond partners.
    @type universe: a MMTK.Universe object
    
    @param heavy_hydrogen_cutoff: the cutoff distance between a heavy and hydrogen atoms under which
                                  they are considered to be bonded. Used only for AtomCluster object
                                  for which no connectivity is available.
    @type heavy_hydrogen_cutoff: float
    
    @params objects: the MMTK chemical object list to search for hydrogen bonds. If None, the whole universe 
                     object is taken by default.
    @type objects: a list of MMTK chemical objects.

    """
    
    if objects is None: objects = universe.objectList()

    # This dictionnary will store resp. the indexes of the acceptors, the donors and hydrogens
    # atoms found to be involved in a H-bond. If a donor has two hydrogens to give, its index is 
    # repeated twice so as the 'don' and 'hyd' lists lengths are the same and the index of one list 
    # corresponds to the other one.
    hb_atoms_indexes = {'acc' : [], 'don' : [], 'hyd' : []}

    # This dictionnary store some information about the heavy atom (name, symbol, initial position ...).
    hb_atoms_info = {}

    # Loop over the chemical objects of the universe.
    for obj in objects:

        # Case of an object for which the connectivity is available.
        if isinstance(obj, (Molecule,PeptideChain,Protein,NucleotideChain)):
            # List of the nitrogens and oxygens of the object.
            heavy_atoms = [at for at in obj.atomList() if at.type.name.lower() in ['oxygen','nitrogen']]
            # Loop over the nitrogen and oxygen atoms of the object. 
            for heavy in heavy_atoms:
                # The index of the atoms.
                ind = heavy.index
                # Its initial position.
                pos = [str(round(v*10.0,3)) for v in universe.configuration().array[ind,:]]
                # An entry for the hb_atoms_info dictionnary is created for that heavy atom.
                hb_atoms_info[ind] = {'mol' : obj.name, 'symb' : heavy.symbol, 'name' : heavy.name, 'pos' : pos}
                
                # List of the atoms connected to that atom.
                neighs = heavy.bondedTo()
                # List of the hydrogens connected to that atoms (can be empty).
                hydrogens = [at for at in neighs if at.type.name.strip().lower() == 'hydrogen']
                
                # Loop over the hydrogens connected to that atom.
                for hyd in hydrogens:
                    # The hydrogen and donor lists are appended with resp. the indexes of the current 
                    # hydrogen and the heavy atom the hydrogen is connected to. 
                    hb_atoms_indexes['don'].append(ind)
                    hb_atoms_indexes['hyd'].append(hyd.index)

                # If the heavy atom is a nitrogen and it is connected to more than 3 atoms, then it is
                # a 'ammonium'-like nitrogen that cannot be an acceptor.
                if heavy.type.name.lower() == 'nitrogen' and len(neighs) > 3: continue

                # The acceptor list is updated with the index of the heavy atom.
                hb_atoms_indexes['acc'].append(ind)
                            
        # Case of an object for which the connectivity is available.
        elif isinstance(obj, AtomCluster):
            # List of the nitrogens and oxygens of the object.
            heavy_atoms = [at for at in obj.atomList() if at.type.name.lower() in ['oxygen','nitrogen']]
            # List of the hydrogens of the object.
            hydrogens_atoms = [at for at in obj.atomList() if at.type.name.lower() == 'hydrogen']

            # Loop over the nitrogen and oxygen atoms of the object. 
            for heavy in heavy_atoms:
                # The index of the atoms.
                ind = heavy.index
                # Its initial position.
                pos = [str(round(v*10.0,3)) for v in universe.configuration().array[ind,:]]
                # An entry for the hb_atoms_info dictionnary is created for that heavy atom.
                hb_atoms_info[ind] = {'mol' : obj.name, 'symb' : heavy.symbol, 'name' : heavy.name, 'pos' : pos}
                
                n_bonds = 0
                
                # Loop over the hydrogens of the object.
                for hyd in hydrogens_atoms:
                    
                    # The heavy - hydrogen atoms distance is more than the cutoff, it can not be a bond.
                    if universe.distance(heavy,hyd) > heavy_hydrogen_cutoff: continue

                    # The hydrogen and donor lists are appended with resp. the indexes of the current 
                    # hydrogen and the heavy atom the hydrogen is close to.                     
                    hb_atoms_indexes['don'].append(ind)
                    hb_atoms_indexes['hyd'].append(hyd.index)

                    n_bonds += 1

                # If the heavy atom is a nitrogen and it is close to more than 3 hydrogens, then it is
                # a 'ammonium'-like nitrogen that cannot be an acceptor.
                if heavy.type.name.lower() == 'nitrogen' and n_bonds >= 3: continue

                # The acceptor list is updated with the index of the heavy atom.
                hb_atoms_indexes['acc'].append(ind)
                
        else:
            raise Error('Can not detect H-bond partners for an object of %s type' % obj.__class__.__name__)

    for k in hb_atoms_indexes.keys():
       hb_atoms_indexes[k] = N.array(hb_atoms_indexes[k], typecode = N.Int32)

    # Returns as output the hb_atoms_indexes and hb_atoms_info dictionnaries.
    return hb_atoms_indexes, hb_atoms_info

def detect_hb_atoms_subset(universe, heavy_hydrogen_cutoff = 0.11, subset = None):
    """Crude Detection of the putative H-Bond acceptor, donor and hydrogens in a subset of atoms.
    For atoms belonging to object for which the connectivity is available (Protein, PeptideChain, NucleotideChain, 
    Molecule), the algorithm uses it otherwise (AtomCluster) it uses distance criteria.
    
    @param universe: the universe in which to search for H-bond partners.
    @type universe: a MMTK.Universe object
    
    @param heavy_hydrogen_cutoff: the cutoff distance between a heavy and hydrogen atoms under which
                                  they are considered to be bonded. Used only for AtomCluster object
                                  for which no connectivity is available.
    @type heavy_hydrogen_cutoff: float
    
    @params subset: the MMTK atom collection to search for hydrogen bonds. If None, the whole universe 
                     is taken by default.
    @type subset: MMTK collection.
    """
    
    if subset is None: subset = universe.atomList()

    # This dictionnary will store resp. the indexes of the acceptors, the donors and hydrogens
    # atoms found to be involved in a H-bond. If a donor has two hydrogens to give, its index is 
    # repeated twice so as the 'don' and 'hyd' lists lengths are the same and the index of one list 
    # corresponds to the other one.
    hb_atoms_indexes = {'acc' : [], 'don' : [], 'hyd' : []}

    # This dictionnary store some information about the heavy atom (name, symbol, initial position ...).
    hb_atoms_info = {}

    # Loop over the chemical objects of the universe.
    for at in subset:
        
        obj = at.topLevelChemicalObject()

        # Case of an object for which the connectivity is available.
        if isinstance(obj, (Molecule,PeptideChain,Protein,NucleotideChain)):
            
            # List of the nitrogens and oxygens of the object.
            if at.type.name.lower() in ['oxygen','nitrogen']:
                            
                # The index of the atoms.
                ind = at.index
                
                # Its initial position.
                pos = [str(round(v*10.0,3)) for v in universe.configuration().array[ind,:]]
                
                # An entry for the hb_atoms_info dictionnary is created for that heavy atom.
                hb_atoms_info[ind] = {'mol' : obj.name, 'symb' : at.symbol, 'name' : at.name, 'pos' : pos}
                
                # List of the atoms connected to that atom.
                neighs = at.bondedTo()
                
                # List of the hydrogens connected to that atoms (can be empty).
                hydrogens = [at for at in neighs if at.type.name.strip().lower() == 'hydrogen']
                
                # Loop over the hydrogens connected to that atom.
                for hyd in hydrogens:
                    # The hydrogen and donor lists are appended with resp. the indexes of the current 
                    # hydrogen and the heavy atom the hydrogen is connected to. 
                    hb_atoms_indexes['don'].append(ind)
                    hb_atoms_indexes['hyd'].append(hyd.index)

                # If the heavy atom is a nitrogen and it is connected to more than 3 atoms, then it is
                # a 'ammonium'-like nitrogen that cannot be an acceptor.
                if at.type.name.lower() == 'nitrogen' and len(neighs) > 3: continue

                # The acceptor list is updated with the index of the heavy atom.
                hb_atoms_indexes['acc'].append(ind)
                            
        # Case of an object for which the connectivity is available.
        elif isinstance(obj, AtomCluster):
            
            # List of the nitrogens and oxygens of the object.
            if at.type.name.lower() in ['oxygen','nitrogen']:

                # List of the hydrogens of the object.
                hydrogens_atoms = [at for at in obj.atomList() if at.type.name.lower() == 'hydrogen']

                # The index of the atoms.
                ind = at.index
                
                # Its initial position.
                pos = [str(round(v*10.0,3)) for v in universe.configuration().array[ind,:]]
                
                # An entry for the hb_atoms_info dictionnary is created for that heavy atom.
                hb_atoms_info[ind] = {'mol' : obj.name, 'symb' : at.symbol, 'name' : at.name, 'pos' : pos}
                
                n_bonds = 0
                
                # Loop over the hydrogens of the object.
                for hyd in hydrogens_atoms:
                    
                    # The heavy - hydrogen atoms distance is more than the cutoff, it can not be a bond.
                    if universe.distance(at,hyd) > heavy_hydrogen_cutoff: continue

                    # The hydrogen and donor lists are appended with resp. the indexes of the current 
                    # hydrogen and the heavy atom the hydrogen is close to.                     
                    hb_atoms_indexes['don'].append(ind)
                    hb_atoms_indexes['hyd'].append(hyd.index)

                    n_bonds += 1

                # If the heavy atom is a nitrogen and it is close to more than 3 hydrogens, then it is
                # a 'ammonium'-like nitrogen that cannot be an acceptor.
                if at.type.name.lower() == 'nitrogen' and n_bonds >= 3: continue

                # The acceptor list is updated with the index of the heavy atom.
                hb_atoms_indexes['acc'].append(ind)
                
        else:
            raise Error('Can not detect H-bond partners for an object of %s type' % obj.__class__.__name__)

    for k in hb_atoms_indexes.keys():
       hb_atoms_indexes[k] = N.array(hb_atoms_indexes[k], typecode = N.Int32)

    # Returns as output the hb_atoms_indexes and hb_atoms_info dictionnaries.
    return hb_atoms_indexes, hb_atoms_info

def belongToAMethyl(atom):
    """Determines whether an atom is part of a methyl group.

    @param atom: the atom to test.
    @type atom: a MMTK.Atom object

    @return: The carbon atom of the methyl group if |atom| belongs to a methyl, None otherwise.
    @rtype: MMTK.Atom or None
    
    @note: this is a recursive function.
    """
    
    # The output variable is set to None.
    belong = None

    # Case where the atom to test is a carbon.
    if atom.type.name.strip().lower() == 'carbon':
        # List of the hydrogen atoms connected to |atom|.
        hydrogens = [neigh for neigh in atom.bondedTo() if neigh.type.name.strip().lower() == 'hydrogen']
        # If |atom| is a carbon and the number of connected hydrogen is 3, then it is a methyl.
        if len(hydrogens) >= 3:
            # In that case, |belong| is set to the instance of the carbon atom of the methyl group.
            belong = atom
    
    # Case where the atom to test is a hydrogen.
    elif atom.type.name.strip().lower() == 'hydrogen':
        # The atom connected to |atom|.
        c = atom.bondedTo()[0]
        # Checks that the connected atom is a carbon of a methyl group.
        if c != atom:
            # Recursive call to the function.
            belong = belongToAMethyl(c)
                
    return belong

def belongToAnAmine(atom):
    """Determines whether an atom is part of an amine group.

    @param atom: the atom to test.
    @type atom: a MMTK.Atom object

    @return: The nitrogen atom of the amine group if |atom| belongs to an amine, None otherwise.
    @rtype: MMTK.Atom or None

    @note: this is a recursive function.
    """
    
    # The output variable is set to None.
    belong = None

    # Case where the atom to test is a nitrogen.
    if atom.type.name.strip().lower() == 'nitrogen':

        # List of the hydrogen atoms connected to |atom|.
        hydrogens = [neigh for neigh in atom.bondedTo() if neigh.type.name.strip().lower() == 'hydrogen']

        # If |atom| is a nitrogen and the number of connected hydrogen is 2, then it is an amine.
        if len(hydrogens) >= 2:
            # In that case, |belong| is set to the instance of the nitrogen atom of the amine group.
            belong = atom
    
    # Case where the atom to test is a hydrogen.
    elif atom.type.name.strip().lower() == 'hydrogen':
        # The atom connected to |atom|.
        n = atom.bondedTo()[0]
        # Checks that the connected atom is a carbon of a methyl group.
        if n != atom:
            # Recursive call to the function.
            belong = belongToAnAmine(n)
                
    return belong

def belongToAThiol(atom):
    """Determines whether an atom is part of a thiol group.

    @param atom: the atom to test.
    @type atom: a MMTK.Atom object

    @return: The sulfur atom of the thiol group if |atom| belongs to a thiol, None otherwise.
    @rtype: MMTK.Atom or None

    @note: this is a recursive function.
    """

    # The output variable is set to None.
    belong = None

    # Case where the atom to test is a sulfur.
    if atom.type.name.strip().lower() in ['sulphur', 'sulfur']:

        # List of the hydrogen atoms connected to |atom|.
        hydrogens = [neigh for neigh in atom.bondedTo() if neigh.type.name.strip().lower() == 'hydrogen']
        
        # If |atom| is a sulfur and the number of connected hydrogen is 1, then it is a thiol.
        if len(hydrogens) >= 1:
            # In that case, |belong| is set to the instance of the sulfur atom of the thiol group.
            belong = atom
    
    # Case where the atom to test is a hydrogen.
    elif atom.type.name.strip().lower() == 'hydrogen':

        # The atom connected to |atom|.
        s = atom.bondedTo()[0]

        # Checks that the connected atom is a sulfur of a thiol group.
        if s != atom:
            # Recursive call to the function.
            belong = belongToAThiol(s)
                
    return belong

def belongToAHydroxy(atom):
    """Determines whether an atom is part of a hydroxy group.

    @param atom: the atom to test.
    @type atom: a MMTK.Atom object

    @return: The oxygen atom of the hydroxy group if |atom| belongs to a hydroxy, None otherwise.
    @rtype: MMTK.Atom or None

    @note: this is a recursive function.
    """

    # The output variable is set to None.
    belong = None

    # Case where the atom to test is an oxygen.
    if atom.type.name.strip().lower() == 'oxygen':
        
        # List of the hydrogen atoms connected to |atom|.
        hydrogens = [neigh for neigh in atom.bondedTo() if neigh.type.name.strip().lower() == 'hydrogen']

        # If |atom| is a oxygen and the number of connected hydrogen is 1, then it is a hydroxy.
        if len(hydrogens) >= 1:
            # In that case, |belong| is set to the instance of the sulfur atom of the thiol group.
            belong = atom
    
    # Case where the atom to test is a hydrogen.
    elif atom.type.name.strip().lower() == 'hydrogen':

        # The atom connected to |atom|.
        o = atom.bondedTo()[0]

        # Checks that the connected atom is an oxygen of a hydroxy group.
        if o != atom:
            # Recursive call to the function.
            belong = belongToAHydroxy(o)
                
    return belong

def setUniqueChemicalObjectName(obj):
    """
    """
        
    # In case of atoms, the nMOLDYN object name is just the atom type.
    if obj.numberOfAtoms() == 1:
        try:
            nMOLDYNObjectName = obj.type.name
        except:
            nMOLDYNObjectName = obj.name
        
    else:
    
        try:
            # If the obj has no 'name' attribute or if is empty throw an error.
            if not obj.name:
                raise AttributeError

        # Case of an atom with no/empty 'name' attribute.
        except AttributeError:
            # Check if the object is referenced in the MMTK database.
            try:
                if not obj.type.name:
                    raise AttributeError                
                
            # Otherwise, sets the |obj.nmoldynname| attribute to an arbitrary name.
            except AttributeError:
                
                # If the object is a nucleotide chain, the name will start with 'NC'.
                if isinstance(obj, NucleotideChain):
                    nMOLDYNObjectName = 'NC'
            
                # If the object is a protein, the name will start with 'P'.
                elif isinstance(obj, Protein):
                    nMOLDYNObjectName = 'P'

                # If the object is a peptide chain, the name will start with 'PC'.
                elif isinstance(obj, PeptideChain):
                    nMOLDYNObjectName = 'PC'
            
                # If the object is an atom cluster, the name will start with 'AC'.
                elif isinstance(obj, AtomCluster):
                    nMOLDYNObjectName = 'AC'                                
                
                # If the object is a molecule, the name will start with 'M'.
                elif isinstance(obj, Molecule):
                    nMOLDYNObjectName = 'M'

                # Then appends the number of atoms of the object to finalize |obj.nmoldynname| attribute.
                nMOLDYNObjectName += str(obj.numberOfAtoms())
                
            else:
                # In that case, sets the |obj.nmoldynname| attribute to the name under which it is referenced in the database.
                nMOLDYNObjectName = obj.type.name
        else:
            nMOLDYNObjectName = obj.name
                
    return nMOLDYNObjectName
    
def hierarchizeUniverse(universe):
    """Hierarchizes the molecules of the universe.
    
    @param universe: the MMTK universe to hierarchize.
    @type universe: an instance of MMTK.Universe.
    
    @return: A list whose 1st element is a dictionnary storing the universe chemical hierarchy and 2nd
             element is a list of all the unique nMOLDYN name defined for each obj of the universe.
    @rtype: List
    
    """
    
    chemObjInfo = {}
    chemObjNamesList = []
    
    # Loop over the object found in the universe.
    for obj in universe.objectList():

        # If the object is not a chemical object then skip it.
        if not isinstance(obj,(Atom, AtomCluster, Molecule, NucleotideChain, PeptideChain, Protein)):
            continue

        # A "unique" name (the nMOLDYN name) is defined for the current object.        
        nMOLDYNObjectName = setUniqueChemicalObjectName(obj)
        
        chemObjNamesList.append(nMOLDYNObjectName)

        if isinstance(obj, NucleotideChain):
            if chemObjInfo.has_key(nMOLDYNObjectName):
                chemObjInfo[nMOLDYNObjectName]['number'] += 1
                
            else:
                chemObjInfo[nMOLDYNObjectName] = {}
                
                chemObjInfo[nMOLDYNObjectName]['number'] = 1
                
                # The objectclass name for a MMTK-NucleotideChain object.
                chemObjInfo[nMOLDYNObjectName]['objectclass']   = 'NucleotideChain'                
                chemObjInfo[nMOLDYNObjectName]['groupinglevel'] = ['atom', 'amine', 'residue', 'nucleicacid']
                
                # A list of the MMTK full names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomname'] = ['*']

                # A set of the MMTK types of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomtype'] = set(['*'])

                # A set of the MMTK element names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = set(['*'])
                
                # A list of the MMTK names of the nucleotides of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['nuclname'] = ['*']

                # A list of the MMTK type of the nucleotides of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['nucltype'] = set(['*'])

                # Loop over the residues of the chain.
                for nucl in obj.residues():
                    chemObjInfo[nMOLDYNObjectName]['nuclname'].append(nucl.fullName())
                    chemObjInfo[nMOLDYNObjectName]['nucltype'].add(nucl.symbol)

                    # Loop over the atoms of the residue.
                    for atom in nucl.atomList():
                        chemObjInfo[nMOLDYNObjectName]['atomname'].append(atom.fullName())
                        chemObjInfo[nMOLDYNObjectName]['atomtype'].add(atom.name)
                        chemObjInfo[nMOLDYNObjectName]['atomelement'].add(atom.type.name)

                chemObjInfo[nMOLDYNObjectName]['atomtype'] = sorted(chemObjInfo[nMOLDYNObjectName]['atomtype'])
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = sorted(chemObjInfo[nMOLDYNObjectName]['atomelement'])

                chemObjInfo[nMOLDYNObjectName]['nucltype'] = sorted(chemObjInfo[nMOLDYNObjectName]['nucltype'])

                chemObjInfo[nMOLDYNObjectName]['misc'] = ['backbone', 'bases']
                        
        elif isinstance(obj, PeptideChain):
            if chemObjInfo.has_key(nMOLDYNObjectName):
                chemObjInfo[nMOLDYNObjectName]['number'] += 1
                
            else:
                chemObjInfo[nMOLDYNObjectName] = {}
                
                chemObjInfo[nMOLDYNObjectName]['number'] = 1
                
                # The objectclass name for a MMTK-PeptideChain object.
                chemObjInfo[nMOLDYNObjectName]['objectclass']   = 'PeptideChain'
                
                chemObjInfo[nMOLDYNObjectName]['groupinglevel'] = ['atom', 'amine', 'hydroxy', 'methyl', 'thiol', 'residue', 'chain']
                
                # A list of the MMTK full names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomname'] = ['*']
                
                # A set of the MMTK types of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomtype'] = set(['*'])
                
                # A set of the MMTK element names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = set(['*'])
                
                # A list of the MMTK residue names of the residues of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['resname'] = ['*']
                
                # A list of the MMTK type of the residues of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['restype'] = set(['*'])
                
                # A list of the classes of the residues of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['resclass'] = sorted(residusChemFamily.keys())

                # Loop over the residues of the chain.
                for residue in obj.residues():
                    chemObjInfo[nMOLDYNObjectName]['resname'].append(residue.fullName())
                    chemObjInfo[nMOLDYNObjectName]['restype'].add(residue.symbol)

                    # Loop over the atoms of the residue.
                    for atom in residue.atomList():
                        chemObjInfo[nMOLDYNObjectName]['atomname'].append(atom.fullName())
                        chemObjInfo[nMOLDYNObjectName]['atomtype'].add(atom.name)
                        chemObjInfo[nMOLDYNObjectName]['atomelement'].add(atom.type.name)

                chemObjInfo[nMOLDYNObjectName]['atomtype']    = sorted(chemObjInfo[nMOLDYNObjectName]['atomtype'])
                
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = sorted(chemObjInfo[nMOLDYNObjectName]['atomelement'])

                chemObjInfo[nMOLDYNObjectName]['restype']     = sorted(chemObjInfo[nMOLDYNObjectName]['restype'])

                chemObjInfo[nMOLDYNObjectName]['chemfragment'] = ['amine', 'c_alphas', 'hydroxy', 'methyl', 'thiol']

                chemObjInfo[nMOLDYNObjectName]['misc'] = ['backbone', 'sidechains']
        
        elif isinstance(obj, Protein):
            if chemObjInfo.has_key(nMOLDYNObjectName):
                chemObjInfo[nMOLDYNObjectName]['number'] += 1
                
            else:
                chemObjInfo[nMOLDYNObjectName] = {}
                
                chemObjInfo[nMOLDYNObjectName]['number'] = 1
                
                # The objectclass name for a MMTK-Protein object.
                chemObjInfo[nMOLDYNObjectName]['objectclass'] = 'Protein'
                
                chemObjInfo[nMOLDYNObjectName]['groupinglevel'] = ['atom', 'amine', 'hydroxy', 'methyl', 'thiol', 'residue','chain', 'protein']
                
                # A list of the MMTK full names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomname'] = ['*']
                
                # A set of the MMTK types of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomtype'] = set(['*'])
                
                # A set of the MMTK element names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = set(['*'])

                # A list of the MMTK residue names of the residues of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['resname'] = ['*']
                
                # A set of the MMTK residue types of the residues of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['restype'] = set(['*'])
                
                # A list of the chemical class of the residues of the object.
                chemObjInfo[nMOLDYNObjectName]['resclass'] = sorted(residusChemFamily.keys())

                # A list of the MMTK chain names of the chains of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['chainname'] = ['*']

                # Loop over the chains of the protein.
                for chain in obj:
                    
                    # Append the chain name.
                    chemObjInfo[nMOLDYNObjectName]['chainname'].append(chain.fullName())

                    # Loop over the residues of the chain.
                    for residue in chain.residues():
                        chemObjInfo[nMOLDYNObjectName]['resname'].append(residue.fullName())
                        chemObjInfo[nMOLDYNObjectName]['restype'].add(residue.symbol)

                        # Loop over the atoms of the residue.
                        for atom in residue.atomList():
                            chemObjInfo[nMOLDYNObjectName]['atomname'].append(atom.fullName())
                            chemObjInfo[nMOLDYNObjectName]['atomtype'].add(atom.name)
                            chemObjInfo[nMOLDYNObjectName]['atomelement'].add(atom.type.name)

                # The atomtype set is converted to a sorted list.
                chemObjInfo[nMOLDYNObjectName]['atomtype']    = sorted(chemObjInfo[nMOLDYNObjectName]['atomtype'])
                
                # The atomelement set is converted to a sorted list.
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = sorted(chemObjInfo[nMOLDYNObjectName]['atomelement'])

                # The restype set is converted to a sorted list.
                chemObjInfo[nMOLDYNObjectName]['restype'] = sorted(chemObjInfo[nMOLDYNObjectName]['restype'])

                # 'chemfragment'/sorted list of the fragment names that can be selected in the object.
                chemObjInfo[nMOLDYNObjectName]['chemfragment'] = ['amine', 'c_alphas', 'hydroxy', 'methyl', 'thiol']

                # 'misc'/list of miscellaneous selectable subparts of the object.
                chemObjInfo[nMOLDYNObjectName]['misc'] = ['backbone', 'sidechains']

        elif isinstance(obj, Atom):
            if chemObjInfo.has_key(nMOLDYNObjectName):
                chemObjInfo[nMOLDYNObjectName]['number'] += 1
                
            else:
                chemObjInfo[nMOLDYNObjectName] = {}
                
                chemObjInfo[nMOLDYNObjectName]['number'] = 1
                
                # The objectclass name for a MMTK-Atom object.
                chemObjInfo[nMOLDYNObjectName]['objectclass']   = 'Atom'
                
                chemObjInfo[nMOLDYNObjectName]['groupinglevel'] = ['atom']

                # A list of the MMTK full names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomname'] = set(['*'])

            chemObjInfo[nMOLDYNObjectName]['atomname'].add(obj.fullName())
                
        elif isinstance(obj, AtomCluster):
            if chemObjInfo.has_key(nMOLDYNObjectName):
                chemObjInfo[nMOLDYNObjectName]['number'] += 1
                
            else:
                chemObjInfo[nMOLDYNObjectName] = {}
                
                chemObjInfo[nMOLDYNObjectName]['number']        = 1
                
                # The objectclass name for a MMTK-AtomCluster object.
                chemObjInfo[nMOLDYNObjectName]['objectclass']   = 'AtomCluster'
                
                chemObjInfo[nMOLDYNObjectName]['groupinglevel'] = ['atom','cluster']

                # A set of the MMTK full names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomname'] = set(['*'])
                
                # A set of the MMTK element names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = set(['*'])

                for atom in obj.atomList():
                    chemObjInfo[nMOLDYNObjectName]['atomname'].add(atom.fullName())
                    chemObjInfo[nMOLDYNObjectName]['atomelement'].add(atom.type.name)
                    
                chemObjInfo[nMOLDYNObjectName]['atomname']    = sorted(chemObjInfo[nMOLDYNObjectName]['atomname'])
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = sorted(chemObjInfo[nMOLDYNObjectName]['atomelement'])
        
        elif isinstance(obj, Molecule):
            if chemObjInfo.has_key(nMOLDYNObjectName):
                chemObjInfo[nMOLDYNObjectName]['number'] += 1
                
            else:
                chemObjInfo[nMOLDYNObjectName] = {}
                
                chemObjInfo[nMOLDYNObjectName]['number'] = 1
                
                # The objectclass name for a MMTK-Molecule object.
                chemObjInfo[nMOLDYNObjectName]['objectclass']   = 'Molecule'
                
                chemObjInfo[nMOLDYNObjectName]['groupinglevel'] = ['atom','amine', 'hydroxy', 'methyl', 'thiol', 'molecule']
                
                chemObjInfo[nMOLDYNObjectName]['atomname']    = set(['*'])
                
                # A set of the MMTK element names of the atoms of the object + the wildcard '*'.
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = set(['*'])
                
                for atom in obj.atomList():
                    chemObjInfo[nMOLDYNObjectName]['atomname'].add(atom.fullName())
                    chemObjInfo[nMOLDYNObjectName]['atomelement'].add(atom.type.name)
                    
                chemObjInfo[nMOLDYNObjectName]['atomname']    = sorted(chemObjInfo[nMOLDYNObjectName]['atomname'])
                chemObjInfo[nMOLDYNObjectName]['atomelement'] = sorted(chemObjInfo[nMOLDYNObjectName]['atomelement'])

                # The list of the fragment names that can be selected for a Molecule object.
                chemObjInfo[nMOLDYNObjectName]['chemfragment'] = ['amine', 'hydroxy', 'methyl', 'thiol']

    # For the ojects of Atom class, the atomname list is sorted in the alphabetical order (e.g. H1,h10,H11,H12)
    for oName in chemObjInfo.keys():
        oClass = chemObjInfo[oName]['objectclass']
        if oClass == 'Atom':
            chemObjInfo[oName]['atomname'] = sorted(chemObjInfo[oName]['atomname'])

                
    # An additional entry is added to |universe.nmoldyncontents| dictionnary to account for the case when one deals
    # with all the atoms of the universe as if they were a single object.
    chemObjInfo['*'] = {'atomelement' : set(['*']),\
                        'number' : universe.numberOfAtoms(),\
                        'objectclass' : 'AllClass',\
                        'groupinglevel' : ['default',]}
    
    for atom in universe.atomList():
        chemObjInfo['*']['atomelement'].add(atom.type.name)
        
    chemObjInfo['*']['atomelement'] = sorted(chemObjInfo['*']['atomelement'])
    
    return chemObjInfo, chemObjNamesList
    