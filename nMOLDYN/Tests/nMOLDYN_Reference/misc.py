import getpass
import re
import sys
from os import *
from string import *
from tempfile import mktemp, gettempdir
from time import asctime,localtime,time

from Scientific import N
from Scientific.IO.NetCDF import NetCDFFile
from Scientific.IO.PDB import PDBFile
from Scientific.IO.TextFile import TextFile

from MMTK import *
from MMTK import Units
from MMTK.Collections import Collection
from MMTK.Database import AtomReference
from MMTK.ParticleProperties import ParticleScalar
from MMTK.PDB import PDBConfiguration
from MMTK.Trajectory import Trajectory, TrajectorySet

from nMOLDYN.Tests.nMOLDYN_Reference.calc import qVectorGenerator, qTrajectory, Quidam

def BincohList(universe, atoms, deuter=None):

    bincoh = universe.getAtomScalarArray('b_incoherent')
    if deuter is not None:
        for atom in deuter:
            bincoh[atom] = atom.b_incoherent_deut
    bincoh = bincoh*bincoh
    if atoms != universe:
        bincoh = bincoh*atoms.booleanMask()
    return bincoh/bincoh.sumOverParticles()


def BcohList(universe, atoms, deuter=None):

    bcoh = universe.getAtomScalarArray('b_coherent')
    if deuter is not None:
        for atom in deuter:
            bcoh[atom] = atom.b_coherent_deut
    if atoms != universe:
        bcoh = bcoh*atoms.booleanMask()
        
    return bcoh/N.sqrt((bcoh*bcoh).sumOverParticles())


## def BincohList(universe, atoms, deuter=None):
 
##     bincoh = universe.getAtomScalarArray('b_incoherent')
##     if deuter is not None:
##         bincoh_d = Atom('d').b_incoherent
##         for atom in deuter:
##             if atom.symbol == 'H':
##                 bincoh[atom] = bincoh_d
##     bincoh = bincoh*bincoh
##     if atoms != universe:
##         bincoh = bincoh*atoms.booleanMask()
##     return bincoh/bincoh.sumOverParticles()

## def BcohList(universe, atoms, deuter=None):

##     bcoh = universe.getAtomScalarArray('b_coherent')
##     if deuter is not None:
##         bcoh_d = Atom('d').b_coherent
##         for atom in deuter:
##             if atom.symbol == 'H':
##                 bcoh[atom] = bcoh_d
##     if atoms != universe:
##         bcoh = bcoh*atoms.booleanMask()
##     return bcoh/sqrt((bcoh*bcoh).sumOverParticles())

def MassList(universe,atoms):
 
    mass = universe.masses()
    if atoms != universe:
        mass = mass*atoms.booleanMask()
    return mass/mass.sumOverParticles()

def OnesList(universe,atoms):
 
    list = ParticleScalar(universe,N.ones(universe.numberOfAtoms(),N.Float))
    list = list/len(atoms)
    return list

def timePrepare(traj,timeInfo=None):
    """
    just preparing a time axis basing on timeInfo given
    or the size of a trajectory
    """
    if timeInfo is None:
        timeInfo = (0,len(traj)-1,1)
    time = traj.time[timeInfo[0]:timeInfo[1]:timeInfo[2]]
    dt = time[1]-time[0]
    return dt*N.arange(len(time))

def getData(filename):
    """
    reads columns of data stored in an ASCII file and
    return an multidimensional array with it
    """ 
    try: f = TextFile(filename,"r")
    except: return None
    lines = f.readlines()
    i = 0; data = []; ole = 0
    for ia in lines:
        data.append([])
        if ia[0] != '#':
            try: ie = index(ia,'#')
            except ValueError: ie = len(ia)
            for ib in split(ia[:ie]): data[i].append(float(ib))
            i = i+1; ile = len(split(ia[:ie]))
            if ile > ole: ole = ile
    res = zeros((i,ole),N.Float)
    for item in range(len(data)):
        if len(data[item]) != 0:
            for num in range(len(data[item])): res[item,num] = data[item][num]
    return res


def inputFileRead(filename):
    """
    read and process an input file
    """
    keywords = ['trajectory','output_files','title','time_info','time_steps',
                'frequency_points','q_vector_set','deuter','projection_vector',
                'reference','rotation_coefficients','ft_window',
                'groups','weights','atoms','units_length','units_frequency',
                'log_file','groups_code','atoms_code','filter_window',
                'results_file','atoms_pdb','differentiation','verbose',
                'symbols','ar_order','ar_precision']
    newvars = {}
    file_text = Utility.readURL(filename)
    exec file_text in vars(sys.modules['__builtin__']), newvars
    input = Quidam()
    for key in keywords:
        if newvars.has_key(key): setattr(input,key,newvars[key])
        else: setattr(input,key,None)

    import os
    print os.getcwd()
    print input.trajectory        
    #
    # general settings
    #
    if input.trajectory:
        if len(input.trajectory)==1:
            traj = Trajectory(None,input.trajectory[0],'r')
            if traj.variables().count('quaternion') > 0:
                traj = qTrajectory(None,input.trajectory[0],'r')
        elif len(input.trajectory)>1:
            traj = TrajectorySet(None,input.trajectory)

        if not input.units_length: input.units_length = Units.nm
        if not input.units_frequency:
            input.units_frequency = Units.tera*Units.Hz

        types = getTypes(traj.universe)
        if input.q_vector_set is None:
            input.q_vector_set = (N.arange(0., 100., 2.), 1., 50) # default
        if input.time_info is None: input.time_info = (0,len(traj),1)
        qVectors = qVectorGenerator(input.q_vector_set,traj)
        #
        # Substitute Hydrogen atoms with Deuter?
        # 
        if input.deuter:
            collection = Collection()
            for i in input.deuter.keys():
                for ia in input.deuter[i]:
                    gj = getChemicalObjects({i: types[i]}, {i: ia})
                    collection.addObject(gj)
            h2d = collection.atomList()
            print 'number of Deuter atoms: ',len(h2d)
        #
        # Atom selection related keywords
        # atoms selected in a different way are stored together
        # and filtered at the end (if there're repetitions only
        # atoms which occur many times are stored and passed to
        # further calculations, else all atoms are passed)
        #
        if input.atoms:
            print 'processing atom selection:\n\t',input.atoms
            collection = []
            for i in input.atoms.keys():
                for ia in input.atoms[i]:
                    typs = {}; typs[i] = types[i]
                    vlst = {}; vlst[i] = ia
                    gj = getChemicalObjects(typs,vlst)
                    collection = collection + gj.atomList()
            input.atoms = collection
            print '\t...done\n\tstored ',len(input.atoms),' atoms\n'
        if input.atoms_pdb:
            print 'processing atom selection from a PDB file\n\t(',
            print input.atoms_pdb,'):'
            atoms_add = parsePDBAtomSelection(input.atoms_pdb,traj)
            print '\t...done\n\tstored ',len(atoms_add.atomList()),' atoms\n'
            if input.atoms: input.atoms = input.atoms + atoms_add.atomList()
            else: input.atoms = atoms_add.atomList()
        if input.atoms_code:
            print 'processing atom selection hardcoded in Python'
            # syntax:
            # def atoms_code(traj,nothing,dummy_a='gj')
            #     # a python code here
            #     return Collection(atom_list)
            #
            # atoms_code is a function object whose first argument is
            # a Trajectory object and which returns a Collection object.
            print '...done\n\tstored ',len(input.atoms_code(traj).atomList()),
            print ' atoms'
            if input.atoms:
                input.atoms = input.atoms + input.atoms_code(traj).atomList()
            else: input.atoms = input.atoms_code(traj).atomList()
        if not input.atoms:
            print ' No atom selection found, taking everything... just in case'
            input.atoms = traj.universe
        else: input.atoms = Collection(input.atoms)
        input.atoms = ghostBusters(input.atoms)
        print ' # atoms in selection: ',len(input.atoms.atomList())
        #
        # Group selection
        #
        if input.groups:
            input.groups, input.reference = parseGroupSelection(types,
                          input.groups,input.reference,verbose=input.verbose)
        if input.groups_code:
            # previous def (if any) overwritten
            # the result returned by groups_code should be consistent
            # with that one from misc.paresGroupSelection
            input.groups, input.reference = input.groups_code(traj)
        #
        # Another piece of general settings
        #
        if input.weights == 'mass':
            weightsList = MassList(traj.universe,input.atoms)
        elif input.weights == 'incoherent':
            if input.deuter:
                weightsList = BincohList(traj.universe, input.atoms, h2d)
            else:
                weightsList = BincohList(traj.universe, input.atoms)
        elif input.weights == 'coherent':
            if input.deuter:
                weightsList = BcohList(traj.universe, input.atoms, h2d)
            else:
                weightsList = BcohList(traj.universe, input.atoms)
        else: weightsList = None

        # input.trajectory = (traj, input.trajectory)
        input.trajectory = traj
        input.q_vector_set = qVectors
        input.weights = weightsList

    return input


def saveText(master=None,data=None,filename=None,title=None):
 
    if not filename:
        if master:
            from nMoldyn.gui import FileDialog
            fd = FileDialog.SaveFileDialog(master)
            filename = fd.go(key='SaveText', pattern='*.txt')
        else:
            raise "Undetectable_System_Error",\
                  "No master for a slave"
    if data is None:
        raise ValueError, "no data to be saved"
    if title is None: title = "pMoldyn data"
    if filename:
        file = TextFile(filename,'w')
        line = '#\n# '+title + '\n#\n'
        file.write(line)
        for i in range(len(data)):
            line = ''
            for ia in range(len(data[i])):
                line = line + str(data[i][ia]) + ' '
            line = line + '\n'
            file.write(line)
        file.close()
    return filename

def saveNetCDF(data,filename,title):
 
    file       = NetCDFFile(filename, 'w', 'Created ')
    file.title = title
    file.createDimension('TIME',len(data[2]))
    file.createDimension('LENGTH',len(data[1]))
    SF      = file.createVariable('SF', N.Float, ('LENGTH','TIME'))
    time    = file.createVariable('time', N.Float, ('TIME',))
    qlenght = file.createVariable('qlenght', N.Float, ('LENGTH',))

    for i in range(len(data[0])):
        for j in range(len(data[1])):
            SF[j,i] = data[2][i][j]
    for i in range(len(data[0])):
        time[i] = data[0][i]
    for i in range(len(data[1])):
        qlenght[i] = data[1][i]
    file.close()

## def PrepBcoh(traj,h2d):
 
##     bcoh  = traj.universe.getAtomScalarArray('b_coherent')
##     for atom in h2d:
##         if atom.symbol == 'H': bcoh[atom.index] = 6.674e-06 # deuter
##     if len(h2d) == len(traj.universe.atomList()):
##         bcoh = bcoh/sqrt((bcoh*bcoh).sumOverParticles())
##     else:
##         a = 0.0
##         for at in h2d: a = a+sqrt(bcoh[at]*bcoh[at])
##         bcoh = bcoh/a
##     return bcoh

## def PrepBincoh(traj,h2d):
 
##     bincoh = traj.universe.getAtomScalarArray('b_incoherent')
##     for atom in h2d:
##         if atom.symbol=='H': bincoh[atom.index] = 4.022e-06 # deuter
##     if len(h2d) == len(traj.universe.atomList()):
##         bincoh =  bincoh/bincoh.sumOverParticles()
##     else:
##         a = 0.0
##         for at in h2d: a = a+bincoh[at.index]
##         bincoh = bincoh/a
##     return bincoh

def tup2dict(tup):
    result = {}
    for i in range(len(tup)):
        result[tup[i][0]] = tup[i][1]
    return result

def getTypes(universe):

    types = {}
    for io in range(len(universe)):
         object = universe[io]
         if object.__class__.__name__ == 'Protein':
              name = 'Protein.'+str(io)
         elif object.__class__.__name__ == 'AtomCluster':
              name = object.name
         else:
              name = object.type.name       
         try:               collection = types[name]
         except KeyError:   collection = Collection()
         collection.addObject(object)
         types[name] = collection
    return types

def getResidues(protein):

    result = {}
    for chain in protein.molecules:
        for res in chain.groups:
            try: result[res.type.name].append(res)
            except: result[res.type.name] = [res]
    return result


def getReference(gr,atom_list,verbose=0):
    """ for given group of atoms construct Configuration object
    containing reference posiotion for each atom in the group """

    verbose = 1

    universe = gr.universe()
    # update in group make-up
    reference = Configuration(universe,cell=universe.cellParameters())
    group = Collection()
    for atom in gr.atomList():
        try:
            idx = atom_list.index(atom)
            reference[atom] = atom_list[idx].position()
            group.addObject(atom)
        except ValueError: # not in list
            if verbose > 0:
                print 'Warning: nMoldyn.misc.getReference'
                print 'Matching problem'
                print atom
    if group.numberOfAtoms() == 0:
        print 'Error: No reference atoms found'
        reference = None
        group = gr

    return group,reference


def parseGroupSelection(types,group_dict,refer_dict,verbose=0):

    if len(group_dict) == 2:
        groups = [ {}, {} ]
        reference = [ {}, {} ]
    else:
        groups = [ {} ]
        reference = [ {} ]

    for ngr in range(len(group_dict)):
        for title in group_dict[ngr].keys():
            master = split(title)
            if len(master) == 2: # Proteins
                mol = types[master[0]][0]
                if master[1] == 'All':
                    groups[ngr][title] = [types[master[0]]] 
                    try:
                        filename = refer_dict[ngr][title]['*']
            # possible extension here: define a RE pattern
            # (instead of '*') to make sub-selection (res-numbers, etc.)
                        if filename:
                            info = parsePDBReference(filename,[mol],
                                                     verbose=verbose)
                            reference[ngr][title] = \
                                     Collection(info['PDB']).atomList()
                        else: reference[ngr][title] = None
                    except:   reference[ngr][title] = None
                elif master[1] == 'SideChain':
                    # SideChain selection was intentionaly made more flexible:
                    # we're picking whole residues (not sidechain attribute)
                    # to be able to define either parts of real sidechains
                    # (aromatic rings,...),
                    # or, for example, select also C_alpha atoms
                    # This means that after entering real calculations
                    # atoms in a group has to be filtered out according to
                    # make-up of the reference structure (if the latter one is
                    # set to None, then use sidechain attribute and use first
                    # configuration as a reference)
                    pat = getResidues(mol)
                    groups[ngr][title] = []
                    for ia in group_dict[ngr][title]:
                        if ia == '*':
                            for g in pat.values():
                                groups[ngr][title] = groups[ngr][title] + g
                            break
                        else:
                            groups[ngr][title] = groups[ngr][title] + pat[ia]
                    reference[ngr][title] = []
                    try:
                        for ia in refer_dict[ngr][title].keys():
                            filename = refer_dict[ngr][title][ia]
                            if filename:
                                info = parsePDBReference(filename,pat[ia],
                                                         verbose=verbose)
                                reference[ngr][title] = reference[ngr][title]+\
                                         Collection(info['PDB']).atomList()
                    except: reference[ngr][title] = None 
                elif master[1] == 'BackBone':
                    pat = getProteinBackbone(mol)
                    groups[ngr][title] = pat
                    try:
                        filename = refer_dict[ngr][title]['*']
                        if filename:
                            info = parsePDBReference(filename,pat,
                                                     verbose=verbose)
                            reference[ngr][title] = \
                                     Collection(info['PDB']).atomList()
                        else:
                            reference[ngr][title] = None
                    except: reference[ngr][title] = None
                elif master[1] == 'Methyl':
                    pat = getMethyls(mol)
                    groups[ngr][title] = pat
                    try:
                        filename = refer_dict[ngr][title]['*']
                        if filename:
                            info = parsePDBReference(filename,pat,
                                                     verbose=verbose)
                            reference[ngr][title] = \
                                     Collection(info['PDB']).atomList()
                        else:
                            reference[ngr][title] = None
                    except: reference[ngr][title] = None
            else:                 # Molecules
                groups[ngr][title] = types[master[0]].objects
                try:
                    filename = refer_dict[ngr][title]['*']
                    if filename:
                        info = parsePDBReference(filename,groups[ngr][title],
                                                 verbose=verbose)
                        reference[ngr][title] = \
                                 Collection(info['PDB']).atomList()
                    else:
                        reference[ngr][title] = None
                except: reference[ngr][title] = None
    return groups, reference


def parsePDBReference(filename,pattern,verbose=None):
    """ for a given filename of a PDB file match atoms against
    objects in a list <pattern> """

    pdb = PDBConfiguration(filename)
    if verbose:
        print 'A quick look into an MMTK pattern reveals ',
        print 'the presence of\n\t',len(pattern),
        print ' atom collections ready to be matched'
    tokens = getTokens(filename)
    if len(tokens) > 0 and verbose:
        print 'RE pattern defined in your PDB file: '
        print tokens
    pdb_collection = {}
    for object in range(len(pattern)):
        natom = pattern[object].numberOfAtoms()
        pdb_natom = 0
        for pdb_object in pdb.objects: 
            if hasattr(pdb_object, 'atom_list'):
                gj = 0
                for atom in pdb_object.atom_list:
                    if atom.properties['element'] == '*': gj=gj+1
                if gj > 0:
                    try: pdb_collection[object].append(pdb_object)
                    except KeyError:
                        pdb_collection[object] = [pdb_object]
                        pdb_natom = pdb_natom + gj
            if hasattr(pdb_object, 'residues'):
                for residue in pdb_object.residues:
                    gj = 0
                    for atom in residue.atom_list:
                        if atom.properties['element'] == '*': gj=gj+1
                    if gj > 0:
                        try: pdb_collection[object].append(residue)
                        except KeyError:
                            pdb_collection[object] = [residue]
                            pdb_natom = pdb_natom + gj
        if pdb_natom < natom:
            print 'Warning: fewer atoms in PDB for object:',
            print pattern[object]
        elif pdb_natom > natom:
            print 'ERROR: the number of atoms does not match'
            print object,len(pattern),len(pdb.objects)
            print pattern[object].atomList()
            print pdb_natom,pdb_collection[object]
            return None

    selection = []
    for object in range(len(pattern)):
        coll = []
        pdbmap = []
        type_name = pattern[object].__class__.__name__
        if verbose:
            print '\n\nPROCESSING ',pattern[object],
            print ' of type: ',type_name
            print     '----------\n'
        if type_name == 'Protein':
            for chain in pattern[object].chains:
                for residue in chain[0]: pdbmap.append(residue)
        elif type_name == 'SubChain':
            for residue in pattern[object]: pdbmap.append(residue)
        elif type_name == 'Residue' or type_name == 'Molecule' or \
             type_name == 'Collection':
            pdbmap.append(pattern[object])

        if type_name == 'Protein' or type_name == 'SubChain' or\
           type_name == 'Residue':
            for it in range(len(pdbmap)):
                pdb_item = pdb_collection[object][it]
                atom_list = pdbmap[it].atomList()
                res_list = [] 
                for ia in pdbmap[it].pdbmap: res_list.append(ia[0])
                res_list = map(upper,res_list)
                if res_list.count(upper(pdb_item.name)) == 0:
                    print 'ERROR: problem with matching'
                    print pdb_item.name,res_list
                    return None
                for ia in pdb_item.atoms.keys():
                    atom = pdb_item.atoms[ia]
                    if atom.properties['element'] == '*':
                        anum = pdbmap[it].pdbmap[0][1][ia].number
                        atra = atom_list[anum]
                        oldpos = atra.position()
                        atra.setPosition(atom.position/10.)
                        newpos = atra.position()
                        coll.append(atra)
                        if verbose:
                            quickCheck(pdb_item,atom,atra,oldpos,newpos)
            selection.append(Collection(coll))
        elif type_name == 'Molecule' or type_name == 'Collection':
            for it in range(len(pdb_collection[object])):
                pdb_item = pdb_collection[object][it]
                atom_list = pdbmap[0].atomList()
                res_list = []
                gj = 0
                for ia in pdbmap[0].pdbmap:
                    if evalREPattern(pdb_item.name,upper(ia[0]),tokens):
                        gj = 1
                        break
                if not gj:
                    print 'ERROR: problem with matching'
                    return None
                atname_list = pdb_item.atoms.keys()
                for ia in range(len(atname_list)):
                    atom = pdb_item.atoms[atname_list[ia]]
                    if atom.properties['element'] == '*' and \
                       evalREPattern(pdb_item.name,
                                     upper(pdbmap[0].pdbmap[it][0]),tokens):
                        mmtk_map = pdbmap[0].pdbmap[it][1]
                        for ib in mmtk_map.keys():
                            #if evalREPattern(atname_list[ia],ib,tokens):
                            if atname_list[ia] == ib:
                                aname = ib
                                break
                        anum = pdbmap[0].pdbmap[it][1][aname].number
                        atra = atom_list[anum]
                        oldpos = atra.position()
                        atra.setPosition(atom.position/10.)
                        newpos = atra.position()
                        coll.append(atra)
                        if verbose:
                            quickCheck(pdb_item,atom,atra,oldpos,newpos)
            selection.append(Collection(coll))
    info = {'MMTK': pattern,'PDB': selection}
    if verbose:
        print len(Collection(selection).atomList()),
        print 'MMTK atoms in ',len(pattern),'objects',
        print 'were mapped successfully on PDB atoms in file ',filename
    return info


def evalREPattern(pdb_pattern,name,tokens):

    re_pat = ''
    for letter in pdb_pattern:
        if tokens.has_key(letter): re_pat = re_pat + tokens[letter]
        else: re_pat = re_pat + letter
    return re.match(re_pat,ljust(name,4))


def quickCheck(pdb_item,pdb_atom,mmtk_atom,oldpos,newpos):

    print 'In PDB residue: ',pdb_item.name,
    print ' atom: ',pdb_atom,' was assigned to ',
    print mmtk_atom,' MMTK atom'
    print 'Previous position of MMTK atom was ',
    print oldpos,
    print ', current position taken from ',
    print 'the reference is ',newpos,'\n'


def getChemicalObjects(types,var_list):
    """ processing atom selections """
    
    myCollection = Collection()

    for k in types.keys():
        if var_list[k] == '*':  # 'All' 
            for molecule in types[k]: myCollection.addObject(molecule)
        elif var_list[k] == 'None': pass
        elif var_list[k] == 'BackBone':
            for molecule in types[k]:
                myCollection.addObject(molecule.backbone())
        elif var_list[k] == 'SideChain':
            for molecule in types[k]:
                reslist = molecule.sidechains()
                for res in reslist: myCollection.addObject(res)
        elif var_list[k] == 'Methyl':
            al = types[k].atomList()
            for atom in al:
                if atom.symbol == 'C':
                    h   = 0
                    al2 = atom.bondedTo()
                    for at2 in al2:
                        if at2.symbol=='H': h = h+1
                    if h == 3:
                        meth = Collection()
                        meth.addObject(atom)
                        for at2 in al2: 
                            if at2.symbol=='H': meth.addObject(at2)
                        myCollection.addObject(meth)
        elif var_list[k] == 'C_alpha':
            al = types[k].atomList()
            for atom in al:
                if atom.name == 'C_alpha': myCollection.addObject(atom)
        elif var_list[k] == 'Carbon' or var_list[k] == 'Phosphorus' or \
             var_list[k] == 'Oxygen'  or var_list[k] == 'Nitrogen' or \
             var_list[k] == 'Hydrogen' or var_list[k] == 'Sulfur':
            al = types[k].atomList()
            # This test is tricky: the second line compares only the
            # first letter, so calcium is treated like carbon!
            # Perhaps this should be rewritten to use chemical element
            # symbols only.
            for atom in al:
                if atom.type.name == lower(var_list[k]) or \
                   upper(atom.type.name[0]) == var_list[k][0]:
                    myCollection.addObject(atom) 
   
    return myCollection


def getAtoms(molecule):

    met=carb=hydr=oxy=nitr=sulf=phos = 0 # met?
    for atom in molecule.atomList():
        if atom.symbol == 'C':
            h = 0
            al2 = atom.bondedTo()
            for at2 in al2:
                if at2.symbol=='H': h = h+1
            if h == 3: met = 1
        if atom.symbol == 'C': carb = 1
        if atom.symbol == 'H': hydr = 1
        if atom.symbol == 'O': oxy = 1
        if atom.symbol == 'N': nitr = 1
        if atom.symbol == 'S': sulf = 1
        if atom.symbol == 'P': phos = 1

    return met,carb,hydr,oxy,nitr,sulf,phos


def getMethyls(molecule):

    al = molecule.atomList()
    result = []
    for atom in al:
        if atom.symbol == 'C':
            h   = 0
            al2 = atom.bondedTo()
            for at2 in al2:
                if at2.symbol=='H': h = h+1
            if h == 3:
                # find a complete pdbmap within a one of parent groups
                object = atom.parent
                try: res_name = object.pdbmap[0][0]
                except AttributeError: res_name = ''
                # an atom can belong only to one residue
                while len(res_name) == 0:
                    object = object.parent
                    try: res_name = object.pdbmap[0]
                    except AttributeError: pass
                atom_number = object.atomList().index(atom)
                aindex = 0
                pdbmap = {}
                for key in object.pdbmap[0][1].keys():
                    if object.pdbmap[0][1][key].number == atom_number:
                        pdbmap[key] = AtomReference(aindex)
                        break
                aindex = aindex + 1
                meth = [atom]
                for at2 in al2: 
                    if at2.symbol=='H':
                        meth.append(at2)
                        atom_number = object.atomList().index(at2)
                        for key in object.pdbmap[0][1].keys():
                            if object.pdbmap[0][1][key].number == atom_number:
                                 pdbmap[key] = AtomReference(aindex)
                                 break
                        aindex = aindex + 1
                myCollection = Collection(meth)
                myCollection.pdbmap = [(res_name,pdbmap)]
                result.append(myCollection)
    return result


def getProteinBackbone(protein):

    result = []
    for chain in protein.chains:
        for residue in chain[0]:
            if dir(residue).count('peptide') > 0: # terminal NHE...
                resmap = residue.pdbmap
                res_name = resmap[0][0]
                bb = residue.backbone().atomList()
                pdbmap = {}
                for atom in range(len(bb)):
                    atom_number = residue.atomList().index(bb[atom])
                    ikey = (map(lambda x: x.number, resmap[0][1].values()))\
                           .index(atom_number)
                    pdbmap[resmap[0][1].keys()[ikey]] = AtomReference(atom)
                myCollection = Collection(bb)
                myCollection.pdbmap = [(res_name,pdbmap)]
                result.append(myCollection)
    return result


def getTokens(filename):
    """ read PDB file and extract RE patterns which are defined in REMARK
    fields """
    file = PDBFile(filename)
    line = (None,None)
    tokens = {'!':'[0-9]*','@':'[A-Z]*','.':'.*','#':'.+'} # default settings
    while line[0] != 'END':
        line = file.readLine()
        if line[0] == 'REMARK':
            x = eval(line[1])
            for key in x.keys(): tokens[key] = x[key]
            # in case of multiple definitions of the same token previous
            # definitions will be overwritten
    return tokens


def ghostBusters(atoms):
    """ apply filtering to an atom List and return collection
    having no repetitions of atoms using a histogram-like approach """
    unique = {}
    ghost = 'No'
    for atom in atoms.atomList():
        try:
            unique[atom] = unique[atom]+1
            ghost = 'Yes'
        except: unique[atom] = 1
    if ghost == 'Yes':
        atsel = []
        for atom in unique.keys():
            if unique[atom] > 1: atsel.append(atom)
    else: atsel = unique.keys()
    atomsNew = Collection(atsel)
    return atomsNew
        
def logfileInit(modName):
    """ create a file with information about currently running calculations
    (progress, start date, etc.) """
    t1       = asctime(localtime(time()))+'\n0\n'
    owner    = getpass.getuser()
    pid      = getpid()
    filename = join([mktemp(), owner, str(pid), modName, 'moldyn'], '.')
    file     = TextFile(filename, 'w')
    file.write(t1)
    file.flush()
    return file, filename

def logfileUpdate(file,i,norm):
    i = i+1
    t = int(100*i/norm)
    file.write(str(t)+"\n")
    file.flush()
    return i

def cleanUp(file,filename):
    file.close()
    unlink(filename)
    return
    
def parsePDBAtomSelection(filename, traj):

    univ = traj.universe
    pdb = PDBConfiguration(filename)

    # find objects in UNIV and PDB with the same number of atoms
    total = 0
    pdb_index = 0
    pdb_collection = {}
    for object in range(len(univ)):
        natom = len(univ[object].atomList())
        pdb_natom = 0
        while pdb_natom < natom:
            pdb_object = pdb.objects[pdb_index]
            if dir(pdb_object).count('atom_list') == 1: # no groups
                try: pdb_collection[object].append(pdb_object)
                except KeyError: pdb_collection[object] = [pdb_object]
                pdb_natom = pdb_natom + len(pdb_object.atom_list)
                                                         # chains ?
            elif dir(pdb_object).count('residues') == 1: # biopolymers
                for residue in pdb_object.residues:
                    try: pdb_collection[object].append(residue)
                    except KeyError: pdb_collection[object] = [residue]
                    pdb_natom = pdb_natom + len(residue.atom_list)
            pdb_index = pdb_index + 1
        if pdb_natom != natom: return None ######### ERROR

    # match PDBMAP names from UNIV against PDB atom names
    # and add selected atoms to the collection
    selection = []
    for object in range(len(univ)):
        pdbmap = []
        if univ[object].__class__.__name__ == 'Protein':
            for chain in univ[object].chains:
                for residue in chain[0]: pdbmap.append(residue)
        elif univ[object].__class__.__name__ == 'Molecule':
            pdbmap.append(univ[object])
        for item in range(len(pdbmap)):
            pdb_item = pdb_collection[object][item]
            atom_list = pdbmap[item].atomList()
            if upper(pdbmap[item].pdbmap[0][0]) != upper(pdb_item.name):
                return None                ######### ERROR
            for ia in pdb_item.atoms.keys():
                atom = pdb_item.atoms[ia]
                if atom.properties['element'] == '*':
                    anum = pdbmap[item].pdbmap[0][1][ia].number
                    atra = atom_list[anum]
                    selection.append(atra)

    return Collection(selection)

def getProgress():
    """ Provide info about currently running calculations.
    The name of a job (eg. the module name), a job ID,
    a job owner, the date of starting
    and a progress in percents are shown """

    list = listdir(gettempdir())
    list2 = []
    owner = getpass.getuser()
    for file in list:
            try:
                a = split(file,'.')
                if a[-1] == 'moldyn' and a[-4] == owner:
                    jobID = int(a[-3])
                    try:
                        kill(jobID, 0)
                        list2.append(file)
                    except:
                        unlink(file)
            except:
                pass
        
    if len(list2) == 0:
        print "\n There are no currently running calculations\n"
    else:
        print '\n  Progress indicator'
        print ' --------------------'
        print '    Module . Job ID .   Job Owner .'+\
              '           Started              . Progress'
        print ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'+\
              '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        for file in list2:
            record = split(file,'.')
            mod = record[-2]
            jobID = int(record[-3])
            jobOwner = record[-4]
            temp = gettempdir()
            filename = path.join(temp, file)
            f = TextFile(filename,'r')
            l = []
            for line in f: l.append(line[:-1])
            if len(l) > 1:
                try:
                    t    = l[0]
                    prog = l[len(l)-1] + ' %'
                    print ' ' + rjust(mod,9),'.' + rjust(str(jobID),7),'.' +\
                          rjust(jobOwner,12),'.' + rjust(t,31),'.' + \
                          rjust(prog,8)
                except: pass
        print


## def uniqueGroups(group_A,group_B,done):
##     """ check whether two groups share any atoms """

##     gj = None
##     gal = group_A.atomList()
##     gbl = group_B.atomList()
##     if done.count([gal[0].index,gbl[0].index]) > 0 or \
##        done.count([gbl[0].index,gal[0].index]) > 0:
##         return None
##     for atom in gal:
##         if gbl.count(atom) > 0:
##             gj = 1
##             break
##     if gj: return None
##     else:
##         done.append([gal[0].index,gbl[0].index])
##         return done

def uniqueGroupPairs_old(grA,grB):
    groups = {}
    for igA in grA:
        for igB in grB:
            if groups.has_key(igA):
                if groups[igA].has_key(igB): pass
                elif igA != igB:
                    try: groups[igA][igB] = 0
                    except: groups[igA] = {igB:0}
            elif groups.has_key(igB):
                if groups[igB].has_key(igA): pass
                elif igA != igB:
                    try: groups[igB][igA] = 0
                    except: groups[igB] = {igA:0}
            elif igA != igB: groups[igA] = {igB:0}
    result = {}
    for igA in groups.keys():
        result[igA] = []
        for igB in groups[igA].keys(): result[igA].append(igB)
    return result

def getVariables(traj,key):

    minor = 0
    try:
        res = traj.trajectory.file.variables[key]
        if traj.trajectory.file.dimensions.keys().count('minor_step_number'):
            minor = traj.trajectory.file.dimensions['minor_step_number']
    except None:
        res=[]
        for t in traj.trajectories:
            res.append(t.trajectory.file.variables[key])
        res = concatenate(res,0)
        if traj.trajectories[0].file.\
           dimensions.keys().count('minor_step_number'):
            minor = traj.trajectory.file.dimensions['minor_step_number']
    if minor > 0: # i know it's ugly
        if len(res.shape) == 3: # particleTrajectory, box
            res = reshape(swapaxes(res,1,2),
                         (res.shape[0]*res.shape[2],res.shape[1]))
        elif len(res.shape) == 2: # time, step
            res = reshape(res,(res.shape[0]*res.shape[1],))
    return res



def uniqueGroupPairs(groups_A,groups_B,refGroup_A,refGroup_B):
    """ Find possible pairs formed by groups in grA and grB.
    Make a list of [ (group_1,reference_1), (group_2,reference_2) ] """
    groups = {}
    gal = []; gbl = []
    itm = 0
    for key_A in groups_A.keys():
        for igA in groups_A[key_A]:
            if refGroup_A.has_key(key_A):
                gal.append(tuple(getReference(igA,refGroup_A[key_A])))
            else: gal.append((igA, None))
    if groups_A != groups_B and refGroup_A != refGroup_B:
        for key_B in groups_B.keys():
            for igB in groups_B[key_B]:
                if refGroup_B.has_key(key_B):
                    gbl.append(tuple(getReference(igB,refGroup_B[key_B])))
                else: gbl.append((igB,None))
    else: gbl = gal
    for itA in gal:
        for itB in gbl:
            if itA != itB:
                if groups.has_key(itA):
                    try: groups[itA][itB] = 0
                    except: groups[itA] = {itB:0}
                elif groups.has_key(itB):
                    try: groups[itB][itA] = 0
                    except: groups[itB] = {itA:0}
                else:
                    groups[itA] = {itB:0}
    result = {}
    for itA in groups.keys():
        result[itA] = groups[itA].keys()
        result[itA].sort()
    return result
