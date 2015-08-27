"""Collections of classes for the determination of NMR-related properties.

Classes:

    * OrderParameter             : sets up an order parameter analysis.
    * OrderParameterContactModel : sets up an order parameter analysis using the contact model approach.
"""

# The python distribution modules
import copy
import math
from numpy import linalg
import operator
import os
import sys
from time import asctime
from timeit import default_timer

# The ScientificPython modules
from Scientific import N 
from Scientific.Geometry import Vector
from Scientific.IO.NetCDF import NetCDFFile

# The MMTK distribution modules
from MMTK.Collections import Collection
from MMTK import Units
from MMTK.NucleicAcids import NucleotideChain
from MMTK.Proteins import PeptideChain, Protein

# The nMOLDYN modules
from order_parameter import order_parameter
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Analysis.Analysis import Analysis
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import convertNetCDFToASCII
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Mathematics.Analysis import correlation

#####################################################################################
# ORDER PARAMETER ANALYSIS
#####################################################################################
class OrderParameter(Analysis):
    """Sets up an order parameter analysis.

    A Subclass of nMOLDYN.Analysis.Analysis. 

    Constructor: OrderParameter(|parameters| = None)

    Arguments:

        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory         -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo           -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                    number to consider, 'last' is an integer specifying the last frame number to consider and 
                                    'step' is an integer specifying the step number between two frames.
            * group              -- a selection string specifying the groups of atoms that will define the vectors on which the 
                                    analysis will be computed. Each group must contain two and only two atoms.
            * atomorder          -- a string of the form 'atom1,atom2,atom3' where 'atom1', 'atom2' and 'atom3' are 
                                    respectively the MMTK atom names of the atoms in the way they should be ordered.
            * pyroserver         -- a string specifying if Pyro will be used and how to run the analysis.
            * output             -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                                    instead of the '.nc' extension.

    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.


    Comments:
    
        - This code is based on a first implementation made by Vania Calandrini.
    """

    def __init__(self, parameters = None, statusBar = None):
        """The constructor.
        
        @param parameters: if not None, a dictionnary storing the input parameters names and their corresponding values.
        @type parameters: dict
        
        @param statusBar: if not None, an instance of nMOLDYN.GUI.Widgets.StatusBar. Will attach a status bar to the 
            selected analysis.
        @type statusBar: instance of nMOLDYN.GUI.Widgets.StatusBar
        """
        
        # The inheritance.
        Analysis.__init__(self, parameters, statusBar)

        path = os.path.join(GVAR['nmoldyn_analysis'],'OP')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
        
        self.bondNames = {}
        
        for aIndexes in self.group:

            for at in self.universe.atomList():
                
                if at.index == aIndexes[0]:
                                        
                    if isinstance(at.topLevelChemicalObject(), (Protein, PeptideChain, NucleotideChain)):
                        self.bondNames[tuple(aIndexes)] = at.parent.parent.sequence_number
            
                    else:
                        self.bondNames[tuple(aIndexes)] = at.index
                        
                    break

        self.referenceDirection = Vector(0.0,0.0,1.0)            
            
        # The results are stored in dictionnary because the order in which the bonds are treated does not
        # always follow the structure.
        self.P2 = {}
        self.S2 = {}
                        
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        self.buildTimeInfo()

        self.group = self.selectGroups(self.groupDefinition, self.atomOrder)
        
        self.nGroups = len(self.group)

        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                
    def calc(self, atomIndexes, trajectory):
        """Calculates the contribution for one group.
        
        @param bondIndex: the index of the group in |self.group| list.
        @type bondIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        atoms = [None]*2

        for at in trajectory.universe.atomList():
            if at.index in atomIndexes:
                atoms[atomIndexes.index(at.index)] = at
        
        # The atoms forming the bond.
        at1, at2 = atoms
        
        at1Traj = trajectory.readParticleTrajectory(at1, first = self.first, last = self.last, skip = self.skip)
        at2Traj = trajectory.readParticleTrajectory(at2, first = self.first, last = self.last, skip = self.skip)
                
        costheta = N.zeros((self.nFrames,), typecode = N.Float)
        sinphi   = N.zeros((self.nFrames,), typecode = N.Float)
        cosphi   = N.zeros((self.nFrames,), typecode = N.Float)
        sintheta = N.zeros((self.nFrames,), typecode = N.Float)

        for comp in range(self.nFrames): 

            pVect = Vector(trajectory.universe.distanceVector(at1Traj[comp], at2Traj[comp]))
            pVect = pVect.normal()
            
            costheta[comp] = pVect * self.referenceDirection
            sintheta[comp] = pVect.cross(self.referenceDirection).length()
            cosphi[comp] = (pVect[0]/sintheta[comp])
            sinphi[comp] = (pVect[1]/sintheta[comp])

        tr2         = 3.0*costheta**2 - 1.
        cos2phi     = 2.0*cosphi**2-1.
        sin2phi     = 2.0*sinphi*cosphi
        cossintheta = costheta*sintheta
        sintheta_sq = sintheta**2
        
        # calcul de <P2(cos(theta))> en termes de somme de fonctions de correlation
        # d'harmoniques spheriques (theoreme d'addition des harmoniques spheriques).
        P2 = (0.25*correlation(tr2) + 3.00*correlation(cosphi*cossintheta) + \
              3.00*correlation(sinphi*cossintheta) + 0.75*correlation(cos2phi*sintheta_sq) + \
              0.75*correlation(sin2phi*sintheta_sq))

        # calcul du parametre d'ordre S^2 (limite pour t->infini de <P2(cos(theta))>).
        S2 = (0.75 * (N.sum(cos2phi*sintheta_sq)**2 + N.sum(sin2phi*sintheta_sq)**2) + \
              3.00 * (N.sum(cosphi*cossintheta)**2 + N.sum(sinphi*cossintheta)**2) +
              0.25 * N.sum(tr2)**2) / self.nFrames**2

        return atomIndexes, (P2, S2)
    
    def combine(self, atomIndexes, x):
        """
        """
        
        bondName = self.bondNames[tuple(atomIndexes)]
        
        # calcul de <P2(cos(theta))> en termes de somme de fonctions de correlation
        # d'harmoniques spheriques (theoreme d'addition des harmoniques spheriques).
        self.P2[bondName] = x[0]

        # calcul du parametre d'ordre S^2 (limite pour t->infini de <P2(cos(theta))>).
        self.S2[bondName] = x[1]
                        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
    
        outputFile = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        outputFile.createDimension('NGROUPS', self.nGroups)
        outputFile.createDimension('NFRAMES', self.nFrames)

        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times
        TIMES.units = 'ps'

        GROUPNUMBER = outputFile.createVariable('group_number', N.Int32, ('NGROUPS',))
        P2 = outputFile.createVariable('p2', N.Float, ('NGROUPS','NFRAMES'))
        P2AVG = outputFile.createVariable('p2-groupavg', N.Float, ('NFRAMES',))
        S2   = outputFile.createVariable('s2', N.Float, ('NGROUPS',))

        p2Avg = N.zeros((self.nFrames), typecode = N.Float)

        comp = 0
        for bKey in sorted(self.bondNames.keys()):
            
            bName = self.bondNames[bKey]
            
            GROUPNUMBER[comp] = bName
            S2[comp] = self.S2[bName]
            P2[comp,:] = self.P2[bName]
            N.add(p2Avg, self.P2[bName], p2Avg)
            comp += 1

        P2AVG[:] = p2Avg/float(self.nGroups)

        asciiVar = sorted(outputFile.variables.keys())            

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'pair', 'yVar' : 'S2'}

        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)

#####################################################################################
# ORDER PARAMETER ANALYSIS USING CONTACT MODEL
#####################################################################################
class OrderParameterContactModel(Analysis):
    """Sets up an order parameter analysis using the contact model .

    A Subclass of nMOLDYN.Analysis.Analysis. 

    Constructor: OrderParameterContactModel(|parameters| = None)

    Arguments:

        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo   -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                            number to consider, 'last' is an integer specifying the last frame number to consider and 
                            'step' is an integer specifying the step number between two frames.
            * pyroserver -- a string specifying if Pyro will be used and how to run the analysis.
            * output     -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                            instead of the '.nc' extension.

    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.

    Comments:        
        
        - This code is adapted from the s2predict code developped by F. Zhang and R. Bruschweiler and available in:
          http://nmr.clarku.edu/software/S2/s2predict.html

        - For more details about the method: Zhang, F., Bruschweiler, R. J. AM. Chem. Soc. 2002, 124, 12654-12655.
    """

    def __init__(self, parameters = None, statusBar = None):
        """The constructor.
        
        @param parameters: if not None, a dictionnary storing the input parameters names and their corresponding values.
        @type parameters: dict
        
        @param statusBar: if not None, an instance of nMOLDYN.GUI.Widgets.StatusBar. Will attach a status bar to the 
            selected analysis.
        @type statusBar: instance of nMOLDYN.GUI.Widgets.StatusBar
        """
        
        # The inheritance.
        Analysis.__init__(self, parameters, statusBar)

        path = os.path.join(GVAR['nmoldyn_analysis'],'OPCM')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
        
        if self.universe.basisVectors() is None:
            raise Error('The universe must be periodic for this kind of calculation.')
        
        self.S2 = {}
        self.sequence = {}
        self.hnLookup = {}

        self.nChemicalObjects = len(self.universe.objectList())        
        chemicalObjects = self.universe.objectList()

        for comp in range(self.nChemicalObjects):
                    
            obj = chemicalObjects[comp]
            objName = self.chemicalObjectNames[comp]
            
            if self.chemicalObjectInfo[objName]['objectclass'] in ['PeptideChain', 'Protein']:

                self.sequence[objName] = []
                self.hnLookup[objName] = []
                
                comp = 1
                for r in obj.residues()[1:]:
                    if r.symbol.lower() != 'pro':
                        self.sequence[objName].append(r.sequence_number)
                        self.hnLookup[objName].append(comp)
                    comp += 1
                    
                self.S2[objName] = N.zeros((len(self.sequence[objName]),self.nFrames), typecode = N.Float)

        # Case where no protein or peptide chain was found in the universe.
        if not self.sequence:
            raise Error('The universe must contains at least one peptide chain to perform the analysis.')

        self.scaleHydrogenPos = N.zeros((3,), typecode = N.Float)
        self.scaleOxygenPos = N.zeros((3,), typecode = N.Float)
                                            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        self.buildTimeInfo()
        
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')

    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one group.
        
        @param frameIndex: the index of the frame in |self.frameIndexes| array.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
        
        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
                    
        directCell = N.ravel(N.array([v for v in trajectory.universe.basisVectors()], typecode = N.Float))
        reverseCell = N.ravel(N.transpose(N.array([v for v in trajectory.universe.reciprocalBasisVectors()], typecode = N.Float)))
         
        S2 = {}

        chemicalObjects = trajectory.universe.objectList()

        for comp in range(self.nChemicalObjects):
                    
            obj = chemicalObjects[comp]
            objName = self.chemicalObjectNames[comp]
            
            if self.chemicalObjectInfo[objName]['objectclass'] in ['PeptideChain', 'Protein']:
                        
                S2[objName] = N.zeros((53,), typecode = N.Float)
                S2[objName] = N.zeros((len(self.sequence[objName]),), typecode = N.Float)
                
                for v in range(len(self.hnLookup[objName])):
                    
                    aa_no = self.hnLookup[objName][v]
                    aa_no_minus = aa_no - 1

                    H = obj.backbone()[aa_no].H
                    O = obj.backbone()[aa_no_minus].O

                    heavyAtoms = Collection()
                    for m in range(len(obj.residues())):
                        if (m != aa_no) and (m != aa_no_minus):
                            for atom in obj.residues()[m].atomList():
                                if atom.symbol != 'H': 
                                    heavyAtoms.addObject(atom)

                    # The indexes of the selected atoms.
                    indexes = N.zeros((heavyAtoms.numberOfAtoms(),), typecode = N.Int32)
                    comp = 0
                    for at in heavyAtoms.atomList():
                        indexes[comp] = at.index
                        comp += 1

                    scaleconfig = N.zeros((3*heavyAtoms.numberOfAtoms(),), typecode = N.Float)

                    val = order_parameter(trajectory.universe.contiguousObjectConfiguration().array,\
                                          H.index, O.index, directCell, reverseCell, indexes,\
                                          scaleconfig, self.scaleHydrogenPos, self.scaleOxygenPos,\
                                          Units.Ang)

                    S2[objName][v] = val
            
        return frameIndex, S2
    
    def combine(self, frameIndex, x):
        """
        """
        
        comp = self.frameIndexes.index(frameIndex)
        
        for oName in x.keys():
            self.S2[oName][:,comp] = x[oName]

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    

        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        outputFile.createDimension('NFRAMES', self.nFrames)
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times[:]
        TIMES.units = 'ps'

        for oName, seq in self.sequence.items():
            outputFile.createDimension('SEQ%s' % oName, len(seq))

            SEQUENCE = outputFile.createVariable('%s_sequence' % oName, N.Int32, ('SEQ%s' % oName,))
            SEQUENCE[:] = N.array(seq)
            SEQUENCE.units = 'unitless'

            S2 = outputFile.createVariable('%s_s2' % oName, N.Float, ('SEQ%s' % oName,'NFRAMES'))
            S2[:] = self.S2[oName][:,:]
            S2.units = 'unitless'

            S2AVG = outputFile.createVariable('%s_s2_timeavg' % oName, N.Float, ('SEQ%s' % oName,))
            S2AVG[:] = self.S2[oName][:,:].sum(1)/float(self.nFrames)
            S2AVG.units = 'unitless'

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = None

        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
