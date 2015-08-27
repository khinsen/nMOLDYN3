"""Collections of classes for the determination of structure-related properties.

Classes:

    * PairDistributionFunction     : sets up a Pair Distribution Function Analysis.
    * CoordinationNumber           : sets up a Coordination Number Analysis.
    * DensityProfile               : sets up a Density Profile Analysis.
    * SpatialDensity               : sets up a Spatial Density Analysis.
    * ScrewFit                     : sets up a Screw Fit Analysis.
    * HydrogenBondSurvivalAnalysis : sets up a Hydrogen Bond Survival Analaysis.
"""

# The python distribution modules
import copy
from numpy import linalg, linspace
import operator
import os
import sys
from time import asctime
from timeit import default_timer

# The ScientificPython modules
from Scientific import N 
from Scientific.Geometry import Quaternion, Transformation, Vector
from Scientific.Geometry.Objects3D import Line 
from Scientific.IO.NetCDF import NetCDFFile

# The MMTK distribution modules
from MMTK import Atom
from MMTK.Collections import Collection
from MMTK import Units
from MMTK.NucleicAcids import NucleotideChain
from MMTK.Proteins import PeptideChain, Protein
from MMTK.Trajectory import Trajectory
from MMTK.Universe import InfiniteUniverse

# The nMOLDYN modules
from coordination_number import coordination_number
from distance_histogram import distance_histogram
from hbond_detection import hbond_detection
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Analysis.Analysis import Analysis
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import convertNetCDFToASCII
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Chemistry.Chemistry import detect_hb_atoms_subset
from nMOLDYN.Mathematics.Analysis import correlation
from nMOLDYN.Mathematics.Geometry import changeBasis, sphericalCoordinates
from nMOLDYN.Core.Misc import parseInterval

#####################################################################################
# PAIR DISTRIBUTION FUNCTION ANALYSIS
#####################################################################################
class PairDistributionFunction(Analysis):
    """Sets up a Pair Distribution Function analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: PairDistributionFunction(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory    -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo      -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                               number to consider, 'last' is an integer specifying the last frame number to consider and 
                               'step' is an integer specifying the step number between two frames.
            * rvalues       -- a string of the form 'rmin:rmax:dr' where 'rmin' is a float specifying the minimum distance to 
                               consider, 'rmax' is a float specifying the maximum distance value to consider and 'dr' is a float
                               specifying the distance increment. 
            * subset        -- a selection string specifying the atoms to consider for the analysis.
            * deuteration   -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights       -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                               scheme to use.
            * distanceunits -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.                             
            * pyroserver    -- a string specifying if Pyro will be used and how to run the analysis.
            * output        -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                               instead of the '.nc' extension.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.

    Comments:        
        
        - This code contains a pyrex function for the distance histogram calculation that is based on a FORTRAN code 
          written by Miguel Gonzalez, Insitut Laue Langevin, Grenoble, France.
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'PDF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                    
        # The histogram of the intramolecular distances.
        self.hIntra = N.zeros((self.nElements,self.nElements,self.nRBins), typecode = N.Float32)
        # The histogram of the intermolecular distances.
        self.hInter = N.zeros((self.nElements,self.nElements,self.nRBins), typecode = N.Float32)
        
        # ...
        self.scaleconfig = N.zeros((3*self.nSelectedAtoms,), typecode = N.Float)

        self.averageDensity = 0.0
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        # Some additional checkings.
        if self.universe.basisVectors() is None:
            raise Error('The universe must be periodic for this kind of calculation.')
                
        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()

        # The list of the species name found in the universe.
        self.elementsList = sorted(self.elementInformation.keys())

        # The number of elements found in the universe.
        self.nElements = len(self.elementsList)

        # The id of the specie to which belong each selected atom.
        self.elements = N.zeros((self.nSelectedAtoms,), typecode = N.Int32)

        # The id of the molecule to which belongs each selected atoms.
        self.molecules = N.zeros((self.nSelectedAtoms,), typecode = N.Int32)

        molList = []
        for comp in range(self.nSelectedAtoms):
            aIndex = self.subset[comp]
            
            element = self.atomInformation[aIndex]['element']
            self.elements[comp] = self.elementsList.index(element)
            
            molecule = self.atomInformation[aIndex]['molecule']

            if molecule not in molList:
                molList.append(molecule)
            self.molecules[comp] = molList.index(molecule)
                        
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                            
    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame in |self.frameIndexes| array.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
                
        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
                        
        directCell = N.ravel(N.array([v for v in trajectory.universe.basisVectors()], typecode = N.Float))
        reverseCell = N.ravel(N.transpose(N.array([v for v in trajectory.universe.reciprocalBasisVectors()], typecode = N.Float)))
        
        cellVolume = trajectory.universe.cellVolume()
            
        hIntraTemp = N.zeros(self.hIntra.shape, typecode = N.Float32)
        hInterTemp = N.zeros(self.hInter.shape, typecode = N.Float32)
        
        distance_histogram(trajectory.universe.contiguousObjectConfiguration().array,\
                             directCell, reverseCell,\
                             N.array(self.subset, typecode = N.Int32),\
                             self.molecules, self.elements, hIntraTemp, hInterTemp,\
                             self.scaleconfig, self.rMin, self.dR)

        N.multiply(hIntraTemp, cellVolume, hIntraTemp)
        N.multiply(hInterTemp, cellVolume, hInterTemp)
                                                                                                
        return frameIndex, (cellVolume, hIntraTemp, hInterTemp)
    
    def combine(self, frameIndex, x):
        """
        """
        
        self.averageDensity += self.nSelectedAtoms/x[0]

        # The temporary distance histograms are normalized by the volume. This is done for each step because the
        # volume can variate during the MD (e.g. NPT conditions). This volume is the one that intervene in the density
        # calculation.
        N.add(self.hIntra, x[1], self.hIntra)
        N.add(self.hInter, x[2], self.hInter)
                                    
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """
        
        self.scaleconfig = None
        
        self.rValues = self.rValues[:-1] + self.dR/2.0

        self.averageDensity /= self.nFrames
        densityFactor = 4.0*N.pi*self.rValues
        shellSurfaces = densityFactor*self.rValues
        shellVolumes  = shellSurfaces*self.dR

        # The intramolecular TCF.
        tcfIntra = {}
        # The intermolecular TCF.
        tcfInter = {}
        # The total TCF.
        tcfTotal = {'total' : N.zeros((self.nRBins), typecode = N.Float)}

        # The intramolecular RDF.
        rdfIntra = {}
        # The intermolecular RDF.
        rdfInter = {}
        # The total RDF.
        rdfTotal = {'total' : N.zeros((self.nRBins), typecode = N.Float)}
            
        # The intramolecular RDF.
        pdfIntra = {}
        # The intermolecular RDF.
        pdfInter = {}
        # The total RDF.
        pdfTotal = {'total' : N.zeros((self.nRBins), typecode = N.Float)}
        
        test = 0.0
        
        for i in range(self.nElements):
            
            element1 = self.elementsList[i]

            nA = self.elementInformation[element1]['number']
            wA = self.elementInformation[element1]['weight']            

            for j in range(i, self.nElements):
                
                element2 = self.elementsList[j]

                nB = self.elementInformation[element2]['number']
                wB = self.elementInformation[element2]['weight']
                
                pair = element1 + '-' + element2

                if i == j:
                    
                    pdfIntra[pair] = self.hIntra[i,j,:].copy()
                    pdfInter[pair] = self.hInter[i,j,:].copy()
                    
                    wAB = wA*wA*nA*nA
                    nAB = nA*(nA-1)/2.0

                else:
                    pdfIntra[pair] = self.hIntra[i,j,:] + self.hIntra[j,i,:]
                    pdfInter[pair] = self.hInter[i,j,:] + self.hInter[j,i,:]
                    
                    wAB = 2.0*wA*wB*nA*nB
                    nAB = nA*nB

                pdfIntra[pair] /= (float(self.nFrames)*shellVolumes*nAB)                                    
                pdfInter[pair] /= (float(self.nFrames)*shellVolumes*nAB)
                pdfTotal[pair] = pdfIntra[pair] + pdfInter[pair]

                rdfIntra[pair] = shellSurfaces*self.averageDensity*pdfIntra[pair]
                rdfInter[pair] = shellSurfaces*self.averageDensity*pdfInter[pair]
                rdfTotal[pair] = rdfIntra[pair] + rdfInter[pair]

                test += wAB

                tcfIntra[pair] = densityFactor*self.averageDensity*(pdfIntra[pair] - 1.0)
                tcfInter[pair] = densityFactor*self.averageDensity*(pdfInter[pair] - 1.0)
                tcfTotal[pair] = tcfIntra[pair] + tcfInter[pair]
                                                                
                N.add(rdfTotal['total'],wAB*rdfTotal[pair],rdfTotal['total'])
                N.add(pdfTotal['total'],wAB*pdfTotal[pair],pdfTotal['total'])
                N.add(tcfTotal['total'],wAB*tcfTotal[pair],tcfTotal['total'])
                                    
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NRVALUES', self.nRBins)

        # Creation of the NetCDF output variables.
        # The Q.
        RVALUES = outputFile.createVariable('r', N.Float, ('NRVALUES',))
        RVALUES[:] = self.distanceConv*self.rValues
        RVALUES.units = self.distanceUnits

        for k in pdfIntra.keys():

            PDFINTRA = outputFile.createVariable('pdf-%s-intra' % k, N.Float, ('NRVALUES',))
            PDFINTRA[:] = pdfIntra[k]
            PDFINTRA.units = 'nm^3'

            PDFINTER = outputFile.createVariable('pdf-%s-inter' % k, N.Float, ('NRVALUES',))
            PDFINTER[:] = pdfInter[k]
            PDFINTER.units = 'nm^3'

            PDFTOTAL = outputFile.createVariable('pdf-%s' % k, N.Float, ('NRVALUES',))
            PDFTOTAL[:] = pdfTotal[k]
            PDFTOTAL.units = 'nm^3'

            RDFINTRA = outputFile.createVariable('rdf-%s-intra' % k, N.Float, ('NRVALUES',))
            RDFINTRA[:] = rdfIntra[k]
            RDFINTRA.units = 'unitless'

            RDFINTER = outputFile.createVariable('rdf-%s-inter' % k, N.Float, ('NRVALUES',))
            RDFINTER[:] = rdfInter[k]
            RDFINTER.units = 'unitless'

            RDFTOTAL = outputFile.createVariable('rdf-%s' % k, N.Float, ('NRVALUES',))
            RDFTOTAL[:] = rdfTotal[k]
            RDFTOTAL.units = 'unitless'

            TCFINTRA = outputFile.createVariable('tcf-%s-intra' % k, N.Float, ('NRVALUES',))
            TCFINTRA[:] = tcfIntra[k]
            TCFINTRA.units = 'unitless'

            TCFINTER = outputFile.createVariable('tcf-%s-inter' % k, N.Float, ('NRVALUES',))
            TCFINTER[:] = tcfInter[k]
            TCFINTER.units = 'unitless'

            TCFTOTAL = outputFile.createVariable('tcf-%s' % k, N.Float, ('NRVALUES',))
            TCFTOTAL[:] = tcfTotal[k]
            TCFTOTAL.units = 'unitless'            
            
        PDF = outputFile.createVariable('pdf-total', N.Float, ('NRVALUES',))
        PDF[:] = pdfTotal['total']
        PDF.units = 'nm^3'
        
        RDF = outputFile.createVariable('rdf-total', N.Float, ('NRVALUES',))
        RDF[:] = rdfTotal['total']
        RDF.units = 'unitless'

        TCF = outputFile.createVariable('tcf-total', N.Float, ('NRVALUES',))
        TCF[:] = tcfTotal['total']
        TCF.units = 'unitless'

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'r', 'yVar' : 'pdf-total'} 
        
        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)

#####################################################################################
# COORDINATION NUMBER ANALYSIS
#####################################################################################
class CoordinationNumber(Analysis):
    """Sets up a Coordination Number analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: CoordinationNumber(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory    -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo      -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                               number to consider, 'last' is an integer specifying the last frame number to consider and 
                               'step' is an integer specifying the step number between two frames.
            * rvalues       -- a string of the form 'rmin:rmax:dr' where 'rmin' is a float specifying the minimum distance to 
                               consider, 'rmax' is a float specifying the maximum distance value to consider and 'dr' is a float
                               specifying the distance increment. 
            * group         -- a selection string specifying the groups of atoms that will be used to define the points around which 
                               the coordination number will be computed. For each group, there is one point defined as the center of 
                               gravity of the group.
            * subset        -- a selection string specifying the atoms to consider for the analysis.
            * deuteration   -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * timeunits     -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * distanceunits -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.                             
            * pyroserver    -- a string specifying if Pyro will be used and how to run the analysis.
            * output        -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                               instead of the '.nc' extension.
    
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.

    Comments:        
        
        - This code contains a pyrex function for the distance histogram calculation than enhances significantly its 
          performance.
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'CN')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                            
        # The histogram of the intramolecular distances.
        self.hIntra = N.zeros((self.nElements,self.nRBins,self.nFrames), typecode = N.Float32)
        # The histogram of the intermolecular distances.
        self.hInter = N.zeros((self.nElements,self.nRBins,self.nFrames), typecode = N.Float32)

        # ...
        self.scaleconfig = N.zeros((3*self.nSelectedAtoms,), typecode = N.Float)
        self.groupcenter = N.zeros((3,), typecode = N.Float)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        if self.universe.basisVectors() is None:
            raise Error('The universe must be periodic for this kind of calculation.')
                
        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)

        self.setElementInformation()

        # The list of the species name found in the universe.
        self.elementsList = sorted(self.elementInformation.keys())

        # The number of elements found in the universe.
        self.nElements = len(self.elementsList)

        # The id of the specie to which belong each selected atom.
        self.elements = N.zeros((self.nSelectedAtoms,), typecode = N.Int32)
                
        # The id of the molecule to which belongs each selected atoms.
        self.molecules = N.zeros((self.nSelectedAtoms,), typecode = N.Int32)

        molList = list(set([self.atomInformation[at.index]['molecule'] for at in self.trajectory.universe.atomList()]))
        
        for comp in range(self.nSelectedAtoms):
            aIndex = self.subset[comp]
            
            element = self.atomInformation[aIndex]['element']
            self.elements[comp] = self.elementsList.index(element)
            
            molecule = self.atomInformation[aIndex]['molecule']
                
            self.molecules[comp] = molList.index(molecule)
                                    
        self.group = self.selectGroups(self.groupDefinition)
        self.nGroups = len(self.group)
                                        
        self.gMolId = []
        for g in self.group:
            self.gMolId.append(molList.index(self.atomInformation[g[0]]['molecule']))
            
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
        
    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame in |self.frameIndexes| array.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
        
        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
            
        directCell = N.ravel(N.array([v for v in trajectory.universe.basisVectors()], typecode = N.Float))
        reverseCell = N.ravel(N.transpose(N.array([v for v in trajectory.universe.reciprocalBasisVectors()], typecode = N.Float)))
                
        hIntraTemp = N.zeros((self.nElements,self.nRBins), typecode = N.Float32)
        hInterTemp = N.zeros((self.nElements,self.nRBins), typecode = N.Float32)
                        
        for comp in range(self.nGroups):
            coordination_number(trajectory.universe.contiguousObjectConfiguration().array,\
                                directCell, reverseCell, N.array(self.group[comp], typecode = N.Int32),\
                                self.gMolId[comp], N.array(self.subset, typecode = N.Int32),\
                                self.molecules, self.elements, hIntraTemp, hInterTemp,\
                                self.scaleconfig, self.groupcenter, self.rMin, self.dR)
                                
        return frameIndex, (hIntraTemp, hInterTemp)

    def combine(self, frameIndex, x):
        """
        """

        comp = self.frameIndexes.index(frameIndex)

        N.add(self.hIntra[:,:,comp], x[0], self.hIntra[:,:,comp])
        N.add(self.hInter[:,:,comp], x[1], self.hInter[:,:,comp])

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """
                                      
        self.scaleconfig = None
        
        coordNumberIntra = {}
        coordNumberInter = {}
        coordNumberTotal = {'total' : N.zeros((self.nRBins, self.nFrames), typecode = N.Float)}
        for i in range(self.nElements):
            element = self.elementsList[i]

            coordNumberIntra[element] = self.hIntra[i,:,:]/float(self.nGroups)
            coordNumberInter[element] = self.hInter[i,:,:]/float(self.nGroups)

            coordNumberTotal[element] = coordNumberIntra[element] + coordNumberInter[element]
                
            N.add(coordNumberTotal['total'], coordNumberTotal[element], coordNumberTotal['total'])
                                            
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NRVALUES', self.nRBins)
        outputFile.createDimension('NTIMES', self.nFrames)
        
        self.rValues = self.rValues[:-1] + self.dR/2.0

        # Creation of the NetCDF output variables.
        # The Q.
        RVALUES = outputFile.createVariable('r', N.Float, ('NRVALUES',))
        RVALUES[:] = self.distanceConv*self.rValues
        RVALUES.units = self.distanceUnits

        TIMES = outputFile.createVariable('time', N.Float, ('NTIMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        for k in coordNumberIntra.keys():

            CNINTRA = outputFile.createVariable('cn-%s-intra' % k, N.Float, ('NRVALUES','NTIMES'))
            CNINTRA[:] = coordNumberIntra[k]
            CNINTRA.units = 'unitless'

            CNINTER = outputFile.createVariable('cn-%s-inter' % k, N.Float, ('NRVALUES','NTIMES'))
            CNINTER[:] = coordNumberInter[k]
            CNINTER.units = 'unitless'

            CNTOTAL = outputFile.createVariable('cn-%s' % k, N.Float, ('NRVALUES','NTIMES'))
            CNTOTAL[:] = coordNumberTotal[k]
            CNTOTAL.units = 'unitless'

            CNCUMULINTRA = outputFile.createVariable('cn-cumul-%s-intra' % k, N.Float, ('NRVALUES','NTIMES'))
            CNCUMULINTRA[:] = N.cumsum(coordNumberIntra[k])
            CNCUMULINTRA.units = 'unitless'

            CNCUMULINTER = outputFile.createVariable('cn-cumul-%s-inter' % k, N.Float, ('NRVALUES','NTIMES'))
            CNCUMULINTER[:] = N.cumsum(coordNumberInter[k])
            CNCUMULINTER.units = 'unitless'

            CNCUMULTOTAL = outputFile.createVariable('cn-cumul-%s' % k, N.Float, ('NRVALUES','NTIMES'))
            CNCUMULTOTAL[:] = N.cumsum(coordNumberTotal[k])
            CNCUMULTOTAL.units = 'unitless'
            
            CNAVGINTRA = outputFile.createVariable('cn-timeavg-%s-intra' % k, N.Float, ('NRVALUES',))
            CNAVGINTRA[:] = coordNumberIntra[k].sum(1)/float(self.nFrames)
            CNAVGINTRA.units = 'unitless'

            CNAVGINTER = outputFile.createVariable('cn-timeavg-%s-inter' % k, N.Float, ('NRVALUES',))
            CNAVGINTER[:] = coordNumberInter[k].sum(1)/float(self.nFrames)
            CNAVGINTER.units = 'unitless'

            CNAVGTOTAL = outputFile.createVariable('cn-timeavg-%s' % k, N.Float, ('NRVALUES',))
            CNAVGTOTAL[:] = coordNumberTotal[k].sum(1)/float(self.nFrames)
            CNAVGTOTAL.units = 'unitless'

            CNAVGCUMULINTRA = outputFile.createVariable('cn-timeavg-cumul-%s-intra' % k, N.Float, ('NRVALUES',))
            CNAVGCUMULINTRA[:] = N.cumsum(coordNumberIntra[k].sum(1))/float(self.nFrames)
            CNAVGCUMULINTRA.units = 'unitless'
            
            CNAVGCUMULINTER = outputFile.createVariable('cn-timeavg-cumul-%s-inter' % k, N.Float, ('NRVALUES',))
            CNAVGCUMULINTER[:] = N.cumsum(coordNumberInter[k].sum(1))/float(self.nFrames)
            CNAVGCUMULINTER.units = 'unitless'
                        
            CNAVGCUMULTOTAL = outputFile.createVariable('cn-timeavg-cumul-%s' % k, N.Float, ('NRVALUES',))
            CNAVGCUMULTOTAL[:] = N.cumsum(coordNumberTotal[k].sum(1))/float(self.nFrames)
            CNAVGCUMULTOTAL.units = 'unitless'
                        
        CNAVG = outputFile.createVariable('cn-timeavg', N.Float, ('NRVALUES',))
        CNAVG[:] = coordNumberTotal['total'].sum(1)/float(self.nFrames)
        CNAVG.units = 'unitless'

        CNAVGCUMUL = outputFile.createVariable('cn-timeavg-cumul', N.Float, ('NRVALUES',))
        CNAVGCUMUL[:] = N.cumsum(coordNumberTotal['total'].sum(1))/float(self.nFrames)
        CNAVGCUMUL.units = 'unitless'
        
        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'r', 'yVar' : 'cn-total-tavg-cumul'} 
        
        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
        
#####################################################################################
# DENSITY PROFILE ANALYSIS
#####################################################################################
class DensityProfile(Analysis):
    """Sets up a Density Profile analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: DensityProfile(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory    -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo      -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                               number to consider, 'last' is an integer specifying the last frame number to consider and 
                               'step' is an integer specifying the step number between two frames.
            * thickness     -- a float specifying the thickness of the bins
            * direction     -- one of X, Y or Z to specify along which direction the profile should be determined.
            * subset        -- a selection string specifying the atoms to consider for the analysis.
            * deuteration   -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * distanceunits -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.                             
            * pyroserver    -- a string specifying if Pyro will be used and how to run the analysis.
            * output        -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                               instead of the '.nc' extension.
    
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.

    Comments:        
        
        - This code contains a pyrex function for the distance histogram calculation than enhances significantly its 
          performance.
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'DP')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
        
        self.average_basis_vector_length = 0.0
        
        self.density_profile = {}

        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            self.density_profile[element] = N.zeros((self.n_bins), typecode = N.Float)
                    
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        if self.universe.basisVectors() is None:
            raise Error('The universe must be periodic for this kind of calculation.')
        
        # The 'thickness' input parameter. It must be a float or a value convertible to a float > 0.                
        try:
            self.thickness = float(self.parameters['thickness'])

        # The value can not be converted to a float.
        except ValueError:
            raise Error('Error when parsing "thickness" parameter: must be a float or a value convertible to a float.')
                
        else:
            # Must > 0.
            if self.thickness <= 0.0:
                raise Error('Error when parsing "thickness" parameter: must be > 0.')
            
        try:
            self.direction = ['x','y','z'].index(self.parameters['direction'].lower())
        except:
            raise Error('Error when parsing "direction" parameter: must be "X", "Y" or "Z"')
                                   
        self.buildTimeInfo()
        
        firstFrame = self.frameIndexes[0]
        self.trajectory.universe.setFromTrajectory(self.trajectory, firstFrame)
        
        self.n_bins = int(self.universe.basisVectors()[self.direction].length()/self.thickness) + 1
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()
                                                                                                    
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
        
    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame in |self.frameIndexes| array.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
        
        dp_per_frame = {}
                
        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])
                
        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
        box_coord = trajectory.universe._realToBoxPointArray(trajectory.universe.configuration().array) + 0.5
        basis_vector_length = trajectory.universe.basisVectors()[self.direction].length()
        
        for at in selectedAtoms:
            element = self.atomInformation[at.index]['element']
            if not dp_per_frame.has_key(element):
                dp_per_frame[element] = N.zeros((self.n_bins), typecode = N.Float)
                
            bin = int(self.n_bins*box_coord[at.index,self.direction])
                        
            if bin < 0: bin = 0
            if bin >= self.n_bins: bin = self.n_bins - 1            
                
            dp_per_frame[element][bin] += 1.0
                                        
        return frameIndex, (basis_vector_length, dp_per_frame)

    def combine(self, frameIndex, x):
        """
        """
        
        basis_vector_length, dp_per_frame = x
        
        self.average_basis_vector_length += basis_vector_length
        
        for element in dp_per_frame:
            N.add(self.density_profile[element], dp_per_frame[element], self.density_profile[element])

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """
        
        self.density_profile['total'] = N.zeros((self.n_bins), typecode = N.Float)
        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']
            self.density_profile[element] /= n
            self.density_profile['total'] += n*w*self.density_profile[element]
                                    
        for element in self.density_profile:
            self.density_profile[element] /= self.nFrames
            
        self.average_basis_vector_length /= self.nFrames
                                                                                  
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NBINS', self.n_bins)
        
        dr = self.average_basis_vector_length/self.n_bins
        
        r_values = dr*N.arange(self.n_bins) + dr/2.0

        # Creation of the NetCDF output variables.
        RVALUES = outputFile.createVariable('r', N.Float, ('NBINS',))
        RVALUES[:] = self.distanceConv*r_values
        RVALUES.units = self.distanceUnits

        for element in self.density_profile:

            DP = outputFile.createVariable('density-%s' % element, N.Float, ('NBINS',))
            DP[:] = self.density_profile[element]
            DP.units = 'unitless'
                                
        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'r', 'yVar' : 'density-total'} 
        
        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
        
#####################################################################################
# SCREW FIT ANALYSIS
#####################################################################################
class ScrewFitAnalysis(Analysis):
    """Set up a Screw Fit analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: ScrewFit(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory  -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo    -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                             number to consider, 'last' is an integer specifying the last frame number to consider and 
                             'step' is an integer specifying the step number between two frames.
            * sfa         -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                             instead of the '.nc' extension.
            * pyroserver  -- a string specifying if Pyro will be used and how to run the analysis.

    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
    Comments:
                                              
        - This code is based on a first implementation made by Paolo Calligari.
        
        - For more details: Kneller, G.R., Calligari, P. Acta Crystallographica , D62, 302-311
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'SFA')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                
        # The dictionnaries that will store respectively the orientational distance, the helix radius and the 
        # straightness for each chain found in the universe.
        self.orientDist = {}
        self.helixRadius = {}
        self.straightness = {}
        self.chainContents = {}

        self.nChemicalObjects = len(self.universe.objectList())        
        chemicalObjects = self.universe.objectList()

        for comp in range(self.nChemicalObjects):
                    
            obj = chemicalObjects[comp]
            objName = self.chemicalObjectNames[comp]
            
            if self.chemicalObjectInfo[objName]['objectclass'] in ['PeptideChain', 'Protein']:
    
                for chain in obj:
                    
                    self.chainContents[chain.name] = [r.sequence_number for r in chain.residues()]
                    self.straightness[chain.name] = N.zeros((self.nFrames,len(chain) - 4), typecode = N.Float)
                    self.helixRadius[chain.name] = N.zeros((self.nFrames,len(chain) - 3), typecode = N.Float)
                    self.orientDist[chain.name] = N.zeros((self.nFrames,len(chain) - 1), typecode = N.Float)
        
        # Case where no protein or peptide chain was found in the universe.
        if not self.orientDist:
            raise Error('The universe must contains at least one peptide chain to perform a screw fit analysis.')
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        # Some additional checkings.
        self.buildTimeInfo()

        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                    
    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame in |self.frameIndexes| array.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
        
        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
            
        od = {}
        hr = {}
        s = {}
        
        chemicalObjects = trajectory.universe.objectList()

        for comp in range(self.nChemicalObjects):
                    
            obj = chemicalObjects[comp]
            objName = self.chemicalObjectNames[comp]
            
            if self.chemicalObjectInfo[objName]['objectclass'] in ['PeptideChain', 'Protein']:
        
                for chain in obj:
                    
                    od[chain.name] = self.angularDistance(chain)
                    hr[chain.name], s[chain.name] = self.screwMotionAnalysis(chain)

        return frameIndex, (od, hr, s)
    
    def combine(self, frameIndex, x):
        """
        """
        
        comp = self.frameIndexes.index(frameIndex)
        for chainName in x[0].keys():
            self.orientDist[chainName][comp,:] = x[0][chainName].copy()
            self.helixRadius[chainName][comp,:] = x[1][chainName].copy()
            self.straightness[chainName][comp,:] = x[2][chainName].copy()
                
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        outputFile = NetCDFFile(self.output, 'w', self.information)
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        outputFile.createDimension('NFRAMES', self.nFrames)

        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times
        TIMES.units = 'ps'

        for chainName in self.straightness.keys():

            straightnessY = chainName + '_STRAIGHTNESS_Y'
            helixRadiusY  = chainName + '_HELIXRADIUS_Y'
            orientDistY   = chainName + '_ORIENTDIST_Y'

            outputFile.createDimension(straightnessY, self.straightness[chainName].shape[1])
            outputFile.createDimension(helixRadiusY, self.helixRadius[chainName].shape[1])
            outputFile.createDimension(orientDistY, self.orientDist[chainName].shape[1])

            STRAIGHTNESS_Y = outputFile.createVariable(straightnessY.lower(), N.Int32, (straightnessY,))
            HELIXRADIUS_Y = outputFile.createVariable(helixRadiusY.lower(), N.Int32, (helixRadiusY,))
            ORIENTDIST_Y = outputFile.createVariable(orientDistY.lower(), N.Int32, (orientDistY,))

            for r in range(self.straightness[chainName].shape[1]):
                res = self.chainContents[chainName][r]
                STRAIGHTNESS_Y[r] = HELIXRADIUS_Y[r] = ORIENTDIST_Y[r] = res

            HELIXRADIUS_Y[r+1] = ORIENTDIST_Y[r+1] = self.chainContents[chainName][r+1]
            
            ORIENTDIST_Y[r+2:r+3] = self.chainContents[chainName][r+2:r+3]
            
            STRAIGHTNESS = outputFile.createVariable(chainName + '_straightness', N.Float, ('NFRAMES',straightnessY))
            STRAIGHTNESS[:,:] = self.straightness[chainName][:,:]

            HELIXRADIUS  = outputFile.createVariable(chainName + '_helixradius' , N.Float, ('NFRAMES',helixRadiusY))
            HELIXRADIUS.units = 'nm'
            HELIXRADIUS[:,:]  = self.helixRadius[chainName][:,:]

            ORIENTDIST   = outputFile.createVariable(chainName + '_orientdist'  , N.Float, ('NFRAMES',orientDistY))
            ORIENTDIST[:,:]   = self.orientDist[chainName][:,:]

        asciiVar = sorted(outputFile.variables.keys())            

        outputFile.close()

        self.toPlot = None

        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)

    def findQuaternionMatrix(self, peptide, point_ref, conf1, conf2 = None, matrix = True):
        """ 
        Returns the complete matrix of quaternions compatibles with linear trasformation.|conf1| is
        the reference configuration. |point_ref| is the reference point about which the fit is calculated
        """

        if conf1.universe != peptide.universe():
            raise Error("conformation is for a different universe")

        if conf2 is None:
            conf1, conf2 = conf2, conf1

        else:
            if conf2.universe != peptide.universe():
                raise Error('conformation is for a different universe')

        ref = conf1
        conf = conf2
        weights = peptide.universe().masses()
        weights = weights/peptide.mass()
        ref_cms = point_ref.position().array
        pos = N.zeros((3,), typecode = N.Float)
        pos = point_ref.position(conf).array
        possq = 0.
        cross = N.zeros((3, 3), typecode = N.Float)

        for a in peptide.atomList():
            r = a.position(conf).array - pos
            r_ref = a.position(ref).array-ref_cms
            w = weights[a]
            possq = possq + w*N.add.reduce(r*r) + w*N.add.reduce(r_ref*r_ref)
            cross = cross + w*r[:, N.NewAxis]*r_ref[N.NewAxis, :]

        k = N.zeros((4, 4), typecode = N.Float)
        k[0, 0] = -cross[0, 0] - cross[1, 1] - cross[2, 2]
        k[0, 1] =  cross[1, 2] - cross[2, 1]
        k[0, 2] =  cross[2, 0] - cross[0, 2]
        k[0, 3] =  cross[0, 1] - cross[1, 0]
        k[1, 1] = -cross[0, 0] + cross[1, 1] + cross[2, 2]
        k[1, 2] = -cross[0, 1] - cross[1, 0]
        k[1, 3] = -cross[0, 2] - cross[2, 0]
        k[2, 2] =  cross[0, 0] - cross[1, 1] + cross[2, 2]
        k[2, 3] = -cross[1, 2] - cross[2, 1]
        k[3, 3] =  cross[0, 0] + cross[1, 1] - cross[2, 2]

        for i in range(1, 4):
            for j in range(i):
                k[i, j] = k[j, i]

        k = 2.*k
        for i in range(4):
            k[i, i] = k[i, i] + possq - N.add.reduce(pos*pos)

        e, v = linalg.eigh(k)
        emin = e.argmin()

        v = v[emin]
        if v[0] < 0: v = -v
        if e[emin] <= 0.:
            rms = 0.
        else:
            rms = N.sqrt(e[emin])

        if matrix:
            emax = e.argmax()
            QuatMatrix = v

            return Quaternion.Quaternion(QuatMatrix),v, e, e[emin],e[emax], rms

        else:
            return Quaternion.Quaternion(v), Vector(ref_cms), Vector(pos), rms      

    def findGenericTransformation(self,peptide, point_ref,conf1, conf2 = None):

        q, cm1, cm2, rms = self.findQuaternionMatrix(peptide,point_ref,conf1,conf2,False)

        return Transformation.Translation(cm2), q.asRotation(), Transformation.Translation(-cm1), rms

    def angularDistance(self, chain):

        universe = InfiniteUniverse()

        Crefpos = Vector(0.00,0.00,0.00)
        Orefpos = Vector(0.00,0.00,0.123)
        Nrefpos = Vector(0.111,0.00,-0.0728)

        Catom = Atom('C', position = Crefpos)
        Oatom = Atom('O', position = Orefpos)
        Natom = Atom('N', position = Nrefpos)

        plane = Collection()
        plane.addObject(Catom)
        plane.addObject(Oatom)
        plane.addObject(Natom)

        universe.addObject(plane)

        refconf = copy.copy(universe.configuration())

        chainLength = len(chain) - 1
        orientDist = N.zeros((chainLength,),typecode = N.Float)

        for resInd in range(chainLength):

            Cpos = chain[resInd].peptide.C.position()
            Opos = chain[resInd].peptide.O.position()
            Npos = chain[resInd+1].peptide.N.position()

            Catom.setPosition(Cpos)
            Oatom.setPosition(Opos)
            Natom.setPosition(Npos)

            distV = universe.distanceVector(Cpos,Crefpos)
            Catom.translateBy(distV)
            Oatom.translateBy(distV)
            Natom.translateBy(distV)

            Qm, v, e, emin, emax, rms = self.findQuaternionMatrix(plane, Catom, refconf, matrix = True)

            rmsD = plane.rmsDifference(refconf)

            orientDist[resInd] = rmsD/ N.sqrt(emax)

            refconf = copy.copy(universe.configuration())

        return orientDist

    def screwMotionAnalysis(self, chain):

        Crefpos = chain[0].peptide.C.position()
        Orefpos = chain[0].peptide.O.position()
        Nrefpos = chain[1].peptide.N.position()

        Catom = Atom('C', position=Crefpos)
        Oatom = Atom('O', position=Orefpos)
        Natom = Atom('N', position=Nrefpos)

        plane = Collection()
        plane.addObject(Catom)
        plane.addObject(Oatom)
        plane.addObject(Natom)

        universe = InfiniteUniverse()
        universe.addObject(plane)

        refconf = copy.copy(universe.configuration())

        chainLength = len(chain)

        DatiAsse=[]
        for i in range(chainLength - 2):

            Cpos = chain[i+1].peptide.C.position()
            Opos = chain[i+1].peptide.O.position()
            Npos = chain[i+2].peptide.N.position()

            Catom.setPosition(Cpos)
            Oatom.setPosition(Opos)
            Natom.setPosition(Npos)

            newconf = copy.copy(universe.configuration())

            TL1, ROT, TL2, rms = self.findGenericTransformation(plane, Catom, refconf, newconf)
            Tr = TL1*ROT*TL2
            DatiAsse.append(Tr.screwMotion())
            refconf = copy.copy(universe.configuration())

        asseIP = []
        for i in range(len(DatiAsse)):
            if DatiAsse[i][1].length() != 0.:
                asseIP.append(Line(DatiAsse[i][0],DatiAsse[i][1]))
            else:
                asseIP.append(Line(DatiAsse[i][0],DatiAsse[i+1][1]))

        helixRadius = N.zeros((chainLength-3,), typecode = N.Float)

        straightness = N.zeros((chainLength-4,), typecode = N.Float)

        for resInd in range(chainLength-3):

            Ca = chain[resInd].peptide.C.position()
            Ca1 = chain[resInd+1].peptide.C.position()
            Ca2 = chain[resInd+2].peptide.C.position()

            helixRadius[resInd]= asseIP[resInd].distanceFrom(Ca)                        

            if resInd <= len(DatiAsse)-3:
                projCa  = asseIP[resInd].projectionOf(Ca)
                projCa1 = asseIP[resInd+1].projectionOf(Ca1)
                projCa2 = asseIP[resInd+2].projectionOf(Ca2)

                pscal = (projCa - projCa1)*(projCa1 - projCa2)
                mod1 = (projCa - projCa1).length()
                mod2 = (projCa1 - projCa2).length()

                straightness[resInd] = pscal/(mod1*mod2)

        return helixRadius, straightness

#####################################################################################
# SPATIAL DENSITY ANALYSIS
#####################################################################################
class SpatialDensity(Analysis):
    """Sets up a Spatial Density analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: SpatialDensity(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo   -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                            number to consider, 'last' is an integer specifying the last frame number to consider and 
                            'step' is an integer specifying the step number between two frames.
            * rvalues    -- a string of the form 'rmin:rmax:dr' where 'rmin' is a float specifying the minimum distance to 
                            consider, 'rmax' is a float specifying the maximum distance value to consider and 'dr' is a float
                            specifying the distance increment. 
            * group      -- a selection string specifying the groups of atoms that will be used to define the points around which 
                            the coordination number will be computed. For each group, there is one point defined as the center of 
                            gravity of the group.
            * atomorder  -- a string of the form 'atom1,atom2,atom3' where 'atom1', 'atom2' and 'atom3' are 
                            respectively the MMTK atom names of the atoms in the way they should be ordered.
            * target     -- a selection string specifying the groups of atoms that will be used to define the points around which 
                            the coordination number will be computed. For each group, there is one point defined as the center of 
                            gravity of the group.
            * pyroserver -- a string specifying if Pyro will be used and how to run the analysis.
            * output     -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                            instead of the '.nc' extension.
    
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.

    Comments:        
        
        - This code contains a pyrex function for the distance histogram calculation than enhances significantly its 
          performance.
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'SD')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()

        if self.universe.basisVectors() is None:
            raise Error('The universe must be periodic for this kind of calculation.')
                                                                                        
        self.SD = N.zeros((self.nRBins,self.nThetaBins,self.nPhiBins), N.Float)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        self.buildTimeInfo()

        self.group = self.selectGroups(self.groupDefinition, self.atomOrder)
                        
        # The self.nGroups is needed when dealing with the template defined for analysis whose mainloop is over groups
        self.nGroups = len(self.group)

        self.target = self.selectGroups(self.targetDefinition)

        self.nTargets = len(self.target)

        self.thetaValues = parseInterval(self.parameters['thetavalues'])
                                        
        if len(self.thetaValues) <= 1:
            raise Error('Error when parsing thetavalues parameter: the theta interval %s could not be parsed properly.' % self.parameters['thetavalues'])

        self.thetaMin = self.thetaValues[0]
                
        self.dTheta = self.thetaValues[1] - self.thetaValues[0]
                
        self.nThetaBins = len(self.thetaValues) - 1
                    
        # The thetamin value must be >= 0.
        if self.thetaMin < 0.0:
            raise Error('Error when parsing "thetavalues" parameter: the thetamin value must be >= 0 deg.')

        # The thetamax value must be <= 180.
        if self.thetaValues[-1] > 180.0:
            raise Error('Error when parsing "thetavalues" parameter: the thetamax value must be <= 180 deg.')

        # The |self.thetaMin| and |self.dTheta| attributes are converted in radians.
        self.thetaMin *= Units.deg
        self.dTheta *= Units.deg
                        
        self.phiValues = parseInterval(self.parameters['phivalues'])
                                        
        if len(self.phiValues) <= 1:
            raise Error('Error when parsing phivalues parameter: the phi interval %s could not be parsed properly.' % self.parameters['phivalues'])

        self.phiMin = self.phiValues[0]
                
        self.dPhi = self.phiValues[1] - self.phiValues[0]
                
        self.nPhiBins = len(self.phiValues) - 1
            
        # The min phi value must be >= 180.
        if self.phiMin < -180.0:
            raise Error('Error when parsing "phivalues" parameter: the min value for phi must be >= -180 deg.')
        
        # The max phi value must be <= 180.
        if self.phiValues[-1] > 180.0:
            raise Error('Error when parsing "phivalues" parameter: the max value for phi must be <= 180 deg.')

        # The |self.phiMin| and |self.dPhi| attributes are converted in radians.
        self.phiMin *= Units.deg
        self.dPhi *= Units.deg
                
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                            
    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame in |self.frameIndexes| array.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
        
        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))

        groups = [Collection([orderedAtoms[ind] for ind in aIndexes]) for aIndexes in self.group]
                            
        conf = trajectory.configuration[frameIndex]
        trajectory.universe.setConfiguration(trajectory.universe.contiguousObjectConfiguration(Collection(groups), conf))
            
        targets = [Collection([orderedAtoms[ind] for ind in aIndexes]) for aIndexes in self.target]
        targetCOMsList = [t.centerOfMass() for t in targets]
            
        sphericalBinList = []

        for g in groups:
            
            gCOM = g.centerOfMass()

            a1, a2, a3 = g
            v1 = trajectory.universe.distanceVector(a1.position(), a2.position()).normal()
            v2 = trajectory.universe.distanceVector(a1.position(), a3.position()).normal()
            vi, vj, vk = self.constructOrthonormalBasis(a1.position(), v1, v2)

            for tCOM in targetCOMsList:
                
                if tCOM == gCOM:
                    continue

                r = (gCOM - tCOM).length()                
                rBin = int((r - self.rMin)/self.dR)
                
                # The group is not considered if it falls outside the desired r interval.
                if (rBin < 0) or (rBin >= self.nRBins):
                    continue                
                
                x, y, z = changeBasis(gCOM,tCOM,vi,vj,vk)

                # The spherical theta is computed
                theta = N.arccos(z/r)
                thetaBin = int((theta - self.thetaMin)/self.dTheta)
                
                # The group is not considered if it falls outside the desired theta interval.
                if (thetaBin < 0) or (thetaBin >= self.nThetaBins):
                    continue                

                # The spherical phi is computed
                phi = N.arctan2(y,x)                
                phiBin = int((phi - self.phiMin)/self.dPhi)
                
                # The group is not considered if it falls outside the desired phi interval.
                if (phiBin < 0) or (phiBin >= self.nPhiBins):
                    continue

                sphericalBinList.append((rBin, thetaBin, phiBin))
                                
        return frameIndex, sphericalBinList
    
    def combine(self, frameIndex, x):
        """
        """
        for rBin, thetaBin, phiBin in x:
            # The Spatial density is incremented by one unit.                
            self.SD[rBin, thetaBin, phiBin] += 1.0
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        self.SD /= float(self.nFrames*self.nGroups*self.nTargets)
        
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NTHETAVALUES', self.nThetaBins)
        outputFile.createDimension('NPHIVALUES', self.nPhiBins)

        self.thetaValues = self.thetaValues[:-1] + self.dTheta/(2.0*Units.deg)
        
        # Creation of the NetCDF output variables.
        THETAS = outputFile.createVariable('theta', N.Float, ('NTHETAVALUES',))
        THETAS[:] = self.thetaValues[:]
        THETAS.units = 'deg'
        
        self.phiValues = self.phiValues[:-1] + self.dPhi/(2.0*Units.deg)
        
        PHIS = outputFile.createVariable('phi', N.Float, ('NPHIVALUES',))
        PHIS[:] = self.phiValues[:]
        PHIS.units = 'deg'

        self.rValues = self.rValues[:-1] + self.dR/2.0

        for comp in range(self.nRBins):
            
            r = self.rValues[comp]
            
            SD = outputFile.createVariable('sd_r%snm' % round(r,3), N.Float, ('NTHETAVALUES','NPHIVALUES'))
            SD[:,:] = self.SD[comp,:,:]
            SD.units = 'unitless'
        
        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = None
        
        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)

    def constructOrthonormalBasis(self, a1, v1, v2):
        """This method construct a set of three oriented orthonormal axes i, j, k from a triplet of atoms
        such as (i,j,k) forms a clockwise orthonormal basis.
        If a1, a2 and a3 stand respectively for the three atoms of the triplet then:
            vector1 = (vector(a1,a2)_normalized + vector(a1,a3)_normalized)_normalized
            vector3 = (vector1 ^ vector(a1,a3))_normalized and correclty oriented
            vector2 = (vector3 ^ vector1)_normalized
        
        @param triplet: the triplet of atoms.
        @type triplet: a list of three MMTK Atoms

        @return: the three axis.
        @rtype: a list of three Scientific Vector   
        
        """

        n1 = (v1 + v2).normal()
        n3 = v1.cross(n1).normal()
        n2 = n3.cross(n1).normal()
                
        return n1, n2, n3

#####################################################################################
# HYDROGEN BOND SURVIVAL ANALYSIS ANALYSIS
#####################################################################################
class HydrogenBondSurvivalAnalysis(Analysis):
    """Sets up a Hydrogen Bond Survival Analysis analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: HydrogenBondSurvivalAnalaysis(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory    -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo      -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                               number to consider, 'last' is an integer specifying the last frame number to consider and 
                               'step' is an integer specifying the step number between two frames.
            * disminmax     -- a string of the form 'min:max' where 'min' and 'max' are floats specifying respectively 
                               the lower and the upper distance bounds to detect HB.
            * angminmax     -- a string of the form 'min:max' where 'min' and 'max' are floats specifying respectively 
                               the lower and the upper angle bounds to detect HB.
            * subset1       -- a selection string specifying the first set of atom to find HB donors and acceptors.
            * subset2       -- a selection string specifying the second set of atom to find HB donors and acceptors.
            * timeunits     -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * distanceunits -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.                             
            * pyroserver    -- a string specifying if Pyro will be used and how to run the analysis.
            * output        -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                               instead of the '.nc' extension.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.

    Comments:        
        
        - This code contains a pyrex function for the distance histogram calculation that is based on a FORTRAN code 
          written by Miguel Gonzalez, Insitut Laue Langevin, Grenoble, France.
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

        path = os.path.join(GVAR['nmoldyn_analysis'], 'HBSA')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                    
        self.HBPerFrame = {}        
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        # Some additional checkings.
        if self.universe.basisVectors() is None:
            raise Error('The universe must be periodic for this kind of calculation.')
                
        self.buildTimeInfo()
        
        orderedAtoms = sorted(self.trajectory.universe.atomList(), key = operator.attrgetter('index'))
        
        self.subset1 = self.selectAtoms(self.subset1Definition)
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset1])
        self.hb_indexes1, self.hb_info1 = detect_hb_atoms_subset(self.trajectory.universe, subset = selectedAtoms)

        self.subset2 = self.selectAtoms(self.subset2Definition)
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset2])
        self.hb_indexes2, self.hb_info2 = detect_hb_atoms_subset(self.trajectory.universe, subset = selectedAtoms)
        
        # The min and max distances are converted in the selected units.
        self.disMin *= self.distanceConv
        self.disMax *= self.distanceConv

        # The angles are converted in radians.
        self.angMin *= Units.deg
        self.angMax *= Units.deg
        
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                            
    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame in |self.frameIndexes| array.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
                        
        hits = []
                
        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
                        
        directCell = N.ravel(N.array([v for v in trajectory.universe.basisVectors()], typecode = N.Float))
        reverseCell = N.ravel(N.transpose(N.array([v for v in trajectory.universe.reciprocalBasisVectors()], typecode = N.Float)))
        
        cellVolume = trajectory.universe.cellVolume()
            
        hbond_matrix = N.zeros((len(self.hb_indexes1['acc']),len(self.hb_indexes2['don'])), typecode = N.Float32)

        hbond_detection(trajectory.universe.contiguousObjectConfiguration().array, directCell, reverseCell,\
                        self.hb_indexes1['acc'], self.hb_indexes2['don'], self.hb_indexes2['hyd'],\
                        hbond_matrix, self.disMin, self.disMax, self.angMin, self.angMax)
    
        hits.extend([tuple([self.hb_indexes1['acc'][hit[0]],self.hb_indexes2['don'][hit[1]]]) for hit in N.transpose(hbond_matrix.nonzero())])
        
        hbond_matrix = N.zeros((len(self.hb_indexes2['acc']),len(self.hb_indexes1['don'])), typecode = N.Float32)

        hbond_detection(trajectory.universe.contiguousObjectConfiguration().array, directCell, reverseCell,\
                        self.hb_indexes2['acc'], self.hb_indexes1['don'], self.hb_indexes1['hyd'],\
                        hbond_matrix, self.disMin, self.disMax, self.angMin, self.angMax)

        hits.extend([tuple([self.hb_indexes2['acc'][hit[0]],self.hb_indexes1['don'][hit[1]]]) for hit in N.transpose(hbond_matrix.nonzero())])
                                                                                                        
        return frameIndex, hits
    
    def combine(self, frameIndex, x):
        """
        """
        
        self.HBPerFrame[frameIndex] = x
                                            
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        HBondsLife = {}
        
        for comp in range(self.nFrames):
            
            frameIndex = self.frameIndexes[comp]
            
            for hit in self.HBPerFrame[frameIndex]:
                
                if HBondsLife.has_key(hit):
                    
                    if comp - HBondsLife[hit][-1][-1] > self.trueBreakStep:
                        HBondsLife[hit].append([comp, comp])
                        
                    else:
                        HBondsLife[hit][-1][1] = comp
                        
                else:
                    HBondsLife[hit] = [[comp, comp]]
                    
        HBondsSurv = {'symb' : {}, 'name' : {}, 'all' : N.zeros((self.nFrames,), typecode = N.Float)}
        nHBonds = {'symb' : {}, 'name' : {}, 'all' : 0}
                    
        for k in HBondsLife.keys():
            
            symb1 = self.atomInformation[k[0]]['element']
            symb2 = self.atomInformation[k[1]]['element']
            symb_pair = tuple(sorted([symb1, symb2]))
            if not HBondsSurv['symb'].has_key(symb_pair):
                HBondsSurv['symb'][symb_pair] = N.zeros((self.nFrames,), typecode = N.Float)
                nHBonds['symb'][symb_pair] = 0
                            
            name1 = self.atomInformation[k[0]]['name']
            name2 = self.atomInformation[k[1]]['name']
            name_pair = tuple(sorted([name1, name2]))
            if not HBondsSurv['name'].has_key(name_pair):
                HBondsSurv['name'][name_pair] = N.zeros((self.nFrames,), typecode = N.Float)
                nHBonds['name'][name_pair] = 0
                
            nBonds = len(HBondsLife[k])
            nHBonds['symb'][symb_pair] += nBonds
            nHBonds['name'][name_pair] += nBonds
            nHBonds['all'] += nBonds

            for hb in HBondsLife[k]:
                duration = hb[1] - hb[0] + 1 
                HBondsSurv['symb'][symb_pair][:duration] += 1.0
                HBondsSurv['name'][name_pair][:duration] += 1.0                                                        
                HBondsSurv['all'][:duration] += 1.0
        
        # The normalisation step.
        for k in nHBonds.keys():

            if k == 'all':
                HBondsSurv[k] /= float(nHBonds[k])
            
            else:                              
                for kk in nHBonds[k].keys():
                    HBondsSurv[k][kk] /= float(nHBonds[k][kk])        
                                    
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NTIMES', self.nFrames)

        TIMES = outputFile.createVariable('time', N.Float, ('NTIMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        for k in HBondsSurv.keys():
                                
            if k == 'all':
                HBSURV = outputFile.createVariable('hbsurv-all', N.Float, ('NTIMES',))
                HBSURV[:] = HBondsSurv[k]
                
            else:
                for kk in nHBonds[k].keys():
                    HBSURV = outputFile.createVariable('hbsurv_%s_%s' % (k, ''.join(kk)), N.Float, ('NTIMES',))
                    HBSURV[:] = HBondsSurv[k][kk]
                    
            HBSURV.units = 'unitless'

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'hbsurv-all'} 
        
        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
