"""Collections of classes for the determination of scattering-related properties.

Classes:
    * DynamicCoherentStructureFactor           : sets up a Dynamic Coherent Structure Factor analysis.
    * DynamicCoherentStructureFactorARModel    : sets up a Dynamic Coherent Structure Factor analysis using an Auto Regressive model.
    * DynamicIncoherentStructureFactor         : sets up an Dynamic Incoherent Structure Factor analysis.
    * DynamicIncoherentStructureFactorGaussian : sets up an Dynamic Incoherent Structure Factor analysis using a Gaussian approximation.
    * IncoherentStructureFactorARModel         : sets up an Dynamic Incoherent Structure Factor analysis using an Auto Regressive model.
    * ElasticIncoherentStructureFactor         : sets up an Elastic Incoherent Structure Factor analysis.
    * StaticCoherentStructureFactor            : sets up a Static Coherent Structure Factor analysis.
    * SmoothedStaticCoherentStructureFactor    : sets up a Smoothed Static Coherent Structure Factor analysis.
    
Procedures:
    * DynamicStructureFactor : returns the Dynamic Structure Factor.
"""

# The python distribution modules
import os
import re
from time import asctime
from timeit import default_timer

import numpy

# The ScientificPython modules
from Scientific import N 
from Scientific.IO.NetCDF import _NetCDFFile, NetCDFFile
from Scientific.Signals.Models import AutoRegressiveModel, AveragedAutoRegressiveModel

# The MMTK distribution modules
from MMTK.Collections import Collection
from MMTK_forcefield import NonbondedList
from MMTK.Trajectory import Trajectory

# The nMOLDYN modules
from smoothed_static_coherent_structure_factor import smoothed_static_coherent_structure_factor
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Analysis.Analysis import Analysis
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import convertNetCDFToASCII
from nMOLDYN.Mathematics.Analysis import correlation, FFT, gaussian, gaussianWindow
from nMOLDYN.Mathematics.ReciprocalSpace import QVectors

#####################################################################################
# DYNAMIC STRUCTURE FACTOR ANALYSIS
#####################################################################################
def DynamicStructureFactor(time, f_qt):
    """Computes the dynamic structure factor from an intermediate scattering function.
    
    @param netcdf: the intermediate scattering function from which the dynamic structure factor will be computed..
    @type netcdf: string or instance of _NetCDFFile
        
    @param alpha: the width, in percentage of the trajectory length, of the gaussian used in the smoothing procedure.
    @type alpha: float
    """
    
#    if isinstance(netcdf, str):
#        try:
#            outputFile = NetCDFFile(netcdf, 'a')
#        except IOError:
#            raise Error('The file %s could not be loaded.' % netcdf)
#    else:
#        return

#    output = {}
    
    # The tim NetCDF variable is extracted of the NC file.
#    time = outputFile.variables['time'].getValue()
#    dt = time[1] - time[0]
    
#    nQValues = len(outputFile.variables['q'][:])
    
    # frequencies = 1D Numeric array. Contains the frequencies at which the dynamic structure is calculated
#    frequencies = N.arange(len(time))/(2.0*len(time)*dt)

#    timeResolution = gaussian(time, sigma)
#    frequencyFunction = gaussian(frequencies, 1.0/(2.0*N.pi*sigma))

    # Creation of the NetCDF dimension FREQUENCY. The number of frequency point present in the NC output file.  
#    outputFile.createDimension('NFREQUENCIES', len(frequencies[:]))

    # Creation of the NetCDF variable |frequency| that will store the frequency at which the
    # dynamic structure factor is calculated
#    FREQUENCIES  = outputFile.createVariable('frequency', N.Float, ('NFREQUENCIES',))
#    FREQUENCIES[:] = frequencies[:]
#    FREQUENCIES.units = 'THz'

#    OMEGAS = outputFile.createVariable('angular_frequency', N.Float, ('NFREQUENCIES',))
#    OMEGAS[:] = 2.0*N.pi*frequencies[:]
#    OMEGAS.units = 'rad ps-1'

    # The resolution function.
    RESOLUTIONFUNCTION = outputFile.createVariable('time_resolution', N.Float, ('NFREQUENCIES',))
    RESOLUTIONFUNCTION[:] = timeResolution[:]
    RESOLUTIONFUNCTION.units = 'unitless'

    RESOLUTIONFUNCTION = outputFile.createVariable('frequency_resolution', N.Float, ('NFREQUENCIES',))
    RESOLUTIONFUNCTION[:] = frequencyFunction[:]
    RESOLUTIONFUNCTION.units = 'unitless'

    # Is there some psf terms in the intermediate scattering function input file ?
    partialTerms = [re.findall('Fqt-(.*)',v)[0] for v in outputFile.variables if 'Fqt-' in v]
    for pTerm in partialTerms:
        
        # Creation of the NetCDF variable |dsf| that will store the dynamic structure factor    
        DCSF = outputFile.createVariable('Sqw-'+pTerm, N.Float, ('NQVALUES','NFREQUENCIES'))
        
        # A loop is done for each q where the scattering function was calculated.
        for qVal in range(nQValues):
            DCSF[qVal] = 0.5 * dt * FFT(gaussianWindow(outputFile.variables['Fqt-'+pTerm][qVal], timeResolution)).real[:len(frequencies)]
            
        DCSF.units = 'unitless'
            
    asciiVar = sorted(outputFile.variables.keys())

    outputFile.close()

    # Creates an ASCII version of the NetCDF output file.
    convertNetCDFToASCII(inputFile = netcdf,\
                         outputFile = os.path.splitext(netcdf)[0] + '.cdl',\
                         variables = asciiVar)

#####################################################################################
# COHERENT STRUCTURE FACTOR ANALYSIS
#####################################################################################
class DynamicCoherentStructureFactor(Analysis):
    """Sets up a Dynamic Coherent Structure Factor analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: DynamicCoherentStructureFactor(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory        -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo          -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                   number to consider, 'last' is an integer specifying the last frame number to consider and 
                                   'step' is an integer specifying the step number between two frames.
            * qshellvalues      -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                                   'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                                   the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * qshellwidth       -- a float specifying the width of the q shells.
            * qvectorspershell  -- a float specifying the number of q vectors to generate per q shell.
            * qvectorsgenerator -- a string being one of 'isotropic', 'anisotropic' or 'explicit' specifying the way the q vectors
                                   will be generated.
            * qvectorsdirection -- a string of the form 'v1x,v1y,v1z;v2x,v2y,v2z...' where 'v1x', 'v2x' ..., 'v1y', 'v2y' ... and
                                   'v1z', 'v2z' ... are floats that represents respectively the x, y and z values of the vectord along 
                                   which the q vectors should be generated.
            * resolution        -- a float specifying the width of the gaussian, that will be used to mimics the experimental resolution.
            * subset            -- a selection string specifying the atoms to consider for the analysis.
            * deuteration       -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights           -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                   scheme to use.
            * timeunits         -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * frequencyunits    -- a string equal to 'THz', 'rad s^-1', "cm^-1", "meV" or 'ueV' that specifies the frequency units.
            * qunits            -- a string equal to 'nm^-1' or 'ang^-1' that specifies the q units.
            * pyroserver        -- a string specifying if Pyro will be used and how to run the analysis.
            * output            -- the output NetCDF file name for the intermediate scattering function.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'DCSF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
                
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                                                                        
        self.FQT = {'total' : N.zeros((self.nQValues,self.nFrames), typecode = N.Float)}

        # Loop over all the atom symbol.      
        for element1 in self.elementInformation.keys():
            # Loop over all the atom symbol.      
            for element2 in self.elementInformation.keys():
                # The symbol pair tuple that will be used as the key for the histogram dictionnary.
                pairName = tuple(sorted((element1,element2)))
                # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
                self.FQT[pairName] = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()
                        
        qv = QVectors(self.universe, self.qVectorsDict)
                    
        self.qRadii = qv.qRadii
        self.qVectors = qv.qVectors
        self.qVectorsStatistics = qv.statistics
        
        self.nQValues = len(self.qRadii)

        # Loop over the user-defined q values
        for comp in range(self.nQValues):
                            
            # For each q vector length, the corresponding set of q vector is transform into an array
            qarray = N.array(self.qVectors[comp], typecode = N.Float)
            
            # TAKE CARE, THIS CONVERTS THE Q VECTORS IN BOX COORDINATES.
            qarray = self.universe._boxToRealPointArray(qarray)
            
            # qVect[1][j] is now of shape (3,Qcount)
            self.qVectors[comp] = N.transpose(qarray)
                    
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
            
    def calc(self, qIndex, trajectory):
        """Calculates the contribution for one Q-shell.
        
        @param qIndex: the index of the Q-shell in |self.qRadii| list.
        @type qIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        rho = {}
        for element in self.elementInformation.keys():
            rho[element] = N.zeros((self.nFrames, self.qVectors[qIndex].shape[1]), N.Complex)

        # loop over the trajectory time steps
        for comp in range(self.nFrames):
            
            frameIndex = self.frameIndexes[comp]
            
            # conf contains the positions of all the atoms of the system for time i.
            conf = trajectory.configuration[frameIndex]
            conf.convertToBoxCoordinates()

            for element in self.elementInformation.keys():
                selAtoms = N.compress(self.elementInformation[element]['mask'], conf.array, 0)
                rho[element][comp,:] = N.add.reduce(N.exp(1j*N.dot(selAtoms, self.qVectors[qIndex])))

        return qIndex, rho
    
    def combine(self, qIndex, x):
        """
        """
        
        for element1 in self.elementInformation.keys():
            for element2 in self.elementInformation.keys():
                pairName = tuple(sorted((element1,element2)))
                corr = correlation(x[element1],x[element2])[:]
                self.FQT[pairName][qIndex,:] += corr/self.qVectors[qIndex].shape[1]

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                             
        
        for pairName in self.FQT.keys():
            
            if pairName == 'total':
                continue

            nA = self.elementInformation[pairName[0]]['number']
            nB = self.elementInformation[pairName[1]]['number']
            nAB = N.sqrt(nA*nB)

            wA = self.elementInformation[pairName[0]]['weight']
            wB = self.elementInformation[pairName[1]]['weight']
            wAB = wA*wB
                        
            self.FQT[pairName] = self.FQT[pairName]/nAB
                        
            N.add(self.FQT['total'], nAB*wAB*self.FQT[pairName], self.FQT['total'])

        # The frequencies are computed.
        freqs = N.arange(len(self.times))/(2.0*len(self.times)*self.dt)
        # The time resolution function is computed.
        timeResolution = gaussian(self.times, self.timeSigma)
        # The frequency resolution function is computed.
        frequencyFunction = gaussian(freqs, 1.0/(2.0*N.pi*self.timeSigma))
            
        # The NetCDF output file intermediate scattering function.
        outputFile       = NetCDFFile(self.output, 'w')
        # The title for the NetCDF file.
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()
        
        # Some NetCDF dimensions.
        outputFile.createDimension('NQVALUES', self.nQValues)
        outputFile.createDimension('NTIMES', self.nFrames)
        outputFile.createDimension('NFREQUENCIES', self.nFrames)
        outputFile.createDimension('NOCTANS', 8)
        outputFile.createDimension('OCTANNAME', 8)
                
        QLENGTH = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        QLENGTH[:] = self.qConv*self.qRadii
        QLENGTH.units = self.qUnits

        # Creation of the NetCDF variable |sf| that will store the structure factor    
        TIMES    = outputFile.createVariable('time', N.Float, ('NTIMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        # The frequencies are computed and stored inthe NetCDF file.
        FREQUENCIES  = outputFile.createVariable('frequency', N.Float, ('NFREQUENCIES',))
        FREQUENCIES[:] = self.frequencyConv*freqs[:]
        FREQUENCIES.units = self.frequencyUnits

        OMEGAS = outputFile.createVariable('angular_frequency', N.Float, ('NFREQUENCIES',))
        OMEGAS[:] = 2.0*N.pi*freqs[:]
        OMEGAS.units = 'rad ps-1'

        # Write the q vector generation statistics.
        OCTANS = outputFile.createVariable('octan', N.Character, ('NOCTANS','OCTANNAME'))
        OCTANS[:,:] = N.array([list('X-.Y-.Z-'),list('X-.Y-.Z+'),list('X-.Y+.Z-'),list('X-.Y+.Z+'),\
                                list('X+.Y-.Z-'),list('X+.Y-.Z+'),list('X+.Y+.Z-'),list('X+.Y+.Z+')])
        OCTANS.units = 'unitless'
                
        STAT = outputFile.createVariable('qvectors_statistics', N.Int32, ('NQVALUES', 'NOCTANS'))
        STAT[:,:] = self.qVectorsStatistics

        # The time resolution function is written.
        TIME_RES = outputFile.createVariable('time_resolution', N.Float, ('NTIMES',))
        TIME_RES[:] = timeResolution[:]
        TIME_RES.units = 'unitless'

        # The time resolution function is written.
        FREQ_RES = outputFile.createVariable('frequency_resolution', N.Float, ('NFREQUENCIES',))
        FREQ_RES[:] = frequencyFunction[:]
        FREQ_RES.units = 'unitless'
                
        # The Fqt total is first written.
        DCSF = outputFile.createVariable('Fqt-total', N.Float, ('NQVALUES','NTIMES'))
        DCSF[:] = self.FQT['total']
        DCSF.units = 'unitless'
        
        # Then the Sqw total is computed and written.     
        DCSF = outputFile.createVariable('Sqw-total', N.Float, ('NQVALUES','NFREQUENCIES'))
                                
        # A loop is done for each q where the scattering function was calculated.
        for qVal in range(self.nQValues):
            DCSF[qVal] = 0.5 * self.dt * FFT(gaussianWindow(self.FQT['total'][qVal], timeResolution)).real[:self.nFrames]            
        DCSF.units = 'unitless'
                                        
        for k in sorted(self.FQT.keys()):
            
            if not isinstance(k,tuple): continue
            
            partial = ''.join(k)
            
            # The partial Fqts are then written.
            DCSF = outputFile.createVariable('Fqt-%s' % partial, N.Float, ('NQVALUES','NTIMES'))                
            DCSF[:] = self.FQT[k]
            DCSF.units = 'unitless'
            
            # The partial Sqws are computed and written.     
            DCSF = outputFile.createVariable('Sqw-%s' % partial, N.Float, ('NQVALUES','NFREQUENCIES'))
                
            # A loop is done for each q where the scattering function was calculated.
            for qVal in range(self.nQValues):
                DCSF[qVal] = 0.5 * self.dt * FFT(gaussianWindow(self.FQT[k][qVal], timeResolution)).real[:self.nFrames]            
            DCSF.units = 'unitless'                                    
                                
        outputFile.close()
        
        self.toPlot = {'netcdf' : self.output, 'xVar' : 'q', 'yVar' : 'time', 'zVar' : 'Fqt-total'}
            
#####################################################################################
# STATIC COHERENT STRUCTURE FACTOR ANALYSIS
#####################################################################################
class StaticCoherentStructureFactor(Analysis):
    """Sets up a Coherent Structure Factor analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: StaticCoherentStructureFactor(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory        -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo          -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                   number to consider, 'last' is an integer specifying the last frame number to consider and 
                                   'step' is an integer specifying the step number between two frames.
            * qshellvalues      -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                                   'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                                   the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * qshellwidth       -- a float specifying the width of the q shells.
            * qvectorspershell  -- a float specifying the number of q vectors to generate per q shell.
            * qvectorsgenerator -- a string being one of 'isotropic', 'anisotropic' or 'explicit' specifying the way the q vectors
                                   will be generated.
            * qvectorsdirection -- a string of the form 'v1x,v1y,v1z;v2x,v2y,v2z...' where 'v1x', 'v2x' ..., 'v1y', 'v2y' ... and
                                   'v1z', 'v2z' ... are floats that represents respectively the x, y and z values of the vectord along 
                                   which the q vectors should be generated.
            * resolution      -- a float specifying the width of the gaussian, that will be used to mimics the experimental resolution.
            * subset            -- a selection string specifying the atoms to consider for the analysis.
            * deuteration       -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * qunits            -- a string equal to 'nm^-1' or 'ang^-1' that specifies the q units.
            * weights           -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                   scheme to use.
            * pyroserver        -- a string specifying if Pyro will be used and how to run the analysis.
            * output            -- the output NetCDF file name for the intermediate scattering function.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'SCSF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                                                        
        self.SQ = {'total' : N.zeros((self.nQValues,), typecode = N.Float)}

        # Loop over all the atom symbol.      
        for element1 in self.elementInformation.keys():
            # Loop over all the atom symbol.      
            for element2 in self.elementInformation.keys():
                # The symbol pair tuple that will be used as the key for the histogram dictionnary.
                pairName = tuple(sorted((element1,element2)))
                # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
                self.SQ[pairName] = N.zeros((self.nQValues,), typecode = N.Float)

        # Loop over the user-defined q values
        for qIndex in range(self.nQValues):
                            
            # For each q vector length, the corresponding set of q vector is transform into an array
            qarray = N.array(self.qVectors[qIndex], typecode = N.Float)
            qarray = self.universe._boxToRealPointArray(qarray)
            # qVect[1][j] is now of shape (3,Qcount)
            self.qVectors[qIndex] = N.transpose(qarray)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)

        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()

        qv = QVectors(self.universe, self.qVectorsDict)
                    
        self.qRadii = qv.qRadii
        self.qVectors = qv.qVectors
        self.qVectorsStatistics = qv.statistics

        self.nQValues = len(self.qRadii)

        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
            
    def calc(self, qIndex, trajectory):
        """Calculates the contribution for one Q-shell.
        
        @param qIndex: the index of the Q-shell in |self.qRadii| list.
        @type qIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
                                                    
        rho = {}
        for element in self.elementInformation.keys():
            rho[element] = N.zeros((self.nFrames, self.qVectors[qIndex].shape[1]), typecode = N.Complex)
        
        comp = 0
        # loop over the trajectory time steps
        for frameIndex in self.frameIndexes:
            
            # conf contains the positions of all the atoms of the system for time i.
            conf = trajectory.configuration[frameIndex]
            conf.convertToBoxCoordinates()

            for element in self.elementInformation.keys():
                selAtoms = N.compress(self.elementInformation[element]['mask'], conf.array, 0)
                rho[element][comp,:] = N.add.reduce(N.exp(1j*N.dot(selAtoms, self.qVectors[qIndex])))
                                        
            comp += 1
                
        return qIndex, rho
    
    def combine(self, qIndex, x):
        """
        """

        for element1 in self.elementInformation.keys():
            
            for element2 in self.elementInformation.keys():
                pairName = tuple(sorted((element1,element2)))
                
                rhoirhoj = N.conjugate(x[element1]) * x[element2]
                corrt0 = N.add.reduce(N.add.reduce(rhoirhoj,1)).real/self.nFrames
                                
                self.SQ[pairName][qIndex] += corrt0/self.qVectors[qIndex].shape[1]

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                

        for pairName in self.SQ.keys():
            if pairName == 'total':
                continue
            
            nA = self.elementInformation[pairName[0]]['number']
            nB = self.elementInformation[pairName[1]]['number']
            nAB = N.sqrt(nA*nB)

            wA = self.elementInformation[pairName[0]]['weight']
            wB = self.elementInformation[pairName[1]]['weight']
            wAB = wA*wB
                                        
            self.SQ[pairName] = self.SQ[pairName]/nAB
            
            N.add(self.SQ['total'], nAB*wAB*self.SQ[pairName], self.SQ['total'])
            
        # The NetCDF output file intermediate scattering function.
        outputFile       = NetCDFFile(self.output, 'w')
        # The title for the NetCDF file.
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some NetCDF dimensions.
        outputFile.createDimension('NQVALUES', None)
        outputFile.createDimension('NOCTANS', 8)
        outputFile.createDimension('OCTANNAME', 8)
        
        QLENGTH = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        QLENGTH[:] = self.qConv*self.qRadii
        QLENGTH.units = self.qUnits

        OCTANS = outputFile.createVariable('octan', N.Character, ('NOCTANS','OCTANNAME'))
        OCTANS[:,:] = N.array([list('X-.Y-.Z-'),list('X-.Y-.Z+'),list('X-.Y+.Z-'),list('X-.Y+.Z+'),\
                                list('X+.Y-.Z-'),list('X+.Y-.Z+'),list('X+.Y+.Z-'),list('X+.Y+.Z+')])
        OCTANS.units = 'unitless'
                
        STAT = outputFile.createVariable('qvectors_statistics', N.Int32, ('NQVALUES', 'NOCTANS'))
        STAT[:,:] = self.qVectorsStatistics

        SCSF = outputFile.createVariable('Sq-total', N.Float, ('NQVALUES',))
        SCSF[:] = self.SQ['total']
        SCSF.units = 'unitless'

        for k in sorted(self.SQ.keys()):
            
            if not isinstance(k,tuple):
                continue
            
            SCSF = outputFile.createVariable('Sq-%s' % ''.join(k), N.Float, ('NQVALUES',))            
            SCSF[:] = self.SQ[k]
            SCSF.units = 'unitless'
            
        outputFile.close()
        
        self.toPlot = {'netcdf' : self.output, 'xVar' : 'q', 'yVar' : 'Sq-total'}
                
#####################################################################################
# COHERENT STRUCTURE FACTOR ANALYSIS WITHIN AR MODEL FRAMEWORK
#####################################################################################
class AutoRegressiveDynamicCoherentStructureFactor(Analysis):
    """Sets up a Dynamic Coherent Structure Factor analysis using an Auto Regressive model.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: DynamicCoherentStructureFactorARModel(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory        -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo          -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                   number to consider, 'last' is an integer specifying the last frame number to consider and 
                                   'step' is an integer specifying the step number between two frames.
            * armodelorder      -- an integer in [1, len(trajectory)[ specifying the order of the model
            * qshellvalues      -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                                   'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                                   the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * qshellwidth       -- a float specifying the width of the q shells.
            * qvectorspershell  -- a float specifying the number of q vectors to generate per q shell.
            * qvectorsgenerator -- a string being one of 'isotropic', 'anisotropic' or 'explicit' specifying the way the q vectors
                                   will be generated.
            * qvectorsdirection -- a string of the form 'v1x,v1y,v1z;v2x,v2y,v2z...' where 'v1x', 'v2x' ..., 'v1y', 'v2y' ... and
                                   'v1z', 'v2z' ... are floats that represents respectively the x, y and z values of the vectord along 
                                   which the q vectors should be generated.
            * subset            -- a selection string specifying the atoms to consider for the analysis.
            * deuteration       -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights           -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                   scheme to use.
            * pyroserver        -- a string specifying if Pyro will be used and how to run the analysis.
            * output            -- the output NetCDF file name.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'ARDCSF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                    
        # The frequency values.
        freqMax = 1.0/(2.0*self.dt)
        df = freqMax/self.nFrames
        self.frequencies = N.arange(self.nFrames + 1) * df        

        self.modelRealPart = N.zeros((self.nQValues, self.arModelOrder + 2), typecode = N.Float)
        self.modelImagPart = N.zeros((self.nQValues, self.arModelOrder + 2), typecode = N.Float)
        self.FQT = N.zeros((self.nQValues, self.nFrames), typecode = N.Float)
        self.FQTMemory = N.zeros((self.nQValues, self.arModelOrder + self.arModelOrder/2), typecode = N.Float)
        self.SQW = N.zeros((self.nQValues, self.nFrames + 1), typecode = N.Float)
        
        temp = [self.atomInformation[ind]['weight'] for ind in sorted(self.atomInformation.keys()) if self.mask[ind]]

        self.weights = N.array(temp, typecode = N.Float)[:, N.NewAxis]
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.mask = N.zeros(self.universe.numberOfAtoms(), typecode = N.Int0)
        for aIndex in self.subset:
            self.mask[aIndex] = 1
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        if (self.arModelOrder <= 0) or (self.arModelOrder >= self.nFrames):
            raise Error('The AR order must be an integer in [1,%d[.' % self.nFrames)

        qv = QVectors(self.universe, self.qVectorsDict)
                    
        self.qRadii = qv.qRadii
        self.qVectors = qv.qVectors
        self.qVectorsStatistics = qv.statistics

        self.nQValues = len(self.qRadii)
                                
        # Loop over the user-defined q values
        for qIndex in range(self.nQValues):
                            
            # For each q vector length, the corresponding set of q vector is transform into an array
            qarray = N.array(self.qVectors[qIndex], typecode = N.Float)
            qarray = self.universe._boxToRealPointArray(qarray)
            # qVect[1][j] is now of shape (3,Qcount)
            self.qVectors[qIndex] = N.transpose(qarray)
            
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                        
    def calc(self, qIndex, trajectory):
        """Calculates the contribution for one Q-shell.
        
        @param qIndex: the index of the Q-shell in |self.qRadii| list.
        @type qIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        rho = N.zeros((self.nFrames, self.qVectors[qIndex].shape[1]), typecode = N.Complex)
        
        comp = 0
        for comp in range(self.nFrames):
            
            frameIndex = self.frameIndexes[comp]
            
            conf = trajectory.configuration[frameIndex]
            conf.convertToBoxCoordinates()
            
            selAtoms = N.compress(self.mask, conf.array, 0)
            
            rho[comp,:] = N.add.reduce(self.weights*N.exp(1j*N.dot(selAtoms, self.qVectors[qIndex])))
                    
        return qIndex, rho

    def combine(self, qIndex, x):
        """
        """

        model = AveragedAutoRegressiveModel(self.arModelOrder, self.dt)

        for i in range(x.shape[1]):
            data = x[:, i]
            data = data - N.add.reduce(data)/len(data)
            m = AutoRegressiveModel(self.arModelOrder, data, self.dt)
            model.add(m)

        parameters = N.concatenate((model.coeff[::-1], N.array([model.sigma, model.variance], typecode = N.Float)))
        
        self.modelRealPart[qIndex, :] = parameters.real
        self.modelImagPart[qIndex, :] = parameters.imag
        
        average = N.add.reduce(N.add.reduce(x))/N.multiply.reduce(x.shape)            
        self.FQT[qIndex,:] = model.correlation(self.nFrames).values.real + (average*N.conjugate(average)).real
            
        mem = model.memoryFunction(self.arModelOrder + self.arModelOrder/2).values.real
        self.FQTMemory[qIndex,:] = mem[:]
        
        spectrum = model.spectrum(2.0*N.pi*self.frequencies)
        self.SQW[qIndex,:] = spectrum.values
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
        
        outputFile = NetCDFFile(self.output, 'w')
        outputFile.title = '%s (order %d)' % (self.__class__.__name__, self.arModelOrder)
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        outputFile.createDimension('NQVALUES', self.nQValues)
        outputFile.createDimension('NTIMES', self.nFrames)
        outputFile.createDimension('NTIMES_MEMORY', self.arModelOrder + self.arModelOrder/2)
        outputFile.createDimension('NPOLES', self.arModelOrder + 2)
        outputFile.createDimension('NFREQUENCIES', self.nFrames + 1)
        outputFile.createDimension('NOCTANS', 8)
        outputFile.createDimension('OCTANNAME', 8)
        
        OCTANS = outputFile.createVariable('octan', N.Character, ('NOCTANS','OCTANNAME'))
        OCTANS[:,:] = N.array([list('X-.Y-.Z-'),list('X-.Y-.Z+'),list('X-.Y+.Z-'),list('X-.Y+.Z+'),\
                                list('X+.Y-.Z-'),list('X+.Y-.Z+'),list('X+.Y+.Z-'),list('X+.Y+.Z+')])
        OCTANS.units = 'unitless'
                
        STAT = outputFile.createVariable('qvectors_statistics', N.Int32, ('NQVALUES', 'NOCTANS'))
        STAT[:,:] = self.qVectorsStatistics        
        
        # Creation of the NetCDF variable |sf| that will store the structure factor    
        TIMES = outputFile.createVariable('time', N.Float, ('NTIMES',))
        TIMES[:] = self.times[:]
        TIMES.units = 'ps'

        QLENGTH = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        QLENGTH[:] = self.qRadii[:]
        QLENGTH.units = 'nm^-1'

        TIMEMEMORY = outputFile.createVariable('time_memory', N.Float, ('NTIMES_MEMORY',))
        TIMEMEMORY[:] = self.dt * N.arange(self.arModelOrder + self.arModelOrder/2)
        TIMEMEMORY.units = 'ps'

        FREQUENCY = outputFile.createVariable('frequency', N.Float, ('NFREQUENCIES',))
        FREQUENCY[:] = self.frequencies[:]
        FREQUENCY.units = 'THz'
        
        MODELORDER = outputFile.createVariable('n', N.Int32, ('NPOLES',))
        MODELORDER[:] = N.arange(0, self.arModelOrder + 2)
        MODELORDER.units = 'unitless'

        MODELREAL = outputFile.createVariable('ar_coefficients_real', N.Float, ('NQVALUES', 'NPOLES'))
        MODELREAL[:,:] = self.modelRealPart[:,:]
        MODELREAL.units = 'unitless'

        MODELIMAG = outputFile.createVariable('ar_coefficients_imag', N.Float, ('NQVALUES', 'NPOLES'))
        MODELIMAG[:,:] = self.modelImagPart[:,:]
        MODELIMAG.units = 'unitless'
        
        FQT     = outputFile.createVariable('Fqt', N.Float, ('NQVALUES','NTIMES'))
        FQT[:,:] = self.FQT[:,:]
        FQT.units = 'unitless'
        
        FQTMEM  = outputFile.createVariable('Fqt_memory', N.Float, ('NQVALUES','NTIMES_MEMORY'))
        FQTMEM[:,:] = self.FQTMemory[:,:]
        FQTMEM.units = 'unitless'

        SQW = outputFile.createVariable('Sqw', N.Float, ('NQVALUES','NFREQUENCIES'))
        SQW[:,:] = self.SQW[:,:]
        SQW.units = 'unitless'
            
        outputFile.close()

        self.toPlot = None

#####################################################################################
# DYNAMIC INCOHERENT STRUCTURE FACTOR ANALYSIS
#####################################################################################
class DynamicIncoherentStructureFactor(Analysis):
    """Sets up an Dynamic Incoherent Structure Factor analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: DynamicIncoherentStructureFactorARModel(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory        -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo          -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                   number to consider, 'last' is an integer specifying the last frame number to consider and 
                                   'step' is an integer specifying the step number between two frames.
            * qshellvalues      -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                                   'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                                   the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * qshellwidth       -- a float specifying the width of the q shells.
            * qvectorspershell  -- a float specifying the number of q vectors to generate per q shell.
            * qvectorsgenerator -- a string being one of 'isotropic', 'anisotropic' or 'explicit' specifying the way the q vectors
                                   will be generated.
            * qvectorsdirection -- a string of the form 'v1x,v1y,v1z;v2x,v2y,v2z...' where 'v1x', 'v2x' ..., 'v1y', 'v2y' ... and
                                   'v1z', 'v2z' ... are floats that represents respectively the x, y and z values of the vectord along 
                                   which the q vectors should be generated.
            * resolution      -- a float specifying the width of the gaussian, that will be used to mimics the experimental resolution.
            * subset            -- a selection string specifying the atoms to consider for the analysis.
            * deuteration       -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights           -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                   scheme to use.
            * timeunits         -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * frequencyunits    -- a string equal to 'THz', 'rad s^-1', "cm^-1", "meV" or 'ueV' that specifies the frequency units.
            * qunits            -- a string equal to 'nm^-1' or 'ang^-1' that specifies the q units.
            * pyroserver        -- a string specifying if Pyro will be used and how to run the analysis.
            * output            -- the output NetCDF file name for the intermediate scattering function.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'DISF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                                                    
        self.FQT = {}

        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
            self.FQT[element] = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)
                            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()

        qv = QVectors(self.universe, self.qVectorsDict)
                    
        self.qRadii = qv.qRadii
        self.qVectors = qv.qVectors
        self.qVectorsStatistics = qv.statistics
            
        self.nQValues = len(self.qRadii)

        for qIndex in range(self.nQValues):
            # For each q vector length, the corresponding set of q vector is transform into an array
            qArray = N.array(self.qVectors[qIndex], typecode = N.Float)
            self.qVectors[qIndex] = N.transpose(qArray)
            
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                    
    def calc(self, atomIndex, trajectory):
        """Calculates the atomic term.
        
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer

        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
        # last step with the selected step increment.
        series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array
                                    
        atomicISF = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)

        for qIndex in range(self.nQValues):

            expSeries = N.exp(1j*N.dot(series, self.qVectors[qIndex]))

            # The ensemble average of (exp(-iqR(0))*exp(iqR(t))) is replaced by an autocorrelation.
            res = correlation(expSeries, expSeries)[:self.nFrames]/self.qVectors[qIndex].shape[1]
            
            N.add(atomicISF[qIndex,:], res, atomicISF[qIndex,:])
                
        return atomIndex, atomicISF
    
    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
        
        element = self.atomInformation[atomIndex]['element']
        N.add(self.FQT[element], x, self.FQT[element])
                    
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """
                
        self.FQT['total'] = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)

        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']
            self.FQT[element] /= n
            N.add(self.FQT['total'], n*w*self.FQT[element], self.FQT['total'])

        # The frequencies are computed.
        freqs = N.arange(len(self.times))/(2.0*len(self.times)*self.dt)
        # The time resolution function is computed.
        timeResolution = gaussian(self.times, self.timeSigma)
        # The frequency resolution function is computed.
        frequencyFunction = gaussian(freqs, 1.0/(2.0*N.pi*self.timeSigma))
                    
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some NetCDF dimensions.
        outputFile.createDimension('NQVALUES', self.nQValues)
        outputFile.createDimension('NTIMES', self.nFrames)
        outputFile.createDimension('NFREQUENCIES', self.nFrames)
        outputFile.createDimension('NOCTANS', 8)
        outputFile.createDimension('OCTANNAME', 8)        

        QLENGTH = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        QLENGTH[:] = self.qConv*self.qRadii
        QLENGTH.units = self.qUnits

        TIMES = outputFile.createVariable('time', N.Float, ('NTIMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits
        
        # The frequencies are computed and stored inthe NetCDF file.
        FREQUENCIES  = outputFile.createVariable('frequency', N.Float, ('NFREQUENCIES',))
        FREQUENCIES[:] = self.frequencyConv*freqs[:]
        FREQUENCIES.units = self.frequencyUnits

        OMEGAS = outputFile.createVariable('angular_frequency', N.Float, ('NFREQUENCIES',))
        OMEGAS[:] = 2.0*N.pi*freqs[:]
        OMEGAS.units = 'rad ps-1'

        OCTANS = outputFile.createVariable('octan', N.Character, ('NOCTANS','OCTANNAME'))
        OCTANS[:,:] = N.array([list('X-.Y-.Z-'),list('X-.Y-.Z+'),list('X-.Y+.Z-'),list('X-.Y+.Z+'),\
                                list('X+.Y-.Z-'),list('X+.Y-.Z+'),list('X+.Y+.Z-'),list('X+.Y+.Z+')])
        OCTANS.units = 'unitless'
                
        STAT = outputFile.createVariable('qvectors_statistics', N.Int32, ('NQVALUES', 'NOCTANS'))
        STAT[:,:] = self.qVectorsStatistics

        # The time resolution function is written.
        TIME_RES = outputFile.createVariable('time_resolution', N.Float, ('NTIMES',))
        TIME_RES[:] = timeResolution[:]
        TIME_RES.units = 'unitless'

        # The time resolution function is written.
        FREQ_RES = outputFile.createVariable('frequency_resolution', N.Float, ('NFREQUENCIES',))
        FREQ_RES[:] = frequencyFunction[:]
        FREQ_RES.units = 'unitless'
        
        # The Fqt total is written.     
        DISF = outputFile.createVariable('Fqt-total', N.Float, ('NQVALUES','NTIMES'))
        DISF[:,:] = self.FQT['total']
        DISF.units = 'unitless'

        # Then the Sqw total is computed and written.     
        DISF = outputFile.createVariable('Sqw-total', N.Float, ('NQVALUES','NFREQUENCIES'))
                
        # A loop is done for each q where the scattering function was calculated.
        for qVal in range(self.nQValues):
            DISF[qVal] = 0.5 * self.dt * FFT(gaussianWindow(self.FQT['total'][qVal], timeResolution)).real[:self.nFrames]            
        DISF.units = 'unitless'
                
        for k in sorted(self.FQT.keys()):
            
            if k == 'total': continue
            
            # The partial Fqt is written. 
            DISF = outputFile.createVariable('Fqt-%s' % k, N.Float, ('NQVALUES','NTIMES'))
            DISF[:,:] = self.FQT[k]
            DISF.units = 'unitless'
            
            # The partial Sqw is computed and written. 
            DISF = outputFile.createVariable('Sqw-%s' % k, N.Float, ('NQVALUES','NTIMES'))
            
            # A loop is done for each q where the scattering function was calculated.
            for qVal in range(self.nQValues):
                DISF[qVal] = 0.5 * self.dt * FFT(gaussianWindow(self.FQT[k][qVal], timeResolution)).real[:self.nFrames]            
            DISF.units = 'unitless'
            
        # The NetCDF output file is closed.
        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'q', 'yVar' : 'time', 'zVar' : 'Fqt-total'}
        
#####################################################################################
# DYNAMIC INCOHERENT STRUCTURE FACTOR WITHIN AR MODEL FRAMEWORK
#####################################################################################
class AutoRegressiveDynamicIncoherentStructureFactor(Analysis):
    """Sets up an Dynamic Incoherent Structure Factor analysis using an Auto Regressive model.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: DynamicIncoherentStructureFactorARModel(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory        -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo          -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                   number to consider, 'last' is an integer specifying the last frame number to consider and 
                                   'step' is an integer specifying the step number between two frames.
            * armodelorder      -- an integer in [1, len(trajectory)[ specifying the order of the model
            * qshellvalues      -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                                   'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                                   the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * qshellwidth       -- a float specifying the width of the q shells.
            * qvectorspershell  -- a float specifying the number of q vectors to generate per q shell.
            * qvectorsgenerator -- a string being one of 'isotropic', 'anisotropic' or 'explicit' specifying the way the q vectors
                                   will be generated.
            * qvectorsdirection -- a string of the form 'v1x,v1y,v1z;v2x,v2y,v2z...' where 'v1x', 'v2x' ..., 'v1y', 'v2y' ... and
                                   'v1z', 'v2z' ... are floats that represents respectively the x, y and z values of the vectord along 
                                   which the q vectors should be generated.
            * subset            -- a selection string specifying the atoms to consider for the analysis.
            * deuteration       -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights           -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                   scheme to use.
            * pyroserver        -- a string specifying if Pyro will be used and how to run the analysis.
            * output            -- the output NetCDF file name for the intermediate scattering function.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'ARDISF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                
        # Array that will store the AR model for each Q-shell.
        self.ARModel = []
        
        # Loop over the user-defined q values
        [self.ARModel.append(AveragedAutoRegressiveModel(self.arModelOrder, self.dt)) for qRadius in self.qRadii]
        
        # The frequency values.
        freqMax = 1.0/(2.0*self.dt)
        df = freqMax/self.nFrames
        self.frequencies = N.arange(self.nFrames + 1) * df        

        self.modelRealPart = N.zeros((self.nQValues, self.arModelOrder + 2), typecode = N.Float)
        self.modelImagPart = N.zeros((self.nQValues, self.arModelOrder + 2), typecode = N.Float)
        self.FQT = N.zeros((self.nQValues, self.nFrames), typecode = N.Float)
        self.FQTMemory = N.zeros((self.nQValues, self.arModelOrder + self.arModelOrder/2), typecode = N.Float)
        self.SQW = N.zeros((self.nQValues, self.nFrames + 1), typecode = N.Float)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        if (self.arModelOrder <= 0) or (self.arModelOrder >= self.nFrames):
            raise Error('The order of the AR model must be in [1,%s[' % self.nFrames)

        qv = QVectors(self.universe, self.qVectorsDict)
                    
        self.qRadii = qv.qRadii
        self.qVectors = qv.qVectors
        self.qVectorsStatistics = qv.statistics

        self.nQValues = len(self.qRadii)
                        
        # Loop over the user-defined q values
        for qIndex in range(self.nQValues):
                            
            # For each q vector length, the corresponding set of q vector is transform into an array
            qarray = N.array(self.qVectors[qIndex], typecode = N.Float)
            # qVect[1][j] is now of shape (3,Qcount)
            self.qVectors[qIndex] = N.transpose(qarray)
            
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')

    def calc(self, atomIndex, trajectory):
        """Calculates the atomic term.
        
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
                    
        series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array
        
        atomicFQT = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)
        atomicSQW = N.zeros((self.nQValues,self.nFrames + 1), typecode = N.Float)

        avgARModel = []
        
        if self.atomInformation[atomIndex]['weight'] == 0.0:
                        
            [avgARModel.append([]) for qRadius in self.qRadii]
                    
        else:
            w = self.atomInformation[atomIndex]['weight']

            for qIndex in range(self.nQValues):
            
                avgARModel.append([])
                rho = N.exp(1j*N.dot(series, self.qVectors[qIndex]))
                rho_av = N.add.reduce(rho)/rho.shape[0]
                qModel = AveragedAutoRegressiveModel(self.arModelOrder, self.dt)

                for i in range(rho.shape[1]):
                    data = rho[:, i] - rho_av[i]
                    m = AutoRegressiveModel(self.arModelOrder, data, self.dt)
                    qModel.add(m, w)
                    avgARModel[-1].append([m, w])
        
                rho_av_sq = (rho_av*N.conjugate(rho_av)).real
                average = N.add.reduce(rho_av_sq)/rho.shape[1]

                res = w*(qModel.correlation(self.nFrames).values.real + average)
                N.add(atomicFQT[qIndex,:], res, atomicFQT[qIndex,:])            
            
                res = w*(qModel.spectrum(2.0*N.pi*self.frequencies).values)
                N.add(atomicSQW[qIndex,:], res, atomicSQW[qIndex,:])
                        
        return atomIndex, (atomicFQT, atomicSQW, avgARModel)
                
    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """

        N.add(self.FQT, x[0], self.FQT)
        N.add(self.SQW, x[1], self.SQW)

        # Loop over the user-defined q values
        for qIndex in range(self.nQValues):
            for m, w in x[2][qIndex]:
                self.ARModel[qIndex].add(m,w)
                                            
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    

        for qIndex in range(self.nQValues):
            parameters = N.concatenate((self.ARModel[qIndex].coeff[::-1],\
                                        N.array([self.ARModel[qIndex].sigma, self.ARModel[qIndex].variance],\
                                                typecode = N.Float)))

            self.modelRealPart[qIndex, :] = parameters.real
            self.modelImagPart[qIndex, :] = parameters.imag
            
            mem = self.ARModel[qIndex].memoryFunction(self.arModelOrder + self.arModelOrder/2).values.real
            self.FQTMemory[qIndex,:] = mem

        outputFile = NetCDFFile(self.output, 'w')
        outputFile.title = '%s (order %d)' % (self.__class__.__name__, self.arModelOrder)
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        outputFile.createDimension('NQVALUES', self.nQValues)
        outputFile.createDimension('NTIMES', self.nFrames)
        outputFile.createDimension('NTIMES_MEMORY', self.arModelOrder + self.arModelOrder/2)
        outputFile.createDimension('NPOLES', self.arModelOrder + 2)
        outputFile.createDimension('NFREQUENCIES', self.nFrames + 1)
        outputFile.createDimension('NOCTANS', 8)
        outputFile.createDimension('OCTANNAME', 8)
        
        OCTANS = outputFile.createVariable('octan', N.Character, ('NOCTANS','OCTANNAME'))
        OCTANS[:,:] = N.array([list('X-.Y-.Z-'),list('X-.Y-.Z+'),list('X-.Y+.Z-'),list('X-.Y+.Z+'),\
                                list('X+.Y-.Z-'),list('X+.Y-.Z+'),list('X+.Y+.Z-'),list('X+.Y+.Z+')])
        OCTANS.units = 'unitless'
                
        STAT = outputFile.createVariable('qvectors_statistics', N.Int32, ('NQVALUES', 'NOCTANS'))
        STAT[:,:] = self.qVectorsStatistics        
        
        # Creation of the NetCDF variable |sf| that will store the structure factor    
        TIMES = outputFile.createVariable('time', N.Float, ('NTIMES',))
        TIMES[:] = self.times[:]
        TIMES.units = 'ps'

        QLENGTH = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        QLENGTH[:] = self.qRadii[:]
        QLENGTH.units = 'nm^-1'

        TIMEMEMORY = outputFile.createVariable('time_memory', N.Float, ('NTIMES_MEMORY',))
        TIMEMEMORY[:] = self.dt * N.arange(self.arModelOrder + self.arModelOrder/2)
        TIMEMEMORY.units = 'ps'

        FREQUENCY = outputFile.createVariable('frequency', N.Float, ('NFREQUENCIES',))
        FREQUENCY[:] = self.frequencies[:]
        FREQUENCY.units = 'THz'
        
        MODELORDER = outputFile.createVariable('n', N.Int32, ('NPOLES',))
        MODELORDER[:] = N.arange(0, self.arModelOrder + 2)
        MODELORDER.units = 'unitless'

        MODELREAL = outputFile.createVariable('ar_coefficients_real', N.Float, ('NQVALUES', 'NPOLES'))
        MODELREAL[:,:] = self.modelRealPart[:,:]
        MODELREAL.units = 'unitless'

        MODELIMAG = outputFile.createVariable('ar_coefficients_imag', N.Float, ('NQVALUES', 'NPOLES'))
        MODELIMAG[:,:] = self.modelImagPart[:,:]
        MODELIMAG.units = 'unitless'

        FQT     = outputFile.createVariable('Fqt', N.Float, ('NQVALUES','NTIMES'))
        FQT[:,:] = self.FQT[:,:]
        FQT.units = 'unitless'
        
        FQTMEM  = outputFile.createVariable('Fqt_memory', N.Float, ('NQVALUES','NTIMES_MEMORY'))
        FQTMEM[:,:] = self.FQTMemory[:,:]
        FQTMEM.units = 'unitless'

        SQW = outputFile.createVariable('Sqw', N.Float, ('NQVALUES','NFREQUENCIES'))
        SQW[:,:] = self.SQW[:,:]
        SQW.units = 'unitless'
            
        outputFile.close()
        
        self.toPlot = None
        
#####################################################################################
# DYNAMIC INCOHERENT STRUCTURE FACTOR ANALYSIS USING GAUSSIAN APPROXIMATION
#####################################################################################
class DynamicIncoherentStructureFactorGaussianApproximation(Analysis):
    """Sets up an Dynamic Incoherent Structure Factor analysis within Gaussian approximation.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: DynamicIncoherentStructureFactorGaussian(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory   -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo     -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                              number to consider, 'last' is an integer specifying the last frame number to consider and 
                              'step' is an integer specifying the step number between two frames.
            * qshellvalues -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                              'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                              the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * resolution      -- a float specifying the width of the gaussian, that will be used to mimics the experimental resolution.
            * subset       -- a selection string specifying the atoms to consider for the analysis.
            * deuteration  -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights      -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                              scheme to use.
            * timeunits         -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * frequencyunits    -- a string equal to 'THz', 'rad s^-1', "cm^-1", "meV" or 'ueV' that specifies the frequency units.
            * qunits            -- a string equal to 'nm^-1' or 'ang^-1' that specifies the q units.
            * pyroserver   -- a string specifying if Pyro will be used and how to run the analysis.
            * output       -- the output NetCDF file name for the intermediate scattering function.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'DISFGA')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """

        # The input parameters are parsed.
        self.interpreteInputParameters()
                                        
        self.FQT = {}

        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
            self.FQT[element] = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()

        self.nQValues = len(self.qShellValues)

        # |qSquare| is an array whose values are the squared of the q values list. These are the q^2 values used in
        # equation 3.35 
        self.qSquare = N.array(self.qShellValues, typecode = N.Float)*N.array(self.qShellValues, typecode = N.Float)

        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                    
    def calc(self, atomIndex, trajectory):
        """Calculates the atomic term.
        
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer

        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
        # last step with the selected step increment.
        series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array
                    
        atomicISFG = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)

        res  = self.getMSD(series)[:self.nFrames]
        
        # This is the formula 3.36. The factor 6 is in fact 2*3 where 3 account for the isotropic calculation.
        for qVal in range(self.nQValues):
            gaussian = N.exp(-res*self.qSquare[qVal]/6.)

            N.add(atomicISFG[qVal,:], gaussian, atomicISFG[qVal,:])
                        
        return atomIndex, atomicISFG
    
    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
        
        element = self.atomInformation[atomIndex]['element']
        N.add(self.FQT[element], x, self.FQT[element])
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """
            
        self.FQT['total'] = N.zeros((self.nQValues,self.nFrames), typecode = N.Float)
        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']
            self.FQT[element] /= n
            N.add(self.FQT['total'], n*w*self.FQT[element], self.FQT['total'])

        # The frequencies are computed.
        freqs = N.arange(len(self.times))/(2.0*len(self.times)*self.dt)
        # The time resolution function is computed.
        timeResolution = gaussian(self.times, self.timeSigma)
        # The frequency resolution function is computed.
        frequencyFunction = gaussian(freqs, 1.0/(2.0*N.pi*self.timeSigma))
            
        # The NetCDF output file intermediate scattering function.
        outputFile       = NetCDFFile(self.output, 'w')
        # The title for the NetCDF file.
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some NetCDF dimensions.
        outputFile.createDimension('NQVALUES', self.nQValues)
        outputFile.createDimension('NTIMES', self.nFrames)
        outputFile.createDimension('NFREQUENCIES', self.nFrames)

        # Creation of the NetCDF variable |sf| that will store the structure factor    
        Qlength = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        Qlength[:] = self.qConv*N.sqrt(self.qSquare)
        Qlength.units = self.qUnits
        
        TIME    = outputFile.createVariable('time', N.Float, ('NTIMES',))
        TIME[:] = self.timeConv*self.times[:]
        TIME.units = self.timeUnits

        # The frequencies are computed and stored inthe NetCDF file.
        FREQUENCIES  = outputFile.createVariable('frequency', N.Float, ('NFREQUENCIES',))
        FREQUENCIES[:] = self.frequencyConv*freqs[:]
        FREQUENCIES.units = self.frequencyUnits

        OMEGAS = outputFile.createVariable('angular_frequency', N.Float, ('NFREQUENCIES',))
        OMEGAS[:] = 2.0*N.pi*freqs[:]
        OMEGAS.units = 'rad ps-1'

        # The time resolution function is written.
        TIME_RES = outputFile.createVariable('time_resolution', N.Float, ('NTIMES',))
        TIME_RES[:] = timeResolution[:]
        TIME_RES.units = 'unitless'

        # The time resolution function is written.
        FREQ_RES = outputFile.createVariable('frequency_resolution', N.Float, ('NFREQUENCIES',))
        FREQ_RES[:] = frequencyFunction[:]
        FREQ_RES.units = 'unitless'

        # The Fqt total is written.     
        DISFGA = outputFile.createVariable('Fqt-total', N.Float, ('NQVALUES','NTIMES'))
        DISFGA[:,:] = self.FQT['total']
        DISFGA.units = 'unitless'

        # Then the Sqw total is computed and written.     
        DISFGA = outputFile.createVariable('Sqw-total', N.Float, ('NQVALUES','NFREQUENCIES'))
                
        # A loop is done for each q where the scattering function was calculated.
        for qVal in range(self.nQValues):
            DISFGA[qVal] = 0.5 * self.dt * FFT(gaussianWindow(self.FQT['total'][qVal], timeResolution)).real[:self.nFrames]            
        DISFGA.units = 'unitless'

        for k in sorted(self.FQT.keys()):
            
            if k == 'total': continue
            
            # The partial Fqt is written.     
            DISFGA = outputFile.createVariable('Fqt-%s' % k, N.Float, ('NQVALUES','NTIMES'))
            DISFGA[:,:] = self.FQT[k]
            DISFGA.units = 'unitless'
            
            # The partial Sw is computed and written.     
            DISFGA = outputFile.createVariable('Sqw-%s' % k, N.Float, ('NQVALUES','NTIMES'))
            # A loop is done for each q where the scattering function was calculated.
            for qVal in range(self.nQValues):
                DISFGA[qVal] = 0.5 * self.dt * FFT(gaussianWindow(self.FQT[k][qVal], timeResolution)).real[:self.nFrames]            
            DISFGA.units = 'unitless'
                        
        Q2length = outputFile.createVariable('q2', N.Float, ('NQVALUES',))
        Q2length[:] = self.qSquare
        Q2length.units = 'nm^-2'

        # The NetCDF output file is closed.
        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'q', 'yVar' : 'time', 'zVar' : 'Fqt-total'}

#        DynamicStructureFactor(self.output, self.timeSigma)                

    def getMSD(self, series):
        """
        Computes the atomic component of the Mean-Square-Displacement.
        This is the exact copy of the version written in nMOLDYN.Simulations.Dynamics but rewritten here
        for to keep the module Scattering independant from module Dynamics.  
        """

        # dsq is the squared norm of the position for each time step
        # dsq refers to DSQ(k) in the published algorithm
        dsq = N.add.reduce(series * series,1)

        # sum_dsq1 is the cumulative sum of dsq
        sum_dsq1 = N.add.accumulate(dsq)

        # sum_dsq1 is the reversed cumulative sum of dsq
        sum_dsq2 = N.add.accumulate(dsq[::-1])

        # sumsq refers to SUMSQ in the published algorithm
        sumsq = 2.*sum_dsq1[-1]

        # this line refers to the instruction SUMSQ <-- SUMSQ - DSQ(m-1) - DSQ(N - m) of the published algorithm
        # In this case, msd is an array because the instruction is computed for each m ranging from 0 to len(traj) - 1
        # So, this single instruction is performing the loop in the published algorithm
        Saabb  = sumsq - N.concatenate(([0.], sum_dsq1[:-1])) - N.concatenate(([0.], sum_dsq2[:-1]))

        # Saabb refers to SAA+BB/(N-m) in the published algorithm
        # Sab refers to SAB(m)/(N-m) in the published algorithm
        Saabb = Saabb / (len(dsq) - N.arange(len(dsq)))
        Sab   = 2.*correlation(series)

        # this line refers to the instruction MSD(m) <-- SUMSQ - 2Sab(m) in the published algorithm
        return (Saabb - Sab)
                            
#####################################################################################
# ELASTIC INCOHERENT STRUCTURE FACTOR
#####################################################################################
class ElasticIncoherentStructureFactor(Analysis):
    """Sets up an Elastic Incoherent Structure Factor.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: ElasticIncoherentStructureFactor(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory        -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo          -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                   number to consider, 'last' is an integer specifying the last frame number to consider and 
                                   'step' is an integer specifying the step number between two frames.
            * qshellvalues      -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                                   'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                                   the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * qshellwidth       -- a float specifying the width of the q shells.
            * qvectorspershell  -- a float specifying the number of q vectors to generate per q shell.
            * qvectorsgenerator -- a string being one of 'isotropic', 'anisotropic' or 'explicit' specifying the way the q vectors
                                   will be generated.
            * qvectorsdirection -- a string of the form 'v1x,v1y,v1z;v2x,v2y,v2z...' where 'v1x', 'v2x' ..., 'v1y', 'v2y' ... and
                                   'v1z', 'v2z' ... are floats that represents respectively the x, y and z values of the vectord along 
                                   which the q vectors should be generated.
            * subset            -- a selection string specifying the atoms to consider for the analysis.
            * deuteration       -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights           -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                   scheme to use.
            * qunits            -- a string equal to 'nm^-1' or 'ang^-1' that specifies the q units.
            * pyroserver        -- a string specifying if Pyro will be used and how to run the analysis.
            * output            -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                                   instead of the '.nc' extension.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'EISF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                                        
        self.EISF = {}
        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            self.EISF[element] = N.zeros((self.nQValues,), typecode = N.Float)

        for qIndex in range(self.nQValues):
                
            # For each q vector length, the corresponding set of q vector is transform into an array
            qArray = N.array(self.qVectors[qIndex], typecode = N.Float)

            # For each q vector length, the corresponding set of q vector is transform into an array
            # qVect[1][qVal] is now of shape (3,Qcount)
            self.qVectors[qIndex] = N.transpose(qArray)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)

        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()
            
        qv = QVectors(self.universe, self.qVectorsDict)
                    
        self.qRadii = qv.qRadii
        self.qVectors = qv.qVectors
        self.qVectorsStatistics = qv.statistics
                        
        self.nQValues = len(self.qRadii)

        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                                
    def calc(self, atomIndex, trajectory):
        """Calculates the atomic term.
        
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer

        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
        
        # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
        # last step with the selected step increment.
        series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array
                            
        atomicEISF = N.zeros((self.nQValues,), typecode = N.Float)

        for qIndex in range(self.nQValues):            
            a = N.add.reduce(N.exp(1j * N.dot(series, self.qVectors[qIndex])))/self.nFrames
            a = (a*N.conjugate(a)).real
            atomicEISF[qIndex] = N.add.reduce(a)/self.qVectors[qIndex].shape[1]
            
        return atomIndex, atomicEISF
    
    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
        
        element = self.atomInformation[atomIndex]['element']
        N.add(self.EISF[element], x, self.EISF[element])
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """
                
        self.EISF['total'] = N.zeros((self.nQValues,), typecode = N.Float)
        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']            
            self.EISF[element] /= n
            N.add(self.EISF['total'], n*w*self.EISF[element], self.EISF['total'])

        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NQVALUES', self.nQValues)
        outputFile.createDimension('NOCTANS', 8)
        outputFile.createDimension('OCTANNAME', 8)
        
        OCTANS = outputFile.createVariable('octan', N.Character, ('NOCTANS','OCTANNAME'))
        OCTANS[:,:] = N.array([list('X+.Y+.Z+'),list('X+.Y+.Z-'),list('X+.Y-.Z+'),list('X+.Y-.Z-'),\
                                list('X-.Y+.Z+'),list('X-.Y+.Z-'),list('X-.Y-.Z+'),list('X-.Y-.Z-')])
        OCTANS.units = 'unitless'
                
        STAT = outputFile.createVariable('qvectors_statistics', N.Int32, ('NQVALUES', 'NOCTANS'))
        STAT[:,:] = self.qVectorsStatistics        

        # Creation of the NetCDF output variables.
        # The Q.
        QRADII = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        QRADII[:] = self.qConv*self.qRadii
        QRADII.units = self.qUnits

        for k in sorted(self.EISF.keys()):
            EISF = outputFile.createVariable('eisf-%s' % ''.join(k), N.Float, ('NQVALUES',))
            EISF[:] = self.EISF[k][:]
            EISF.units = 'unitless'

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'q', 'yVar' : 'eisf-total'}

        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
                                
#####################################################################################
# SMOOTHED STATIC COHERENT STRUCTURE FACTOR
#####################################################################################
class SmoothedStaticCoherentStructureFactor(Analysis):
    """Sets up an Smoothed Static Coherent Structure Factor.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: SmoothedStaticCoherentStructureFactor(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory        -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo          -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                   number to consider, 'last' is an integer specifying the last frame number to consider and 
                                   'step' is an integer specifying the step number between two frames.
            * qshellvalues      -- a string of the form 'qmin1:qmax1:dq1;qmin2:qmax2:dq2...' where 'qmin1', 'qmin2' ... , 
                                   'qmax1', 'qmax2' ... and 'dq1', 'dq2' ... are floats that represents respectively 
                                   the q minimum, the q maximum and the q steps for q interval 1, 2 ...
            * subset            -- a selection string specifying the atoms to consider for the analysis.
            * deuteration       -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights           -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                   scheme to use.
            * qunits            -- a string equal to 'nm^-1' or 'ang^-1' that specifies the q units.
            * pyroserver        -- a string specifying if Pyro will be used and how to run the analysis.
            * output            -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                                   instead of the '.nc' extension.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
    Comments:
        
        - The analysis is based on the angular averaged coherent static structure factor formula where the
          summation over the q vectors is replaced by an integral over the q space. The formula used is taken 
          from equation 2.35 of Fischer et al. Rep. Prog. Phys. 69 (2006) 233-299.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'SSCSF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                                                              
        self.sscsfIJ = N.zeros((self.nElements,self.nElements,self.nQValues), typecode = N.Float)

        # ...
        self.scaleconfig = N.zeros((3*self.nSelectedAtoms,), typecode = N.Float)

        self.qShellArray = N.array(self.qShellValues, typecode = N.Float)
            
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
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)
        
        self.setElementInformation()

        self.nQValues = len(self.qShellValues)        

        # The list of the species name found in the universe.
        self.elementsList = sorted(self.elementInformation.keys())

        # The number of elements found in the universe.
        self.nElements = len(self.elementsList)

        # The id of the specie to which belong each selected atom.
        self.elements = N.zeros((self.nSelectedAtoms,), typecode = N.Int32)

        for comp in range(self.nSelectedAtoms):
            aIndex = self.subset[comp]
            element = self.atomInformation[aIndex]['element']
            self.elements[comp] = self.elementsList.index(element)
            
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
        reverseCell = N.ravel(N.transpose(N.array([v for v in trajectory.universe.reciprocalBasisVectors()],\
                                                   typecode = N.Float)))

        sscsfIJTemp = N.zeros((self.nElements,self.nElements,self.nQValues), typecode = N.Float)
                            
        smoothed_static_coherent_structure_factor(trajectory.universe.contiguousObjectConfiguration().array,\
                                                  directCell, reverseCell, self.qShellArray,\
                                                  N.array(self.subset, typecode = N.Int32),\
                                                  self.elements, sscsfIJTemp, self.scaleconfig)
        
        N.multiply(sscsfIJTemp,2.0,sscsfIJTemp)
        
        return frameIndex, sscsfIJTemp
    
    def combine(self, frameIndex, x):
        """
        """
        
        N.add(self.sscsfIJ, x, self.sscsfIJ)
                        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                       

        self.scaleconfig = None
            
        # The SSCSF.
        SSCSF = {'total' : N.zeros((self.nQValues), typecode = N.Float)}
        for i in range(self.nElements):
            
            element1 = self.elementsList[i]

            nA = self.elementInformation[element1]['number']
            wA = self.elementInformation[element1]['weight']            
                    
            for j in range(i, self.nElements):
                
                element2 = self.elementsList[j]
                
                nB = self.elementInformation[element2]['number']
                nAB = N.sqrt(nA*nB)
                
                wB = self.elementInformation[element2]['weight']            
        
                wAB = wA*wB
                deltaAB = 1.0

                pair = element1 + '-' + element2
        
                if i == j:
                    SSCSF[pair] = self.sscsfIJ[i,j,:][:]
                
                else:
                    deltaAB = 0.0
                    SSCSF[pair] = self.sscsfIJ[i,j,:][:] + self.sscsfIJ[j,i,:][:]

                SSCSF[pair] = deltaAB + SSCSF[pair]/(float(self.nFrames)*nAB)

                N.add(SSCSF['total'],nAB*wAB*SSCSF[pair],SSCSF['total'])

        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NQVALUES', self.nQValues)

        # Creation of the NetCDF output variables.
        # The Q.
        QRADII = outputFile.createVariable('q', N.Float, ('NQVALUES',))
        QRADII[:] = self.qConv*self.qShellValues
        QRADII.units = self.qUnits

        for key, value in SSCSF.items():
            SSCSF = outputFile.createVariable('Sq-%s' % ''.join(key), N.Float, ('NQVALUES',))
            SSCSF[:] = value[:]
            SSCSF.units = 'unitless'

        asciiVar = sorted(outputFile.variables.keys())            
            
        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'q', 'yVar' : 'sscsf-total'}
        
        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
