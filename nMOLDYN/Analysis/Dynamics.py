"""Collections of classes for the determination of dynamics-related properties.

Classes:

    * MeanSquareDisplacement                   : sets up a Mean-Square-Displacement analysis.
    * RootMeanSquareDeviation                  : sets up a Root Mean-Square-Deviation analysis.
    * GyrationRadius                           : sets up a Gyration Radius analysis.
    * AngularCorrelation                       : sets up an Angular Correlation analysis.
    * CartesianVelocityAutoCorrelationFunction : sets up a Cartesian Velocity AutoCorrelation analysis. 
    * CartesianPositionAutoCorrelationFunction : sets up a Cartesian Position AutoCorrelation analysis. 
    * DensityOfStates                          : sets up a Density Of States analysis.
    * AutoRegressiveAnalysis                   : sets up an Auto-Regressive analysis.
"""

# The python distribution modules
import os
from time import asctime

# The ScientificPython modules
from Scientific import N 
from Scientific.Functions.Interpolation import InterpolatingFunction
from Scientific.IO.NetCDF import NetCDFFile
from Scientific.Signals.Models import AutoRegressiveModel, AveragedAutoRegressiveModel

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Analysis.Analysis import Analysis
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import convertNetCDFToASCII
from nMOLDYN.Mathematics.Analysis import correlation, differentiate, FFT, gaussian, gaussianWindow

#####################################################################################
# MEAN-SQUARE DISPLACEMENT ANALYSIS
#####################################################################################
class MeanSquareDisplacement(Analysis):
    """Sets up a Mean Square Displacement analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: MeanSquareDisplacement(parameters = None, statusBar = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory    -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo      -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                               number to consider, 'last' is an integer specifying the last frame number to consider and 
                               'step' is an integer specifying the step number between two frames.
            * projection    -- a string of the form 'vx,vy,vz' specifying the vector along which the analysis
                               will be computed. 'vx', 'vy', and 'vz' are floats specifying respectively the x, y and z value 
                               of that vector.
            * subset        -- a selection string specifying the atoms to consider for the analysis.
            * deuteration   -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights       -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                               scheme to use.
            * timeunits     -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * distanceunits -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.
            * output        -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                               instead of the '.nc' extension.
            * pyroserver    -- a string specifying if Pyro will be used and how to run the analysis.

        - statusBar -- an instance of nMOLDYN.Widgets.StatusBar if not None
        
    Comments:
    
        - The algorithm is based on the Fast Correlation Algorithm (FCA) algorithm
    """
    
    def __init__(self, parameters = None, statusBar = None):
        """The constructor.
        
        @param parameters: if not None, a dictionnary storing the input parameters names and their corresponding values.
        @type parameters: dict
        
        @param statusBar: if not None, an instance of nMOLDYN.GUI.Widgets.StatusBar. Will attach a status bar to the 
            selected analysis.
        @type statusBar: instance of nMOLDYN.GUI.Widgets.StatusBar
        """
        
        # The parent class is constructed.
        Analysis.__init__(self, parameters, statusBar)
                
        path = os.path.join(GVAR['nmoldyn_analysis'],'MSD')
        exec file(path) in None, self.__dict__
                
    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """

        # The input parameters are parsed and some class attributes are derived necessary for the anamysis
        # are derived from the parsing.
        self.interpreteInputParameters()
                                            
        self.MSD = {}
        
        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            
            # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
            self.MSD[element] = N.zeros((self.nFrames), typecode = N.Float)
                        
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
                            
        if self.projection is None:
            # dsq is the squared norm of the position for each time step
            # dsq refers to DSQ(k) in the published algorithm
            dsq = N.add.reduce(series * series,1)
        else:
            # if a projection vector is given, the trajectory is then projected onto the projection vector.
            series = N.dot(series,self.projection)
            dsq = series * series

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

        atomicMSD = Saabb - Sab
                                
        return atomIndex, atomicMSD
    
    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """

        element = self.atomInformation[atomIndex]['element']
        N.add(self.MSD[element], x, self.MSD[element])
                            
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        self.MSD['total'] = N.zeros((self.nFrames), typecode = N.Float)

        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']
            self.MSD[element] /= n            
            N.add(self.MSD['total'], n*w*self.MSD[element], self.MSD['total'])
                    
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        # The total and partial MSD.
        for k in self.MSD.keys():
            MSD = outputFile.createVariable('msd-%s' % k, N.Float, ('NFRAMES',))
            MSD[:] = self.distanceConv*self.distanceConv*self.MSD[k]
            MSD.units = self.distanceUnits + "^2"

        asciiVar = sorted(outputFile.variables.keys())
            
        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'msd-total'}

        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
           
#####################################################################################
# ROOT MEAN-SQUARE DEVIATION ANALYSIS
#####################################################################################
class RootMeanSquareDeviation(Analysis):
    """Sets up a Root Mean Square Deviation analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: RootMeanSquareDeviation(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory     -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo       -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                number to consider, 'last' is an integer specifying the last frame number to consider and 
                                'step' is an integer specifying the step number between two frames.
            * referenceframe -- an integer in [1,len(trajectory)] specifying which frame should be the reference.
            * subset         -- a selection string specifying the atoms to consider for the analysis.
            * deuteration    -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights        -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                scheme to use.
            * timeunits      -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * distanceunits  -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.
            * pyroserver     -- a string specifying if Pyro will be used and how to run the analysis.
            * output         -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'RMSD')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
        
        self.RMSD = {}
        
        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            
            # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
            self.RMSD[element] = N.zeros((self.nFrames), typecode = N.Float)

    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        if (self.referenceFrame < 0) or (self.referenceFrame > len(self.trajectory)):
            raise Error('The reference frame must be an integer in [1,%s].' % len(self.trajectory))
        
        self.buildTimeInfo()
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)

        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()
        
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
                    
#    def calc(self, frameIndex, trajname):
#        """Calculates the contribution for one frame.
#        
#        @param frameIndex: the index of the frame in |self.frameIndexes| array.
#        @type frameIndex: integer.
#        
#        @param trajname: the name of the trajectory file name.
#        @type trajname: string
#        """

#        if self.architecture[0] == 'no':
#            t = self.trajectory            
#        else:
#            # Load the whole trajectory set.
#            t = Trajectory(None, trajname, 'r')
            
#        frame = self.frameIndexes[frameIndex]
#        msd = t.configuration[frame] - t.configuration[self.referenceFrame]
#        msd = self.weights * msd * msd
#        self.RMSD[frameIndex] = N.add.reduce(msd)
            
#        return None

    def calc(self, atomIndex, trajectory):
        """Calculates the atomic term.
        
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer

        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        
        @note: an atom-by-atom implementation was prefered than a frame-by-frame
            implementation of the type:
            msd = t.configuration[frame] - t.configuration[self.referenceFrame]
            msd = self.weights * msd * msd
            self.RMSD[frameIndex] = N.sqrt(N.add.reduce(msd))
        """
            
        # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
        # last step with the selected step increment.
        series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array
                    
        N.add(series, -series[self.referenceFrame,:], series)

        atomicRMSD = N.add.reduce(series**2,1)
                
        return atomIndex, atomicRMSD

    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
        
        element = self.atomInformation[atomIndex]['element']
        N.add(self.RMSD[element], x, self.RMSD[element])
            
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """
        

        self.RMSD['total'] = N.zeros((self.nFrames), typecode = N.Float)        
        for element in self.elementInformation.keys():
            w = self.elementInformation[element]['weight']
            N.multiply(self.RMSD[element],w,self.RMSD[element])
            N.add(self.RMSD['total'], self.RMSD[element], self.RMSD['total'])
            self.RMSD[element] = N.sqrt(self.RMSD[element])

        self.RMSD['total'] = N.sqrt(self.RMSD['total'])
        
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        # The total and partial RMSD.
        for k in self.RMSD.keys():
            RMSD = outputFile.createVariable('rmsd-%s' % k, N.Float, ('NFRAMES',))
            RMSD[:] = self.distanceConv*self.RMSD[k]
            RMSD.units = self.distanceUnits

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()
    
        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'rmsd-total'}

        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)

#####################################################################################
# CARTESIAN VELOCITY AUTOCORRELATION FUNCTION ANALYSIS
#####################################################################################
class CartesianVelocityAutoCorrelationFunction(Analysis):
    """Sets up a Cartesian Velocity AutoCorrelation analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: CartesianVelocityAutoCorrelationFunction(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory      -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo        -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                 number to consider, 'last' is an integer specifying the last frame number to consider and 
                                 'step' is an integer specifying the step number between two frames.
            * differentiation -- an integer in [0,5] specifying the order of the differentiation used to get the velocities
                                 out of the coordinates. 0 means that the velocities are already present in the trajectory loaded
                                 for analysis.
            * projection      -- a string of the form 'vx,vy,vz' specifying the vector along which the analysis
                                 will be computed. 'vx', 'vy', and 'vz' are floats specifying respectively the x, y and z value 
                                 of that vector.
            * normalize       -- a string being one of 'Yes' or 'No' specifying whether the analysis should be normalized to 1
                                 at t = 0 ('Yes') or not ('No').
            * subset          -- a selection string specifying the atoms to consider for the analysis.
            * deuteration     -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights         -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                 scheme to use.
            * timeunits       -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * distanceunits   -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.
            * pyroserver      -- a string specifying if Pyro will be used and how to run the analysis.
            * output          -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'CVACF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()

        self.VACF = {}
        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
            self.VACF[element] = N.zeros((self.nFrames), typecode = N.Float)
                            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        # Parses the analysis specific parameters.        

        # The 'normalize' parameter must be a string.
        if isinstance(self.parameters['normalize'],str):
            
            p = self.parameters['normalize'].strip().lower()
            
            # If 'yes', sets the |self.normalize| attribute to True. The signal will be normalized.
            if p == 'yes':
                self.normalize = True
                
            # If 'no', sets the |self.normalize| attribute to False. The signal will not be normalized.
            elif p == 'no':
                self.normalize = False
                
            # Otherwise, raises an error.
            else:
                raise Error('Error when parsing "normalize" parameter: must be "yes" or "no".')
                    
        # Otherwise, raises an error.
        else:
            raise Error('Error when parsing "normalize" parameter: must be a string equal to "yes" or "no".')

        self.buildTimeInfo()
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()
        
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
        
        if self.differentiation != 0:

            # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
            # last step with the selected step increment.
            series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array

            # the x, y and z axis
            for axis in range(3):
                # diff is a function of CalcFunctions.py module of nMoldyn package
                # it generates the velocities out of positions
                # diffScheme is the algorithm used to derive velocities from position and can be one
                # of the keywords 1 (default), 2, 3, 4, 5
                # dt is the differentiation time
                series[:,axis] = differentiate(series[:,axis], self.differentiation, self.dt)
                
        else:
            # series contains the velocities of the selected atom at from timeInfo[O] to timeInfo[1] with a time
            # increment of timeInfo[2]
            series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip, variable = 'velocities').array

        if self.projection is None:
            # if the trajectory is not projected, the calculations is done on the actual dependant variable
            atomicVACF = correlation(series)/3.

        else:
            # if a projection vector is given, the trajectory is then projected onto the projection vector.
            projectedTrajectory = N.dot(series, self.projection)
            # if the trajectory has been projected, the calculations is performed on the projected trajectory
            atomicVACF = correlation(projectedTrajectory)
                            
        return atomIndex, atomicVACF

    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
        
        element = self.atomInformation[atomIndex]['element']
        N.add(self.VACF[element], x, self.VACF[element])

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """

        self.VACF['total'] = N.zeros((self.nFrames), typecode = N.Float)

        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']
            self.VACF[element] /= n            
            N.add(self.VACF['total'], n*w*self.VACF[element], self.VACF['total'])
        
        # It is of common use to normalize the VACF by its value at t=0.
        if self.normalize:
            
            for k in self.VACF.keys():            
                if self.VACF[k][0] == 0:
                    raise Error('The normalization factor is equal to zero !!!')
                else:
                    self.VACF[k] = self.VACF[k]/self.VACF[k][0]

        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        # The total and partial VACF.
        for k in self.VACF.keys():
            VACF = outputFile.createVariable('vacf-%s' % k, N.Float, ('NFRAMES',))
            if self.normalize:
                VACF[:] = self.VACF[k]
                VACF.units = 'Unitless'
            else:
                VACF[:] = self.distanceConv*self.distanceConv*self.VACF[k]
                VACF.units = self.distanceUnits+"^2*ps^-2"

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'vacf-total'}
        
        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)         

#####################################################################################
# CARTESIAN POSITION AUTOCORRELATION FUNCTION ANALYSIS
#####################################################################################
class CartesianPositionAutoCorrelationFunction(Analysis):
    """Sets up a Cartesian Position AutoCorrelation analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: CartesianPositionAutoCorrelationFunction(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory      -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo        -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                 number to consider, 'last' is an integer specifying the last frame number to consider and 
                                 'step' is an integer specifying the step number between two frames.
            * projection      -- a string of the form 'vx,vy,vz' specifying the vector along which the analysis
                                 will be computed. 'vx', 'vy', and 'vz' are floats specifying respectively the x, y and z value 
                                 of that vector.
            * normalize       -- a string being one of 'Yes' or 'No' specifying whether the analysis should be normalized to 1
                                 at t = 0 ('Yes') or not ('No').
            * subset          -- a selection string specifying the atoms to consider for the analysis.
            * deuteration     -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights         -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                 scheme to use.
            * timeunits       -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * distanceunits   -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.
            * pyroserver      -- a string specifying if Pyro will be used and how to run the analysis.
            * output          -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'CPACF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()

        self.PACF = {}
        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
            self.PACF[element] = N.zeros((self.nFrames), typecode = N.Float)
                            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        # Parses the analysis specific parameters.        

        # The 'normalize' parameter must be a string.
        if isinstance(self.parameters['normalize'],str):
            
            p = self.parameters['normalize'].strip().lower()
            
            # If 'yes', sets the |self.normalize| attribute to True. The signal will be normalized.
            if p == 'yes':
                self.normalize = True
                
            # If 'no', sets the |self.normalize| attribute to False. The signal will not be normalized.
            elif p == 'no':
                self.normalize = False
                
            # Otherwise, raises an error.
            else:
                raise Error('Error when parsing "normalize" parameter: must be "yes" or "no".')
                    
        # Otherwise, raises an error.
        else:
            raise Error('Error when parsing "normalize" parameter: must be a string equal to "yes" or "no".')

        self.buildTimeInfo()
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

        self.setElementInformation()
        
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
        
        # series contains the velocities of the selected atom at from timeInfo[O] to timeInfo[1] with a time
        # increment of timeInfo[2]
        series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip, variable = 'configuration').array
        series -= N.add.reduce(series)/len(series)

        if self.projection is None:
            # if the trajectory is not projected, the calculations is done on the actual dependant variable
            atomicPACF = correlation(series)/3.

        else:
            # if a projection vector is given, the trajectory is then projected onto the projection vector.
            projectedTrajectory = N.dot(series, self.projection)
            # if the trajectory has been projected, the calculations is performed on the projected trajectory
            atomicPACF = correlation(projectedTrajectory)
                            
        return atomIndex, atomicPACF

    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
        
        element = self.atomInformation[atomIndex]['element']
        N.add(self.PACF[element], x, self.PACF[element])

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """

        self.PACF['total'] = N.zeros((self.nFrames), typecode = N.Float)

        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']
            self.PACF[element] /= n            
            N.add(self.PACF['total'], n*w*self.PACF[element], self.PACF['total'])
        
        # It is of common use to normalize the PACF by its value at t=0.
        if self.normalize:
            
            for k in self.PACF.keys():            
                if self.PACF[k][0] == 0:
                    raise Error('The normalization factor is equal to zero !!!')
                else:
                    self.PACF[k] = self.PACF[k]/self.PACF[k][0]

        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        # The total and partial PACF.
        for k in self.PACF.keys():
            PACF = outputFile.createVariable('pacf-%s' % k, N.Float, ('NFRAMES',))
            if self.normalize:
                PACF[:] = self.PACF[k]
                PACF.units = 'Unitless'
            else:
                PACF[:] = self.distanceConv*self.distanceConv*self.PACF[k]
                PACF.units = self.distanceUnits+'^2'

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'vacf-total'}
        
        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)         

#####################################################################################
# CARTESIAN DENSITY OF STATES ANALYSIS
#####################################################################################
class CartesianDensityOfStates(Analysis):
    """Sets up a Cartesian Density Of States analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: CartesianDensityOfStates(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory      -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo        -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                 number to consider, 'last' is an integer specifying the last frame number to consider and 
                                 'step' is an integer specifying the step number between two frames.
            * differentiation -- an integer in [0,5] specifying the order of the differentiation used to get the velocities
                                 out of the coordinates. 0 means that the velocities are already present in the trajectory loaded
                                 for analysis.
            * projection      -- a string of the form 'vx,vy,vz' specifying the vector along which the analysis
                                 will be computed. 'vx', 'vy', and 'vz' are floats specifying respectively the x, y and z value 
                                 of that vector.
            * resolution      -- a float in ]0.0,100.0[ specifying the width of the gaussian, in percentage of the trajectory length
                                 that will be used in the smoothing procedure.
            * subset          -- a selection string specifying the atoms to consider for the analysis.
            * deuteration     -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights         -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                 scheme to use.
            * frequencyunits  -- a string equal to 'THz', 'rad s^-1', 'cm^-1', 'meV' or 'ueV' that specifies the frequency units.
            * pyroserver      -- a string specifying if Pyro will be used and how to run the analysis.
            * output          -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'CDOS')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                                        
        self.DOS = {}
        # Loop over all the atom symbol.      
        for element in self.elementInformation.keys():
            # For each |pairName| key, the entry is a subdictionnary that stores the intra and intermolecular distances histograms.
            self.DOS[element] = N.zeros((self.nFrames), typecode = N.Float)
                    
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

        self.resolutionFunction = gaussian(self.times, self.timeSigma)
                        
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
        
        if self.differentiation != 0:

            # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
            # last step with the selected step increment.
            series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array

            # the x, y and z axis
            for axis in range(3):
                # diff is a function of CalcFunctions.py module of nMoldyn package
                # it generates the velocities out of positions
                # diffScheme is the algorithm used to derive velocities from position and can be one
                # of the keywords 1 (default), 2, 3, 4, 5
                # dt is the differentiation time
                series[:,axis] = differentiate(series[:,axis], self.differentiation, self.dt)
                
        else:
            # series contains the velocities of the selected atom at from timeInfo[O] to timeInfo[1] with a time
            # increment of timeInfo[2]
            series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip, variable = 'velocities').array
                
        if self.projection is None:
            # if the trajectory is not projected, the calculations is done on the actual dependant variable
            # The 1/3 factor is used to normalize the VACF
            aVACF = correlation(series)/3.
        else:
            # if a projection vector is given, the trajectory is then projected onto the projection vector.
            projectedTrajectory = N.dot(series, self.projection)            
            # if the trajectory has been projected, the calculations is performed on the projected trajectory
            aVACF = correlation(projectedTrajectory)
            
        # the dos is the fft of the product between the gaussian kernel and Sab
        # The DOS for atom 'at' is added to the dos array
        atomicDOS = FFT(gaussianWindow(aVACF, self.resolutionFunction)).real[:len(aVACF)]

        return atomIndex, atomicDOS

    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
        
        element = self.atomInformation[atomIndex]['element']
        N.add(self.DOS[element], x, self.DOS[element])
                
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """

        # 'freqencies' = 1D Numeric array. Frequencies at which the DOS was computed
        frequencies = N.arange(self.nFrames)/(2.0*self.nFrames*self.dt)

        self.DOS['total'] = N.zeros((self.nFrames), typecode = N.Float)

        for element in self.elementInformation.keys():
            n = self.elementInformation[element]['number']
            w = self.elementInformation[element]['weight']
            self.DOS[element] *= 0.5*self.dt/n
            N.add(self.DOS['total'], n*w*self.DOS[element], self.DOS['total'])

        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The frequencies.
        FREQUENCIES = outputFile.createVariable('frequency', N.Float, ('NFRAMES',))
        FREQUENCIES[:] = self.frequencyConv*frequencies[:]
        FREQUENCIES.units = self.frequencyUnits

        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times[:]
        TIMES.units = 'ps'
        
        # The resolution function.
        RESOLUTIONFUNCTION = outputFile.createVariable('resolution_function', N.Float, ('NFRAMES',))
        RESOLUTIONFUNCTION[:] = self.resolutionFunction[:]
        RESOLUTIONFUNCTION.units = 'unitless'

        # The total and partial DOS.
        for k in self.DOS.keys():
            DOS = outputFile.createVariable('dos-%s' % k, N.Float, ('NFRAMES',))
            DOS[:] = self.DOS[k]
            DOS.units = 'nm^2*ps^-1'

        asciiVar = sorted(outputFile.variables.keys())
        
        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'frequency', 'yVar' : 'dos-total'}

        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)

#####################################################################################
# AUTO-REGRESSIVE ANALYSIS ALYSIS
#####################################################################################
class AutoRegressiveAnalysis(Analysis):
    """Sets up an AutoRegressive Analysis analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: AutoRegressiveAnalysis(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory      -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo        -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                 number to consider, 'last' is an integer specifying the last frame number to consider and 
                                 'step' is an integer specifying the step number between two frames.
            * differentiation -- an integer in [0,5] specifying the order of the differentiation used to get the velocities
                                 out of the coordinates. 0 means that the velocities are already present in the trajectory loaded
                                 for analysis.
            * projection      -- a string of the form 'vx,vy,vz' specifying the vector along which the analysis
                                 will be computed. 'vx', 'vy', and 'vz' are floats specifying respectively the x, y and z value 
                                 of that vector.
            * armodelorder    -- an integer in [1, len(trajectory)[ specifying the order of the model
            * subset          -- a selection string specifying the atoms to consider for the analysis.
            * deuteration     -- a selection string specifying the hydrogen atoms whose atomic parameters will be those of the deuterium.
            * weights         -- a string equal to 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber' that specifies the weighting
                                 scheme to use.
            * pyroserver      -- a string specifying if Pyro will be used and how to run the analysis.
            * output          -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'ARA')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
        
        self.ARModel = AveragedAutoRegressiveModel(self.arModelOrder, self.dt)

    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()

        if (self.arModelOrder <= 0) or (self.arModelOrder >= self.nFrames):
            raise Error('The AR order must be an integer in [1,%d[.' % self.nFrames)
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
        
        self.deuterated = self.selectAtomsForDeuteration(self.deuterationDefinition)
        
        self.defineWeights(self.subset, self.deuterated, self.weightingScheme)

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
            
        if self.differentiation != 0:

            # series = 2D Numeric array. The positions of the selected atom |at| from the first step to the
            # last step with the selected step increment.
            series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip).array

            # the x, y and z axis
            for axis in range(3):
                # diff is a function of CalcFunctions.py module of nMoldyn package
                # it generates the velocities out of positions
                # diffScheme is the algorithm used to derive velocities from position and can be one
                # of the keywords 1 (default), 2, 3, 4, 5
                # dt is the differentiation time
                series[:,axis] = differentiate(series[:,axis], self.differentiation, self.dt)
                
        else:
            
            # series contains the velocities of the selected atom at from timeInfo[O] to timeInfo[1] with a time
            # increment of timeInfo[2]
            series = trajectory.readParticleTrajectory(atomIndex, first = self.first, last = self.last, skip = self.skip, variable = 'velocities').array
                            
        return atomIndex, series

    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """

        w = self.atomInformation[atomIndex]['weight']
        # The x, y and z axis
        for axis in range(3):
            self.ARModel.add(AutoRegressiveModel(self.arModelOrder, x[:, axis], self.dt), w)
                
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """               

        # This is the former fonction "evaluate_model". This fonction was used only twice 
        # in the AutoRegressiveAnalysis and the AutoRegressiveAnalysisXYZ functions. As these
        # functions has been fused in the current implementation no need for a function.

        # The VACF within the ARA framework.        
        vacf = self.ARModel.correlation(self.nFrames)
        c = vacf.values
        c_init = N.fabs(c[0].real)
        for i in range(len(c)):
            if N.fabs(c[i].real)/c_init < 1.e-10:
                break
        if i < len(c):
            vacf = InterpolatingFunction((vacf.axes[0][:i],), c[:i])

        # The DOS within the ARA framework.        
        omega_max = 1.1*N.pi/self.ARModel.delta_t
        omega = omega_max*N.arange(self.nFrames) / float(self.nFrames)
        dos = self.ARModel.spectrum(omega)
        dos = InterpolatingFunction((omega/(2.0*N.pi),), dos.values)
        dos.axes = (dos.axes[0],)

        # The MSD within the ARA framework.
        poles = self.ARModel.poles()
        cpoles = N.conjugate(poles)
        coeff0 = N.conjugate(self.ARModel.coeff[0])
        beta = N.zeros((self.ARModel.order), typecode = N.Complex)

        for i in range(self.ARModel.order):
            pole = poles[i]
            beta[i] = -(self.ARModel.sigsq*pole**(self.ARModel.order-1)/coeff0)/\
                (N.multiply.reduce((pole-poles)[:i])*\
                 N.multiply.reduce((pole-poles)[i+1:])*\
                 N.multiply.reduce(pole-1./cpoles)*\
                 self.ARModel.variance)
                
        beta = beta/sum(beta)
        msd = N.zeros((self.nFrames), typecode = N.Float)
        n = N.arange(self.nFrames)

        for i in range(self.ARModel.order):
            pole = poles[i]
            msd = msd + (beta[i]*((pole**n-1.)*pole/(1-pole)**2 + n/(1-pole))).real
        msd = 6.0 * self.ARModel.delta_t**2 * self.ARModel.variance * msd

        msd = InterpolatingFunction((self.ARModel.delta_t*n,), msd)
        
        # The Memory function
        friction = self.ARModel.frictionConstant()
        try:
            memory = self.ARModel.memoryFunction(self.ARModel.order+self.ARModel.order/2)
        except OverflowError:
            null = N.zeros((0), typecode = N.Float)
            memory = InterpolatingFunction((null,), null)

        memory.friction = friction        
                
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.comment = 'Friction constant: %20.15f' % memory.friction
        outputFile.comment += 'Variance: %.15f\nSigma: %.15f' % (self.ARModel.variance, self.ARModel.sigma)
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()
                    
        # Some dimensions are created.
        outputFile.createDimension('NFRAMES_VACF', len(vacf.axes[0]))
        outputFile.createDimension('NFRAMES_MSD', len(msd.axes[0]))
        outputFile.createDimension('NFRAMES_MEMORY', len(memory.axes[0]))
        outputFile.createDimension('FREQSTEPS', len(dos.axes[0]))
        outputFile.createDimension('ARMODELORDER', self.arModelOrder)
            
        # Creation of the NetCDF output variables.

        # VACF related variables.
        TIMES = outputFile.createVariable('time_vacf', N.Float, ('NFRAMES_VACF',))
        TIMES[:] = vacf.axes[0]
        TIMES.units = 'ps'

        VACF = outputFile.createVariable('vacf', N.Float, ('NFRAMES_VACF',))
        VACF[:] = vacf.values.real
        VACF.units = 'nm^2*ps^-2'

        # DOS related variables.
        FREQUENCIES = outputFile.createVariable('frequency', N.Float, ('FREQSTEPS',))
        FREQUENCIES[:] = dos.axes[0]
        FREQUENCIES.units = 'THz'

        OMEGA = outputFile.createVariable('angular_frequency', N.Float, ('FREQSTEPS',))
        OMEGA[:] = 2.0*N.pi*dos.axes[0]
        OMEGA.units = 'rad ps-1'

        DOS = outputFile.createVariable('dos', N.Float, ('FREQSTEPS',))
        DOS[:] = dos.values.real
        DOS.units = 'nm2*ps-1'

        # MSD related variables.
        TIMES = outputFile.createVariable('time_msd', N.Float, ('NFRAMES_MSD',))
        TIMES[:] = msd.axes[0]
        TIMES.units = 'ps'

        MSD = outputFile.createVariable('msd', N.Float, ('NFRAMES_MSD',))
        MSD[:] = msd.values.real
        MSD.units = 'nm^2'

        # Memory function related variables.
        TIMES = outputFile.createVariable('time_memory', N.Float, ('NFRAMES_MEMORY',))
        TIMES[:] = memory.axes[0]
        TIMES.units = 'ps'

        MEMORY = outputFile.createVariable('memory_function', N.Float, ('NFRAMES_MEMORY',))
        MEMORY[:] = memory.values.real
        MEMORY.units = 'unitless'
        
        # ARA parameters related variables.
        COEFFINDEXES = outputFile.createVariable('n', N.Int32, ('ARMODELORDER',))
        COEFFINDEXES[:] = 1 + N.arange(0, self.arModelOrder)
        COEFFINDEXES.units = 'unitless'

        COEFFS = outputFile.createVariable('ar_coefficients', N.Float, ('ARMODELORDER',))
        COEFFS[:] = self.ARModel.coeff[::-1]
        COEFFS.units = 'unitless'

        asciiVar = sorted(outputFile.variables.keys())

        outputFile.close()

        self.toPlot = None

        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
                        
#####################################################################################
# GYRATION RADIUS ANALYSIS
#####################################################################################
class RadiusOfGyration(Analysis):
    """Sets up a Radius Of Gyration analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: RadiusOfGyration(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory     -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo       -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                number to consider, 'last' is an integer specifying the last frame number to consider and 
                                'step' is an integer specifying the step number between two frames.
            * subset         -- a selection string specifying the atoms to consider for the analysis.
            * timeunits      -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * distanceunits  -- a string equal to 'nm', 'ang' or 'fm' that specifies the distance units.
            * pyroserver     -- a string specifying if Pyro will be used and how to run the analysis.
            * output         -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'ROG')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
        
        self.subsetMass = 0.0
            
        # The center of gravity.
        self.cog = N.zeros((self.nFrames,3), typecode = N.Float)

        # The sum of each position.
        self.sumRi = N.zeros((self.nFrames,3), typecode = N.Float)
            
        # ROG = 1D Numeric array. The array that stores the Gyration Radius over time.
        self.ROG = N.zeros((self.nFrames), typecode = N.Float)

    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)

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

        atom = [at for at in trajectory.universe.atomList() if at.index == atomIndex][0]

        return atomIndex, (atom.mass(),series)
    
    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """

        self.subsetMass += x[0]
        N.add(self.cog, x[0]*x[1], self.cog)
        N.add(self.sumRi, x[1], self.sumRi)
        N.add(self.ROG, N.add.reduce(x[1]**2,1), self.ROG)
                    
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
        
        self.cog = self.cog/self.subsetMass
                        
        N.add(self.ROG, -2.0 * N.add.reduce(self.sumRi*self.cog,1), self.ROG)
        N.add(self.ROG, float(self.nSelectedAtoms)*N.add.reduce(self.cog**2,1), self.ROG)
        
        self.ROG = N.sqrt(self.ROG/float(self.nSelectedAtoms))
                
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.timeConv*self.times[:]
        TIMES.units = self.timeUnits

        # The Gyration radius.
        ROG = outputFile.createVariable('rog', N.Float, ('NFRAMES',))
        ROG[:] = self.distanceConv*self.ROG[:]
        ROG.units = self.distanceUnits

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'rog'}

        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = ['time', 'rog'])  
        
#####################################################################################
# ANGULAR CORRELATION ANALYSIS
#####################################################################################
class AngularCorrelation(Analysis):
    """Sets up an Angular Correlation analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: AngularCorrelation(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory      -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo        -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                 number to consider, 'last' is an integer specifying the last frame number to consider and 
                                 'step' is an integer specifying the step number between two frames.
            * group           -- a selection string specifying the groups of three atoms that will define the plane on 
                                 which the angular correlation will be computed.
            * atomorder       -- a string of the form 'atom1,atom2,atom3' where 'atom1', 'atom2' and 'atom3' are 
                                 respectively the MMTK atom names of the atoms in the way they should be ordered.
            * timeunits       -- a string equal to 'ps', 'ns' or 'fs' that specifies the time units.
            * ac              -- the output NetCDF file name. A CDL version of this file will also be generated with the '.cdl' extension
                                 instead of the '.nc' extension.
            * pyroserver      -- a string specifying if Pyro will be used and how to run the analysis.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'AC')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()

        self.angle1 = N.zeros((self.nGroups,self.nFrames), typecode = N.Float)
        self.angle2 = N.zeros((self.nGroups,self.nFrames), typecode = N.Float)
        self.angle3 = N.zeros((self.nGroups,self.nFrames), typecode = N.Float)
                                
        self.angCorr1 = N.zeros((self.nGroups,self.nFrames), typecode = N.Float)
        self.angCorr2 = N.zeros((self.nGroups,self.nFrames), typecode = N.Float)
        self.angCorr3 = N.zeros((self.nGroups,self.nFrames), typecode = N.Float)
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
            
        self.buildTimeInfo()

        if len(self.atomOrder) != 3:
            raise Error('You must define exactly three atoms to set the correlation vectors v1, v2 and v3.')
                
        self.group = self.selectGroups(self.groupDefinition, self.atomOrder)
                                                            
        # The self.nGroups is needed when dealing with the template defined for analysis whose mainloop is over groups
        self.nGroups = len(self.group)

        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')

    def calc(self, atomIndexes, trajectory):
        """Calculates the contribution for one group.
        
        @param atomIndexes: the index of the atoms of the group.
        @type atomIndexes: list of integers.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        atoms = [None]*3

        for at in trajectory.universe.atomList():
            if at.index in atomIndexes:
                atoms[atomIndexes.index(at.index)] = at
                                    
        at1, at2, at3 = atoms
                
        at1Traj = trajectory.readParticleTrajectory(at1, first = self.first, last = self.last, skip = self.skip)
        at2Traj = trajectory.readParticleTrajectory(at2, first = self.first, last = self.last, skip = self.skip)
        at3Traj = trajectory.readParticleTrajectory(at3, first = self.first, last = self.last, skip = self.skip)
            
        v1 = N.zeros((self.nFrames,3), typecode = N.Float)
        v2 = N.zeros((self.nFrames,3), typecode = N.Float)
        v3 = N.zeros((self.nFrames,3), typecode = N.Float)

        for comp in range(self.nFrames):

            v12 = trajectory.universe.distanceVector(at1Traj[comp], at2Traj[comp]).normal()
            v13 = trajectory.universe.distanceVector(at1Traj[comp], at3Traj[comp]).normal()
            
            n1 = (v12 + v13).normal()
            n3 = v12.cross(n1).normal()                                        
            n2 = n3.cross(n1).normal()
                
            v1[comp,:] = N.array(n1, typecode = N.Float)
            v2[comp,:] = N.array(n2, typecode = N.Float)
            v3[comp,:] = N.array(n3, typecode = N.Float)
        
        return atomIndexes, (v1, v2, v3)
    
    def combine(self, atomIndexes, x):
        """
        """
        
        comp = self.group.index(atomIndexes)
        
        cosAng = N.dot(x[0][0], N.transpose(x[0]))
        sinAng = N.dot(N.cross_product(x[0][0],x[0]),x[2][0])
        ang = N.arccos(cosAng)
        self.angle1[comp,:] = N.where(sinAng >= 0, ang, -ang)*180.0/N.pi
                
        cosAng = N.dot(x[1][0], N.transpose(x[1]))
        sinAng = N.dot(N.cross_product(x[1][0],x[1]),x[2][0])
        ang = N.arccos(cosAng)
        self.angle2[comp,:] = N.where(sinAng >= 0, ang, -ang)*180.0/N.pi
        
        cosAng = N.dot(x[2][0], N.transpose(x[2]))
        sinAng = N.dot(N.cross_product(x[2][0],x[2]),x[2][0])
        ang = N.arccos(cosAng)
        self.angle3[comp,:] = N.where(sinAng >= 0, ang, -ang)*180.0/N.pi
                                        
        # The autocorrelation for v1, v2, v3 are computed. 
        # As they are unit vectors AC1(0) = AC2(0) = AC3(0)= 1
        self.angCorr1[comp,:] = correlation(x[0])
        self.angCorr2[comp,:] = correlation(x[1])
        self.angCorr3[comp,:] = correlation(x[2])
                    
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
        
        outputFile = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Dictionnary whose keys are of the form Gi where i is the group number 
        # and the entries are the list of the index of the atoms building the group.
        comp = 1
        for g in self.group:
            outputFile.jobinfo += 'Group %s: %s\n' % (comp,[index for index in g])
            comp += 1

        outputFile.createDimension('NGROUPS', self.nGroups)
        outputFile.createDimension('NFRAMES', self.nFrames)

        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.timeConv*self.times
        TIMES.units = self.timeUnits

        GROUPNUMBER = outputFile.createVariable('group_number', N.Int32, ('NGROUPS',))
        GROUPNUMBER[:] = 1 + N.arange(self.nGroups)

        AC1 = outputFile.createVariable('ac1_per_group', N.Float, ('NGROUPS','NFRAMES'))
        AC1[:,:] = self.angCorr1[:,:]

        AC2 = outputFile.createVariable('ac2_per_group', N.Float, ('NGROUPS','NFRAMES'))
        AC2[:,:] = self.angCorr2[:,:]

        AC3 = outputFile.createVariable('ac3_per_group', N.Float, ('NGROUPS','NFRAMES'))
        AC3[:,:] = self.angCorr3[:,:]

        ANGLE = outputFile.createVariable('angle1_per_group', N.Float, ('NGROUPS','NFRAMES'))
        ANGLE[:,:] = self.angle1[:,:]
        ANGLE.units = "angle (deg)"

        ANGLE = outputFile.createVariable('angle2_per_group', N.Float, ('NGROUPS','NFRAMES'))
        ANGLE[:,:] = self.angle2[:,:]
        ANGLE.units = "angle (deg)"

        ANGLE = outputFile.createVariable('angle3_per_group', N.Float, ('NGROUPS','NFRAMES'))
        ANGLE[:,:] = self.angle3[:,:]
        ANGLE.units = "angle (deg)"

        AC1AVG = outputFile.createVariable('ac1', N.Float, ('NFRAMES',))
        AC1AVG[:] = self.angCorr1.sum(0)/float(self.nGroups)

        AC2AVG = outputFile.createVariable('ac2', N.Float, ('NFRAMES',))
        AC2AVG[:] = self.angCorr2.sum(0)/float(self.nGroups)

        AC3AVG = outputFile.createVariable('ac3', N.Float, ('NFRAMES',))
        AC3AVG[:] = self.angCorr3.sum(0)/float(self.nGroups)

        asciiVar = sorted(outputFile.variables.keys())            

        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'ac3'}

        # Creates an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
