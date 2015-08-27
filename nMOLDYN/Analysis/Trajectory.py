"""Collections of classes for the determination of dynamics-related properties.

Classes:

    * QuasiHarmonicAnalysis        : sets up a Quasi-Harmonic Analysis analysis.
    * PassBandTrajectoryFilter     : sets up a Pass-Band Trajectory Filter analysis.
    * GlobalMotionTrajectoryFilter : sets up a Global Motion Trajectory Filter analysis.
    * CenterOfMassTrajectory       : sets up a Center Of Mass Trajectory analysis.
    * ReducedTrajectory            : sets up a Center Of Mass Trajectory analysis.
"""

# The python distribution modules
import copy
import operator
import os
from time import asctime
from timeit import default_timer

# The ScientificPython modules
from Scientific import N 
from Scientific.Geometry import Vector
from Scientific.IO.NetCDF import NetCDFFile
from Scientific.LA import Heigenvectors

# The MMTK distribution modules
from MMTK import Atom, AtomCluster
from MMTK import Units
from MMTK.Collections import Collection
from MMTK.ParticleProperties import Configuration, ParticleVector
from MMTK.Trajectory import SnapshotGenerator, Trajectory, TrajectoryOutput

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Analysis.Analysis import Analysis
from nMOLDYN.Core.Error import Error
from nMOLDYN.Mathematics.Analysis import FFT, invFFT

#####################################################################################
# PASSBANDFILTEREDTRAJECTORY ANALYSIS
#####################################################################################
class PassBandFilteredTrajectory(Analysis):
    """Sets up a Pass-Band Trajectory Filter analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: PassBandFilteredTrajectory(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo   -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                            number to consider, 'last' is an integer specifying the last frame number to consider and 
                            'step' is an integer specifying the step number between two frames.
            * filter     -- a string of the form 'low:high' where 'low' and 'high' are floats specifying respectively 
                            the lower and the upper bounds of the pass-band filter.
            * subset     -- a selection string specifying the atoms to consider for the analysis.
            * pyroserver -- a string specifying if Pyro will be used and how to run the analysis.
            * output     -- the output NetCDF file name.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'PBFT')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                
        # The frequencies.
        frequencies = N.arange(self.nFrames)/(2.0*self.nFrames*self.dt)
        # Computation of the index corresponding to the lowest frequency of the pass-band filter.
        indexMin = N.searchsorted(frequencies, self.freqMin)
        # Computation of the index corresponding to the highest frequency of the pass-band filter.
        indexMax = N.searchsorted(frequencies, self.freqMax)
            
        halfFilter = N.zeros(self.nFrames, N.Float)
        halfFilter[indexMin:indexMax] = 1

        self.pbFilter = N.zeros(2*self.nFrames - 2, N.Float)
        self.pbFilter[:self.nFrames] = halfFilter
        self.pbFilter[self.nFrames:] = halfFilter[-2:0:-1]

        # This dictionnary works a little bit like a Particle Tensor object. Its keys are the index of the atoms of
        # the selection. Its value are the filtered atomic trajectories.
        self.filteredTrajectory = {}

    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        # Parses the analysis specific parameters.

        # The 'filter' input parameter. Must be a string of the form 'fmin:fmax' where fmin and fmax are respectively 
        # the min and max frequencies of the filter.
        if isinstance(self.parameters['filter'],str):
            try:
                # Sets the |self.freqMin| and |self.freqMax| attributes to respectively the min and max frequencies.
                self.freqMin, self.freqMax = [float(v) for v in self.parameters['filter'].split(':')]
            except:
                raise Error('Error when parsing "filter" parameter: must be a string of the form fmin:fmax.')
                    
            # The min frequency must be >= 0.
            if self.freqMin < 0.0:
                raise Error('Error when parsing "filter" parameter: the min frequency must be >= 0.')

            # The max frequency must be > min frequency.
            if self.freqMin >= self.freqMax:
                raise Error('Error when parsing "filter" parameter: the max frequency must be > min frequency.')
                    
        # Otherwise raises an error.
        else:
            raise Error('Error when parsing "filter" parameter: must be a string.')

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

        unfiltered = N.zeros(2*self.nFrames - 2, 'd')
        atomicPBFT = N.zeros((self.nFrames,3), 'd')

        # Loop over the x, y and z coordinates.
        for coord in range(3):
            unfiltered[:self.nFrames] = series[:,coord]
            unfiltered[self.nFrames:] = series[:,coord][-2:0:-1]
            temp = FFT(unfiltered)
            atomicPBFT[:,coord] = invFFT(self.pbFilter*temp)[:self.nFrames].real

        return atomIndex, atomicPBFT
    
    def combine(self, atomIndex, x):
        """
        @param atomIndex: the index of the selected atom.
        @type atomIndex: integer
        """
            
        self.filteredTrajectory[atomIndex] = x
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """               
        
        if self.architecture == 'monoprocessor':
            t = self.trajectory
            
        else:
            # Load the whole trajectory set.
            t = Trajectory(None, self.trajectoryFilename, 'r')

        orderedAtoms = sorted(t.universe.atomList(), key = operator.attrgetter('index'))
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])
                
        # traj_new is the filtered trajectory
        outputFile = Trajectory(selectedAtoms, self.output, "w") 
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Create the snapshot generator
        snapshot = SnapshotGenerator(t.universe, actions = [TrajectoryOutput(outputFile, ["configuration","time"], 0, None, 1)])
                        
        # Loop over the output frames.
        for comp in range(self.nFrames):
            
            frameIndex = self.frameIndexes[comp]                        
            t.universe.setFromTrajectory(t, frameIndex)
            
            for comp1 in range(self.nSelectedAtoms):
                at = selectedAtoms[comp1]
                at.setPosition(self.filteredTrajectory[at.index][comp,:])
                
            snapshot(data = {'time': self.times[comp]})
                                    
        outputFile.close()

        t.close()
        
        self.toPlot = None

#####################################################################################
# REDUCED TRAJECCTORY ANALYSIS
#####################################################################################
class ReducedTrajectory(Analysis):
    """Sets up a Reduced Trajectory analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: ReducedTrajectory(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo   -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                            number to consider, 'last' is an integer specifying the last frame number to consider and 
                            'step' is an integer specifying the step number between two frames.
            * subset     -- a selection string specifying the atoms to consider for the analysis.
            * pyroserver -- a string specifying if Pyro will be used and how to run the analysis.
            * output     -- the output NetCDF file name.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'RT')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                
        self.reducedTrajectory = {}

    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
            
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)

        orderedAtoms = sorted(self.trajectory.universe.atomList(), key = operator.attrgetter('index'))
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])

        # traj_new is the filtered trajectory
        self.outputFile = Trajectory(selectedAtoms, self.output, "w")
                 
        # Create the snapshot generator
        self.snapshot = SnapshotGenerator(self.trajectory.universe, actions = [TrajectoryOutput(self.outputFile, ["all"], 0, None, 1)])

    def calc(self, frameIndex, trajectory):
        """Calculates the atomic term.
        
        @param frameIndex: the index of the frame.
        @type frameIndex: integer.

        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
        
        idx = self.frameIndexes.index(frameIndex)
        
        self.snapshot(data = {'time': self.times[idx]})
                                                    
        return frameIndex, None

    def combine(self, frameIndex, x):
        """
        @param frameIndex: the index of the frame.
        @type frameIndex: integer.
        """
            
        return

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """               
 
        self.outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()
                                                                                                   
        self.outputFile.close()

        self.trajectory.close()
        
        self.toPlot = None        

#    def interpreteInputParameters(self):
#        """Parse the input parameters for the analysis.
#        """
#
#        # Parses the parameters that are common to different analysis.
#        Analysis.interpreteInputParameters(self)
#        
#        self.buildTimeInfo()
#            
#        self.subset = self.selectAtoms(self.subsetDefinition)
#        self.nSelectedAtoms = len(self.subset)
#
#        if self.architecture != 'monoprocessor':
#            # The attribute trajectory is removed because it can not be pickled by Pyro.
#            delattr(self, 'trajectory')
                        
#    def calc(self, frameIndex, trajectory):
#        """Calculates the atomic term.
#        
#        @param frameIndex: the index of the frame.
#        @type frameIndex: integer.
#
#        @param trajectory: the trajectory.
#        @type trajectory: MMTK.Trajectory.Trajectory object
#        """
#
#        trajectory.universe.setFromTrajectory(trajectory, frameIndex)
#
#        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
#        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])
#        
#        frame = {}        
#        for at in selectedAtoms:
#            frame[at.index] = at.position()
#                                                    
#        return frameIndex, frame
            
#    def combine(self, frameIndex, x):
#        """
#        @param frameIndex: the index of the frame.
#        @type frameIndex: integer.
#        """
#            
#        self.reducedTrajectory[frameIndex] = x
#        
#    def finalize(self):
#        """Finalizes the calculations (e.g. averaging the total term, output files creations ...)
#        """               
#        
#        if self.architecture == 'monoprocessor':
#            t = self.trajectory
#            
#        else:
#            # Load the whole trajectory set.
#            t = Trajectory(None, self.trajectoryFilename, 'r')
#
#        orderedAtoms = sorted(t.universe.atomList(), key = operator.attrgetter('index'))
#        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])
#                
#        # traj_new is the filtered trajectory
#        outputFile = Trajectory(selectedAtoms, self.output, "w") 
#        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()
#
#        # Create the snapshot generator
#        snapshot = SnapshotGenerator(t.universe, actions = [TrajectoryOutput(outputFile, ["all"], 0, None, 1)])
#
#        # Loop over the output frames.
#        for comp in range(self.nFrames):
#            
#            frameIndex = self.frameIndexes[comp]
#            t.universe.setFromTrajectory(t, frameIndex)
#            
#            for at in selectedAtoms:
#                at.setPosition(Vector(self.reducedTrajectory[frameIndex][at.index]))
#                
#            snapshot(data = {'time': self.times[comp]})
#                                    
#        outputFile.close()
#
#        t.close()
#        
#        self.toPlot = None        
                
#####################################################################################
# GlobalMotionTrajectoryFilter ANALYSIS
#####################################################################################
class GlobalMotionFilteredTrajectory(Analysis):
    """Sets up a Global Motion Trajectory Filter analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: GlobalMotionFilteredTrajectory(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo   -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                            number to consider, 'last' is an integer specifying the last frame number to consider and 
                            'step' is an integer specifying the step number between two frames.
            * subset     -- a selection string specifying the atoms used to perform the global motion filtering.
            * target     -- a selection string specifying the atoms to include in the output trajectory.
            * pyroserver -- a string specifying if Pyro will be used and how to run the analysis.
            * output     -- the output NetCDF file name.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'GMFT')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
        
        self.firstFrame = self.frameIndexes[0]
                    
        # This dictionnary works a little bit like a Particle Tensor object. Its keys are the index of the atoms of
        # the selection. Its value are the filtered atomic trajectories.
        self.filteredTrajectory = {}
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        self.buildTimeInfo()
                
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)
                
        self.target = self.selectAtoms(self.targetDefinition)
                        
        self.trajectory.universe.setFromTrajectory(self.trajectory, self.firstFrame)

        orderedAtoms = sorted(self.trajectory.universe.atomList(), key = operator.attrgetter('index'))
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])
                
        # For the fist frame, defines principal axes transformation.
        tr = selectedAtoms.normalizingTransformation()
        
        # And apply it to the system.
        self.trajectory.universe.applyTransformation(tr)
#        selectedAtoms.applyTransformation(tr)
        
        self.initialConfArray = copy.deepcopy(self.trajectory.universe.configuration().array)
        
        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')
        
    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """
        
        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])
        
        targetAtoms = Collection([orderedAtoms[ind] for ind in self.target])
        
        initialConf = Configuration(trajectory.universe, self.initialConfArray)
        
        # First frame, nothing to do because the initialConf already stores the initial transformation.
        if frameIndex == self.firstFrame:
            trajectory.universe.setConfiguration(initialConf)
            
        else:
            trajectory.universe.setFromTrajectory(trajectory, frameIndex)                    
            # Find and apply the linear transformation that will minimize the RMS with the configuration
            # resulting from the principal axes transformation.
            tr, rms = selectedAtoms.findTransformation(initialConf)
            #selectedAtoms.applyTransformation(tr)
            trajectory.universe.applyTransformation(tr)
        
        globalMotionFilteredFrame = {}
        for at in targetAtoms:
            globalMotionFilteredFrame[at.index] = at.position()
                                                    
        return frameIndex, globalMotionFilteredFrame
    
    def combine(self, frameIndex, x):
        """
        @param frameIndex: the index of the frame.
        @type frameIndex: integer.

        @param x: 
        @type x: 
        """
                    
        self.filteredTrajectory[frameIndex] = x
                
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
        
        if self.architecture == 'monoprocessor':
            t = self.trajectory            
        else:
            # Load the whole trajectory set.
            t = Trajectory(None, self.trajectoryFilename, 'r')

        orderedAtoms = sorted(t.universe.atomList(), key = operator.attrgetter('index'))
        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])

        targetAtoms = Collection([orderedAtoms[ind] for ind in self.target])
        
        # traj_new is the filtered trajectory
        outputFile = Trajectory(targetAtoms, self.output, "w") 
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Create the snapshot generator
        snapshot = SnapshotGenerator(t.universe, actions = [TrajectoryOutput(outputFile, ["configuration"], 0, None, 1)])

        # Loop over the output frames.
        for comp in range(self.nFrames):
            
            frameIndex = self.frameIndexes[comp]
            t.universe.setFromTrajectory(t, frameIndex)
            
            for at in targetAtoms:
                at.setPosition(Vector(self.filteredTrajectory[frameIndex][at.index]))
                
            snapshot(data = {'time': self.times[comp]})

        outputFile.close()

        t.close()
        
        self.toPlot = None

#####################################################################################
# CENTER OF MASS TRAJECTORY ANALYSIS
#####################################################################################
class CenterOfMassTrajectory(Analysis):
    """Sets up a Center Of Mass Trajectory analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: CenterOfMassTrajectory(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo   -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                            number to consider, 'last' is an integer specifying the last frame number to consider and 
                            'step' is an integer specifying the step number between two frames.
            * group      -- a selection string specifying the groups of atoms on which the center of mass will be defined
                            (one center of mass per group).
            * pyroserver -- a string specifying if Pyro will be used and how to run the analysis.
            * output     -- the output NetCDF file name.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'COMT')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
                                                            
        self.comsTrajectory = {}

    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        self.buildTimeInfo()

        self.group = self.selectGroups(self.groupDefinition)
        self.nGroups = len(self.group)

        if self.architecture != 'monoprocessor':
            # The attribute trajectory is removed because it can not be pickled by Pyro.
            delattr(self, 'trajectory')

    def calc(self, frameIndex, trajectory):
        """Calculates the contribution for one frame.
        
        @param frameIndex: the index of the frame.
        @type frameIndex: integer.
        
        @param trajectory: the trajectory.
        @type trajectory: MMTK.Trajectory.Trajectory object
        """

        trajectory.universe.setFromTrajectory(trajectory, frameIndex)

        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
        groups = [Collection([orderedAtoms[ind] for ind in g]) for g in self.group]
                                        
        centersOfMass = [g.centerOfMass() for g in groups]
                
        return frameIndex, centersOfMass
    
    def combine(self, frameIndex, x):
        """
        @param frameIndex: the index of the frame.
        @type frameIndex: integer
        """

        self.comsTrajectory[frameIndex] = x
                
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """     

        if self.architecture == 'monoprocessor':
            t = self.trajectory            
        else:
            # Load the whole trajectory set.
            t = Trajectory(None, self.trajectoryFilename, 'r')
            
        comsUniverse = t.universe.__copy__()
        
        comsUniverse.removeObject(comsUniverse.objectList()[:])

        orderedAtoms = sorted(t.universe.atomList(), key = operator.attrgetter('index'))
        groups = [Collection([orderedAtoms[ind] for ind in g]) for g in self.group]

        comp = 1
        for g in groups:
            
            comAtom = Atom('H', name = 'COM'+str(comp))
            comAtom._mass = g.mass()
            comsUniverse.addObject(comAtom)
            comp += 1
                                                            
        # traj_new is the filtered trajectory
        outputFile = Trajectory(comsUniverse, self.output, "w") 
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()
        
        # Each time |snapshot| is called, the universe contents i flushed into the output file.
        snapshot = SnapshotGenerator(comsUniverse,\
                                     actions = [TrajectoryOutput(outputFile, ["configuration","time"], 0, None, 1)])

        # Loop over the output frames.
        for comp in range(self.nFrames):
            
            frameIndex = self.frameIndexes[comp]
            
            t.universe.setFromTrajectory(t, frameIndex)
            
            comsUniverse.setCellParameters(t.universe.cellParameters())
            
            aComp = 0
            for at in comsUniverse.atomList():
                
                at.setPosition(self.comsTrajectory[frameIndex][aComp])
                aComp += 1
                
            snapshot(data = {'time': self.times[comp]})
        
        # The output COM trajectory is closed.
        outputFile.close()
        
        self.toPlot = None
                
#####################################################################################
# QUASI-HARMONIC ANALYSIS ANALYSIS
#####################################################################################
class QuasiHarmonicAnalysis(Analysis):
    """Sets up a Quasi Harmonic Analysis analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: QuasiHarmonicAnalysis(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory  -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo    -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                             number to consider, 'last' is an integer specifying the last frame number to consider and 
                             'step' is an integer specifying the step number between two frames.
            * temperature -- the temperature at which the MD was performed.
            * subset      -- a selection string specifying the atoms to consider for the analysis.
            * qha         -- the output NetCDF file name.
        
    Running modes:
    
        - To run the analysis do: a.runAnalysis() where a is the analysis object.
        - To estimate the analysis do: a.estimateAnalysis() where a is the analysis object.
        - To save the analysis to 'file' file name do: a.saveAnalysis(file) where a is the analysis object.
        
    Comments:
    
        - This analysis is used to get effective modes of vibration from fluctuations calculated by an MD simulation. 
          The results of such an analysis can be seen by generating pseudo-trajectories reproducing the vibrations along
          a vibration mode.
        - For more details: Brooks et al., J. Comp. Chem. 1995, 16, 1522-1542.
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'QHA')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        self.interpreteInputParameters()
        
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        # The 'temperature' input parameter. It must be a float or a value convertible to a float > 0.                
        # Must be a float or convertible to a float.
        try:
            self.temperature = float(self.parameters['temperature'])

        # The value can not be converted to a float.
        except ValueError:
            raise Error('Error when parsing "temperature" parameter: must be a float or a value convertible to a float.')
                
        else:
            # Must > 0.
            if self.temperature <= 0.0:
                raise Error('Error when parsing "temperature" parameter: must be > 0.')
            
        self.buildTimeInfo()
        
        self.subset = self.selectAtoms(self.subsetDefinition)
        self.nSelectedAtoms = len(self.subset)

        self.mask = N.zeros(self.universe.numberOfAtoms(), typecode = N.Int0)
        for aIndex in self.subset:
            self.mask[aIndex] = 1
                            
    def internalRun(self):
        """Runs the analysis."""
        
        self.chrono = default_timer()

        orderedAtoms = sorted(self.universe.atomList(), key = operator.attrgetter('index'))

        selectedAtoms = Collection([orderedAtoms[ind] for ind in self.subset])
                
        M1_2 = N.zeros((3*self.nSelectedAtoms,), N.Float)
        
        weightList = [N.sqrt(el[1]) for el in sorted([(at.index, at.mass()) for at in selectedAtoms])]
        
        for i in range(len(weightList)):
            M1_2[3*i:3*(i+1)] = weightList[i]
            
        invM1_2 = (1.0/M1_2)*N.identity(3*self.nSelectedAtoms, typecode = N.Float)
            
        # The initial structure configuration.
        initialStructure = self.trajectory.configuration[self.first]

        averageStructure = ParticleVector(self.universe)

        for conf in self.trajectory.configuration[self.first:self.last:self.skip]:
            averageStructure += self.universe.configurationDifference(initialStructure, conf)

        averageStructure = averageStructure/self.nFrames
        averageStructure = initialStructure + averageStructure        
        
        mdr  = N.zeros((self.nFrames, 3*self.nSelectedAtoms), N.Float)
        
        # Calculate the fluctuation matrix.
        sigmaPrim = N.zeros((3*self.nSelectedAtoms, 3*self.nSelectedAtoms), N.Float)
        comp = 0
        for conf in self.trajectory.configuration[self.first:self.last:self.skip]:
            mdr[comp,:] = M1_2*N.ravel(N.compress(self.mask,(conf-averageStructure).array,0)) 
            sigmaPrim += mdr[comp,:, N.NewAxis] * mdr[N.NewAxis, comp, :]            
            comp += 1
                        
        sigmaPrim = sigmaPrim/float(self.nFrames)

        try:
            # Calculate the quasiharmonic modes
            omega, dx = Heigenvectors(sigmaPrim)
            
        except MemoryError:
            raise Error('Not enough memory to diagonalize the %sx%s fluctuation matrix.' % sigmaPrim.shape)
            
        # Due to numerical imprecisions, the result can have imaginary parts.
        # In that case, throw the imaginary parts away.
        # Conversion from uma*nm2 to kg*m2 (SI)
        omega = omega.real/(Units.kg*Units.m**2)        
        omega = (Units.Hz/Units.invcm)*N.sqrt((Units.k_B*Units.K*self.temperature/Units.J)*(1.0/omega))

        dx = dx.real
        dx = N.dot(invM1_2, dx)

        # Sort eigen vectors by decreasing fluctuation amplitude.
        indices = N.argsort(omega)[::-1]
        
        omega = N.take(omega, indices)
        dx = N.take(dx, indices)

        # Eq 66 of the reference paper.
        mdr = N.take(mdr, indices, axis = 1)
        at = N.dot(mdr,N.transpose(dx))
        
        # The NetCDF output file is opened.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # The universe is emptied from its objects keeping just its topology.
        self.universe.removeObject(self.universe.objectList()[:])
        # The atoms of the subset are copied
        atoms = copy.deepcopy(selectedAtoms.atomList())
        
        # And their parent attribute removed to allow their transfer in the empty universe.
        for a in atoms:
            a.parent = None
        ac = AtomCluster(atoms,name='QHACluster')
        self.universe.addObject(ac)

        # Some dimensions are created.
        # NEIVALS = the number eigen values
        outputFile.createDimension('NEIGENVALS', len(omega))

        # UDESCR = the universe description length
        outputFile.createDimension('UDESCR', len(self.universe.description()))

        # NATOMS = the number of atoms of the universe
        outputFile.createDimension('NATOMS', self.nSelectedAtoms)

        # NFRAMES = the number of frames.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # NCOORDS = the number of coordinates (always = 3).
        outputFile.createDimension('NCOORDS', 3)

        # 3N.
        outputFile.createDimension('3N', 3*self.nSelectedAtoms)

        if self.universe.cellParameters() is not None:
            outputFile.createDimension('BOXDIM', len(self.universe.cellParameters()))

        # Creation of the NetCDF output variables.
        # EIVALS = the eigen values.
        OMEGA = outputFile.createVariable('omega', N.Float, ('NEIGENVALS',))
        OMEGA[:] = omega

        # EIVECS = array of eigen vectors.
        DX = outputFile.createVariable('dx', N.Float, ('NEIGENVALS','NEIGENVALS'))
        DX[:,:] = dx

        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times[:]
        TIMES.units = 'ps'

        # MODE = the mode number.
        MODE = outputFile.createVariable('mode', N.Float, ('3N',))
        MODE[:] = 1 + N.arange(3*self.nSelectedAtoms)
                
        # LCI = local character indicator. See eq 56.
        LCI = outputFile.createVariable('local_character_indicator', N.Float, ('3N',))
        LCI[:] = (dx**4).sum(0)
        
        # GCI = global character indicator. See eq 57.
        GCI = outputFile.createVariable('global_character_indicator', N.Float, ('3N',))
        GCI[:] = (N.sqrt(3.0*self.nSelectedAtoms)/(N.absolute(dx)).sum(0))**4

        # Projection of MD traj onto normal modes. See eq 66..
        AT = outputFile.createVariable('at', N.Float, ('NFRAMES','3N'))
        AT[:,:] = at[:,:]
        
        # DESCRIPTION = the universe description.
        DESCRIPTION = outputFile.createVariable('description', N.Character, ('UDESCR',))
        DESCRIPTION[:] = self.universe.description()

        # AVGSTRUCT = the average structure.
        AVGSTRUCT = outputFile.createVariable('avgstruct', N.Float, ('NATOMS','NCOORDS'))
        AVGSTRUCT[:,:] = N.compress(self.mask,averageStructure.array,0)

        # If the universe is periodic, create an extra variable storing the box size.
        if self.universe.cellParameters() is not None:
            BOXSIZE = outputFile.createVariable('box_size', N.Float, ('BOXDIM',))
            BOXSIZE[:] = self.universe.cellParameters()

        outputFile.close()

        self.toPlot = None

        self.chrono = default_timer() - self.chrono
        
        return None
