"""Collections of classes for the determination of dynamics-related properties.

Classes:

    * RigidBodyTrajectory                    : sets up a Rigid-Body Trajectory analysis.
    * ReorientationalCorrelationFunction     : sets up an Reorientational Correlation Function analysis.
    * AngularVelocityAutoCorrelationFunction : sets up an Angular Velocity AutoCorrelation Function analysis.
    * AngularDensityOfStates                 : sets up an Angular Density Of States analysis.
"""

# The python distribution modules
import copy
import operator
import os
from time import asctime

# The ScientificPython modules
from Scientific import N 
from Scientific.Geometry import Vector
from Scientific.Geometry.Quaternion import Quaternion
from Scientific.Geometry.Transformation import Translation
from Scientific.IO.NetCDF import NetCDFFile

# The MMTK distribution modules
from MMTK.Collections import Collection
from MMTK.Trajectory import SnapshotGenerator, Trajectory, TrajectoryOutput

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Analysis.Analysis import Analysis
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import convertNetCDFToASCII
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Mathematics.Analysis import correlation, differentiate, factorial, FFT, gaussian, gaussianWindow
from nMOLDYN.Mathematics.Geometry import getAngularVelocity

#####################################################################################
# RIGID BODY TRAJECTORY ANALYSIS
#####################################################################################
class RigidBodyTrajectory(Analysis):
    """Sets up a Rigid Body Trajectory analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: RigidBodyTrajectory(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory     -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo       -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                number to consider, 'last' is an integer specifying the last frame number to consider and 
                                'step' is an integer specifying the step number between two frames.
            * referenceframe -- an integer in [1,len(trajectory)] specifying which frame should be the reference.
            * stepwiserbt    -- a string being one of 'Yes' or 'No' specifying whether the reference frame for frame i should be 
                                the frame i - 1 ('Yes') or should be a fixed frame defined with |referenceframe| ('No').
            * group          -- a selection string specifying the groups of atoms on which the rigid body trajectory will be defined.
                                (each group being a rigid body).
            * pyroserver     -- a string specifying if Pyro will be used and how to run the analysis.
            * output         -- the output NetCDF file name.
        
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'RBT')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()

        self.firstFrame = self.frameIndexes[0]
        
        self.RBT = {'trajectory' : {}}
                                            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()

        if (self.referenceFrame < 0) or (self.referenceFrame > len(self.trajectory)):
            self.referenceFrame = 0
            LogMessage('warning', 'The reference frame must be an integer in [1,%s].\n\
            It will be set to 1 for the running analysis.' % len(self.trajectory), ['console'])

        self.group = self.selectGroups(self.groupDefinition)
        self.nGroups = len(self.group)
        self.groupsSize = [len(g) for g in self.group]

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

        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
        group = Collection([orderedAtoms[ind] for ind in atomIndexes])

        rbtPerGroup = {}        
        # Those matrix will store the quaternions and the CMS coming from the RBT trajectory.
        rbtPerGroup['quaternions'] = N.zeros((self.nFrames, 4), typecode = N.Float)
        rbtPerGroup['com'] = N.zeros((self.nFrames, 3), typecode = N.Float)
        rbtPerGroup['fit'] = N.zeros((self.nFrames,), typecode = N.Float)
        rbtPerGroup['trajectory'] = {}
        
        # Case of a moving reference.
        if self.stepwiseRBT:
            
            # The reference configuration is always the one of the previous frame excepted for the first frame
            # where it is set by definition to the first frame (could we think about a cyclic alternative way ?).
            for comp in range(self.nFrames):
                
                frameIndex = self.frameIndexes[comp]
                
                if comp == 0:
                    previousFrame = self.firstFrame
                    
                else:
                    previousFrame = self.frameIndexes[comp-1]
                    
                refConfig = trajectory.configuration[previousFrame]

                # The RBT is created just for the current step.
                rbt = trajectory.readRigidBodyTrajectory(group,\
                                                         first = frameIndex,\
                                                         last = frameIndex + 1,\
                                                         skip = 1,\
                                                         reference = refConfig)

                # The corresponding quaternions and cms are stored in their corresponding matrix.
                rbtPerGroup['quaternions'][comp,:] = copy.copy(rbt.quaternions)
                rbtPerGroup['com'][comp,:] = copy.copy(rbt.cms)
                rbtPerGroup['fit'][comp] = copy.copy(rbt.fit)
                                            
        # The simplest case, the reference frame is fixed.
        # A unique RBT is performed from first to last skipping skip steps and using refConfig as the reference.
        else:

            # If a fixed reference has been set. We can already set the reference configuration here.
            refConfig = trajectory.configuration[self.referenceFrame]
            
            # The RBT is created.
            rbt = trajectory.readRigidBodyTrajectory(group,\
                                                     first = self.first,\
                                                     last = self.last,\
                                                     skip = self.skip,\
                                                     reference = refConfig)
            
            # The corresponding quaternions and cms are stored in their corresponding matrix.
            rbtPerGroup['quaternions'] = copy.copy(rbt.quaternions)
            rbtPerGroup['com'] = copy.copy(rbt.cms)
            rbtPerGroup['fit'] = copy.copy(rbt.fit)

        # I can not use the centers of mass defined by rbt.cms because the reference frame
        # selected can be out of the selected frames for the Rigid Body Trajectory.
        centerOfMass = group.centerOfMass(refConfig)

        # Loop over the atoms of the group to set the RBT trajectory.
        for atom in group:
            
            rbtPerGroup['trajectory'][atom.index] = N.zeros((self.nFrames,3), typecode = N.Float)
            
            # The coordinates of the atoms are centered around the center of mass of the group.
            xyz = refConfig[atom] - centerOfMass

            # Loop over the selected frames.
            for comp in range(self.nFrames):
                
                # The rotation matrix corresponding to the selected frame in the RBT.
                transfo = Quaternion(rbtPerGroup['quaternions'][comp,:]).asRotation()

                if self.removeTranslation:
                    # The transformation matrix corresponding to the selected frame in the RBT.
                    transfo = Translation(centerOfMass)*transfo
                        
                # Compose with the CMS translation if the removeTranslation flag is set off.
                else:
                    # The transformation matrix corresponding to the selected frame in the RBT.
                    transfo = Translation(Vector(rbtPerGroup['com'][comp,:]))*transfo

                # The RBT is performed on the CMS centered coordinates of atom at.
                rbtPerGroup['trajectory'][atom.index][comp,:] = transfo(Vector(xyz))
                                        
        return atomIndexes, rbtPerGroup
    
    def combine(self, atomIndexes, x):
        """
        """
        
        comp = self.group.index(atomIndexes)
        self.RBT[comp] = {'quaternions' : x['quaternions'], 'com' : x['com'], 'fit' : x['fit']}
        self.RBT['trajectory'].update(x['trajectory'])
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
                
        if self.architecture == 'monoprocessor':
            t = self.trajectory            
        else:
            # Load the whole trajectory set.
            t = Trajectory(None, self.trajectoryFilename, 'r')
                    
        selectedAtoms = Collection()
        orderedAtoms = sorted(t.universe.atomList(), key = operator.attrgetter('index'))
        groups = [[selectedAtoms.addObject(orderedAtoms[index]) for index in atomIndexes] for atomIndexes in self.group]

        # Create trajectory
        outputFile = Trajectory(selectedAtoms, self.output, 'w')

        # Create the snapshot generator
        snapshot = SnapshotGenerator(t.universe, actions = [TrajectoryOutput(outputFile, ["configuration","time"], 0, None, 1)])
                
        # The output is written
        for comp in range(self.nFrames):

            frameIndex = self.frameIndexes[comp]
            t.universe.setFromTrajectory(t, frameIndex)
            
            for atom in selectedAtoms:
                atom.setPosition(self.RBT['trajectory'][atom.index][comp,:])
            snapshot(data = {'time' : self.times[comp]})
  
        outputFile.close()
        
        outputFile = NetCDFFile(self.output, 'a')

        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        outputFile.jobinfo += 'Input trajectory: %s\n\n' % self.trajectoryFilename
                                                      
        outputFile.createDimension('NFRAMES', self.nFrames)
        outputFile.createDimension('NGROUPS', self.nGroups)
        outputFile.createDimension('QUATERNIONLENGTH',4)

        # The NetCDF variable that stores the quaternions.
        QUATERNIONS = outputFile.createVariable('quaternion', N.Float, ('NGROUPS', 'NFRAMES','QUATERNIONLENGTH'))

        # The NetCDF variable that stores the centers of mass.
        COM = outputFile.createVariable('com', N.Float, ('NGROUPS','NFRAMES','xyz'))
            
        # The NetCDF variable that stores the rigid-body fit.
        FIT = outputFile.createVariable('fit', N.Float, ('NGROUPS','NFRAMES'))

        # Loop over the groups.
        for comp in range(self.nGroups):
            
            aIndexes = self.group[comp]
            
            outputFile.jobinfo += 'Group %s: %s\n' % (comp+1, [index for index in aIndexes])

            QUATERNIONS[comp,:,:] = self.RBT[comp]['quaternions'][:,:]
            COM[comp,:,:] = self.RBT[comp]['com'][:,:]
            FIT[comp,:] = self.RBT[comp]['fit'][:]
                        
        outputFile.close()
        
        self.toPlot = None        

#####################################################################################
# REORIENTATIONAL CORRELATION FUNCTION ANALYSIS
#####################################################################################
class ReorientationalCorrelationFunction(Analysis):
    """Sets up a Reorientational Correlation Function analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: ReorientationalCorrelationFunction(|parameters| = None)
    
    Arguments:
    
        - |parameters| -- a dictionnary of the input parameters, or 'None' to set up the analysis without parameters.
            * trajectory     -- a trajectory file name or an instance of MMTK.Trajectory.Trajectory class.
            * timeinfo       -- a string of the form 'first:last:step' where 'first' is an integer specifying the first frame 
                                number to consider, 'last' is an integer specifying the last frame number to consider and 
                                'step' is an integer specifying the step number between two frames.
            * referenceframe -- an integer in [1,len(trajectory)] specifying which frame should be the reference.
            * stepwiserbt    -- a string being one of 'Yes' or 'No' specifying whether the reference frame for frame i should be 
                                the frame i - 1 ('Yes') or should be a fixed frame defined with |referenceframe| ('No').
            * wignerindexes  -- a string of the form 'j,m,n' where 'j', 'm' and 'n' are respectively the j, m and n indexes of the
                                Wigner function Djmn.
            * group          -- a selection string specifying the groups of atoms on which the rigid body trajectory will be defined.
                                (each group being a rigid body).
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'RCF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
            
        self.firstFrame = self.frameIndexes[0]
                    
        self.RCF = {}
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)

        self.buildTimeInfo()

        if (self.referenceFrame < 0) or (self.referenceFrame > len(self.trajectory)):
            self.referenceFrame = 0
            LogMessage('warning', 'The reference frame must be an integer in [1,%s].\n\
            It will be set to 1 for the running analysis.' % len(self.trajectory), ['console'])

        # The 'wignerindexes' input parameter. It must be  string of the form j,m,n where j, m and n are the Wigner 
        # indexes.            
        if isinstance(self.parameters['wignerindexes'], str):
            # Must be a string of three comma-separated integers.
            try:
                j, m, n = [int(v) for v in self.parameters['wignerindexes'].split(',')]

            except:
                raise Error('Error when parsing "wignerindexes" parameter: must be a string of the form "j,m,n".')
                        
            else:
                if j < 0:
                    raise Error('Error when parsing "wignerindexes" parameter: j must be >= 0.')
                        
                if m > j:
                    raise Error('Error when parsing "wignerindexes" parameter: m must be <= j.')
                        
                if abs(n) > m:
                    raise Error('Error when parsing "wignerindexes" parameter: n must be <= n.')
                    
                self.wignerIndexes = (j,m,n)
        
        self.group = self.selectGroups(self.groupDefinition)
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
        
        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
        group = Collection([orderedAtoms[ind] for ind in atomIndexes])

        j, m, n = self.wignerIndexes

        # Those matrix will store the quaternions and the CMS coming from the RBT trajectory.
        quaternions = N.zeros((self.nFrames, 4), typecode = N.Float)
                
        # Case of a moving reference.
        if self.stepwiseRBT:
            
            # The reference configuration is always the one of the previous frame excepted for the first frame
            # where it is set by definition to the first frame (could we think about a cyclic alternative way ?).
            for comp in range(self.nFrames):
                
                frameIndex = self.frameIndexes[comp]

                if comp == 0:
                    previousFrame = self.firstFrame
                    
                else:
                    previousFrame = self.frameIndexes[comp-1]
                    
                refConfig = trajectory.configuration[previousFrame]

                # The RBT is created just for the current step.
                rbt = trajectory.readRigidBodyTrajectory(group,\
                                                         first = frameIndex,\
                                                         last = frameIndex + 1,\
                                                         skip = 1,\
                                                         reference = refConfig)
                
                # The corresponding quaternions and cms are stored in their corresponding matrix.
                quaternions[comp,:] = copy.copy(rbt.quaternions)
                
        # The simplest case, the reference frame is fixed.
        # A unique RBT is performed from first to last skipping skip steps and using refConfig as the reference.
        else:
            # If a fixed reference has been set. We can already set the reference configuration here.
            refConfig = trajectory.configuration[self.referenceFrame]
            
            # The RBT is created.
            rbt = trajectory.readRigidBodyTrajectory(group,\
                                                     first = self.first,\
                                                     last = self.last,\
                                                     skip = self.skip,\
                                                     reference = refConfig)

            quaternions = rbt.quaternions
            
        # c1 is the scaling factor converting Wigner function into spherical harmonics. It depends only on j.
        c1  = N.sqrt(((2.0*j+1)/(4.0*N.pi)))

        quat2 = N.zeros(quaternions.shape, typecode = N.Complex)  
        # quat2[:,0] refers to the (q0+iq3) of equation 3.55
        quat2[:,0] = quaternions[:,0] + 1j*quaternions[:,3]
        # quat2[:,2] refers to the (q2+iq1) of equation 3.55
        quat2[:,2] = quaternions[:,2] + 1j*quaternions[:, 1]
        # quat2[:,1] refers to the (q0-iq3) of equation 3.55
        quat2[:,1] = N.conjugate(quat2[:,0])
        # quat2[:,3] refers to the (q2-iq1) of equation 3.55
        quat2[:,3] = N.conjugate(quat2[:,2])

        pp  = self.preparePP(j, m, n)
        Djmn = N.add.reduce(N.multiply.reduce(quat2[:,N.NewAxis,:]**pp[0][N.NewAxis,:,:],-1)*pp[1],1)
        
        if m == n:
            Djnm = Djmn
        else:
            pp = self.preparePP(j, n, m)
            Djnm = N.add.reduce(N.multiply.reduce(quat2[:,N.NewAxis,:]**pp[0][N.NewAxis,:,:],-1)*pp[1],1)

        Djmn = Djmn*c1
        Djnm = Djnm*c1
                        
        return atomIndexes, (Djmn, Djnm)
    
    def combine(self, atomIndexes, x):
        """
        """
        
        comp = self.group.index(atomIndexes) + 1

        self.RCF[comp] = correlation(x[0],x[1])
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
        
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Dictionnary whose keys are of the form Gi where i is the group number 
        # and the entries are the list of the index of the atoms building the group.
        comp = 1
        for g in self.group:
            outputFile.jobinfo += 'Group %s: %s\n' % (comp,[index for index in g])
            comp += 1

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times
        TIMES.units = 'ps'

        rcfTotal = N.zeros((self.nFrames), typecode = N.Float)
        
        for k in self.RCF.keys():

            RCF = outputFile.createVariable('rcf-group%s' % k, N.Float, ('NFRAMES',))
            RCF[:] = self.RCF[k][:]
            RCF.units = 'unitless'

            N.add(rcfTotal, self.RCF[k], rcfTotal)
                        
        rcfTotal *= 4.0*N.pi/self.nGroups
                
        # The RCF.
        RCF = outputFile.createVariable('rcf-total', N.Float, ('NFRAMES',))
        RCF[:] = rcfTotal[:]
        RCF.units = 'unitless'

        asciiVar = sorted(outputFile.variables.keys())
        
        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'rcf-total'}

        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
        
    def preparePP(self, j, m, n):
        """Intermediate function used to setup the calculation of spherical harmonics."""

        # c2 is the sqrt[(j+m)! * (j-m)!] * j! of equation 3.55
        c2 = N.sqrt((factorial(j+m) * factorial(j-m))) * factorial(j) 
        
        c3 = j+m
        pp = []
        aa = []
        for p in range(0,j+1,1):
            if (p <= c3) and (p >= m):
                # p1 is the j+m-p of equation 3.55
                p1 = c3-p
                # p2 is the j-p of equation 3.55
                p2 = j-p
                p3 = p-m

                # a1 is the first part of equation 3.55 before the x
                a1 = ((-1)**p)*c2/(factorial(p1)*factorial(p2)*factorial(p)*factorial(p3))

                pp.append((p1,p2,p3,p))

                aa.append(a1)

        return N.array(pp, typecode = N.Float), N.array(aa, typecode = N.Float)
    
#####################################################################################
# ANGULAR VELOCITY AUTOCORRELATION FUNCTION ANALYSIS
#####################################################################################
class AngularVelocityAutoCorrelationFunction(Analysis):
    """Sets up an Angular Velocity AutoCorrelation Function analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: AngularVelocityAutoCorrelationFunction(|parameters| = None)
    
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
            * referenceframe  -- an integer in [1,len(trajectory)] specifying which frame should be the reference.
            * stepwiserbt     -- a string being one of 'Yes' or 'No' specifying whether the reference frame for frame i should be 
                                 the frame i - 1 ('Yes') or should be a fixed frame defined with |referenceframe| ('No').
            * group           -- a selection string specifying the groups of atoms on which the rigid body trajectory will be defined.
                                 (each group being a rigid body).
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'AVACF')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
            
        self.AVACF = {}
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()

        if (self.referenceFrame < 0) or (self.referenceFrame > len(self.trajectory)):
            self.referenceFrame = 0
            LogMessage('warning', 'The reference frame must be an integer in [1,%s].\n\
            It will be set to 1 for the running analysis.' % len(self.trajectory), ['console'])

        self.group = self.selectGroups(self.groupDefinition)
        self.nGroups = len(self.group)
        self.groupsSize = [len(g) for g in self.group]

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

        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
        group = Collection([orderedAtoms[ind] for ind in atomIndexes])
        
        angvel = getAngularVelocity(trajectory,\
                                    group,\
                                    self.frameIndexes,\
                                    self.dt,\
                                    self.stepwiseRBT,\
                                    self.referenceFrame,\
                                    self.differentiation)

        return atomIndexes, angvel
    
    def combine(self, atomIndexes, x):
        """
        """
        
        comp = self.group.index(atomIndexes) + 1

        # the AVACF for group g
        if self.projection is None:
            self.AVACF[comp] = correlation(x)/3.0
            
        else:
            projectedAngVel = N.dot(x, self.projection)
            self.AVACF[comp] = correlation(projectedAngVel)

        self.AVACF[comp] /= self.AVACF[comp][0]
        
    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
        
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()

        # Dictionnary whose keys are of the form Gi where i is the group number 
        # and the entries are the list of the index of the atoms building the group.
        comp = 1
        for g in self.group:
            outputFile.jobinfo += 'Group %d: %s\n' % (comp,[index for index in g])
            comp += 1

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)

        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times[:]
        TIMES.units = 'ps'

        avacfTotal = N.zeros((self.nFrames), typecode = N.Float)
        
        for k in self.AVACF.keys():

            AVACF = outputFile.createVariable('avacf-group%s' % k, N.Float, ('NFRAMES',))
            AVACF[:] = self.AVACF[k][:]
            AVACF.units = 'rad^2*ps^-2'

            N.add(avacfTotal, self.AVACF[k], avacfTotal)
            
        avacfTotal /= self.nGroups
                        
        AVACF = outputFile.createVariable('avacf-total', N.Float, ('NFRAMES',))
        AVACF[:] = avacfTotal[:]
        AVACF.units = 'rad^2*ps^-2'

        asciiVar = sorted(outputFile.variables.keys())
        
        outputFile.close()

        self.toPlot = {'netcdf' : self.output, 'xVar' : 'time', 'yVar' : 'avacf-total'}        

        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
        
#####################################################################################
# ANGULAR DENSITY OF STATES ANALYSIS
#####################################################################################
class AngularDensityOfStates(Analysis):
    """Sets up an Angular Density Of States analysis.
    
    A Subclass of nMOLDYN.Analysis.Analysis. 
    
    Constructor: AngularDensityOfStates(|parameters| = None)
    
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
            * referenceframe  -- an integer in [1,len(trajectory)] specifying which frame should be the reference.
            * stepwiserbt     -- a string being one of 'Yes' or 'No' specifying whether the reference frame for frame i should be 
                                 the frame i - 1 ('Yes') or should be a fixed frame defined with |referenceframe| ('No').
            * resolution      -- a float specifying the width of the gaussian, that will be used to mimics the experimental resolution.
            * group           -- a selection string specifying the groups of atoms on which the rigid body trajectory will be defined.
                                 (each group being a rigid body).
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

        path = os.path.join(GVAR['nmoldyn_analysis'],'ADOS')
        exec file(path) in None, self.__dict__

    def initialize(self):
        """Initializes the analysis (e.g. parses and checks input parameters, set some variables ...).
        """
        
        # The input parameters are parsed.
        self.interpreteInputParameters()
            
        self.AVACF = {}
        self.ADOS  = {}
            
    def interpreteInputParameters(self):
        """Parse the input parameters for the analysis.
        """

        # Parses the parameters that are common to different analysis.
        Analysis.interpreteInputParameters(self)
        
        self.buildTimeInfo()
        
        if (self.referenceFrame < 0) or (self.referenceFrame > len(self.trajectory)):
            self.referenceFrame = 0
            LogMessage('warning', 'The reference frame must be an integer in [1,%s].\n\
            It will be set to 1 for the running analysis.' % len(self.trajectory), ['console'])

        self.group = self.selectGroups(self.groupDefinition)
        self.nGroups = len(self.group)
        self.groupsSize = [len(g) for g in self.group]

        # The resolution is converted from ueV to rad/ns (1.519).
        self.resolutionFunction = gaussian(self.times, self.timeSigma)
        
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

        orderedAtoms = sorted(trajectory.universe.atomList(), key = operator.attrgetter('index'))
        group = Collection([orderedAtoms[ind] for ind in atomIndexes])
        
        angvel = getAngularVelocity(trajectory,\
                                    group,\
                                    self.frameIndexes,\
                                    self.dt,\
                                    self.stepwiseRBT,\
                                    self.referenceFrame,\
                                    self.differentiation)
                            
        return atomIndexes, angvel
    
    def combine(self, atomIndexes, x):
        """
        """
        
        comp = self.group.index(atomIndexes) + 1
                
        # the AVACF for group g
        if self.projection is None:
            self.AVACF[comp] = correlation(x)/3.0
            
        else:
            projectedAngVel = N.dot(x, self.projection)
            self.AVACF[comp] = correlation(projectedAngVel)
                
        self.ADOS[comp] = FFT(gaussianWindow(self.AVACF[comp], self.resolutionFunction)).real[:self.nFrames]

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """                                    
        
        if self.architecture == 'monoprocessor':
            t = self.trajectory            
        else:
            # Load the whole trajectory set.
            t = Trajectory(None, self.trajectoryFilename, 'r')
            
        orderedAtoms = sorted(t.universe.atomList(), key = operator.attrgetter('index'))
        groups = [Collection([orderedAtoms[ind] for ind in g]) for g in self.group]
        
        # 'freqencies' = 1D Numeric array. Frequencies at which the DOS was computed
        frequencies = N.arange(self.nFrames)/(2.0*self.nFrames*self.dt)
        
        # The NetCDF output file is opened for writing.
        outputFile       = NetCDFFile(self.output, 'w')
        outputFile.title = self.__class__.__name__
        outputFile.jobinfo = self.information + '\nOutput file written on: %s\n\n' % asctime()
        
        # Dictionnary whose keys are of the form Gi where i is the group number 
        # and the entries are the list of the index of the atoms building the group.
        comp = 1
        for g in self.group:
            outputFile.jobinfo += 'Group %d: %s\n' % (comp,[index for index in g])
            comp += 1

        # Some dimensions are created.
        outputFile.createDimension('NFRAMES', self.nFrames)
        
        # Creation of the NetCDF output variables.
        # The time.
        TIMES = outputFile.createVariable('time', N.Float, ('NFRAMES',))
        TIMES[:] = self.times[:]
        TIMES.units = 'ps'

        # The resolution function.
        RESOLUTIONFUNCTION = outputFile.createVariable('resolution_function', N.Float, ('NFRAMES',))
        RESOLUTIONFUNCTION[:] = self.resolutionFunction[:]
        RESOLUTIONFUNCTION.units = 'unitless'
            
        # Creation of the NetCDF output variables.
        # The frequencies.
        FREQUENCIES = outputFile.createVariable('frequency', N.Float, ('NFRAMES',))
        FREQUENCIES[:] = frequencies[:]
        FREQUENCIES.units = 'THz'

        OMEGAS = outputFile.createVariable('angular_frequency', N.Float, ('NFRAMES',))
        OMEGAS[:] = 2.0*N.pi*frequencies[:]
        OMEGAS.units = 'rad ps-1'

        avacfTotal = N.zeros((self.nFrames), typecode = N.Float)
        adosTotal  = N.zeros((self.nFrames), typecode = N.Float)
                
        comp = 1
        totalMass = 0.0
        for g in groups:

            AVACF = outputFile.createVariable('avacf-group%s' % comp, N.Float, ('NFRAMES',))
            AVACF[:] = self.AVACF[comp][:]
            AVACF.units = 'rad^2*ps^-2'

            N.add(avacfTotal, self.AVACF[comp], avacfTotal)

            ADOS = outputFile.createVariable('ados-group%s' % comp, N.Float, ('NFRAMES',))
            ADOS[:] = self.ADOS[comp][:]
            ADOS.units = 'rad^2*ps^-1'

            N.add(adosTotal, g.mass()*self.ADOS[comp], adosTotal)
            
            comp += 1
            totalMass += g.mass()
            
        adosTotal *= 0.5*self.dt/(self.nGroups*totalMass)

        AVACF = outputFile.createVariable('avacf-total', N.Float, ('NFRAMES',))
        AVACF[:] = avacfTotal
        AVACF.units = 'rad^2*ps^-2'
            
        ADOS = outputFile.createVariable('ados-total', N.Float, ('NFRAMES',))
        ADOS[:] = adosTotal
        ADOS.units = 'rad^2*ps^-1'
                
        asciiVar = sorted(outputFile.variables.keys())
        
        outputFile.close()            
        
        self.toPlot = {'netcdf' : self.output, 'xVar' : 'angular_frequency', 'yVar' : 'ados-total'}

        # Create an ASCII version of the NetCDF output file.
        convertNetCDFToASCII(inputFile = self.output,\
                             outputFile = os.path.splitext(self.output)[0] + '.cdl',\
                             variables = asciiVar)
