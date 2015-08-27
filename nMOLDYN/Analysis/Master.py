"""
This modules implements the estimate, serial and parrallel templates for all analysis.
"""

# The Python distribution modules.
import getpass
import os
import subprocess
import sys
from time import asctime, sleep
from timeit import default_timer

# The Scientific module.

# The nMOLDYN modules.
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Analysis import Analysis

class ParallelAnalysis(object):
                    
    runningMode = 'parallel'
    
    taskName = None
    
def setupPyroServer():
        
    import Pyro.errors
    import Pyro.naming

    # Gets a PyRO proxy for the name server.
    locator = Pyro.naming.NameServerLocator()

    # Try to get an existing name server.        
    try:
        LogMessage('info', 'Searching for a Pyro Name server. Please wait ...', ['console'])
        ns = locator.getNS()
        LogMessage('info', 'A Pyro Name server was found on your machine. nMOLDYN will use it.', ['console'])
        return None
            
    # Otherwise, start a new one.        
    except Pyro.errors.NamingError:
        
        LogMessage('info', 'No Pyro Name server found on your machine.', ['console'])
        LogMessage('info', 'Setting a new Pyro Name server. Please wait ...', ['console'])
        
        s = subprocess.Popen([sys.executable, '-O', '-c', "import  Pyro.naming; Pyro.naming.main([])"], stdout = subprocess.PIPE)
        ns = None
        while ns is None:
            try:
                ns = locator.getNS()
            except Pyro.errors.NamingError:
                pass
            
        LogMessage('info', 'Pyro Name server ready.', ['console'])
        return s
            
def terminate(process):
    """
    Kills a process, useful on 2.5 where subprocess.Popens don't have a 
    terminate method.


    Used here because we're stuck on 2.5 and don't have Popen.terminate 
    goodness.
    """

    def terminate_win(process):
        import win32process
        return win32process.TerminateProcess(process._handle, -1)

    def terminate_nix(process):
        import os
        import signal
        return os.kill(process.pid, signal.SIGTERM)

    terminate_default = terminate_nix

    handlers = {
        "win32": terminate_win, 
        "linux2": terminate_nix
    }

    return handlers.get(sys.platform, terminate_default)(process)

def startSlaves(taskName, architecture, numberOfProcs):
    """Starts the slaves.

    @param architecture: the type of pyro server. One of 'multiprocessor' or 'cluster'.
    @type architecture: string.

    @param numberOfProcs: the number of procs allocated for the analysis.
    @type numberOfProcs: int.
    """

    # Case of an analysis launched only on localhost.
    if architecture == 'multiprocessor':
        
        script = os.path.abspath(os.path.join(GVAR['nmoldyn_path'], 'Analysis', 'Slave.py'))

        for comp in range(numberOfProcs):
            subprocess.Popen([sys.executable, script, taskName])
        
    elif architecture == 'cluster':
        LogMessage('warning', "\nThe analysis will be run in cluster mode. \n\n\
The analysis will not be run right now. To do so, you have to:\n\n\
1) log in on one or several remote machine on which the analysis should be dispatched\n\
2) enter:\n\
\ttask_manager slave %s \n\n\
once per process to be run on the remote machine." % taskName, ['gui', 'console'])

def setUpPyroServer():
    """
    """

    import Pyro.errors
    import Pyro.naming

    # Gets a PyRO proxy for the name server.
    locator = Pyro.naming.NameServerLocator()

    # Try to get an existing name server.        
    try:
        LogMessage('info', 'Searching for a Pyro Name server. Please wait ...', ['console'])
        ns = locator.getNS()
        LogMessage('info', 'A Pyro Name server was found on your machine. nMOLDYN will use it.', ['console'])
        return None
            
    # Otherwise, start a new one.        
    except Pyro.errors.NamingError:
        
        LogMessage('info', 'No Pyro Name server found on your machine.', ['console'])
        LogMessage('info', 'Setting a new Pyro Name server. Please wait ...', ['console'])
        
        s = subprocess.Popen([sys.executable, '-O', '-c', "import  Pyro.naming; Pyro.naming.main([])"], stdout = subprocess.PIPE)
        ns = None
        while ns is None:
            try:
                ns = locator.getNS()
            except Pyro.errors.NamingError:
                pass
            
        LogMessage('info', 'Pyro Name server ready.', ['console'])
        return s

# #################################
# atom-by-atom abstract classes
# #################################                
class SerialPerAtom(object):
    """Template class for an analysis atom-by-atom ran in serial mode.
    """
    
    runningMode = 'serial'
    
    def internalRun(self):
        """Performs the analysis in serial mode.
        """
       
        # Estimate mode. First evaluate the time taken to parse the trajectory file for all selected atoms.
        if self.estimate:
            self.offset = default_timer()
                
            first = True
            for aIndex in self.subset:
                junk = self.trajectory.readParticleTrajectory(aIndex, first = self.first, last = self.last, skip = self.skip).array
                if first:
                    timeToReadFirstAtom = default_timer() - self.offset
                    first = False
                
            self.offset = default_timer() - self.offset - self.nSelectedAtoms*timeToReadFirstAtom

            self.chrono = default_timer()                    

            for aIndex in self.subset:
                atomIndex, x = self.calc(aIndex, self.trajectory)
                self.combine(atomIndex, x)
                self.updateJobProgress(1)
                break
                
            # The estimation is computed as:
            # Tau = Sum(ri - r1) + N*(r1 + c1) with
            # N  = number of atoms.
            # ri = time to read trajectory of atom i
            # r1 = time to read the trajectory of atom 1
            # c1 = time to process atom 1 (by hypothesis c1 ~ c2 ~ c3 ... ~ cN)
            self.chrono = int(self.offset + self.nSelectedAtoms*(default_timer() - self.chrono))
            
        else:
            # The analysis actual starting time.
            self.chrono = default_timer()
            
            for aIndex in self.subset:
                atomIndex, x = self.calc(aIndex, self.trajectory)
                assert aIndex == atomIndex
                self.combine(atomIndex, x)
                self.updateJobProgress(self.nSelectedAtoms)
                
            self.finalize()
                        
            # The actual time taken for the analysis.
            self.chrono = int(default_timer() - self.chrono)
                              
class ParallelPerAtom(ParallelAnalysis):
    """Template class for an analysis atom-by-atom ran in parallel mode.
    """
        
    def internalRun(self):
        """Performs the analysis in parallel mode.
        """

        from Scientific.DistributedComputing.MasterSlave import initializeMasterProcess, TaskRaisedException, GlobalStateValue

        # The Pyro server is setup.
        pyroServer = setUpPyroServer()
        
        # If no task name was assigned to the job, build one.
        if self.taskName is None:                                
            self.taskName = '%s_%s_%s' % (self.db_shortname, getpass.getuser(), '_'.join(asctime().split()))

        if self.architecture == 'multiprocessor':
            tasks = initializeMasterProcess(self.taskName)
            
        elif self.architecture == 'cluster':
            tasks = initializeMasterProcess(self.taskName, slave_module='nMOLDYN.Analysis.Slave')
            
        else:
            raise Error('Illegal parallel mode %s' % self.architecture)

        tasks.setGlobalState(analysis=self)
        
        for aIndex in self.subset:
            task_id = tasks.requestTask("analysisPerElement",
                                        GlobalStateValue(1, 'analysis'),
                                        aIndex,
                                        self.trajectoryFilename)

        startSlaves(self.taskName, self.architecture, self.numberOfProcs)
        
        # The analysis actual starting time.
        self.chrono = default_timer()

        for aIndex in self.subset:
            try:                
                task_id, tag, (atomIndex, x) = tasks.retrieveResult("analysisPerElement")
                self.combine(atomIndex, x)
                self.updateJobProgress(self.nSelectedAtoms)
                
            except TaskRaisedException, e:
                LogMessage('error', e.traceback, ['console'])
                raise
                
            except:
                raise Error('Error when retrieving the results over the pyroserver.')
                
        self.finalize()

        # The actual time taken for the analysis.
        self.chrono = int(default_timer() - self.chrono)

        tasks.shutdown()
        if pyroServer is not None:
            terminate(pyroServer)
                
# #################################
# frame-by-frame abstract classes
# #################################
class SerialPerFrame(object):
    """Template class for an analysis frame-by-frame ran in serial mode.
    """

    runningMode = 'serial'
    
    def internalRun(self):
        """Performs the analysis in serial mode.
        """
        
        if self.estimate:
            self.offset = default_timer()
        
            first = True
            for frameIndex in self.frameIndexes:
                self.universe.setFromTrajectory(self.trajectory, frameIndex)
                if first:
                    timeToReadFirstFrame = default_timer() - self.offset
                    first = False
                
            self.offset = default_timer() - self.offset - self.nFrames*timeToReadFirstFrame

            self.chrono = default_timer()
                
            # Treat only the first frame.
            for fIndex in self.frameIndexes:
                frameIndex, x = self.calc(fIndex, self.trajectory)
                assert fIndex == frameIndex                
                self.combine(frameIndex,x)
                self.updateJobProgress(1)
                break
            
            # The estimation is computed as:
            # Tau = Sum(ri - r1) + N*(r1 + c1) with
            # N  = number of frames.
            # ri = time to read the frame i
            # r1 = time to read the frame 1
            # c1 = time to process frame 1 (by hypothesis c1 ~ c2 ~ c3 ... ~ cN)
            self.chrono = int(self.offset + self.nFrames*(default_timer() - self.chrono))

        else:
            # The analysis actual starting time.
            self.chrono = default_timer()

            for fIndex in self.frameIndexes:
                frameIndex, x = self.calc(fIndex, self.trajectory)
                assert fIndex == frameIndex
                self.combine(frameIndex,x)
                self.updateJobProgress(self.nFrames)
                                            
            self.finalize()
            
            # The actual time taken for the analysis.
            self.chrono = int(default_timer() - self.chrono)

class ParallelPerFrame(ParallelAnalysis):
    """Template class for an analysis frame-by-frame ran in parallel mode.
    """
    
    def internalRun(self):
        """Performs the analysis in parallel mode.
        """

        from Scientific.DistributedComputing.MasterSlave import initializeMasterProcess, TaskRaisedException, GlobalStateValue

        pyroServer = setUpPyroServer()

        if self.taskName is None:                                
            # This should be enough to create a unique task name.
            self.taskName = '%s_%s_%s' % (self.db_shortname, getpass.getuser(), '_'.join(asctime().split()))

        if self.architecture == 'multiprocessor':
            tasks = initializeMasterProcess(self.taskName)
            
        elif self.architecture == 'cluster':
            tasks = initializeMasterProcess(self.taskName, slave_module='nMOLDYN.Analysis.Slave')
            
        else:
            raise Error('Illegal parallel mode %s' % self.architecture)

        # The instance of the running analysis is made global for all the process.
        tasks.setGlobalState(analysis=self)
        
        # The master registers the tasks to be done by the slave.
        for fIndex in self.frameIndexes:
            task_id = tasks.requestTask("analysisPerElement",
                                        GlobalStateValue(1, 'analysis'),
                                        fIndex,
                                        self.trajectoryFilename)

        # The slaves are started but the calculation is not started actually.
        # In case of a cluster run, wait for calls to the task_manager.
        startSlaves(self.taskName, self.architecture, self.numberOfProcs)

        # The analysis actual starting time.
        self.chrono = default_timer()

        for fIndex in self.frameIndexes:
            try:
                task_id, tag, (frameIndex, x) = tasks.retrieveResult("analysisPerElement")
                self.combine(frameIndex,x)
                self.updateJobProgress(self.nFrames)
                
            except TaskRaisedException, e:
                LogMessage('error', e.traceback, ['console'])
                raise
            
            except:
                raise Error('Error when retrieving the results over the pyroserver.')
                
        self.finalize()

        # The actual time taken for the analysis.
        self.chrono = int(default_timer() - self.chrono)

        tasks.shutdown()
        if pyroServer is not None:
            terminate(pyroServer)

# #################################
# group-by-group abstract classes
# #################################
class SerialPerGroup(object):
    """Template class for an analysis group-by-group ran in serial mode.
    """

    runningMode = 'serial'
    
    def internalRun(self):
        """Performs the analysis in serial mode.
        """
       
        # The analysis actual starting time.
        self.chrono = default_timer()

        for aIndexes in self.group:
            
            atomIndexes, x = self.calc(aIndexes, self.trajectory)
            assert atomIndexes == aIndexes
            self.combine(atomIndexes, x)
            self.updateJobProgress(self.nGroups)

        self.finalize()
            
        # The actual time taken for the analysis.
        self.chrono = int(default_timer() - self.chrono)
                        
class ParallelPerGroup(ParallelAnalysis):
    """Template class for an analysis group-by-group ran in parallel mode.
    """
    
    def internalRun(self):
        """Performs the analysis in parallel mode.
        """

        from Scientific.DistributedComputing.MasterSlave import initializeMasterProcess, TaskRaisedException, GlobalStateValue

        pyroServer = setUpPyroServer()
        
        if self.taskName is None:                                
            self.taskName = '%s_%s_%s' % (self.db_shortname, getpass.getuser(), '_'.join(asctime().split()))

        if self.architecture == 'multiprocessor':
            tasks = initializeMasterProcess(self.taskName)
            
        elif self.architecture == 'cluster':
            tasks = initializeMasterProcess(self.taskName, slave_module='nMOLDYN.Analysis.Slave')
            
        else:
            raise Error('Illegal parallel mode %s' % self.architecture)

        tasks.setGlobalState(analysis=self)
        for aIndexes in self.group:
            task_id = tasks.requestTask("analysisPerElement",
                                        GlobalStateValue(1, 'analysis'),
                                        aIndexes,
                                        self.trajectoryFilename)

        startSlaves(self.taskName, self.architecture, self.numberOfProcs)

        # The analysis actual starting time.
        self.chrono = default_timer()

        for aIndexes in self.group:
            try:                
                task_id, tag, (atomIndexes, x) = tasks.retrieveResult("analysisPerElement")
                self.combine(atomIndexes, x)
                self.updateJobProgress(self.nGroups)
                
            except TaskRaisedException, e:
                LogMessage('error', e.traceback, ['console'])
                raise
            
            except:
                raise Error('Error when retrieving the results over the pyroserver.')
                
        self.finalize()

        # The actual time taken for the analysis.
        self.chrono = int(default_timer() - self.chrono)

        tasks.shutdown()
        if pyroServer is not None:
            terminate(pyroServer)

class SerialPerQShell(object):
    """Template class for an analysis qshell-by-qshell ran in serial mode.
    """

    runningMode = 'serial'
    
    def internalRun(self):
        """Performs the analysis in serial mode.
        """
        # The analysis actual starting time.
        self.chrono = default_timer()
       
        for qInd in range(self.nQValues):
            qIndex, x = self.calc(qInd, self.trajectory)
            assert qInd == qIndex
            self.combine(qIndex, x)
            self.updateJobProgress(self.nQValues)
                                            
        self.finalize()

        # The actual time taken for the analysis.
        self.chrono = int(default_timer() - self.chrono)
        
class ParallelPerQShell(ParallelAnalysis):
    """Template class for an analysis qshell-by-qshell ran in parallel mode.
    """
    
    def internalRun(self):
        """Performs the analysis in parallel mode.
        """

        from Scientific.DistributedComputing.MasterSlave import initializeMasterProcess, TaskRaisedException, GlobalStateValue
    
        pyroServer = setUpPyroServer()

        if self.taskName is None:                                
            # This should be enough to create a unique task name.
            self.taskName = '%s_%s_%s' % (self.db_shortname, getpass.getuser(), '_'.join(asctime().split()))

        if self.architecture == 'multiprocessor':
            tasks = initializeMasterProcess(self.taskName)
            
        elif self.architecture == 'cluster':
            tasks = initializeMasterProcess(self.taskName,\
                                            slave_module='nMOLDYN.Analysis.Slave')
            
        else:
            raise Error('Illegal parallel mode %s' % self.architecture)

        tasks.setGlobalState(analysis=self)
        for qIndex in range(self.nQValues):
            task_id = tasks.requestTask("analysisPerElement",
                                        GlobalStateValue(1, 'analysis'),
                                        qIndex,
                                        self.trajectoryFilename)
            
        startSlaves(self.taskName, self.architecture, self.numberOfProcs)

        # The analysis actual starting time.
        self.chrono = default_timer()

        for qInd in range(self.nQValues):
            try:
                task_id, tag, (qIndex, x) = tasks.retrieveResult("analysisPerElement")
                self.combine(qIndex, x)
                self.updateJobProgress(self.nQValues)
                
            except TaskRaisedException, e:
                LogMessage('error', e.traceback, ['console'])
                raise
            
            except:
                raise Error('Error when retrieving the results over the pyroserver.')
                
        self.finalize()
        
        # The actual time taken for the analysis.
        self.chrono = int(default_timer() - self.chrono) 

        tasks.shutdown()
        if pyroServer is not None:
            terminate(pyroServer)
