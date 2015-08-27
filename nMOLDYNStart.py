#!python
"""The nMOLDYN launcher.
"""

# The python distribution modules.
from optparse import OptionGroup, OptionParser, OptionValueError
import subprocess
import sys
import unittest
    
# The MMTK distribution modules.
from MMTK.Trajectory import Trajectory

# The nMOLDYN modules.
from nMOLDYN.__pkginfo__ import __version__ as NMOLDYN_VERSION
from nMOLDYN.Analysis.Templates import *
from nMOLDYN.Chemistry.Chemistry import hierarchizeUniverse
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Core.Monitoring import checkForNewVersion, downloadVersion, getCPUInfo, searchProgram
from nMOLDYN.Preferences import findPreferenceHelp, PREFERENCES

# Try to import the nMOLDYN stability tests.
try:
    import nMOLDYN.Tests.StabilityTests
except ImportError:
    LogMessage("warning", "No tests available.", ["console"])

def trajectoryContents(option, opt_str, value, parser):
    """
    """

    if len(parser.rargs) == 0:
        raise Error('No arguments for %s option.' % option)
    
    trajName = parser.rargs[0]
    
    try:
        t = Trajectory(None, trajName, 'r')
        
    except IOError:
        raise Error('The trajectory could not be opened.')
        
    contents, dummy = hierarchizeUniverse(t.universe)
    
    LogMessage('info', '#'*110, ['console'])
    LogMessage('info', 'Contents of the universe contained in %s trajectory:\n' % trajName, ['console'])

    if len(parser.rargs) == 1:
        
        LogMessage('info', 'Objects found in %s:' % trajName, ['console'])        
        
        for k, v in sorted(contents.items()):
            
            LogMessage('info', '    * %s (%s) with selection keywords:' % (k, v['objectclass']), ['console'])
            
            for kk in sorted(contents[k].keys()):
                
                if kk in ['number', 'objectclass']:
                    continue
                
                LogMessage('info', '        - %s' % kk, ['console'])

    elif len(parser.rargs) == 2:
        
        objectName = parser.rargs[1]
        
        if contents.has_key(objectName):
            
            LogMessage('info', 'Selection keyword found in "%s" nMOLDYN name:' % objectName, ['console'])
            
            for kk in sorted(contents[objectName].keys()):
                
                if kk in ['number', 'objectclass']:
                    continue
                
                LogMessage('info', '    * %s' % kk, ['console'])

        else:
            LogMessage('info', 'The trajectory %s does not contain object %s.' % (trajName, objectName), ['console'])
    
    elif len(parser.rargs) == 3:
        
        objectName, selectionKeyword = parser.rargs[1:]
        
        if contents.has_key(objectName):
            
            if contents[objectName].has_key(selectionKeyword):
                
                LogMessage('info', 'Selection values associated to %s objet and %s selection keyword:' % (objectName, selectionKeyword), ['console'])
                
                for kk in sorted(contents[objectName][selectionKeyword]):
                    
                    if kk in ['number', 'objectclass']:
                        continue
                    
                    LogMessage('info', '    * %s' % kk, ['console'])
            else:
                LogMessage('info', '%s is not a valid selection keyword for %s object.' % (selectionKeyword, objectName), ['console'])

        else:
            LogMessage('info', 'The trajectory %s does not contain object %s.' % (value, objectName), ['console'])
        
    else:
        raise OptionValueError('Too many arguments for %s option.' % option)

    LogMessage('info', '#'*110, ['console'])
    
    sys.exit(0)
    
def testAnalysis(option, opt_str, value, parser):
    """
    """

    if len(parser.rargs) == 0:
        raise Error('No arguments for %s option.' % option)
        
    if 'ALL' in parser.rargs:
        selectedTests = availableTests
        
    else:
        selectedTests = parser.rargs

    for n in selectedTests:
        try:
            t = eval(n.strip().upper() + 'Tests()')
            
        except NameError:
            LogMessage('error', '%s test is not a valid analysis benchmark name.' % n, ['console'])
            
        else:
                        
            try:
                t.run()
                                
            except:
                raise Error('An error occured when running %s test.' % n)
            
            else:
                LogMessage('info', t.summary, ['console'])

    sys.exit(0)
    
def checkForUpdate(option, opt_str, value, parser):
    
    newVersion = checkForNewVersion()
    
    if newVersion is not None:
        
        LogMessage('info', 'nMOLDYN %s is available.' % newVersion, ['console'])
            
        downloadIt = ''
        while not downloadIt in ['y', 'yes', 'n', 'no']:
            LogMessage('info', 'Do you want to download it ? [y]es|[n]o', ['console'])
            downloadIt = raw_input().strip().lower()
                
        if downloadIt in ['y', 'yes']:                    
            downloadVersion(newVersion)        
            
    sys.exit(0)

def runSlaveJobs(option, opt_str, value, parser):
    """
    """
    
    if len(parser.rargs) == 0:
        raise Error('You must provide at least one arguments for %s option.' % option)
    
    if not PREFERENCES['task_manager_path']:
        task_manager_path = searchProgram('task_manager')
    
    else:
        if os.path.isfile(PREFERENCES['task_manager_path']):
            task_manager_path = PREFERENCES['task_manager_path']
            
        else:
            task_manager_path = searchProgram('task_manager')
            
    if task_manager_path is None:
        raise Error('The path for task manager script of Scientific module could not be found.')
    
    nTotalProcs, nLoadedProcs, nFreeProcs = getCPUInfo()
    
    if len(parser.rargs) == 1:
        
        slaveJobName = parser.rargs[0]        
                
        if nFreeProcs == 0:
            LogMessage('warning', 'There is no free processor on this node. It would not be nice to force the analysis \
to run on this node. However, you can still override this by restarting the script with the number of processors to allocate \
as the second argument.', ['console', 'gui'])
            sys.exit(0)
        
        nAllocatedProcs = nFreeProcs
        
    elif len(parser.rargs) == 2:
        slaveJobName, nAllocatedProcs = parser.rargs[0:2]
        
        try:
            nAllocatedProcs = int(nAllocatedProcs)
            
            if nAllocatedProcs < 1:
                raise ValueError
            
        except ValueError:
            raise Error('The number of processors to allocate for the analysis must be an integer >= 1')
                        
        if nAllocatedProcs > nTotalProcs:
            LogMessage('warning', 'The number of procs to allocate (%d) is higher than the number of processors available \
on the node (%d). The time taken for the analysis might be much slower than you expect.' % (nAllocatedProcs, nTotalProcs), ['console', 'gui'])
            
    try:
        for slave in range(nAllocatedProcs):
            subprocess.Popen([task_manager_path, 'slave', slaveJobName])
            
    except:
        raise Error('The analysis could not be run properly on cluster mode.')
    
    else:
        sys.exit(0)
                        
def runFromScript(option, opt_str, value, parser):
    """
    """

    try:
        f = open(value, 'r')
        
        parameters = {}
        exec f in None, parameters
        
        analysis = None
        estimate = False
        for k in parameters.keys():
            
            if k.lower() == 'analysis':
                analysis = parameters[k]
                del parameters[k]
                break
                
            elif k.lower() == 'estimate':
                estimate = True
                analysis = parameters[k]
                del parameters[k]
                break

        if analysis is None:
            raise
        
        a = analysis(parameters)
        if estimate:            
            a.estimateAnalysis()
        else:
            a.runAnalysis()
        
        f.close()
        
    except:
        Error('Error when running nMOLDYN from %s input file.' % value)
        
    sys.exit(0)
    
def set_pref_from_command_line(option, opt_str, value, parser):
    """This function will change the selected PREFERENCES variables with the value given in command-line.
    """
    
    PREFERENCES[option.dest] = value
    LogMessage('info', 'Preferences variable "%s" changed to %s' % (option.dest, value), ["console", "file"])        
    
def set_nmoldyn_args_parser():
        
    parser = OptionParser(version='nMOLDYN %s' % NMOLDYN_VERSION)

    parser.add_option('-c', \
                      '--contents', \
                      action='callback', \
                      callback=trajectoryContents, \
                      help='If FILE provided, displays the contents of the trajectory file FILE. If FILE and MOLNAME \
                      provided displays the selection keywords associated to MOLNAME. If FILE and MOLNAME and SELKWD \
                      provided displays the selection values associated to selection keyword SELKWD.', \
                      metavar='FILE MOLNAME SELKWD')

    parser.add_option('-i', \
                      '--input', \
                      action='callback', \
                      callback=runFromScript, \
                      dest='input', \
                      type='string', \
                      help='Start nMOLDYN in CONSOLE mode using the input file FILE', \
                      metavar='FILE')

    parser.add_option('-n', \
                      '--netcdf', \
                      action='store', \
                      dest='netcdf', \
                      type='string', \
                      help='Start nMOLDYN GUI using the netCDF file FILE', \
                      metavar='FILE')

    parser.add_option('-s', \
                      '--slave', \
                      action='callback', \
                      callback=runSlaveJobs, \
                      help='.', \
                      metavar='JOBNAME NPROCS')

    parser.add_option('-t', \
                      '--test', \
                      action='callback', \
                      callback=testAnalysis, \
                      help='Run  nMOLDYN in test mode for analysis ANALYSIS1, ANALYSIS2 ...', \
                      metavar='ANALYSIS1 ANALYSIS2 ...')

    parser.add_option('-u', '--update', action='callback', callback=checkForUpdate, \
                      help='Check for update.', metavar='ANALYSIS1 ANALYSIS2 ...')

    # Adds the option allowing for settings the preferences variables on the fly.
    group = OptionGroup(parser, 'PREFERENCES options')
    for k in PREFERENCES.keys():
    
        group.add_option('--' + k, \
                         action="callback", \
                         callback=set_pref_from_command_line, \
                         type='string', \
                         nargs=1, \
                         dest=k, \
                         help=findPreferenceHelp(k))
    
    parser.add_option_group(group)
    
    return parser

if __name__ == "__main__":
    
    parser = set_nmoldyn_args_parser()
    
    (options, args) = parser.parse_args()
    
    try:
                        
        from nMOLDYN.GUI.MainDialog import MainDialog
                        
        application = MainDialog(options.netcdf)                
                
    except:
           
        raise Error("Could not set the nMOLDYN GUI. Only the command-line mode will be available.") 

