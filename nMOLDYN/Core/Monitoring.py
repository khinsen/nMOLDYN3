"""This modules implements the functions used when monitoring the system to get information for nMOLDYN.

Functions:
    * findNestedDirectories : parses recursively a list of directories and returns the list of nested directories.
    * findExecutable        : searches for an executable and returns its abosulte path if it was found.
    * determineNumberOfCPUs : returns the number of CPUs available on the local machine.
    * cpuInfo               : returns the total numbers of CPUs and the number of loaded and free CPUs on a machine.
"""

# The python distribution modules
from distutils.version import LooseVersion
import platform
import os
import re
import stat
import subprocess
import shutil
import sys
import urllib2

from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Preferences import PREFERENCES

# The platform name.
PLATFORM = platform.system().upper()

def findNestedDirectories(whereToScan):
    """Parses recursively a list of directories and returns the list of their nested directories.

    @param whereToScan: a list specifying the directories to search in.
    @type whereToScan: list

    @return : the nested directories list.
    @rtype : list
    
    @note: the base directories are included in the output list of the nested directories.
    """

    # The list that will contains the nested directories.
    dirList = []    
    
    # Loop over the directories to search in.
    for d in whereToScan:
        
        # Loop of the directory tree generator values.
        # |dirPath| is the path to the current directory
        # |dirNames| is a list of the names of the subdirectories that are in |dirPath|
        # |fileNames| is a list of the names of the files in |dirPath|.
        for dirPath, dirNames, fileNames in os.walk(d):
            
            # If the directory name is not in |dirList| then append it.
            if dirPath not in dirList:
                dirList.append(dirPath)
                
    return dirList
        
def findExecutable(whereToSearch, name):
    """Searches for an executable and returns its absolute path if found.

    @param whereToSearch: a list of directories where to search the executable in.
    @type whereToSearch: string

    @param name: the short name of the executable.
    @type name: string
    
    @return: the absolute path of the executable if found, an empty string otherwise.
    @rtype: string    
    """                   

    # Loop over the directories where to search the executable in.
    for directory in whereToSearch:
        
        # The absolute path of the executable is built.
        full_execname = os.path.join(directory, name)
        
        # Checks if this file exists.
        if os.path.isfile(full_execname):
            
            # If so, checks that it is an executable.
            mode = stat.S_IMODE(os.stat(full_execname)[stat.ST_MODE])            
            if mode & 0111 > 0:
                # If so, OK the executable is found, returns its absolute path.
                return full_execname
            
    # If at the end of the loop, the executable could not be found then return an empty string.
    else:
        return ''

# The list that will store the paths to explore when searching for an external program. 
# This is not a preferences variable.
paths = []
for d, l ,f in os.walk(sys.prefix):
    paths.append(d)
paths.append(os.getcwd())
            
if PLATFORM == 'WINDOWS':
    for k in ['PROGRAMFILES', 'SYSTEMDRIVE', 'HOMEDRIVE']:
        if os.environ.has_key(k):
            paths.append(os.environ[k]+ 'University of Illinois' + 'VMD')
                
else:
    if os.environ.has_key('PATH'):
        paths.extend(os.environ['PATH'].split(os.pathsep))
            
def searchProgram(progName):
    """
    """
                 
    # File processing part        
    # Loop over the directories to search in.
    for d in paths:
                                
        # The absolute path of the executable is built.
        full_execname = os.path.join(d, progName)
                        
        # Checks if this file exists.
        if os.path.isfile(full_execname):
            
            # If so, checks that it is an executable.
            mode = stat.S_IMODE(os.stat(full_execname)[stat.ST_MODE])
                    
            if mode & 0111 > 0:
                # If so, OK the executable is found, returns its absolute path.
                return full_execname
                                                   
    return None                

def getCPUInfo():
    """Sets the total numbers of processors, the number of loaded and free processors on the host 
    machine or on the different nodes of a cluster.
    """

    cpuInfo = None

    # Pyro server in cluster mode implemented only for linux platform.
    if PLATFORM == 'LINUX':
        nProcs = file('/proc/cpuinfo','r').read().count('processor\t:')
            
        nLoadedProcs = min(nProcs,int(os.getloadavg()[1] + 0.5))
        nFreeProcs = max(0,nProcs - nLoadedProcs)
            
    elif PLATFORM == 'DARWIN':
        try:
            nProcs = int(subprocess.Popen(['/usr/sbin/sysctl', '-n', 'hw.ncpu'], stdout = subprocess.PIPE).stdout.read())
        except:
            nProcs = 1
                
        nLoadedProcs = min(nProcs,int(os.getloadavg()[1] + 0.5))
        nFreeProcs = max(0,nProcs - nLoadedProcs)
        hostname = subprocess.Popen(['/bin/hostname'], stdout = subprocess.PIPE).stdout.read()
                
    elif PLATFORM == 'WINDOWS':
        try:
            from win32com.client import GetObject

            nProcs = 0
            loadAvg = 0.0

            wmi = GetObject('winmgmts:')
            cpu = wmi.InstancesOf('Win32_Processor')
            for c in cpu:
                nProcs += c.Properties_('NumberOfCores').Value
                loadAvg += c.Properties_('LoadPercentage').Value
            loadAvg /= nProcs

            nLoadedProcs = int(nProcs * loadAvg / 100.0)
            nFreeProcs = max(0,nProcs - nLoadedProcs)
            
        except:
            LogMessage('warning', 'You do not have the priviledge to perform win32 calls.\n\
nMOLDYN can not define the actual number of loaded processors.', ['gui','console'])
            nProcs = nFreeProcs = int(os.environ['NUMBER_OF_PROCESSORS'])
            nLoadedProcs = 0
            
    cpuInfo = (nProcs,nLoadedProcs,nFreeProcs)

    return cpuInfo

def checkForNewVersion():
                        
    # Find automatically the proxies if some are defined.
    httpProxy = urllib2.getproxies()        

    # The requests will go through the http proxy.
    proxy = urllib2.ProxyHandler(httpProxy)

    # Open the connection possibly through the proxy..
    opener = urllib2.build_opener(proxy)

    # The url for the file storing the last nMOLDYN version.    
    url = 'http://dirac.cnrs-orleans.fr/~nmoldyn/last_version'

    # Build the url request for the file storing the last nMOLDYN version.
    req = urllib2.Request(url)
        
    # Open the url.
    try:
        f = opener.open(req)
        
    # The url could not be opened. Raises an error.
    except urllib2.URLError:
        LogMessage('warning', 'Could not open the url %s.\nPerhaps a problem with a proxy.\n\n\
Can not check whether a new version of nMOLDYN was released.' % url, ['console', 'gui'])
        return

    # The url could be opened. Its contents will be extracted.
    else:
        # The name of the last nMOLDYN version.
        lastVersion = f.read().strip()
        f.close()
        opener.close()
        
    if LooseVersion(vstring = lastVersion) > LooseVersion(vstring = Cfg.nmoldyn_version):
        return lastVersion
    else:
        return None
                                                                
def downloadVersion(version):
    """Fetch the requested version of nMOLDYN from the nMOLDYN server.
    """    

    # Find automatically the proxies if some are defined.
    httpProxy = urllib2.getproxies()        

    # The requests will go through the http proxy.
    proxy = urllib2.ProxyHandler(httpProxy)

    # Open the connection possibly through the proxy..
    opener = urllib2.build_opener(proxy)
    
    if PLATFORM == 'WINDOWS':
        filename = 'nMOLDYN-%s.exe' % version
                                
    # Case of Linux.
    else:
        filename = 'nMOLDYN-%s.zip' % version
        
    url = 'http://dirac.cnrs-orleans.fr/~nmoldyn/' + filename
    
    dest_filename = os.path.join(PREFERENCES['outputfile_path'], filename)
    
    try:
        fileReq = urllib2.Request(url)            
        src = opener.open(fileReq)
                      
    except urllib2.URLError:
        LogMessage('warning', 'Could not open the url %s.' % url, ['console', 'gui'])
      
    else:
        try:
            dst = open(dest_filename, 'w');
            
        except IOError:
            LogMessage('warning', 'Can not open the file %s for writing. Maybe a permission problem.' % dest_filename, ['console', 'gui'])
            return
        
        else:
            LogMessage('info', 'Downloading %s file. Please wait ...' % filename, ['console'])
            shutil.copyfileobj(src, dst)
            dst.close()
            LogMessage('info', '%s file successfully downloaded in %s' % (filename, dest_filename), ['console'])

    opener.close()
    