"""This modules implements the nMOLDYN Preferences and the procedures to load or save Preferences.

The preferences IO are handled in nMOLDYN via ConfigParser mechanism.

Procedures:
    * savePreferencesFile : save Preferences settings to a Preferences file.
    * loadPreferencesFile : load Preferences settings from a Preferences file.
"""

# The python distribution modules
from ConfigParser import ConfigParser
import copy
import os
import platform
import sys


# The nMOLDYN Atoms Database is defined as the default database to search in.
from MMTK import Database, Proteins

from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage

PLATFORM = platform.system().upper()

# This dictionnary will store the global preferences for nMOLDYN. It will be used as a global variable 
# throughout all the code.
PREFERENCES = {}

# the |progress_rate| variable is used to define the rate at which the status of a running job 
# will be displayed on the console, file and gui loggers.
PREFERENCES['progress_rate'] = '10'
            
# The |outputfile_path| variable is used to set the path for the nMOLDYN output file path.
if PLATFORM == 'WINDOWS':
    PREFERENCES['outputfile_path'] = os.environ['USERPROFILE']
                
elif PLATFORM == 'DARWIN':
    PREFERENCES['outputfile_path'] = os.environ['HOME']
                
else:
    PREFERENCES['outputfile_path'] = os.environ['HOME']

# The |trajfile_path| variable is used to set the default path where to search the input trajectories in.
PREFERENCES['trajfile_path'] = PREFERENCES['outputfile_path']

# The |vmd_path| variable is used to set the path for VMD molecular viewer.
PREFERENCES['vmd_path'] = ''
            
# The |ncdump_path| variable is used to set the path for ncdump program.
PREFERENCES['ncdump_path'] = ''

# The |ncgen_path| variable is used to set the path for ncgen program.
PREFERENCES['ncgen_path'] = ''
            
# The |acroread_path| variable is used to set the path for acrobat reader program.
PREFERENCES['acroread_path'] = ''

# The |taksmanager_path| variable is used to set the path for Scientific task_manager script.
PREFERENCES['task_manager_path'] = ''

# The |documentation_style| variable is used to set the format to visualize the documentation in nMOLDYN.
# It can be HTML or PDF.
PREFERENCES['documentation_style'] = 'html'

# Automatically check or not whether a new version of nMOLDYN is installed on the Dirac server.
PREFERENCES['warning_check_for_new_version'] = 'no'
            
# The version of the database that will be used.            
PREFERENCES['database_version'] = 'nMOLDYN'

def savePreferencesFile(filename = None, prefToSave = None):
    """Saves the Preferences settings to a Preferences file.
    
    @param filename: if not None, the name of the file where to save the preferences otherwise the file name
        will be set to an OS dependant location:
            -$USERPROFILE/Application Data/nMOLDYN/nMOLDYN.ini on Windows
            -$HOME/Library/Preferences/nMLDYN/nMOLDYN.pref on MacOS
            -$HOME/.nMOLDYN on Linux
    @type filename: string.

    @param prefToSave: the preferences to save. If None, the default preferences will be saved.
    @type prefToLoad: instance of a dummy class whose attributes must be the nMOLDYN preferences variable names.
    """
    
    prefs = copy.copy(PREFERENCES)
    prefs.update(prefToSave)

    # If no file name was given as input, constructs an arbitrary OS-dependant one.
    if filename is None:
        
        try:
            # Case of Win32.
            if PLATFORM == 'WINDOWS':
                directory = os.path.join(os.environ['USERPROFILE'],'Application Data','nMOLDYN')
                filename = os.path.abspath(os.path.join(directory, 'nMOLDYN.ini'))
                
            # Case of Mac OS.
            elif PLATFORM == 'DARWIN':
                directory = os.path.join(os.environ['HOME'],'Library','Preferences','nMOLDYN')
                filename = os.path.abspath(os.path.join(directory, 'nMOLDYN.pref'))
                
            # Case of Linux.
            else:
                directory = os.environ['HOME']
                filename = os.path.abspath(os.path.join(directory, '.nMOLDYN'))
         
        # The file name could not be built. Raises an error.
        except KeyError:
            raise Error('Invalid default location for the nMOLDYN configuration file.')
            
        # The file name could be built.
        else:
            # Checks that the directory where the preferences file will be saved exists. If not, creates it. 
            if not os.path.exists(directory):
                os.makedirs(directory)
                        
    # Tries to open the preferences file for writing.
    try:
        prefFile = open(filename,'w')
        
    # The file could not be opened. Raises an error.
    except IOError:
        raise Error('Impossible to open the file %s for writing. May be a permission problem.' % filename)
    
    # The file could be opened.
    else:

        # Sets up an instance of configuration parser.
        cfg = ConfigParser()
        
        # Creates a section called 'nmoldyn'.
        cfg.add_section('nmoldyn')
        
        # Loops over the nMOLDYN preferences variables.
        for k, v in prefs.items():
                                            
            if not v:                
                cfg.set('nmoldyn', k, 'undefined')
                
            else:                
                cfg.set('nmoldyn', k, v)
                
        # Write the configuration parser into the preferences file.        
        cfg.write(prefFile)
        
        # Close the preferences file.
        prefFile.close()

def loadPreferencesFile(filename = None):
    """Loads preferences from a preferences file.
    
    @param filename: if not None, the name of the file where to load the preferences otherwise the file name
        will be set to an OS dependant location:
            -$USERPROFILE/Application Data/nMOLDYN/nMOLDYN.ini on Windows
            -$HOME/Library/Preferences/nMLDYN/nMOLDYN.pref on MacOS
            -$HOME/.nMOLDYN on Linux
    @type filename: string.

    """
       
    loadedPrefs = copy.copy(PREFERENCES)
    
    # If no file name was given as input, constructs an arbitrary OS-dependant one.
    if filename is None:
        try:
            # Case of Win32.
            if PLATFORM == 'WINDOWS':
                filename = os.path.join(os.environ['USERPROFILE'],'Application Data','nMOLDYN','nMOLDYN.ini')
                
            # Case of Mac OS.
            elif PLATFORM == 'DARWIN':
                filename = os.path.join('Library','Preferences','nMOLDYN','nMOLDYN.pref')
                
            # Case of Linux.
            else:
                filename = os.path.join(os.environ['HOME'],'.nMOLDYN')

        except KeyError:
            return

    # Case where the preferences file to load exists.
    if os.path.exists(filename):
        try:
            # Sets up an instance of configuration parser.
            prefFile = ConfigParser()
            
            # And load in the preferences file.
            prefFile.read(filename)
            
            if not prefFile.has_section('nmoldyn'):
                raise

        # Something went wrong when loading the preferences file. Raises an error.
        except:                        
            LogMessage('warning','The file %s is not a proper configuration file.' % filename,['console'])
        
        # The preferences file could be loaded.
        else:            
            LogMessage('info','Read the preferences file %s' % filename,['console'])
            
            # Loop over the variables stored in the preferences file as nMOLDYN preferences variable.
            for key, value in prefFile.items('nmoldyn'):
                                
                if value.lower() in ['none','undefined']:
                    continue
                                    
                # If so, sets its value to the corresponding attribute of the target preferences class. 
                loadedPrefs[key] = value
    
    return loadedPrefs

def setPreferences(prefs):
            
    for k, v in prefs.items():
                
        if k in PREFERENCES.keys():
                        
            if k == "database_version":

                PREFERENCES[k] = v.lower()
                                
                if PREFERENCES[k] == "nmoldyn":
                    path = os.path.join(GVAR['nmoldyn_path'], 'Database')
                else:
                    path = os.path.join(GVAR['mmtk_path'], 'Database')
                        
                if path not in Database.path:                      
                    Database.path.insert(0,path)
                else:
                    idx = Database.path.index(path)
                    if idx != 0:
                        Database.path.insert(0,Database.path.pop(idx))                    
                        Database.atom_types = Database.Database('Atoms', Database.AtomType)
                        Database.group_types = Database.Database('Groups', Database.GroupType)
                        Database.molecule_types = Database.Database('Molecules', Database.MoleculeType)
                        Database.crystal_types = Database.Database('Crystals', Database.CrystalType)
                        Database.complex_types = Database.Database('Complexes', Database.ComplexType)
                        Database.protein_types = Database.Database('Proteins', Database.ProteinType)
                        Proteins._residue_blueprints = {}
                                             
            elif k in ["documentation_style", "warning_check_for_new_version"]:
                PREFERENCES[k] = v.lower()
                
            else:
                PREFERENCES[k] = v
                                                
def findPreferenceHelp(name):

    if name == 'progress_rate':
        message = 'The step in percentage at which the job progress will be displayed on the console and/or on the logfile.'
                                                
    elif name == 'outputfile_path':
        message = 'Directory. The path for output files.'
                
    elif name == 'trajfile_path':
        message = 'Directory. The NetCDF trajectory default path.'
                
    elif name == 'vmd_path':
        message  = 'Filename. The path for VMD molecular viewer executable.'
                
    elif name == 'ncdump_path':
        message = 'Filename. The path for NetCDF ncdump program.'
                
    elif name == 'ncgen_path':
        message = 'Filename. The path for NetCDF ncgen program.'
                
    elif name == 'acroread_path':
        message = 'Filename. The path for Acrobat Reader program.'
                                
    elif name == 'documentation_style':
        message = 'html|pdf. The format for any kind of documentation (users guide, API, contextual).'
                                
    elif name == 'warning_check_for_new_version':
        message = 'yes|no. Should nMOLDYN search for new version when it starts ?'
                
    elif name == 'database_version':
        message = 'nMOLDYN|MMTK. Which version of the atoms database to use ?'
                
    # The preferences variables is unknown or still undocumented.
    else:
        message = 'No help available for %s preference variable.' % name
                
    return message
                    
#class nMOLDYNPreferences(object):
#    """This class defines the nMOLDYN preferences variables.
#    
#    This class is built on the Singleton principle. That means that one and just one instance of
#    that class will be created.
#    """
#    
#    class __Singleton:
#                
#        def __init__(self):
#            """The constructor.
#            """
#            
#            self.changed = False
#                                
#            # the |progress_rate| variable is used to define the rate at which the status of a running job 
#            # will be displayed on the console, file and gui loggers.
#            self.progress_rate = '10'
#            
#            # The |outputfile_path| variable is used to set the path for the nMOLDYN output file path.
#            if PLATFORM == 'WINDOWS':
#                self.outputfile_path = os.environ['USERPROFILE']
#                
#            elif PLATFORM == 'DARWIN':
#                self.outputfile_path = os.environ['HOME']
#                
#            else:
#                self.outputfile_path = os.environ['HOME']
#
#            # The |trajfile_path| variable is used to set the default path where to search the input trajectories in.
#            self.trajfile_path = self.outputfile_path
#
#            # The |vmd_path| variable is used to set the path for VMD molecular viewer.
#            self.vmd_path = ''
#            
#            # The |ncdump_path| variable is used to set the path for ncdump program.
#            self.ncdump_path = ''
#
#            # The |ncgen_path| variable is used to set the path for ncgen program.
#            self.ncgen_path = ''
#            
#            # The |acroread_path| variable is used to set the path for acrobat reader program.
#            self.acroread_path = ''
#
#            # The |taksmanager_path| variable is used to set the path for Scientific task_manager script.
#            self.task_manager_path = ''
#
#            # The |documentation_style| variable is used to set the format to visualize the documentation in nMOLDYN.
#            # It can be HTML or PDF.
#            self.documentation_style = 'html'
#
#            self.warning_check_for_new_version = 'no'
#            
#            self.database_version = 'nMOLDYN'
#            
#        def help(self, name):
#            """Returns some information about a given preferences variables.
#            
#            @param name: the name of the preferences variables about which some information is required.
#            @type name: string
#            
#            @return: the information about the selected preferences variables.
#            @rtype: string
#            """
#            
#            if name == 'progress_rate':
#                message = 'The step in percentage at which the job progress will be displayed on the console and/or on the logfile.'
#                                                
#            elif name == 'outputfile_path':
#                message = 'Directory. The path for output files.'
#                
#            elif name == 'trajfile_path':
#                message = 'Directory. The NetCDF trajectory default path.'
#                
#            elif name == 'vmd_path':
#                message  = 'Filename. The path for VMD molecular viewer executable.'
#                
#            elif name == 'ncdump_path':
#                message = 'Filename. The path for NetCDF ncdump program.'
#                
#            elif name == 'ncgen_path':
#                message = 'Filename. The path for NetCDF ncgen program.'
#                
#            elif name == 'acroread_path':
#                message = 'Filename. The path for Acrobat Reader program.'
#                                
#            elif name == 'documentation_style':
#                message = 'html|pdf. The format for any kind of documentation (users guide, API, contextual).'
#                                
#            elif name == 'warning_check_for_new_version':
#                message = 'yes|no. Should nMOLDYN search for new version when it starts ?'
#                
#            elif name == 'database_version':
#                message = 'nMOLDYN|MMTK. Which version of the atoms database to use ?'
#                
#            # The preferences variables is unknown or still undocumented.
#            else:
#                message = 'No help available for %s preference variable.' % name
#                
#            return message
#            
#        def __setattr__(self, name, val):
#            if val is None:
#                val = ''
#                
#            if not hasattr(self, name):
#                self.__dict__['changed'] = True
#            else:
#                if getattr(self, name) != val:
#                    self.__dict__['changed'] = True
#                    
#            self.__dict__[name] = val
#                
#    # This is the singleton mechanism.
#    instance = None
#   
#    # The method __new__ is the cornerstone of a Singleton class.
#    def __new__(c):
#        if not nMOLDYNPreferences.instance:
#            nMOLDYNPreferences.instance = nMOLDYNPreferences.__Singleton()
#        return nMOLDYNPreferences.instance
#
#    def __getattr__(self, name):
#        return getattr(self.instance, name)
#    
#
## The singleton instance of the |nMOLDYNPreferences| class is created.
#PREFERENCES = nMOLDYNPreferences()
