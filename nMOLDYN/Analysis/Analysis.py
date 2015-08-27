"""
This modules implements the base class for all the analysis available in nMOLDYN.

Classes:

    * Analysis: the base class from which all nMOLDYN analysis classes inherits from.
"""

# The python distribution modules
from distutils.version import LooseVersion
import getpass
import glob
import operator
import os
import platform
import re
import sys
from time import asctime
from timeit import default_timer
from unittest import TestSuite

# The Scientific modules
from Scientific import N
from Scientific.Geometry import Vector

# The MMTK modules
from MMTK import Atom, AtomCluster, Molecule
from MMTK.Collections import Collection
from MMTK.NucleicAcids import NucleotideChain
from MMTK.ParticleProperties import ParticleScalar
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain, Protein
from MMTK.Trajectory import Trajectory, TrajectorySet
from MMTK import Units

# The nMOLDYN modules
from nMOLDYN.__pkginfo__ import __version__ as NMOLDYN_VERSION
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Chemistry.Chemistry import belongToAnAmine, belongToAHydroxy, belongToAMethyl, belongToAThiol, hierarchizeUniverse, residusChemFamily
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Core.Misc import convertTime, parseInterval
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Tests.StabilityTests import StabilityTests, StabilityTestsResults

class Analysis(object):
    """Base class for all analysis defined in nMOLDYN. All analysis inherits from this class.

    The class |Analysis| defines attributes and methods common to all the analysis available in nMOLDYN. 
    """

    def __init__(self, parameters = None, statusBar = None):
        """The constructor.

        @param parameters: if not None, a dictionnary storing the input parameters names and their corresponding values.
        @type parameters: dict

        @param statusBar: if not None, an instance of nMOLDYN.GUI.Widgets.StatusBar. Will attach a status bar to the 
            selected analysis.
        @type statusBar: instance of nMOLDYN.GUI.Widgets.StatusBar
        """

        # Will store the time taken for the analysis in seconds.        
        self.chrono = None

        # Check that the preferences value for the progress rate is an integer in ]0,100[.
        try:    
            self.displayProgressRate = int(PREFERENCES['progress_rate'])
            if (self.displayProgressRate <= 1) or (self.displayProgressRate >= 100):
                self.displayProgressRate = 10
            else:
                raise ValueError
        # Otherwise set the progress rate to 10 %.
        except ValueError:
            self.displayProgressRate = 10

        # The statusbar instance displaying the job progress.
        self.statusBar = statusBar

        # Checks that |parameters| argument is a dictionnary and copy it to |self.parameters| attribute.
        # This core attribute will store internally the input parameters from which the analysis will be set and run.
        if isinstance(parameters, dict):
            self.parameters = parameters
        # Otherwise, sets it to None.
        else:        
            self.parameters = None
            
    def setInputParameters(self, parameters):
        """Sets the input parameters for the analysis.

        @param parameters: if not None, a dictionnary associating the input parameters names of an analysis (key) to
            their value (value).
        @type parameters: dict
        """

        # Checks that |parameters| argument is a dictionnary and copy it to |self.parameters| attribute.
        # This core attribute will store internally the input parameters from which the analysis will be set and run.
        if isinstance(parameters, dict):
            self.parameters = parameters

        # Otherwise, throws an error.
        else:        
            raise Error('Bad input for parameters constructor argument. It must be a dictionnary')        

    def interpreteInputParameters(self):
        """Parses and checks the input parameters for the analysis.

        The input parameters are stored in |self.parameters| dictionnary.
        """

        # Checks for the occurence of some optional input parameters.
        # The version number of nMOLDYN. If not stored in |parameters| then sets it to 'unknown' string.
        if not self.parameters.has_key('version'):
            self.parameters['version'] = 'unknown'

        # Case of an estimable analysis.
        if self.db_estimable:
            # If 'estimate' is not a key of |self.parameters| dictionnary, sets it to 'no'.
            if not self.parameters.has_key('estimate'):
                self.parameters['estimate'] = 'no'

        # Case of an unestimable analysis.
        else:
            # If 'estimate' is a key of |self.parameters| dictionnary and its value is 'yes', throws an error because
            # the analysis is not estimable.
            if self.parameters.has_key('estimate'):
                if self.parameters['estimate'] == 'yes':
                    raise Error('%s analysis can not be estimated.' % self.__class__.__name__.split('_')[0])

            # If 'estimate' is not a key of |self.parameters| dictionnary, sets it to 'no'.
            else:
                self.parameters['estimate'] = 'no'

        # Checks if some input parameters are not useful for the analysis. If so, deletes them and sends a warning.
        for p in self.parameters.keys():
            # By construction 'version' and 'estimate' are not specific any analysis. No needs to check for them.
            if p.lower() in ['version','estimate']:
                continue

            if p.lower() not in [v.lower() for v in self.db_parameternames]:
                LogMessage('warning','%s is not a useful parameter for %s analysis.' % (p, self.__class__.__name__),['file','console'])
                del self.parameters[p]

        # Checks if some input parameters are missing. If so, raises an error.
        for p in self.db_parameternames:
            if p.lower() not in [k.lower() for k in self.parameters.keys()]:
                raise Error('%s input parameter is missing for %s analysis.' % (p, self.__class__.__name__))

        # Now loops over the updated |self.parameters| dictionnary keys to check in details all the parameters.
        # Once checked, new class attributes will be defined from the processed value of the corresponding input parameters.
        for pName in self.parameters.keys():

            # The value corresponding to |pName| key.
            pValue = self.parameters[pName]
            # The input parameter name is lowerized to avoid any case-sensitive typo.
            pName = pName.lower()

            # The 'angminmax' input parameter.
            if pName == 'angminmax':                    
                try:
                    # Sets the |self.dMin| and |self.dMax| attributes to respectively the min and max distance.
                    self.angMin, self.angMax = [float(v) for v in pValue.split(':')]
                except:
                    raise Error('Error when parsing "angminmax" parameter: must be a string of the form angmin:angmax.')
                    
                # The max distance must be > min distance.
                if self.angMin >= self.angMax:
                    raise Error('Error when parsing "angminmax" parameter: the max angle must be > min angle.')

            # The 'armodelorder' input parameter. Must be an integer or a float or string that can be converted to an 
            # integer.
            elif pName == 'armodelorder':

                try:
                    # Sets the |self.arModelOrder| attribute to the corresponding input parameter value.
                    self.arModelOrder = int(pValue)

                # Otherwise, raises an error.
                except:
                    raise Error('Error when parsing % parameter: must be an integer.' % pName)

            # The 'atomorder' input parameter. Must be None or a string equal to 'no' or of the form 'atomname1,atomnam2,...'.
            elif pName == 'atomorder':

                # Case of a string.
                if isinstance(pValue, str):
                    # Sets the |self.atomOrder| attribute to the ','-splitted list.
                    self.atomOrder = [p.strip() for p in pValue.split(',')]

                    # If the list is empty, raises an error beacuse there might be a typo in the input parameter value..
                    if not self.atomOrder:
                        raise Error('Error when parsing % parameter: must be a string of the form "atomname1,atomname2...".' % pName)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)

            # The 'truebreakstep' input parameter. Must be an integer.
            elif pName == 'truebreakstep':

                try:
                    # Sets the |self.arModelOrder| attribute to the corresponding input parameter value.
                    self.trueBreakStep = int(pValue)
                    
                # Otherwise, raises an error.
                except:
                    raise Error('Error when parsing % parameter: must be an integer.' % pName)

            # The 'comselection' input parameter. Not many controls here. They will be done further.
            elif pName == 'comselection':

                # If None, sets the |self.comselection| attribute to 'all'.
                if pValue is None:                
                    self.comDefinition = 'all'

                # Otherwise, sets the |self.comselection| attribute to the corresponding input parameter value.
                else:
                    self.comDefinition = pValue

            # The 'deuteration' input parameter. Not many controls here. They will be done further.
            elif pName == 'deuteration':

                # If None, sets the |self.deuteration| attribute to 'no'. There will be no deuteration.
                if pValue is None:
                    self.deuterationDefinition = 'no'

                # Otherwise, sets the |self.deuteration| attribute to the corresponding input parameter value.
                else:
                    self.deuterationDefinition = pValue

            # The 'disminmax' input parameter.
            elif pName == 'disminmax':                    
                try:
                    # Sets the |self.dMin| and |self.dMax| attributes to respectively the min and max distance.
                    self.disMin, self.disMax = [float(v) for v in pValue.split(':')]
                except:
                    raise Error('Error when parsing "disminmax" parameter: must be a string of the form dismin:dismax.')
                    
                # The max distance must be > min distance.
                if self.disMin >= self.disMax:
                    raise Error('Error when parsing "disminmax" parameter: the max distance must be > min distance.')

            # The 'distanceunits' input parameter.
            elif pName == 'distanceunits':
                p = pValue.strip().lower()
                self.distanceUnits = pValue
                if p == 'nm':
                    self.distanceConv = 1.0
                    
                elif p == 'ang':
                    self.distanceConv = 10.0
                    
                elif p == 'fm':
                    self.distanceConv = 1000000.0

            # The 'differentiation' input parameter.
            elif pName == 'differentiation':

                # Must be an integer or a value that can be converted to an integer.
                try:
                    # Sets the |self.differentiation| attribute to the corresponding input parameter value.
                    self.differentiation = int(pValue)

                # Otherwise, raises an error.
                except:
                    raise Error('Error when parsing %s parameter: must be an integer.' % pName)

            # The 'estimate' input parameter. Must be a string equal to 'yes' or 'no'.
            elif pName == 'estimate':

                # Must be a string.
                if isinstance(pValue,str):                
                    p = pValue.strip().lower()
                    # If 'yes', sets the |self.estimate| attribute to True. The analysis will be estimated.
                    if p == 'yes':
                        self.estimate = True
                    # If 'no', sets the |self.estimate| attribute to False. The analysis will not be estimated.
                    elif p == 'no':
                        self.estimate = False

                    # Otherwise raises an error.
                    else:
                        raise Error('Error when parsing %s parameter: must be "yes" or "no".' % pName)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string equal to "yes" or "no".' % pName)

            # The 'frequencyunits' input parameter.
            elif pName == 'frequencyunits':
                p = pValue.strip().lower()
                self.frequencyUnits = pValue
                if p == 'thz':
                    self.frequencyConv = 1.0
                    
                elif p == 'rad s^-1':
                    self.frequencyConv = 2.0*N.pi
                    
                elif p == 'cm^-1':
                    self.frequencyConv = 1.0/Units.invcm
                    
                elif p == 'mev':
                    self.frequencyConv = Units.C*Units.tera*Units.h/(Units.J*Units.s*Units.milli)
                    
                elif p == 'uev':
                    self.frequencyConv = Units.C*Units.tera*Units.h/(Units.J*Units.s*Units.micro)
                    
            # The 'group' input parameter. Not many controls here. They will done further.
            elif pName == 'group':

                # If None, sets the |self.group| attribute to 'all'.
                if pValue is None:
                    self.groupDefinition = 'all'

                # Otherwise, sets the |self.group| attribute to the corresponding input parameter value.
                else:
                    self.groupDefinition = pValue

            # The 'output' input parameter. The name of the output file name.
            elif pName == 'hkls':
                                
                self.hkls = pValue

            # The 'output' input parameter. The name of the output file name.
            elif pName == 'output':

                self.output = pValue

            # The 'projection' input parameter. Must be None or a string equal to 'no' or of 
            # the form 'x,y,z' where x, y and z are the coordinates of the projection vector.
            elif pName == 'projection':

                # Case of None. Sets the |self.projection| attribute to None. There will be no projection.
                if pValue is None:
                    self.projection = None

                # Case of a string.
                elif isinstance(pValue, str):

                    # if 'no', then sets the |self.projection| attribute to None. There will be no projection.
                    if pValue.strip().lower() == 'no':
                        self.projection = None

                    # Otherwise, it must be a string of the form 'x,y,z'.
                    else:

                        # Sets the |self.projection| attribute to the comma-splitted list.
                        try:
                            self.projection = Vector([float(v) for v in pValue.split(',')]).normal()

                        except:
                            raise Error('Error when parsing % parameter: must be a string of the form "x,y,z".' % pName)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)

            elif pName == 'pyroserver':

                if self.runningMode == 'serial':
                    self.architecture = 'monoprocessor'
                    self.numberOfProcs = None

                else:
                    
                    try:
                        import Pyro
                    except:                        
                        raise Error('The PyRO module is not properly installed (if installed).')
                        
                    if pValue is None:
                        self.architecture = 'monoprocessor'
                        self.numberOfProcs = None

                    elif isinstance(pValue, str):

                        try:
                            if pValue.strip().lower() in ['monoprocessor', 'no']:
                                self.architecture = 'monoprocessor'
                                self.numberOfProcs = None

                            elif pValue.strip()[:14].lower() == 'multiprocessor':
                                self.architecture = 'multiprocessor'
                                self.numberOfProcs = int(pValue.split(':')[1])

                            elif pValue.strip().lower() == 'cluster':
                                self.architecture = 'cluster'
                                self.numberOfProcs = None

                        except:
                            raise Error('Error when parsing %s parameter string in %s case.' % (pName,pValue.strip().lower()))

                        else:
                            if self.architecture != 'monoprocessor':
                                if self.statusBar is not None:
                                    LogMessage('warning', 'The job status bar will be desactivated for compatibility with Pyro module.',['console'])
                                    self.statusBar = None

            # The 'qshellvalues' input parameter It must be a list/tuple of floats or a string of the form 
            # "qmin1:qmax1:dq1,qmin2:qmax2:dq2..." where qmin1,qmin2 ... qmax1,qmax2 ... and dq1, dq2 ... are
            # respectively the minimum, maximum and step q values for interval 1, 2 ...
            elif pName == 'qshellvalues':

                # Case of a list or a tuple.
                if isinstance(pValue,(list,tuple)):

                    # It must be a list/tuple of floats.
                    try:                                        
                        self.qShellValues = [float(v) for v in pValue]

                    # Otherwise, raises an error.
                    except:
                        raise Error('Error when parsing % parameter: must be a list/tuple of floats.' % pName)

                # Case of a string.
                elif isinstance(pValue, str):

                    # The |self.qShellValues| attribute will store the generated q values.
                    self.qShellValues = []

                    # Loop over all the defined semi colon-separated q-intervals.
                    for v in pValue.split(';'):

                        qValues = parseInterval(v)

                        if len(qValues) == 0:
                            raise Error('Error when parsing %s parameter: the q interval %s could not be parsed properly.' % (pName,v))

                        # The min q value must be >= 0.
                        if qValues[0] < 0:
                            raise Error('Error when parsing % parameter: the minimum value for q interval %s must be >= 0.' % (pName,v))

                        # The q values are sorted.                
                        self.qShellValues.extend(qValues)

                    self.qShellValues = sorted(set(self.qShellValues))
                    
                else:
                    self.qShellValues = None
                    
            # The 'qshellwidth' input parameter. Must be a float or a value convertible to a float > 0. 
            elif pName == 'qshellwidth':
                # The value must be a float or convertible to a float.
                try:
                    self.qShellWidth = float(pValue)                    

                # Otherwise, raises an error.
                except:
                    self.qShellWidth = None

                else:
                    # If the value is <= 0, raises an error.
                    if self.qShellWidth <= 0.0:
                        raise Error('Error when parsing %s parameter: must be a float >= 0.' % pName)

            # The 'qunits' input parameter.
            elif pName == 'qunits':
                p = pValue.strip().lower()
                self.qUnits = pValue
                if p == 'nm^-1':
                    self.qConv = 1.0
                    
                elif p == 'ang^-1':
                    self.qConv = 0.1
                    
                elif p == 'fm^-1':
                    self.qConv = 0.000001
                    
            # The 'qvectorsdirection' input parameter. Must be None, or a list/tuples of vectors or a string of the form 
            # 'qvx1,qvy1,qvz1;qvx2,qvy2,qvz2 ...' where qvx1, qvx2 ..., qvy1, qvy2 ..., qvz1, qvz2 ... are respectively
            # the x, y and z values of q vector 1, 2 ...
            elif pName == 'qvectorsdirection':

                # If None, sets the |self.qVectorsDirection| attribute to None. The generation will be isotropic.
                if pValue is None:
                    self.qVectorsDirection = None

                # Case of a string.
                elif isinstance(pValue, str):

                    # If 'no', sets the |self.qVectorsDirection| attribute to None. The generation will be isotropic.
                    if pValue.strip().lower() == 'no':
                        self.qVectorsDirection = None

                    # Otherwise must be a string of the form 'qvx1,qvy1,qvz1;qvx2,qvy2,qvz2 ...' otherwise raises an error.
                    else:
                        try:
                            self.qVectorsDirection = [Vector([float(vv.strip()) for vv in v.strip().split(',')]).normal() for v in pValue.split(';')]
                        except:
                            raise Error('Error when parsing %s parameter: must be a string of the form "qvx1,qvy1,qvz1;qvx2,qvy2,qvz2 ...".' % pName)

                # Case of a list or tuple. Must be a list/tuple of Vectors.
                elif isinstance(pValue, (list,tuple)):
                    try:
                        self.qVectorsDirection = [v.normal() for v in pValue]
                    except:
                        raise Error('Error when parsing %s parameter: must be a list/tuple of non-null vectors".' % pName)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)

            # The 'qvectorsgenerator' input parameter. Must be a string equal to '3d isotropic', '2d isotropic' or 'anisotropic'. 
            elif pName == 'qgenerator':

                # Case of a string.
                if isinstance(pValue,str):

                    self.qGenerator = pValue.strip().lower()                    
                    if self.qGenerator not in ['automatic', 'userdefined']:
                        raise Error('Error when parsing %s parameter: must be a string equal to "automatic" or "userdefined".' % pName)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)

            # The 'qgeometry' input parameter. Must be a string equal to 'spatial', 'planar' or 'axial'. 
            elif pName == 'qgeometry':

                # Case of a string.
                if isinstance(pValue,str):

                    self.qGeometry = pValue.strip().lower()                    
                    if self.qGeometry not in ['spatial', 'planar', 'axial']:
                        raise Error('Error when parsing %s parameter: must be a string equal to "spatial", "planar" or "axial".' % pName)

                # Otherwise, raises an error.
                else:
                    self.qGeometry = None

            # The 'qvectorsgenerator' input parameter. Must be a string equal to '3d isotropic', '2d isotropic' or 'anisotropic'. 
            elif pName == 'qvectors':
                self.qVectorsDict = pValue

            # The 'qvectorspershell' input parameter. Must be an integer or a value convertible to an integer > 0.
            elif pName == 'qvectorspershell':

                # The value must be an int or convertible to an int.
                try:
                    self.qVectorsPerShell = int(pValue)

                # Otherwise, raises an error.
                except:
                    self.qVectorsPerShell = None

                else:
                    # If the value is <= 0, raises an error.
                    if self.qVectorsPerShell <= 0:
                        raise Error('Error when parsing %s parameter: must be an integer > 0.' % pName)

            # The 'referenceframe' input parameter. Must be an integer or a value convertible to an integer > 0.
            elif pName == 'referenceframe':

                # The value must be an int or convertible to an int.
                try:
                    self.referenceFrame = int(pValue) - 1

                # Otherwise, raises an error.
                except:
                    raise Error('Error when parsing %s parameter: must be an integer.' % pName)

                else:
                    # If the value is < 0, raises an error.
                    if self.referenceFrame < 0:
                        raise Error('Error when parsing %s parameter: must be an integer > 0.' % pName)

            # The 'removetranslation' input parameter. Must be a string equal to 'yes' or 'no'.
            elif pName == 'removetranslation':

                # Case of a string.
                if isinstance(pValue,str):
                    p = pValue.strip().lower()
                    # If 'yes', sets the |self.removeTranslation| attribute to True. The translation will be removed.
                    if p == 'yes':
                        self.removeTranslation = True
                    # If 'no', sets the |self.removeTranslation| attribute to False. The translation will be removed.
                    elif p == 'no':
                        self.removeTranslation = False
                    # Otherwise, raises an error.
                    else:
                        raise Error('Error when parsing %s parameter: must be a string equal to "yes" or "no".' % pName)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)

            # The 'resolution' input parameter. Must be a float in ]0,100[.
            elif pName == 'resolution':

                # Must be a float or a value that can be converted to an integer.
                try:
                    # Sets the |self.energyFullWidthHalfMaxf| attribute to the corresponding input parameter value.
                    # This must be the FWHM of the energy in meV.
                    self.energyFullWidthHalfMax = float(pValue)

                    # The float must be in > 0.
                    if self.energyFullWidthHalfMax < 0.0:
                        raise Error('Error when parsing %s parameter: must be a float > 0.' % pName)

                    # The sigma of the energy gaussian (in meV).
                    self.energySigma = self.energyFullWidthHalfMax/(2.0*N.sqrt(2.0*N.log(2.0)))

                    # The corresponding sigma of the time gaussian (in ps).
                    self.timeSigma = 1.0/(1.5192669*self.energySigma)

                # Otherwise raises an error.
                except:
                    raise Error('Error when parsing %s parameter: must be a float.' % pName)


            # The 'rvalues' input parameter. Must be a string of the form 'rmin:rmax:dr' where rmin, rmax and dr are respectively 
            # the minimum, maximum and step distances.
            elif pName == 'rvalues':

                self.rValues = parseInterval(pValue)

                if len(self.rValues) <= 1:
                    raise Error('Error when parsing %s parameter: the r interval %s could not be parsed properly.' % (pName,v))

                self.rMin = self.rValues[0]

                self.dR = self.rValues[1] - self.rValues[0]

                self.nRBins = len(self.rValues) - 1

                # The rmin value must be >= 0.
                if self.rMin < 0.0:
                    raise Error('Error when parsing %s parameter: rmin must be >= 0.' % pName)

            # The 'stepwiserbt' input parameter. Must be a string equal to 'yes' or 'no'.
            elif pName == 'stepwiserbt':

                # Case of a string.
                if isinstance(pValue,str):
                    p = pValue.strip().lower()
                    # If 'yes', sets the |self.stepwiseRBT| attribute to True.
                    if p == 'yes':
                        self.stepwiseRBT = True
                    # If 'no', sets the |self.stepwiseRBT| attribute to False.
                    elif p == 'no':
                        self.stepwiseRBT = False
                    # Otherwise, raises an error.
                    else:
                        raise Error('Error when parsing %s parameter: must be a string equal to "yes" or "no".' % pName)

                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)

            # The 'subset' input parameter. Not many controls here. They will done further.
            elif pName == 'subset':

                # If None, sets the |self.subset| attribute to 'all'.
                if pValue is None:                
                    self.subsetDefinition = 'all'

                # Otherwise, sets the |self.subset| attribute to the corresponding input parameter value.
                else:
                    self.subsetDefinition = pValue

            # The 'subset' input parameter. Not many controls here. They will done further.
            elif pName == 'subset1':

                # If None, sets the |self.subset| attribute to 'all'.
                if pValue is None:                
                    self.subset1Definition = 'all'

                # Otherwise, sets the |self.subset| attribute to the corresponding input parameter value.
                else:
                    self.subset1Definition = pValue

            # The 'subset' input parameter. Not many controls here. They will done further.
            elif pName == 'subset2':

                # If None, sets the |self.subset| attribute to 'all'.
                if pValue is None:                
                    self.subset2Definition = 'all'

                # Otherwise, sets the |self.subset| attribute to the corresponding input parameter value.
                else:
                    self.subset2Definition = pValue

            # The 'group' input parameter. Not many controls here. They will done further.
            elif pName == 'target':

                # If None, sets the |self.target| attribute to 'all'.
                if pValue is None:
                    self.targetDefinition = 'all'

                # Otherwise, sets the |self.target| attribute to the corresponding input parameter value.
                else:
                    self.targetDefinition = pValue

            # The 'timeinfo' input parameter. Must be a string of the form 'framemin:framemax:framestep' where framemin
            # framemax, framestep are respectively the minimum, maximum and step values of the frames to consider in the 
            # analysis.
            elif pName == 'timeinfo':

                # Case of a string.
                if isinstance(pValue,str):

                    # The value must be three semi colon separated integer.
                    try:
                        self.first, self.last, self.skip = [int(v) for v in pValue.split(':')]
                    # Otherwise, throws an error.
                    except:
                        raise Error('Error when parsing %s parameter: must be a string of the form "framemin:framemax:framestep".' % pName)
                    else:

                        # framemin  must be > 0.
                        if self.first <= 0:
                            raise Error('Error when parsing %s parameter: framemin must be >= 1".' % pName)

                        # framemax must be >= framemin
                        if self.last < self.first:
                            raise Error('Error when parsing %s parameter: framemax must be >= framemin".' % pName)

                        self.first = self.first - 1

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)
                
            # The 'timeunits' input parameter.
            elif pName == 'timeunits':
                p = pValue.strip().lower()
                self.timeUnits = pValue
                if p == 'ps':
                    self.timeConv = 1.0
                    
                elif p == 'ns':
                    self.timeConv = 0.001
                    
                elif p == 'fs':
                    self.timeConv = 1000.0

            # 'trajectory' input parameter. It must be a string or an instance of MMTK.Trajectory.
            elif pName == 'trajectory':

                # If not a string or a Trajectory, raises an error.
                if not isinstance(pValue, (str,Trajectory)):
                    raise Error('Error when parsing %s parameter: it must be a string or an instance of MMTK.Trajectory class.' % pName)

                # If the value is a string, it must be a filename of the trajectory that will be loaded.
                # Depending on the extension, it will be loaded as MMTK.TrajectorySet or as a MMTK.Trajectory.
                if isinstance(pValue, str):

                    # The extension of the filename of the trajectory to load.
                    ext = os.path.splitext(pValue)[1]

                    # Some info displayed on the console.
                    LogMessage('info', 'Reading trajectory. Please wait ...', ['console'])

                    # Case of a Trajectory set. Loads a MMTK.TrajectorySet object.
                    if ext == '.ncs':

                        try:

                            # The trajectory set is opened for reading and its contents sent to |trajSet| list.
                            trajSetFile = open(pValue, 'r')
                            trajSet = [t.strip() for t in trajSetFile.readlines()]
                            trajSetFile.close()

                            # Load the whole trajectory set.
                            self.trajectory = TrajectorySet(None, trajSet)

                        except:
                            raise Error('Error when parsing %s parameter: could not load the trajectory set file %s.' % (pName,pValue))

                    # Otherwise try to load it as a trajectory..
                    else:

                        try:
                            # Load the whole trajectory set.
                            self.trajectory = Trajectory(None, pValue, 'r')
                            
                        except:
                            raise Error('Error when parsing %s parameter: could not load the trajectory file %s.' % (pName,pValue))

                # If the value is already a Trajectory MMTK object, then sets it directly to |self.trajectory| attribute.
                elif isinstance(pValue, Trajectory):
                    self.trajectory = pValue

                # Some info to display on the console.
                LogMessage('info', 'Trajectory read %s' % str(pValue), ['console'])

                # The |self.trajectoryFilename| attribute is set to the name of the trajectory file.
                self.trajectoryFilename = self.trajectory.filename

                # The |self.universe| attribute is set to the universe contained in the loaded trajectory.
                self.universe = self.trajectory.universe

                # Some info to display on the console.
                LogMessage('info', 'Setting universe contents. Please wait ...', ['console'])

                # This procedure will set up the chemical contents of the universe. This will be used when performing atom selections.
                self.chemicalObjectInfo, self.chemicalObjectNames = hierarchizeUniverse(self.universe)

                self.atomInformation = {}
                for at in self.universe.atomList():
                    self.atomInformation[at.index] = {'type' : at.type.name,\
                                                      'element' : at.symbol,\
                                                      'name' : at.name,\
                                                      'molecule' : id(at.topLevelChemicalObject()),\
                                                      'weight' : None,\
                                                      }

                # Some info to display on the console.                
                LogMessage('info', 'Done', ['console'])
                LogMessage('info', '', ['console'])

            # The 'version' input parameter. Must be a strong of the 'x.y.z' where x, y and z define the version number
            # of nMOLDYN.
            elif pName == 'version':

                # Case of a string.
                if isinstance(pValue,str):

                    # If 'unknown' sets |self.version| attribute to a dummy version number equal to '0.0.0'.
                    if pValue == 'unknown':
                        self.version = LooseVersion(vstring = '0.0')

                    # Otherwise, sets |self.version| attribute to the parameter value.
                    else:
                        self.version = LooseVersion(vstring = pValue)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string of the form "x.y.z".' % pName)

            # The 'weights' input parameter. Must be a string equal to 'equal', 'mass', 'coherent', 'incoherent' or 'atomicnumber'.
            elif pName == 'weights':

                # Case of a string.
                if isinstance(pValue, str):

                    # It is lowerized for better flexibility.
                    p = pValue.strip().lower()

                    # If 'no', the |self.weightingScheme| attribute is set to 'equal'.
                    if p == 'no':
                        self.weightingScheme = 'equal'
                    # In the other cases, the |self.weightingScheme| attribute is set to the parameter value.
                    elif p in ['equal', 'mass', 'coherent', 'incoherent', 'atomicnumber']:
                        self.weightingScheme = p
                    # Otherwise, raises an error.
                    else:
                        raise Error('Error when parsing %s parameter: must be a string equal to "equal", "mass", "coherent", "incoherent" or "atomicNumber".' % pName)

                # Otherwise, raises an error.
                else:
                    raise Error('Error when parsing %s parameter: must be a string.' % pName)

    def buildTimeInfo(self):
        """Builds some attributes related to the frame selection string entered using 'timeinfo' input parameter. 

        They will be used to define at which times a given analysis should be run.
        """

        # The |self.frameIndexes| attribute is set to a list storing the indexes of the selected frames for the analysis.
        self.frameIndexes = range(self.first, self.last, self.skip)

        self.firstFrame = self.frameIndexes[0]

        # If the trajectory has 'time' attribute, then the trajectory contains the actual times of the MD.
        # In that case, sets |t| to the actual times corresponding to the selected frames. 
        if hasattr(self.trajectory, 'time'):
            t = self.trajectory.time[self.first:self.last:self.skip]

        # Otherwise, sets |t| to virtual times.
        else:
            t = range(len(self.trajectory))[self.first:self.last:self.skip]

        # |self.dt| attribute will contain the analysis time step.
        # If there was just one selected frame, sets the |self.dt| attribute to a virtual time step of 1 ps.
        if len(t) == 1:
            self.dt = 1
        else:
            self.dt = t[1] - t[0]

        if self.dt == 0.0:
            raise Error('The time step is null.')

        # The |self.times| attribute is set as t - t initial.
        self.times = self.dt * N.arange(len(t))

        # The |self.nFrames| attribute is set to the number of selected frames of the analysis.        
        self.nFrames = len(self.times)

    def saveAnalysis(self, filename):
        """Saves the settings of an analysis to an input file that can be further reused without the GUI.

        There are two possible formats. The first one corresponds to nMOLDYN autostart files 
        (extension '.py') that are python scripts that can be launched directly from a terminal. The second 
        corresponds to nMOLDYN input files (extension '.nmi') that corresponds to the historic pMoldyn input 
        files of nMOLDYN 2. These one must be launched using 'nMOLDYNStart.py -i file.nmi'.

        @param filename: the name of the file. If the extension is '.nmi' the file will be a
            nMOLDYN input file otherwise the file will be a nMOLDYN autostart file.
        @type: filename: string
        """
        
        orderedParameterNames = sorted(self.db_parameternames)

        if 'trajectory' in orderedParameterNames:
            orderedParameterNames.insert(0,orderedParameterNames.pop(orderedParameterNames.index('trajectory')))

        if 'output' in orderedParameterNames:
            orderedParameterNames.append(orderedParameterNames.pop(orderedParameterNames.index('output')))

        # Case of an file name with '.nmi' extension. The file will be saved as a nMOLDYN input file.
        if filename[-4:].lower() == '.nmi':

            # Opens a unit for writing.
            pScript = open(filename, 'w')

            # Write the input file header.
            pScript.write('################################################################\n')
            pScript.write('# This is an automatically generated python-nMOLDYN run script #\n')
            pScript.write('################################################################\n\n')

            # Sets |analysis| variable to the name of the analysis to save. 
            pScript.write(('\nanalysis = %s\n\n') % (self.__class__.__name__))

            # Sets first the general parameters |version| and |estimate| that will respectively specify the
            # nMOLDYN current version and if the analysis should be run in estimate mode or not.
            pScript.write('################################################################\n')
            pScript.write('# General parameters                                           #\n')
            pScript.write('################################################################\n\n')
            pScript.write(('version = "%s"\n') % NMOLDYN_VERSION)
            pScript.write(('estimate = "%s"\n\n') % self.parameters['estimate'])

            # Sets the analysis-specific parameters.
            pScript.write('################################################################\n')
            pScript.write('# Analysis-specific parameters                                 #\n')
            pScript.write('################################################################\n\n')

            # Loops over the input parameters names. For each of them write 'name = value' where name and value 
            # are respectively the name and value of the current input paramter.
            for pName in orderedParameterNames:
                if pName == 'trajectory':
                    prefix = '\n# The input trajectory\n'
                    suffix = '\n'

                elif pName == 'output':
                    prefix = '\n# The output trajectory\n'
                    suffix = '\n'

                else:
                    prefix = ''
                    suffix = ''

                # The value of the current input parameter.
                pValue = self.parameters[pName]
                # If the parameter is of string type, quotes its value explicitely.
                if isinstance(pValue, str):
                    pScript.write(('%s%s = "%s"%s\n') % (prefix, pName, pValue, suffix))
                else:
                    pScript.write(('%s%s = %s%s\n') % (prefix, pName, pValue, suffix))

            # The unit is closed.
            pScript.close()            

        # Others kind of extension extension. The file will be saved as a nMOLDYN autostart file.
        else:
            
            # Opens a unit for writing.
            pScript = open(filename, 'w')

            # The first line contains the call to the python executable. This is necessary for the file to
            # be autostartable.
            pScript.write('#!%s\n' % sys.executable)

            # Writes the input file header.
            pScript.write('################################################################\n')
            pScript.write('# This is an automatically generated python-nMOLDYN run script #\n')
            pScript.write('################################################################\n\n')
            
            # Some additional entries if the run is parallel.
            if self.runningMode == "parallel":
                pScript.write('from optparse import OptionParser\n')        
                pScript.write('parser = OptionParser()\n')        
                pScript.write('parser.add_option("-t", "--task", dest = "task", help = "The taskname under which the job will be distributed by PyRO server")\n')
                pScript.write('(options, args) = parser.parse_args()\n\n')                                                        
            
            # Writes the 'import' lines for MMTk and for the selected analysis.
            pScript.write("from MMTK import *\n")
                    
            pScript.write("from nMOLDYN.Core.Logger import LogCreate\n")

            pScript.write(("from %s import %s\n\n") % (self.__module__, self.__class__.__name__))

            pScript.write("LogCreate(['console', 'file'])\n\n")
                        
            # Writes the line that will initialize the |parameters| dictionnary.
            pScript.write("parameters = {}\n\n")

            # Sets first the general parameters |version| and |estimate| that will respectively specify the
            # nMOLDYN current version and if the analysis should be run in estimate mode or not.
            pScript.write('################################################################\n')
            pScript.write('# General parameters                                           #\n')
            pScript.write('################################################################\n\n')
            pScript.write('parameters[\'version\'] = "%s"\n' % NMOLDYN_VERSION)
            pScript.write('parameters[\'estimate\'] = "%s"\n\n' % self.parameters['estimate'])

            # Sets the analysis-specific parameters.
            pScript.write('################################################################\n')
            pScript.write('# Analysis-specific parameters                                 #\n')
            pScript.write('################################################################\n\n')

            # Loops over the input parameters names. For each of them write 'parameters[name] = value' where name 
            # and value are respectively the name and value of the current input paramter.
            for pName in orderedParameterNames:
                if pName == 'trajectory':
                    prefix = '# The input trajectory\n'
                    suffix = '\n'

                elif pName == 'output':
                    prefix = '\n# The output trajectory\n'
                    suffix = '\n'

                else:
                    prefix = ''
                    suffix = ''

                # The value of the current input parameter.
                pValue = self.parameters[pName]
                # If the parameter is of string type, quotes its value explicitely.
                if isinstance(pValue, str):
                    pScript.write(('%sparameters[\'%s\'] = "%s"%s\n') % (prefix, pName, pValue, suffix))
                else:
                    pScript.write(('%sparameters[\'%s\'] = %s%s\n') % (prefix, pName, pValue, suffix))

            pScript.write('################################################################\n')
            pScript.write('# Setup and run the analysis                                   #\n')
            pScript.write('################################################################\n\n')
            # Sets |analysis| variable to an instance analysis to save. 
            pScript.write(('analysis = %s(parameters)\n\n') % (self.__class__.__name__))
            
            # Some additional entries if the run is parallel.
            if self.runningMode == "parallel":
                pScript.write('analysis.taskName = options.task\n\n')
                
            # The line that will actually launch the analysis.
            pScript.write('analysis.runAnalysis()')
            
            # The unit is closed.
            pScript.close()

            # Make the script executable on posix systems.
            if os.name == 'posix':
                os.system('chmod u+x %s' % filename)

    def runAnalysis(self):
        """Runs an analysis.

        @return: a dictionnary of the form {'days' : d, 'hours' : h, 'minutes' : m, 'seconds' : s} specifying the 
            time the analysis took in days, hours, minutes and seconds.
        @rtype: dict
        """

        # Stores the number of steps of the analysis loop that has been processed so far.
        self.jobCounter = 0

        try:
            # Initialize the analysis.
            self.initialize()

        except:
            raise Error('An error occured while initializing the analysis.')

        # Checks the version number of the nMOLDYN used to run the analysis. 
        # If the version is anterior to the current version (can happen only for an analysis run 
        # in script mode), displays a warning. Some parameters may be deprecated.
        if self.version < LooseVersion(vstring = NMOLDYN_VERSION):
            LogMessage('warning','The nMOLDYN version number of the script is unknown or deprecated.\n\
Some analysis keywords name/values may have changed. Please check in case of troubles.',['file','console'])
            LogMessage('info', '', ['console'])

        # Builds and displays the informations about the running analysis.
        self.buildJobInfo()

        # Actually run the analysis by calling the analysis-specific |internalRun| method.
        try:
            self.internalRun()

        # Raises an error, if something went wrong during the analysis.
        except:
            raise Error('An error occured while running the analysis.')

        # The ending date of the analysis is displayed.
        LogMessage('info','Job finished on '+asctime()+'.\n',['console','file'])

        # The time taken for the analysis is computed.
        timeTakenForAnalysis = convertTime(self.chrono)

        # The log message that displays the time taken for the analysis is a little bit different 
        # depending on the running mode (estimate or not).
        if self.estimate:

            # The elapsed time for the (estimated) analysis is written in the log file.
            LogMessage('info','Estimated time for analysis: %(days)s days %(hours)s hours %(minutes)s minutes %(seconds)s seconds.\n' % timeTakenForAnalysis,['file','console'])

        else:

            # The elapsed time for the (estimated) analysis is written in the log file.
            LogMessage('info','Elapsed time for analysis: %(days)s days %(hours)s hours %(minutes)s minutes %(seconds)s seconds.\n' % timeTakenForAnalysis,['file','console'])

        return timeTakenForAnalysis

    def updateJobProgress(self, norm):
        """Checks the progress of the running analysis and displays periodically on the console and the logfile 
            its status. 

        Called each time a step of an analysis loop is achieved.

        @param norm: the maximum number of steps of the analysis.
        @type: norm: integer        
        """

        # Computes the old percentage and stores it in |oldProgress|.
        t = int(100*self.jobCounter/norm)        
        oldProgress = (t/self.displayProgressRate)*self.displayProgressRate

        # One step of the analysis loop has finished so increments the |self.jobCounter| by one unit.
        self.jobCounter += 1

        # Computes the new percentage.
        t = int(100*self.jobCounter/norm)
        newProgress = (t/self.displayProgressRate)*self.displayProgressRate

        # If the new percentage is different from the previous one, displays it on the console and file loggers. 
        if newProgress != oldProgress:
            LogMessage('info','Progress rate = %3d %%' % newProgress,['file','console'])

        # If the analysis is run from the GUI, updates the status bar.
        if self.statusBar is not None:
            self.statusBar.setValue(t)            

    def buildJobInfo(self):
        """Builds and displays on the console and file loggers the main information about the analysis.
        """

        # The analysis starting time.
        jobCreationTime = asctime()

        self.information = ''

        self.information += '#'*90 +'\n'
        self.information += 'Job information for %s analysis.\n' % self.__class__.__name__
        self.information += '#'*90 +'\n\n'

        # General informations.
        self.information += 'Job launched on: %s\n\n' % jobCreationTime
        self.information += 'General informations\n'
        self.information += '--------------------\n'
        self.information += 'User: %s\n' % getpass.getuser()
        self.information += 'OS: %s-%s\n' % (platform.system(), platform.release())
        self.information += 'Processor: %s\n' % platform.machine()
        self.information += 'nMOLDYN version: %s\n' % self.parameters['version']
        self.information += 'Estimate run: %s\n\n' % self.parameters['estimate']

        # Information about the parameters of the analysis.
        self.information += 'Parameters\n'
        self.information += '----------\n'
        self.information += '\n'.join(['%s = %s' % (k,self.parameters[k]) for k in self.parameters if k not in ['version','estimate']])
        self.information += '\n\n'

        self.information += 'Job status\n'
        self.information += '----------\n'

        # The informations are displayed on the console and file loggers.
        [LogMessage('info', l, ['file', 'console']) for l in self.information.splitlines()]

    def defineWeights(self, subset, deuterated, weightingScheme):
        """Returns the weights of a given subset of atoms. 

        The values of the weights depends on the selected weighting scheme. Some of the weighting schemes 
        are sensitive to deuterium (e.g. coherent, incoherent).

        @param subset: the subset of atoms to consider when defining the weights.
        @type subset: an instance MMTK.Collections.Collection

        @param weightingScheme: the weighting scheme.
        @type weightingScheme: string being one of 'equal', 'mass', 'coherent' , 'incoherent' or 'atomicNumber'.

        @note: the weighting in nMOLDYN use some extra properties of the atom database that are only stored in the
        database provided with nMOLDYN and that might be missing from the one provided by MMTK (e.g. atomic_number, mass_deut).
        """

        orderedAtoms = sorted(self.universe.atomList(), key = operator.attrgetter('index'))

        selectedAtoms = Collection([orderedAtoms[ind] for ind in subset])
        selectedAtomsForDeuteration = Collection([orderedAtoms[ind] for ind in deuterated])

        # The equal weighting scheme.
        if weightingScheme == 'equal':
            try:
                # MMTK ParticleScalar object containing 1 for each atom of the universe
                weights = ParticleScalar(self.universe, N.ones(self.universe.numberOfAtoms(), typecode = N.Float))

            # Something went wrong when setting the weights. Raises an error.
            except:
                raise Error('Error when defining weights based on equal weighting scheme.')

        # The mass weighting scheme.
        elif weightingScheme == 'mass':
            try:
                # The atomic masses of the universe.
                weights = self.universe.masses()

                for at in selectedAtomsForDeuteration:
                    weights[at] = at.mass_deut

            # Something went wrong when setting the weights. Raises an error.
            except:
                raise Error('Error when defining weights based on mass weighting scheme.')

        # The atomic number weighting scheme.
        elif self.weightingScheme == 'atomicnumber':

            try:
                # MMTK ParticleScalar object containing 1 for each atom of the universe
                weights = ParticleScalar(self.universe)

                for at in self.universe.atomList():
                    weights[at] = at.atomic_number

            # Something went wrong when setting the weights. Raises an error.
            except:
                raise Error('Error when defining weights based on atomic number weighting scheme.')

        # The coherent weighting scheme.
        elif self.weightingScheme == 'coherent':
            try:
                # getAtomScalarArray is a wrap for getParticleScalar function of Universe.py module
                # MMTK ParticleScalar object containing the bcoh. 
                weights = self.universe.getAtomScalarArray('b_coherent')
                
                # The b value of the hydrogen to be deuterated are changed to b value of deuterium.
                for at in selectedAtomsForDeuteration:
                    weights[at] = at.b_coherent_deut
            
            # Something went wrong when setting the weights. Raises an error.
            except:
                raise Error('Error when defining weights based on coherent weighting scheme.')

        # The incoherent weighting scheme.
        elif self.weightingScheme == 'incoherent':
            try:
                # getAtomScalarArray is a wrap for getParticleScalar function of Universe.py module
                # MMTK ParticleScalar object containing the bincoh. 
                weights = self.universe.getAtomScalarArray('b_incoherent')
                
                # The b value of the hydrogen to be deuterated are changed to b value of deuterium.
                for at in selectedAtomsForDeuteration:
                    weights[at] = at.b_incoherent_deut

                # b^2
                weights = weights*weights

            # Something went wrong when setting the weights. Raises an error.
            except:
                raise Error('Error when defining weights based on incoherent weighting scheme.')

        # Unknown weighting scheme. Raises an error.
        else:
            raise Error('%s weighting scheme is unknown.' % scheme)

        # If there is an atom selection. A boolean mask is applied on the atoms to select.
        weights = weights*selectedAtoms.booleanMask()

        # The normalization factor is different for the coherent weighting scheme.
        if weightingScheme == 'coherent':
            weights = weights/N.sqrt((weights*weights).sumOverParticles())
        else:
            weights = weights/weights.sumOverParticles()

        for aIndex in subset:
            self.atomInformation[aIndex]['weight'] = weights[aIndex]

    def setElementInformation(self):
        """Returns the number of each element found in a given subset and deuteration selection and their corresponding weight.

        @param subset: a subset of atoms.
        @type subset: an instance MMTK.Collections.Collection

        @param weights: the atomic weights.
        @type weights: an instance of MMTK.ParticleProperties.ParticledScalar

        @return: a list of two dictionnaries. The first dictionnary gives some informations about the different elements
        found in |subset| and |deuteration| collections. Its keys are the elements found in |subset| and 'D' for the ones in
        |deuterium| and its values are dictionnaries storing their corresponding number (key = 'number') and
        weight (key  = 'weight'). The second dictionnary is a llok up table between the atom index and the corresponding
        element.
        @rtype: a list of two dictionnaries
        """

        numberOfAtoms = len(self.atomInformation)        

        # The dictionnary storing the informations about the elements.
        self.elementInformation = {}

        for aIndex in self.subset:
            element = self.atomInformation[aIndex]['element']

            # If the |elements| dictionnary have already the key corresponding to the atom symbol,
            # just increment the 'number' subkey by 1.0.
            if self.elementInformation.has_key(element):

                self.elementInformation[element]['number'] += 1.0
                self.elementInformation[element]['mask'][aIndex] = 1

            # If the |elements| dictionnary does not have the key corresponding to the atom symbol, sets
            # the 'weight' and 'number' subkeys respectively to the atom symbol and 1.0.
            else:

                self.elementInformation[element] = {'weight' : self.atomInformation[aIndex]['weight'],\
                                                    'number' : 1.0,\
                                                    'mask' : N.zeros((numberOfAtoms,), typecode = N.Int0)}
                self.elementInformation[element]['mask'][aIndex] = 1

    def selectAtoms(self, selection, message = None):
        """Returns a subset of atoms matching a given selection.

        @param selection: which atoms to select. It can be:
            - a string: in that case it will be parsed and the results will be returned.
            - a list or tuple of MMTK atoms: in that case it will be converted to a MMTK Collection and returned.
            - a MMTK Collection, in that case it will be directly returned.
        @type selection: string or list or tuple or a MMTK Collection.

        @return: the atoms that matches the selection.
        @rtype: instance of MMTK.Collections.Collection
        """                

        subset = []

        # Case where |selection| is a MMTK Collection.
        if isinstance(selection, Collection):
            subset = [at.index for at in selection]
            
        # Case where |selection| is a function. The function must only have the universe as arguments and
        # return the selected atoms as a list.
        elif callable(selection):
            subset = [at.index for at in selection(self.universe)]

        # Case where |selection| is a list or a tuple.
        elif isinstance(selection, (list, tuple)):            
            # Tries to build a MMTK Collection out of it. That triggers that the list or the tuple must contain only 
            # MMTK Atoms.
            try:
                subset = [at.index for at in selection]

            # Failed to create the MMTK Collection. Raises an error.
            except:
                raise Error('The results of a subset selection must be a MMTK Collection.')

        # Remaining case where |selection| is a selection string.
        elif isinstance(selection, str):

            # If the string is equal to 'all', returns all the atoms of the universe in a MMTK Collection.
            if selection.lower() == 'all':
                subset = [at.index for at in self.universe.atomList()]

            # Otherwise, parses the selection string.
            else:

                # If the selection string starts with "filename name-of-file" then it will be a selection from a file.
                if selection.lower().strip()[0:8] == 'filename':
                    # The format for such a selection is 'filename name-of-the-file'.
                    try:
                        filename = selection[8:].strip()

                    # Bad format. Raises an error.
                    except:
                        raise Error('Invalid format for\n%s\nfile-based selection string.' % selection)

                    # The format is OK. Performs the selection from |filename| file.
                    else:
                        subset = self.__parseFileSelectionString(filename, 'subset')[0]

                # If the selection string starts with "expression" then it will be a selection from an valid Python expression 
                # expression.
                elif selection.lower().strip()[0:10] == 'expression':

                    # Checks that the expression assigns somehwere selection variable to something.
                    if not re.findall('selection\s*=',selection):
                        raise Error('Invalid format for\n%s\n expression-based selection string. It must contain "selection = list-of-atoms".' % selection)

                    subset = self.__parseExpressionSelectionString(selection, 'subset')

                # If the selection string starts with "object" then it will be a selection from the chemical objects 
                # available in the universe.
                elif selection.lower().strip()[0:6] == 'object':

                    subset = self.__parseObjectSelectionString(selection[7:].strip())

                else:
                    raise Error('The object-based selection string\n%s\n does not start with "object", "filename" or "expression".' % selection)

                # Tries to convert the result of the parsing of the selection string to a MMTK Collection.
                try:
                    subset = list(subset)

                # If the conversion fails, raises an error.
                except:
                    raise Error('The parsing of the object-based selection string \n%s\nproduced an invalid result.\
There might be something wrong with it.' % selection)
                    
        else:
            raise Error('Wrong format for the subset selection. Must be a string or a MMTK Collection' % str(selection))

        subsetSelection = sorted(subset)

        if len(subsetSelection) == 0:
            raise Error('Empty subset selection.')
        
        if message is None:
            message = 'Subset selection'

        LogMessage('info', '', ['console', 'file'])
        LogMessage('info', message, ['console', 'file'])
        LogMessage('info', '-'*len(message), ['console', 'file'])
        LogMessage('info', 'Number of selected atoms: %s' % len(subsetSelection), ['console', 'file'])
        LogMessage('info', '', ['console', 'file'])
        LogMessage('info', '', ['console', 'file'])

        return subsetSelection

    def selectAtomsForDeuteration(self, selection, message = None):
        """Returns a subset of hydrogens atoms that will treated as deuterium in the analysis that match a given selection.

        @param universe: the universe on which the selection will be performed.
        @type universe: instance of MMTK.Universe

        @param selection: which atoms to select. It can be:
            - a string: in that case it will be parsed and the results will be returned.
            - a list or tuple of MMTK atoms: in that case it will be converted to a MMTK Collection and returned.
            - a MMTK Collection, in that case it will be directly returned.
        @type selection: string or list or tuple or a MMTK Collection.

        @return: the atoms that matches the selection.
        @rtype: instance of MMTK.Collections.Collection
        """

        # Case where |selection| is a MMTK Collection.
        if isinstance(selection, Collection):
            deuterated = [at.index for at in selection if at.type.name == 'hydrogen']
            
        # Case where |selection| is a function. The function must only have the universe as arguments and
        # return the selected atoms as a list.
        elif callable(selection):
            deuterated = [at.index for at in selection(self.universe) if at.type.name == 'hydrogen'] 

        # Case where |selection| is a list or a tuple.
        elif isinstance(selection, (list, tuple)):
            # Tries to build a MMTK Collection out of it. That triggers that the list or the tuple must contain only 
            # MMTK Atoms.
            try:
                deuterated = [at.index for at in selection if at.type.name == 'hydrogen']

            # Failed to create the MMTK Collection. Raises an error.
            except:
                raise Error('The results of a deuteration selection must be a MMTK Collection.')

        # Remaining case where |selection| is a selection string.
        elif isinstance(selection, str):

            # If the string is equal to 'no', returns an empty collection. No hydrogen will be deuterated.
            if selection.lower() == 'no':
                deuterated = []

            # If the string is equal to 'all', returns all the hydrogen atoms of the universe in a MMTK Collection.
            elif selection.lower() == 'all':
                deuterated = [at.index for at in self.universe.atomList() if at.type.name == 'hydrogen']

            # Otherwise, parses the selection string.
            else:

                # If the selection string starts with "filename name-of-file" then it will be a selection from a file.
                if selection.lower().strip()[0:8] == 'filename':
                    # The format for such a selection is 'filename name-of-the-file'.
                    try:
                        filename = selection[8:].strip()

                    # Bad format. Raises an error.
                    except:
                        raise Error('Invalid format for\n%s\nfile-based selection string.' % selection)

                    # The format is OK. Performs the deuteration selection from |filename| file.
                    else:
                        deuterated = self.__parseFileSelectionString(filename, 'deuteration')[0]

                # If the selection string starts with "expression" then it will be a selection from an valid Python expression 
                # expression.
                elif selection.lower().strip()[0:10] == 'expression':

                    # Checks that the expression assigns somehwere selection variable to something.
                    if not re.findall('selection\s*=',selection):
                        raise Error('Invalid format for\n%s\n expression-based selection string. It must contain "selection = list-of-atoms".' % selection)

                    deuterated = self.__parseExpressionSelectionString(selection, 'deuteration')

                # If the selection string starts with "object" then it will be a selection from the chemical objects 
                # available in the universe.
                elif selection.lower().strip()[0:6] == 'object':
                    deuterated = self.__parseObjectSelectionString(selection[7:].strip())

                # Otherwise, raises an error.
                else:
                    raise Error('The object-based selection string\n%s\n does not start with "object", "filename" or "expression".' % selection)


                # Tries to convert the result of the parsing of the selection string to a MMTK Collection.
                try:
                    # Only keep the hydrogen atoms. 
                    deuterated = [ind for ind in deuterated if self.atomInformation[ind]['type'] == 'hydrogen']

                # If the conversion fails, raises an error.
                except:
                    raise Error('The parsing of the object-based selection string \n%s\nproduced an invalid result.\
There might be something wrong with it.' % selection)
                    
        else:
            raise Error('Wrong format for the deuteration selection. Must be a string or a MMTK Collection')            

        deuterationSelection = N.array(sorted(deuterated), typecode = N.Int32)

        for at in self.universe.atomList():
            if at.index in deuterationSelection:
                self.atomInformation[at.index]['type'] = 'deuterium'
                self.atomInformation[at.index]['element'] = 'D'

        if message is None:
            message = 'Deuteration selection'

        LogMessage('info', '', ['console', 'file'])
        LogMessage('info', 'Deuteration selection', ['console', 'file'])
        LogMessage('info', '-'*len(message), ['console', 'file'])
        LogMessage('info', 'Number of deuterated atoms: %s' % len(deuterationSelection), ['console', 'file'])
        LogMessage('info', '', ['console', 'file'])
        LogMessage('info', '', ['console', 'file'])

        return deuterationSelection

    def selectGroups(self, selection, atomOrder = None, message = None):
        """Returns a list of MMTK collections of atoms where each collection defines a group of atom on which an analysis will be 
        performed collectively.

        @param universe: the universe on which the selection will be performed.
        @type universe: instance of MMTK.Universe

        @param selection: which atoms to select. It can be:
            - a string: in that case it will be parsed and the results will be returned.
            - a list or tuple of MMTK collections: in that case it will be converted to a list of MMTK Collections and returned.
        @type selection: string or list or tuple or a MMTK Collection.

        @return: list of MMTK collection of atoms where each collection defines a group.
        @rtype: list
        """

        # Case where |selection| is a list or a tuple.
        if isinstance(selection, (list, tuple)):

            # Loop over each element of the list/tuple.
            for el in selection:
                # Check that it is actually a MMTK Collection.
                if not isinstance(el, (list,tuple)):
                    # If not, raises an error.
                    raise Error('Wrong format for the group selection. It must be a list or tuple of MMTK Collections or a group selection string.')

                for at in el:
                    if not isinstance(at, Atom):
                        # If not, raises an error.
                        raise Error('Wrong format for the group selection. It must be a nested-list/tuple of MMTK Atoms or a group selection string.')

            # If all the elements are MMTK Collections, then converts the list/tuple to a list and return the result.
            else:
                groupSelection = [[at.index for at in l] for l in selection]
                
        elif callable(selection):
            groupSelection = [[at.index for at in l] for l in selection(self.universe)]

        # Remaining case where |selection| is a selection string.
        elif isinstance(selection, str):

            # If the string is equal to 'all', returns a list of the MMTK Collections corresponding to each chemical
            # object of the universe.
            if selection.lower() == 'all':

                # This list will contain the MMTK Collections corresponding to each chemical object of the universe.
                groupSelection = []

                # Loop over the chemical objects of the universe.
                for obj in self.universe.objectList():

                    # If the object has an internal hierarchy such a protein, a peptide or nucleotide chain,
                    # do not add the MMTK Collections of the whole object but adds separately MMTK Collections 
                    # corresponding to each residue of the object.
                    if isinstance(obj,(PeptideChain, Protein, NucleotideChain)):

                        # Loop over the residues of the object.
                        for res in obj.residues():
                            # And append the corresponding MMTK Collection to |groupSelection| list.
                            groupSelection.append([at.index for at in res.atomList()])

                    # Otherwise, add the MMTK Collection of the whole object to |groupSelection| list.
                    elif isinstance(obj,(Atom, AtomCluster, Molecule)):

                        groupSelection.append([at.index for at in obj.atomList()])

            # Otherwise, parses the selection string.
            else:

                # This list will contain the MMTK Collections corresponding to each chemical object of the universe.
                groupSelection = []

                # If the selection string starts with "filename name-of-file" then it will be a selection from a file.
                if selection.lower().strip()[0:8] == 'filename':
                    # The format for such a selection is 'filename name-of-the-file'.
                    try:
                        filename = selection[8:].strip()

                    # Bad format. Raises an error.
                    except:
                        raise Error('Invalid format for\n%s\nfile-based selection string.' % selection)

                    # The format is OK. Performs the deuteration selection from |filename| file.
                    else:
                        groupSelection = self.__parseFileSelectionString(filename, 'group')

                # If the selection string starts with "expression" then it will be a selection from an valid Python expression 
                # expression.
                elif selection.lower().strip()[0:10] == 'expression':

                    # Checks that the expression assigns somehwere selection variable to something.
                    if not re.findall('selection\s*=',selection):
                        raise Error('Invalid format for\n%s\n expression-based selection string. It must contain "selection = list-of-atoms".' % selection)

                    groupSelection = self.__parseExpressionSelectionString(selection, 'group')

                # If the selection string starts with "object" then it will be a selection from the chemical objects 
                # available in the universe.
                elif selection.lower().strip()[0:6] == 'object':

                    selectionStringsPerGroup = [v.strip() for v in re.split('groupinglevel (\w+)', selection[7:]) if v != '']

                    if len(selectionStringsPerGroup) == 0:
                        raise Error('The object-based selection string\n%s\ncould not be parsed properly.' % selection)

                    if len(selectionStringsPerGroup)%2 != 0:
                        raise Error('The object-based selection string\n%s\ncould not be parsed properly.' % selection)

                    orderedAtoms = sorted(self.universe.atomList(), key = operator.attrgetter('index'))

                    for comp in range(0, len(selectionStringsPerGroup), 2):

                        selectionString = selectionStringsPerGroup[comp]

                        groupingLevel = selectionStringsPerGroup[comp+1]

                        grp = self.__parseObjectSelectionString(selectionString)

                        groupSelection.extend(self.__buildGroup(grp, groupingLevel, orderedAtoms))

                # Otherwise, raises an error.
                else:
                    raise Error('The object-based selection string\n%s\n does not start with "object", "filename" or "expression".' % selection)
                
        else:
            raise Error('Invalid input format for\n%s\ngroup selection. It must be a list/tuple of MMTK Collections of atoms or a selection string.')

        try:
            groupSelection = [sorted(list(g)) for g in groupSelection]

        except:
            raise Error('The parsing of the object-based selection string \n%s\nproduced an invalid result.\
There might be something wrong with it.' % selection)

        if atomOrder is not None:
            orderedGroups = []
            for g in groupSelection:

                if len(g) < len(atomOrder):
                    continue

                temp = []
                for atName in atomOrder:

                    found = [index for index in g if self.atomInformation[index]['name'] == atName]

                    if len(found) != 1:
                        break

                    temp.append(found[0])

                else:
                    orderedGroups.append(temp)

            nSkipped = len(groupSelection) - len(orderedGroups)
            if nSkipped > 0:
                LogMessage('warning', '%d groups were skipped because their contents did not match the atom order.' % nSkipped, ['console','file'])

            groupSelection = orderedGroups

        if message is None:
            message = 'Group selection'
            
        LogMessage('info', '', ['console', 'file'])
        LogMessage('info', message, ['console', 'file'])
        LogMessage('info', '-'*len(message), ['console', 'file'])
        LogMessage('info', 'Number of groups selected: %s' % len(groupSelection), ['console', 'file'])
        LogMessage('info', '', ['console', 'file'])
        LogMessage('info', '', ['console', 'file'])

        return groupSelection

    def __buildGroup(self, g, gLevel, atoms):
        """
        """

        group = Collection([atoms[ind] for ind in g])

        gDict= {}

        if gLevel.lower() == 'chain':
            for at in group:
                chain = at.parent.parent.parent
                gDict[at.index] = chain

        elif gLevel.lower() == 'residue':
            for at in group:
                residue = at.parent.parent
                gDict[at.index] = residue

        elif gLevel.lower() in ['molecule', 'cluster', 'nucleicacid', 'protein']:
            for at in group:
                molecule = at.topLevelChemicalObject()
                gDict[at.index] = molecule

        elif gLevel.lower() == 'atom':
            for at in group:
                gDict[at.index] = at

        elif gLevel.lower() == 'amine':
            for at in group:
                amine = belongToAnAmine(at)
                if amine is not None:
                    gDict[at.index] = amine

        elif gLevel.lower() == 'hydroxy':
            for at in group:
                hydroxy = belongToAHydroxy(at)
                if hydroxy is not None:
                    gDict[at.index] = hydroxy

        elif gLevel.lower() == 'methyl':
            for at in group:
                methyl = belongToAMethyl(at)
                if methyl is not None:
                    gDict[at.index] = methyl

        elif gLevel.lower() == 'thiol':
            for at in group:
                thiol = belongToAnHydroxy(at)
                if thiol is not None:
                    gDict[at.index] = thiol

        elif gLevel.lower() == 'default':
            for at in group:
                obj = at.topLevelChemicalObject()

                if isinstance(obj, (NucleotideChain,PeptideChain,Protein)):
                    residue = at.parent.parent
                    gDict[at.index] = residue

                elif isinstance(obj, (Atom, AtomCluster,Molecule)):
                    gDict[at.index] = obj

        groupSelection = []
        temp = []
        for ind, name in sorted(gDict.items()):
            if name not in temp:
                temp.append(name)
                groupSelection.append([ind])

            else:
                el = temp.index(name)
                groupSelection[el].append(ind)

        return groupSelection

    def __parseExpressionSelectionString(self, selectionString, selectionMode):
        """Parses a nMOLDYN selection string that starts with 'expression' keyword. 

        @param selectionString: the selection string to parse.
        @type selectionString: string

        @param selectionMode: a string being one of 'subset', 'deuteration' or 'group' that will specify which kind of 
            selection should be performed.
        @type selectionMode: string

        @return: will depend on |selectionMode|.
        @rtype: Python set for |selectionMode| = 'subset' or 'deuteration' and dict for |selectionMode| = 'group'
        @return: a set of the atoms that match |selectionString| selection string.
        @rtype: set
        """

        try:
            expression = ' '.join(selectionString.split()[1:])
            exec(expression)

        except:
            raise Error('Could not parse the expression-based selection string.')

        if not isinstance(selection, list):
            raise Error('The final statement of an expression-based selection string must be "selection = list-of-atoms.".')

        if selectionMode in ['subset', 'deuteration']:

            selection = set([at.index for at in selection])

        elif selectionMode == 'group':
            selection = [[at.index for at in g] for g in selection]

        if not selection:
            raise Error('The parsing of \n%s\n expression-based selection string gave an empty selection.\nSomething might be wrong with it.' % selectionString)

        return selection

    def __parseFileSelectionString(self, filename, selectionMode):
        """Parses a nMOLDYN Selection file (nms file) in order to perform a subset, a deuteration or a group 
        selection. 

        @param filename: the selection string to parse.
        @type filename: string

        @param selectionMode: a string being one of 'subset', 'deuteration' or 'group' that will specify which kind of 
            selection should be performed.
        @type selectionMode: string

        @return: will depend on |selectionMode|.
        @rtype: Python set for |selectionMode| = 'subset' or 'deuteration' and dict for |selectionMode| = 'group'
        """
                                        
        # First try to open the nms file that contains the informations about the atom to select.
        try:
            nmsFile = open(filename, 'r')
            exec nmsFile
            nmsFile.close()
        except:
            raise Error('Unable to open the file %s.' % filename)

        try:
            pdb = pdb.strip()
            if not os.path.exists(pdb):
                raise Error('The pdb file %s does not exist.' % pdb)

        except NameError:
            raise Error('The "pdb" field was not set in the %s nms file.' % filename)

        # Sets up the output variable depending on the selected selection type.
        if selectionMode in ['subset', 'deuteration','group']:
            try:
                selectedAtoms = [[int(vv) for vv in v] for v in eval(selectionMode)]
            except:
                raise Error('Wrong format for %s selection field.' % (selectionMode))

            selection = []
            [selection.append(set()) for i in range(len(selectedAtoms))]

        else:
            raise Error('Unknown selection type: %s.' % selectionMode)

        try:
            
            # The PDB file is opened for reading.
            pdb = PDBConfiguration(pdb)
                                
            univCoord = {}
            for at in self.universe.atomList():
                univCoord[tuple([round(v,4) for v in self.trajectory.configuration[frame][at]])] = at.index
                
            atomList = []
            for res in pdb.residues:
                for at in res.atom_list:
                    atomList.append(at.position)
                                                                                                                
            if selectionMode in ['subset', 'deuteration', 'group']:
                
                pdbCoords = []
                for p in selectedAtoms:
                    temp = []
                    for idx in set(p):
                        if idx <= 0 or idx > self.universe.numberOfAtoms(): continue
                        temp.append(tuple([round(v,4) for v in atomList[idx-1]]))
                    pdbCoords.append(temp)

                for i in range(len(pdbCoords)):
                    pdbCoord = pdbCoords[i]
                    coords = set(pdbCoord).intersection(univCoord.keys())
                    for c in coords:                    
                        selection[i].add(univCoord[c])                    
                    if len(selection[i]) != len(selectedAtoms[i]):
                      raise                                                
                
        except:
            raise Error('Could not parse properly for selection the PDB file associated with %s selection file.' % filename)

        if not selection:
            raise Error('The selection from file %s gave an empty selection. Something might be wrong with it.' % filename)

        return selection

    def __parseObjectSelectionString(self, selectionString):
        """Parses a selection string.

        @param selectionString: the selection string to parse.
        @type selectionString: string

        @return: a set of the atoms that match |selectionString| selection string.
        @rtype: set
        """

        selStr = re.sub('\(',' ( ',selectionString)
        selStr = re.sub('\)',' ) ',selStr)

        l = [v.strip() for v in re.split('\s+|,',selStr) if v.strip()]
        
        processedList = []
        selectionValue = []

        try:
            while True:

                if not l:
                    break

                l0 = l[0].lower()
                if l0 == 'objectname':
                    objectName = l[1]
                    del l[0:2]
                    objectClass = self.chemicalObjectInfo[objectName]['objectclass'].lower()
                    continue

                elif l0 in self.chemicalObjectInfo[objectName].keys():
                    selectionKeyword = l0
                    del l[0]

                    while True:

                        ll0 = l[0].lower()
                        if ll0 == 'objectname':
                            raise

                        elif ll0 in self.chemicalObjectInfo[objectName].keys():
                            raise

                        elif ll0 in ['and','or',')']:
                            if not selectionValue:
                                raise
                            indexes = self.__retrieveAtomIndexes(objectClass, objectName, selectionKeyword, selectionValue)
                            processedList.append(str(set(indexes)))
                            selectionValue = []
                            break

                        elif ll0 == '(':
                            raise

                        else:
                            selectionValue.append(ll0)
                            del l[0]
                            if not l:
                                indexes = self.__retrieveAtomIndexes(objectClass, objectName, selectionKeyword, selectionValue)
                                processedList.append(str(set(indexes)))
                                selectionValue = []
                                break

                elif l0 in ['and','or','(']:
                    processedList.append(l0)
                    del l[0]
                    if not l:
                        raise

                elif l0 == ')':
                    processedList.append(l0)
                    del l[0]

                else:
                    raise

        except:
            raise Error('Bad format for an object-based selection string.')

        processedString = ''.join(processedList)
        processedString = processedString.replace('and', '&')
        processedString = processedString.replace('or', '|')

        # The MMTK objects corresponding to the selection indexes are returned.
        selection = eval(processedString)

        if not selection:            
            raise Error('The selection associated with\n%s\nselection string gave an empty selection.\nSomething might be wrong with it.' % selectionString)

        return selection

    def __retrieveAtomIndexes(self, objectClass, objectName, selectionKeyword, selectionValue):
        """Retrieves the MMTK index of the atoms matching |objectClass| MMTK chemical object class, 
           |objectName| nMOLDYN name, |selectionKeyword| selection keyword and |selectionValue| selection value.

           @param objectClass: the MMTK chemical object to match.
           @type objectClass: string

           @param objectName: the nMOLDYN name to match.
           @type objectName: string

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        if objectClass == 'allclass':
            # The atom parser.
            indexes = self.__allClassParser(selectionKeyword, selectionValue)

        elif objectClass == 'atom':
            # The atom parser.
            indexes = self.__atomParser(objectName.lower(), selectionKeyword, selectionValue)

        elif objectClass == 'atomcluster':
            # The atom cluster parser.
            indexes = self.__atomClusterParser(objectName.lower(), selectionKeyword, selectionValue)

        elif objectClass == 'molecule':
            # The molecule parser.
            indexes = self.__moleculeParser(objectName.lower(), selectionKeyword, selectionValue)

        elif objectClass == 'nucleotidechain':
            # The nucleotide chain parser.
            indexes = self.__nucleotideChainParser(objectName.lower(), selectionKeyword, selectionValue)

        elif objectClass == 'peptidechain':
            # The peptide chain parser.
            indexes = self.__peptideChainParser(objectName.lower(), selectionKeyword, selectionValue)

        elif objectClass == 'protein':
            # The protein parser.
            indexes = self.__proteinParser(objectName.lower(), selectionKeyword, selectionValue)

        else:
            raise Error('The chemical class %s is not defined in MMTK.' % objectClass)

        return indexes

    def __allClassParser(self, selectionKeyword, selectionValue):
        """Retrieves the MMTK indexes of the atoms whose nMOLDYN name is '*' matching |selectionKeyword| 
           selection keyword and |selectionValue| selection value.

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        nChemicalObjects = len(self.universe.objectList())
        chemicalObjects = self.universe.objectList()

        selection = []        

        if selectionKeyword == 'atomelement':
            if '*' in selectionValue:
                for comp in range(nChemicalObjects):
                    cObj = chemicalObjects[comp]
                    selection.extend([at.index for at in cObj])
            else:
                selection = []
                for v in selectionValue:
                    for comp in range(nChemicalObjects):
                        cObj = chemicalObjects[comp]
                        selection.extend([at.index for at in cObj.atomList() if at.type.name.strip().lower() == v])

        return selection

    def __atomParser(self, objectName, selectionKeyword, selectionValue):
        """Retrieves the MMTK index of the atoms whose MMTK chemical object class is 'Atom' matching 
           |objectName| nMOLDYN name, |selectionKeyword| selection keyword and |selectionValue| selection value.

           @param objectName: the nMOLDYN name to match.
           @type objectName: string

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        nChemicalObjects = len(self.universe.objectList())
        chemicalObjects = self.universe.objectList()

        selection = []

        selectedObjects = []
        for comp in range(nChemicalObjects):

            cObj = chemicalObjects[comp]
            cObjName = self.chemicalObjectNames[comp].lower()

            if cObjName == objectName:
                selectedObjects.append(cObj)

        if selectionKeyword == 'atomname':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.fullName().strip().lower() == v])

        return selection

    def __atomClusterParser(self, objectName, selectionKeyword, selectionValue):
        """Retrieves the MMTK index of the atoms whose MMTK chemical object class is 'AtomCluster' matching 
           |objectName| nMOLDYN name, |selectionKeyword| selection keyword and |selectionValue| selection value.

           @param objectName: the nMOLDYN name to match.
           @type objectName: string

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        nChemicalObjects = len(self.universe.objectList())
        chemicalObjects = self.universe.objectList()

        selection = []

        selectedObjects = []
        for comp in range(nChemicalObjects):

            cObj = chemicalObjects[comp]
            cObjName = self.chemicalObjectNames[comp].lower()

            if cObjName == objectName:
                selectedObjects.append(cObj)

        if selectionKeyword == 'atomname':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.fullName().strip().lower() == v])

        elif selectionKeyword == 'atomelement':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    if v in ['sulfur','sulphur']:
                        v = ['sulfur','sulphur']
                    else:
                        v = [v]
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.type.name.strip().lower() in v])

        return selection

    def __moleculeParser(self, objectName, selectionKeyword, selectionValue):
        """Retrieves the MMTK index of the atoms whose MMTK chemical object class is 'Molecule' matching 
           |objectName| nMOLDYN name, |selectionKeyword| selection keyword and |selectionValue| selection value.

           @param objectName: the nMOLDYN name to match.
           @type objectName: string

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        nChemicalObjects = len(self.universe.objectList())
        chemicalObjects = self.universe.objectList()

        selection = []

        selectedObjects = []
        for comp in range(nChemicalObjects):

            cObj = chemicalObjects[comp]
            cObjName = self.chemicalObjectNames[comp].lower()

            if cObjName == objectName:
                selectedObjects.append(cObj)

        if selectionKeyword == 'atomname':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.fullName().strip().lower() == v])

        elif selectionKeyword == 'atomelement':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    if v in ['sulfur','sulphur']:
                        v = ['sulfur','sulphur']
                    else:
                        v = [v]
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.type.name.strip().lower() in v])

        elif selectionKeyword == 'chemfragment':
            for v in selectionValue:

                if v == 'amine':
                    for obj in selectedObjects:
                        nitrogens = [at for at in obj.atomList() if at.type.name.strip().lower() == 'nitrogen']
                        for n in nitrogens:
                            neighbours = n.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The amine.
                            if len(hydrogens) == 2:
                                selection.extend([n.index] + [h.index for h in hydrogens])
                            # The ammonium.
                            elif (len(hydrogens) == 3) and (obj.numberOfAtoms() == 4):
                                selection.extend([n.index] + [h.index for h in hydrogens])

                elif v == 'hydroxy':
                    for obj in selectedObjects:
                        oxygens = [at for at in obj.atomList() if at.type.name.strip().lower() == 'oxygen']
                        for o in oxygens:
                            neighbours = o.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The hydroxy.
                            if len(hydrogens) == 1:
                                selection.extend([o.index] + [h.index for h in hydrogens])
                            # The water.
                            elif (len(hydrogens) == 2) and (obj.numberOfAtoms() == 3):
                                selection.extend([o.index] + [h.index for h in hydrogens])

                elif v == 'methyl':                        
                    for obj in selectedObjects:
                        carbons = [at for at in obj.atomList() if at.type.name.strip().lower() == 'carbon']
                        for c in carbons:
                            neighbours = c.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen']
                            # The methyl
                            if len(hydrogens) == 3:
                                selection.extend([c.index] + [h.index for h in hydrogens])
                            # The methane
                            elif (len(hydrogens) == 4) and (obj.numberOfAtoms() == 5):
                                selection.extend([c.index] + [h.index for h in hydrogens])

                elif v == 'thiol':                        
                    for obj in selectedObjects:
                        sulphurs = [at for at in obj.atomList() if at.type.name.strip().lower() in ['sulphur','sulfur']]
                        for s in sulphurs:
                            neighbours = s.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The thiol.
                            if len(hydrogens) == 1:
                                selection.extend([o.index] + [h.index for h in hydrogens])
                            # The SH2. Quite unusual ...
                            elif (len(hydrogens) == 2) and (obj.numberOfAtoms() == 3):
                                selection.extend([s.index] + [h.index for h in hydrogens])

        return selection

    def __nucleotideChainParser(self, objectName, selectionKeyword, selectionValue):
        """Retrieves the MMTK index of the atoms whose MMTK chemical object class is 'NucleotideChain' matching 
           |objectName| nMOLDYN name, |selectionKeyword| selection keyword and |selectionValue| selection value.

           @param objectName: the nMOLDYN name to match.
           @type objectName: string

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        nChemicalObjects = len(self.universe.objectList())
        chemicalObjects = self.universe.objectList()

        selection = []

        selectedObjects = []
        for comp in range(nChemicalObjects):

            cObj = chemicalObjects[comp]
            cObjName = self.chemicalObjectNames[comp].lower()

            if cObjName == objectName:
                selectedObjects.append(cObj)

        if selectionKeyword == 'atomname':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.fullName().strip().lower() == v])

        elif selectionKeyword == 'atomtype':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])

                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.name.strip().lower() == v])

        elif selectionKeyword == 'atomelement':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    if v in ['sulfur','sulphur']:
                        v = ['sulfur','sulphur']
                    else:
                        v = [v]
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.type.name.strip().lower() in v])

        elif selectionKeyword == 'nuclname':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        for nucl in obj.residues():
                            if nucl.fullName().strip().lower() == v:
                                selection.extend([at.index for at in nucl.atomList()])

        elif selectionKeyword == 'nucltype':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        for nucl in obj.residues():
                            if nucl.symbol.strip().lower() == v:
                                selection.extend([at.index for at in nucl.atomList()])

        elif selectionKeyword == 'misc':
            for v in selectionValue:
                if v == 'bases':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.bases().atomList()])

                elif v == 'backbone':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.backbone().atomList()])

        return selection

    def __peptideChainParser(self, objectName, selectionKeyword, selectionValue):
        """Retrieves the MMTK index of the atoms whose MMTK chemical object class is 'PeptideChain' matching 
           |objectName| nMOLDYN name, |selectionKeyword| selection keyword and |selectionValue| selection value.

           @param objectName: the nMOLDYN name to match.
           @type objectName: string

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        nChemicalObjects = len(self.universe.objectList())
        chemicalObjects = self.universe.objectList()

        selection = []

        selectedObjects = []
        for comp in range(nChemicalObjects):

            cObj = chemicalObjects[comp]
            cObjName = self.chemicalObjectNames[comp].lower()

            if cObjName == objectName:
                selectedObjects.append(cObj)

        if selectionKeyword == 'atomname':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.fullName().strip().lower() == v])

        elif selectionKeyword == 'atomtype':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.name.strip().lower() == v])

        elif selectionKeyword == 'atomelement':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    if v in ['sulfur','sulphur']:
                        v = ['sulfur','sulphur']
                    else:
                        v = [v]
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.type.name.strip().lower() in v])

        elif selectionKeyword == 'chemfragment':
            for v in selectionValue:                    
                if v == 'amine':                        
                    for obj in selectedObjects:
                        nitrogens = [at for at in obj.atomList() if at.type.name.strip().lower() == 'nitrogen']
                        for n in nitrogens:
                            neighbours = n.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The amine.
                            if len(hydrogens) == 2:
                                selection.extend([n.index] + [h.index for h in hydrogens])
                            # The ammonium.
                            elif (len(hydrogens) == 3) and (obj.numberOfAtoms() == 4):
                                selection.extend([n.index] + [h.index for h in hydrogens])

                elif v == 'c_alphas':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.name.strip().lower() == 'c_alpha'])

                elif v == 'hydroxy':                        
                    for obj in selectedObjects:
                        oxygens = [at for at in obj.atomList() if at.type.name.strip().lower() == 'oxygen']
                        for o in oxygens:
                            neighbours = o.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The hydroxy.
                            if len(hydrogens) == 1:
                                selection.extend([o.index] + [h.index for h in hydrogens])
                            # The water.
                            elif (len(hydrogens) == 2) and (obj.numberOfAtoms() == 3):
                                selection.extend([o.index] + [h.index for h in hydrogens])

                elif v == 'methyl':                        
                    for obj in selectedObjects:
                        carbons = [at for at in obj.atomList() if at.type.name.strip().lower() == 'carbon']
                        for c in carbons:
                            neighbours = c.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen']
                            # The methyl
                            if len(hydrogens) == 3:
                                selection.extend([c.index] + [h.index for h in hydrogens])
                            # The methane
                            elif (len(hydrogens) == 4) and (obj.numberOfAtoms() == 5):
                                selection.extend([c.index] + [h.index for h in hydrogens])

                elif v == 'thiol':                        
                    for obj in selectedObjects:
                        sulphurs = [at for at in obj.atomList() if at.type.name.strip().lower() in ['sulphur', 'sulfur']]
                        for s in sulphurs:
                            neighbours = s.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The thiol.
                            if len(hydrogens) == 1:
                                selection.extend([o.index] + [h.index for h in hydrogens])
                            # The SH2. Quite unusual ...
                            elif (len(hydrogens) == 2) and (obj.numberOfAtoms() == 3):
                                selection.extend([s.index] + [h.index for h in hydrogens])

        elif selectionKeyword == 'resname':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        for res in obj.residues():
                            if res.fullName().strip().lower() == v:
                                selection.extend([at.index for at in res.atomList()])

        elif selectionKeyword == 'restype':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        for res in obj.residues():
                            if res.symbol.strip().lower() == v:
                                selection.extend([at.index for at in res.atomList()])

        elif selectionKeyword == 'resclass':
            for v in selectionValue:
                for obj in selectedObjects:
                    for res in obj.residues():
                        if res.symbol.strip().lower() in residusChemFamily[v]:
                            selection.extend([at.index for at in res.atomList()])

        elif selectionKeyword == 'misc':
            for v in selectionValue:
                if v == 'sidechains':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.sidechains().atomList()])

                elif v == 'backbone':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.backbone().atomList()])

        return selection                

    def __proteinParser(self, objectName, selectionKeyword, selectionValue):
        """Retrieves the MMTK index of the atoms whose MMTK chemical object class is 'Protein' matching 
           |objectName| nMOLDYN name, |selectionKeyword| selection keyword and |selectionValue| selection value.

           @param objectName: the nMOLDYN name to match.
           @type objectName: string

           @param selectionKeyword: the selection keyword to match.
           @type selectionKeyword: string

           @param selectionValue: the selection value to match.
           @type selectionValue: string

           @return: a list of MMTK indexes of the selected atoms.
           @rtype: list
        """

        nChemicalObjects = len(self.universe.objectList())
        chemicalObjects = self.universe.objectList()

        selection = []

        selectedObjects = []
        for comp in range(nChemicalObjects):

            cObj = chemicalObjects[comp]
            cObjName = self.chemicalObjectNames[comp].lower()

            if cObjName == objectName:
                selectedObjects.append(cObj)

        if selectionKeyword == 'atomname':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.fullName().strip().lower() == v])

        elif selectionKeyword == 'atomtype':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.name.strip().lower() == v])

        elif selectionKeyword == 'atomelement':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    if v in ['sulfur','sulphur']:
                        v = ['sulfur','sulphur']
                    else:
                        v = [v]
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.type.name.strip().lower() in v])

        elif selectionKeyword == 'chemfragment':
            for v in selectionValue:

                if v == 'amine':                        
                    for obj in selectedObjects:
                        nitrogens = [at for at in obj.atomList() if at.type.name.strip().lower() == 'nitrogen']
                        for n in nitrogens:
                            neighbours = n.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The amine.
                            if len(hydrogens) == 2:
                                selection.extend([n.index] + [h.index for h in hydrogens])
                            # The ammonium.
                            elif (len(hydrogens) == 3) and (obj.numberOfAtoms() == 4):
                                selection.extend([n.index] + [h.index for h in hydrogens])

                elif v == 'c_alphas':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList() if at.name.strip().lower() == 'c_alpha'])

                elif v == 'hydroxy':                        
                    for obj in selectedObjects:
                        oxygens = [at for at in obj.atomList() if at.type.name.strip().lower() == 'oxygen']
                        for o in oxygens:
                            neighbours = o.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The hydroxy.
                            if len(hydrogens) == 1:
                                selection.extend([o.index] + [h.index for h in hydrogens])
                            # The water.
                            elif (len(hydrogens) == 2) and (obj.numberOfAtoms() == 3):
                                selection.extend([o.index] + [h.index for h in hydrogens])

                elif v == 'methyl':                        
                    for obj in selectedObjects:
                        carbons = [at for at in obj.atomList() if at.type.name.strip().lower() == 'carbon']
                        for c in carbons:
                            neighbours = c.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen']
                            # The methyl
                            if len(hydrogens) == 3:
                                selection.extend([c.index] + [h.index for h in hydrogens])
                            # The methane
                            elif (len(hydrogens) == 4) and (obj.numberOfAtoms() == 5):
                                selection.extend([c.index] + [h.index for h in hydrogens])

                elif v == 'thiol':                        
                    for obj in selectedObjects:
                        sulphurs = [at for at in obj.atomList() if at.type.name.strip().lower() in ['sulphur','sulfur']]
                        for s in sulphurs:
                            neighbours = s.bondedTo()
                            hydrogens = [neigh for neigh in neighbours if neigh.type.name.strip().lower() == 'hydrogen'] 
                            # The thiol.
                            if len(hydrogens) == 1:
                                selection.extend([o.index] + [h.index for h in hydrogens])
                            # The SH2. Quite unusual ...
                            elif (len(hydrogens) == 2) and (obj.numberOfAtoms() == 3):
                                selection.extend([s.index] + [h.index for h in hydrogens])

        elif selectionKeyword == 'resname':                
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        for res in obj.residues():
                            if res.fullName().strip().lower() == v:
                                selection.extend([at.index for at in res.atomList()])

        elif selectionKeyword == 'restype':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        for res in obj.residues():
                            if res.symbol.strip().lower() == v:
                                selection.extend([at.index for at in res.atomList()])

        elif selectionKeyword == 'resclass':
            for v in selectionValue:
                for obj in selectedObjects:
                    for res in obj.residues():
                        if res.symbol.strip().capitalize() in residusChemFamily[v]:
                            selection.extend([at.index for at in res.atomList()])

        elif selectionKeyword == 'chainname':
            for v in selectionValue:
                if v == '*':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.atomList()])
                else:
                    for obj in selectedObjects:
                        for chain in obj:
                            if chain.name.strip().lower() == v:
                                selection.extend([at.index for at in chain.atomList()])

        elif selectionKeyword == 'misc':
            for v in selectionValue:
                if v == 'sidechains':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.sidechains().atomList()])

                elif v == 'backbone':
                    for obj in selectedObjects:
                        selection.extend([at.index for at in obj.backbone().atomList()])

        return selection

    def test(self, selectedTests = None):
        
        nTests = len(glob.glob(os.path.join(GVAR['nmoldyn_tests'], self.db_shortname,"*_Reference.py")))
        
        if isinstance(selectedTests, (list, tuple)):
            selectedTests = [t for t in sorted(selectedTests) if t > 0 and t <= nTests]
        else:
            selectedTests = range(1, nTests + 1)
                    
        tests = TestSuite([StabilityTests(self.db_shortname, self.db_longname, self.pmoldyn_arg, t) for t in selectedTests])
        
        testResults = StabilityTestsResults()
               
        start = default_timer()
                        
        tests.run(testResults)
        
        testResults.duration = default_timer() - start
        
        return testResults                 
        