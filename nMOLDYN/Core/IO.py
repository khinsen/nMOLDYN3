"""This module implements IO-related classes and procedures.

Classes:
    * EndOfFile                : an empty dummy class used by |DCDReader|.
    * FortranBinaryFile        : sets up a binary file reader.
    * DCDFile                  : sets up a DCD file reader.
    * AmberNetCDFConverter     : converts a trajectory from Amber > 9 to a MMTK NetCDF trajectory.
    * CASTEPConverter          : converts a trajectory from CASTEP to a MMTK NetCDF trajectory.
    * CHARMMConverter          : converts a trajectory from CHARMM to a MMTK NetCDF trajectory.
    * DL_POLYConverter         : converts a trajectory from DL_POLY > 9 to a MMTK NetCDF trajectory.
    * LAMMPSConverter          : converts a trajectory from LAMMPS to a MMTK NetCDF trajectory.
    * MaterialsStudioConverter : converts a trajectory from MaterialsStudio > 9 to a MMTK NetCDF trajectory.
    * NAMDConverter            : converts a trajectory from NAMD to a MMTK NetCDF trajectory.
    * VASP4Converter           : converts a trajectory from VASP4 to a MMTK NetCDF trajectory.
    * VASP5Converter           : converts a trajectory from VASP5 to a MMTK NetCDF trajectory.
    * VASPBackConverter        : converts back a trajectory from MMTK to VASP > 9.
    * PDBConverter             : converts a single PDB frame to a MMTK NetCDF trajectory.

Procedures:
    * convertNetCDFToASCII: converts a NetCDF file into an ASCII file.
    * convertASCIIToNetCDF: converts an ASCII file into a NetCDF file.
"""

# The python distribution modules
from copy import copy
import operator
import os
import re
import struct
from struct import calcsize
import string
import subprocess
import sys
from tempfile import mktemp
from time import asctime, ctime, localtime, strftime, time
import xml.dom.minidom

import numpy

# The ScientificPython modules
from Scientific import N
from Scientific.IO.FortranFormat import FortranFormat, FortranLine
from Scientific.Geometry import Vector
from Scientific.IO.NetCDF import _NetCDFFile, NetCDFFile
from Scientific.IO.TextFile import TextFile

# The MMTK distribution modules
from MMTK import Atom, AtomCluster
from MMTK import Units
from MMTK.ParticleProperties import Configuration, ParticleVector
from MMTK.PDB import PDBConfiguration
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
from MMTK.Universe import InfiniteUniverse, OrthorhombicPeriodicUniverse, ParallelepipedicPeriodicUniverse

# The nMOLDYN modules
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Chemistry.Chemistry import M_to_Symbol, Symbols
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Mathematics.Geometry import basisVectors

def load_trajectory_file(filename):
                    
    # filename must be a path.
    if not isinstance(filename, (str, unicode)):
        return None
    
    # The path is made absolute and normalized.
    filename = os.path.normpath(os.path.abspath(filename))
    
    # If the file does not exists. Return None.
    if not os.path.exists(filename):
        return None
        
    # The extension of the file.
    ext = os.path.splitext(filename)[1].strip().lower()

    netcdf = None
    
    # Case where the extension is ncs. Load the file as a trajectory set.
    if ext == '.ncs':    

        trajSetFile = open(filename, 'r')
        trajSet = [t.strip() for t in trajSetFile.readlines()]
        trajSetFile.close()
        
        # Checks that all the trajectory file of the set exsists.
        if not all([os.path.exists(f) for f in trajSet]):
            return None
        
        # Try to load the trajectory set. Return the loaded trajectory is OK.
        try:
            netcdf = TrajectorySet(None, trajSet)
            
        except:
            netcdf = None
        
    # Case where the extension is nc or cdf. Load the file as trajectory of a NetCDF file.
    else:
        
        try:
            netcdf = Trajectory(None, filename, 'r')
            
        except:
            netcdf = None
        
    return netcdf

def load_netcdf_file(filename):

    # filename must be a path.
    if not isinstance(filename, (str, unicode)):
        return None
    
    # The path is made absolute and normalized.
    filename = os.path.normpath(os.path.abspath(filename))
    
    # If the file does not exists. Return None.
    if not os.path.exists(filename):
        return None
        
    # The extension of the file.
    ext = os.path.splitext(filename)[1].strip().lower()
    
    try:
        netcdf = NetCDFFile(filename, 'r')
        
    except:
        netcdf = None
        
    return netcdf
    
def ASCII_To_NetCDF(inputfile, outputfile, twod = False, varname = None, dimname = None, comments = "#", \
                    delimiter = None, skiprows = 0, bycolumn = False, append = False):

    try:
        datas = numpy.loadtxt(inputfile, comments = comments, skiprows = skiprows, delimiter = delimiter)
    except:
        raise Error("Could not parse the file properly.")

    if len(datas.shape) == 1:
        datas = datas[N.NewAxis,:]

    if bycolumn:
        datas = N.transpose(datas)

    if append:
        output = NetCDFFile(outputfile, "a")        
    else:        
        output = NetCDFFile(outputfile, "w")

    if twod:

        if instance(dimname, (list, tuple)):
            if len(dimname) != 2:
                raise Error("Wrong dimension names")        
        else:        
            dimname = ["dim1", "dim2"]
            
        dimname = tuple(dimname)
              
        if varname is None:
            varname = ["var"]
       
        output.createDimension(dimname[0], datas.shape[0])
        output.createDimension(dimname[1], datas.shape[1])
                    
        VAR = output.createVariable(varname[0], N.Float, dimname)

        VAR[:,:] = datas
        

    else:
        
        try:
            if len(varname) != datas.shape[0]:
                raise
        except:
            varname = ["var"+str(v+1) for v in range(datas.shape[0])]
            
        if not dimname:
            dimname = 'dim'        

        output.createDimension(dimname, datas.shape[1])

        for i in range(len(varname)):
            VAR = output.createVariable(varname[i], N.Float, (dimname,))
            VAR[:] = datas[i,:]

    output.close()

class EndOfFile(Exception):
    """A subclass of Exception.
    """
    pass

class FortranBinaryFile(object):
    """Sets up a Fortran binary file reader. 

    @note: written by Konrad Hinsen.
    """
    def __init__(self, filename, byte_order = '='):
        """The constructor.

        @param filename: the input file.
        @type filename: string.

        @param byte_order: the byte order to read the binary file.
        @type byte_order: string being one '@', '=', '<', '>' or '!'.
        """
        self.file = file(filename, 'rb')
        self.byte_order = byte_order

    def __iter__(self):
        return self

    def next(self):
        data = self.file.read(4)
        if not data:
            raise StopIteration
        reclen = struct.unpack(self.byte_order + 'i', data)[0]
        data = self.file.read(reclen)
        reclen2 = struct.unpack(self.byte_order + 'i', self.file.read(4))[0]
        assert reclen==reclen2
        return data

    def skipRecord(self):
        data = self.file.read(4)
        reclen = struct.unpack(self.byte_order + 'i', data)[0]
        self.file.seek(reclen, 1)
        reclen2 = struct.unpack(self.byte_order + 'i', self.file.read(4))[0]
        assert reclen==reclen2

    def getRecord(self, format, repeat = False):
        """Reads a record of the binary file.

        @param format: the format corresponding to the binray structure to read.
        @type format: string.        

        @param repeat: if True, will repeat the reading.
        @type repeat: bool.        
        """

        try:
            data = self.next()
        except StopIteration:
            raise EndOfFile()
        if repeat:
            unit = struct.calcsize(self.byte_order + format)
            assert len(data) % unit == 0
            format = (len(data)/unit) * format
        try:
            return struct.unpack(self.byte_order + format, data)
        except:
            raise

class DCDFile(object):
    """Sets up a DCD file reader.

    @note: written by Konrad Hinsen.
    """

    def __init__(self, dcd_filename):
        """The constructor.

        @param dcd_filename: the name of the DCD file to read.
        @type dcd_filename: string.
        """
        # Identity the byte order of the file by trial-and-error
        self.byte_order = None
        data = file(dcd_filename, 'rb').read(4)
        for byte_order in ['<', '>']:
            reclen = struct.unpack(byte_order + 'i', data)[0]
            if reclen == 84:
                self.byte_order = byte_order
                break
        if self.byte_order is None:
            raise IOError("%s is not a DCD file" % dcd_filename)
        # Open the file
        self.binary = FortranBinaryFile(dcd_filename, self.byte_order)
        # Read the header information
        header_data = self.binary.next()
        if header_data[:4] != 'CORD':
            raise IOError("%s is not a DCD file" % dcd_filename)
        self.header = struct.unpack(self.byte_order + '9id9i', header_data[4:])
        self.nset = self.header[0]
        self.istart = self.header[1]
        self.nsavc = self.header[2]
        self.namnf = self.header[8]
        self.charmm_version = self.header[-1]
        self.has_pbc_data = False
        self.has_4d = False
        if self.charmm_version != 0:
            self.header = struct.unpack(self.byte_order + '9if10i',
                                        header_data[4:])
            if self.header[10] != 0:
                self.has_pbc_data = True
            if self.header[11] != 0:
                self.has_4d = True
        self.delta = self.header[9]*Units.akma_time
        # Read the title
        title_data = self.binary.next()
        nlines = struct.unpack(self.byte_order + 'i', title_data[:4])[0]
        assert len(title_data) == 80*nlines+4
        title_data = title_data[4:]
        title = []
        for i in range(nlines):
            title.append(title_data[:80].rstrip())
            title_data = title_data[80:]
        self.title = '\n'.join(title)
        # Read the number of atoms.
        self.natoms = self.binary.getRecord('i')[0]
        # Stop if there are fixed atoms.
        if self.namnf > 0:
            raise Error('NAMD converter can not handle fixed atoms yet')

    def readStep(self):
        """Reads a frame of the DCD file.
        """
        if self.has_pbc_data:
            unit_cell = N.array(self.binary.getRecord('6d'), typecode = N.Float)
            a, gamma, b, beta, alpha, c = unit_cell
            if -1. < alpha < 1. and -1. < beta < 1. and -1. < gamma < 1.:
                # assume the angles are stored as cosines
                # (CHARMM, NAMD > 2.5)
                alpha = 0.5*N.pi-N.arcsin(alpha)
                beta = 0.5*N.pi-N.arcsin(beta)
                gamma = 0.5*N.pi-N.arcsin(gamma)
            else:
                # assume the angles are stored in degrees (NAMD <= 2.5)
                alpha *= Units.deg
                beta *= Units.deg
                gamma *= Units.deg
            unit_cell = (a*Units.Ang, b*Units.Ang, c*Units.Ang,
                         alpha, beta, gamma)
        else:
            unit_cell = None
        format = '%df' % self.natoms
        x = N.array(self.binary.getRecord(format), typecode = N.Float32)*Units.Ang
        y = N.array(self.binary.getRecord(format), typecode = N.Float32)*Units.Ang
        z = N.array(self.binary.getRecord(format), typecode = N.Float32)*Units.Ang
        if self.has_4d:
            self.binary.skipRecord()
        return unit_cell, x, y, z

    def skipStep(self):
        """Skips a frame of the DCD file.
        """
        nrecords = 3
        if self.has_pbc_data:
            nrecords += 1
        if self.has_4d:
            nrecords += 1
        for i in range(nrecords):
            self.binary.skipRecord()

    def __iter__(self):
        return self

    def next(self):
        try:
            return self.readStep()
        except EndOfFile:
            raise StopIteration

def convertNetCDFToASCII(inputFile, outputFile, variables, floatPrecision = 9, doublePrecision = 17):
    """Converts a file in NetCDF format to a file in ASCII/CDL format.

    @param inputFile: the name of the NetCDF input file.
    @type inputFile: string

    @param outputFile: the name of the CDL output file.
    @type outputFile: string

    @param variables: the NetCDF variables names (string) to extract when performing the conversion.
    @type variables: list

    @param floatPrecision: the precision on the float numbers.
    @type floatPrecision: integer

    @param doublePrecision: the precision on the double numbers.
    @type doublePrecision: integer

    @note: this is a wrapper for ncdump program of the NetCDF library.
    """

    # If the input file does not exist, raises an error.
    if not os.path.exists(inputFile):
        raise Error('The NetCDF input file name does not exist.')

    # If the path to the ncdump is not set, raises an error.
    if not os.path.exists(PREFERENCES['ncdump_path']):
        LogMessage('warning', 'The path for ncdump program is not correctly set: \
the file %s will not be converted in ASCII.' % inputFile, ['console'])
        return

    # If the output file is not defined, raises an error.
    if not outputFile:
        raise Error('No ASCII output file name.')

    # If the |variables| argument is not a list or a tuple, raises an error.
    if not isinstance(variables,(list,tuple)):
        raise Error('No NetCDF variables list given for conversion.')

    # If the list is empty, raises an error.
    if not variables:
        raise Error('Empty NetCDF variables list given for conversion.')

    # Joins the variables list with ',' separator.
    variablesList = ','.join(variables)

    # Calls the ncdump program for the conversion.
    try:
        output = open(outputFile, 'w')
        s = subprocess.call([PREFERENCES['ncdump_path'], '-b', 'c', '-v',\
                             variablesList, '-p', '%d,%d' % (floatPrecision,doublePrecision),\
                             inputFile], stdout = output)

        output.close()

        # Something went wrong during the conversion so raises an error.
        if s:
            raise

    except:    
        raise Error('Conversion failed.')

def convertASCIIToNetCDF(inputFile, outputFile):
    """Converts a file in ASCII format to a file in NetCDF format.

    If the input file extension is '.cdl', the converter will consider the file to be a CDL file and will use
    the ncgen program of the NetCDF library to perform the conversion.

    @param inputFile: the name of the NetCDF input file.
    @type inputFile: string

    @param outputFile: the name of the CDL output file.
    @type outputFile: string
    """

    # If the input file does not exist, raises an error.
    if not os.path.exists(inputFile):
        raise Error('The ASCII input file name does not exist.')

    # If the path to the ncgen is not set, raises an error.
    if not os.path.exists(PREFERENCES['ncgen_path']):
        LogMessage('warning', 'The path for ncgen program is not correctly set: \
the file %s will not be converted in NetCDF.' % inputFile, ['console'])
        return

    # If the input file extension is '.cdl', uses ncgen to perform the conversion.
    if inputFile[-4:] == '.cdl':
        try:
            subprocess.call([PREFERENCES['ncgen_path'], '-x', '-o', outputFile, inputFile])

        # Something went wrong during the conversion so raises an error.
        except:
            raise Error('Conversion failed.')

    # Otherwise performs the conversion using the nMOLDYN internal converter.
    else:

        # The base name of the input file.
        baseName = os.path.splitext(os.path.basename(inputFile))[0]

        # Open the input file for reading.
        try:
            ascii = open(inputFile, 'r')

        # The file could not be opened for reading.
        except IOError:
            raise Error('Can not open for reading %s ASCII file.' % inputFile)

        else:
            # Loads the files into |asciiContents| list.
            asciiContents = ascii.readlines()
            ascii.close()

            # Remove the empty lines of the ASCII file.
            asciiContents = [d.strip() for d in asciiContents if d.strip() != '']

            # This dictionnary will store the global attributes declared in the input file. Its keys are
            # the names of the global attribute and its values their corresponding values.
            globalAttributes = {}

            # This 'comment' NetCDF global attribute is set up by default.
            globalAttributes['comment'] = ''

            # A list of dictionnary where each dictionnary will define a NetCDF variable.
            variables = []

            # This list will contain nested lists where each nested list will stored the values of each numeric line.
            numLines = []

            # Loop over the non-empty lines of the ascii file.
            for line in asciiContents:

                # If the line starts with '#', it can be either for global variable declaration, either for 
                # a NetCDF variables declaration either for a comment.
                if line[0] == '#':

                    # Check whether this is a NetCDF global variable declaration line.
                    # The format for a global attribute declaration is '# global name_of_the_variable = value_of_the_variable'
                    g = re.findall('#\s*global(.*=.*)',line.lower())
                    # That is the case.
                    if g:
                        try:

                            # The name (|name|) and the value (|value|) of the global attribute.
                            name, value = g[0].split('=')

                            # If so, updates the |globalAttributes| dictionnary with the name/value entry.
                            globalAttributes[name.strip()] = value.strip()

                            # Parse the next line.
                            continue

                        except:
                            raise Error('The line %s can not be parsed.' % line)

                    # Check whether this is a NetCDF variable declaration line.
                    # The format for a NetCDF variable declaration line is 
                    # '# variable name = name_of_the_variable ; name_of_attribute_1 = value_of_attribute_1 ; ...'
                    # where attribute1 ... will be NetCDF attribute for the NetCDF variable.
                    v = re.findall('#\s*variable(.*)',line.lower())
                    # That is the case.
                    if v:
                        try:
                            # The couples 'name = name_of_the_variable', 'name_of_attribute_1 = value_of_attribute_1' ... are extracted.
                            # and stored into a dictionnary whose keys are 'name', 'name_of_attribute_1' ... 
                            # and values are 'name_of_the_variable', 'value_of_attribute_1' ...
                            d = dict(re.findall('(\w+)\s*=\s*(\w+)',v[0]))

                            # The dictionnary must contain the 'name' key.
                            if not d.has_key('name'):
                                raise

                            # The dictionnary for the declared NetCDF variable is appended to |variables| list.
                            variables.append(d)

                            # Parse the next line.
                            continue

                        except:
                            raise Error('The line %s can not be parsed.' % line)

                    # If a line starting with '#' is not a NetCDF global variable declaration line neither 
                    # a NetCDF variable declaration line, then it will be interpredted as a comment and added to 
                    # the 'comment' global attribute.
                    globalAttributes['comment'] += line[1:] + '\n'

                # Otherwise the line must contain ',', ';' or white space(s) separated numbers.
                else:
                    # The numeric values of the line are stored in a list that is appended to |numLines| list. 
                    numLines.append(re.split('[;,\s]+',line.strip()))

            # The transposition of |numLines| is stored in |data| so that now, each nested list of |data| 
            # represent a numeric column of the ascii file.
            data = [[r[col].strip() for r in numLines] for col in range(len(numLines[0]))]

            # If the number of declared NetCDF variable is null or not the same that the number of columns, then 
            # (re)assign an arbitrary variable name to each column.
            if len(variables) != len(data):
                variables = [{'name' : 'column'+str(i)} for i in range(1, len(data)+1)]

            # Retrieves the owner of the job using OS dependent modules.
            try:
                import win32api
                owner = win32api.GetUserName()

            except ImportError:
                from pwd import getpwuid
                from os import getuid
                owner = getpwuid(getuid())[0]

            # pid = the pid number for the job
            pid = os.getpid()

            # This creates a temporary CDL file that will be further converted to a NetCDF file using ncgen
            # with a uniue suffix.
            suffix = '.'.join(['',owner, str(pid), 'moldyn'])            
            cdlFile = mktemp(suffix)

            # The CDL file is opened for writing.
            file  = open(cdlFile, 'w')

            file.write('netcdf %s {\n' % baseName)

            # The NetCDF dimensions declaration section.
            file.write('dimensions:\n')
            file.write('\tnvalues = %d ;\n' % len(data[0]))

            # The NetCDF variables declaration section.
            file.write('\nvariables:\n')
            for var in variables:
                file.write('\tdouble %s(nvalues) ;\n' % var['name'])
                for k, v in var.items():
                    if k == 'name':
                        continue
                    file.write('\t%s:%s = "%s";\n' % (var['name'],k,v))

            # The NetCDF global attribute declaration section.
            file.write('\n// global attributes:\n')
            for k, v in globalAttributes.items():
                file.write('\t:%s = "%s";\n' % (k,v))

            # The NetCDF data declaration section.
            file.write('\ndata:\n\n')

            for i in range(len(variables)):
                file.write(' %s =\n    ' % variables[i]['name'])
                toFlush = '    '
                for val in data[i]:
                    if len(toFlush + val + ', ') > 80:
                        file.write(toFlush.strip() + '\n')
                        toFlush = '    ' + val + ', '
                    else:
                        toFlush += val + ', '

                toFlush = toFlush.strip()
                if toFlush[-1] == ',':
                    toFlush = toFlush[:-1]

                file.write(toFlush + ' ;\n\n')

            file.write('}')            

            # The CDL file is closed.
            file.close()

            # ncgen is called to performed the conversion.
            try:
                subprocess.call([PREFERENCES['ncgen_path'], '-x', '-o', outputFile, cdlFile])

            # Something went wrong during the conversion so raises an error.
            except:
                raise Error('Conversion failed.')

            os.unlink(cdlFile)

class AmberNetCDFConverter(object):
    """Converts an Amber NetCDF Trajectory into a MMTK NetCDFFile. 

    @note: this code is an improved version of the original converter written by Paolo Calligari.
    """

    def __init__(self, pdbFile, amberNetCDFFile, outputFile, timeStep = 1.0):
        """The constructor.

        @param pdbFile: the Amber PDB file name of one frame.
        @type pdbFile: string

        @param amberNetCDFFile: the Amber NetCDF file name of the trajectory to convert.
        @type amberNetCDFFile: string

        @param outputFile: the MMTK NetCDF output file name of the trajectory to convert.
        @type outputFile: string

        @param timeStep: the MD timestep. Default to 1 ps.
        @type timeStep: float.

        @note: calling the constructor will actually perform the conversion.

        @note: the user must provide the MD time step as it is not stored in the Amber NetCDF input file.
        """

        # The arguments are copied to instance attributes.
        self.pdbFile = pdbFile
        self.amberNetCDFFile = amberNetCDFFile
        self.outputFile = outputFile
        self.timeStep = timeStep

        # Do the conversion.
        self.__convert()

    def __convert(self):
        """Performs the actual conversion.
        """

        # The Amber NetCDF input file is opened for reading.
        try:
            netcdfInputFile = NetCDFFile(self.amberNetCDFFile, 'r')

        # Could not be opened. Raises an error.
        except:
            raise Error('The Amber NetCDF file %s could not be opened for reading.' % self.amberNetCDFFile)

        # The Amber NetCDF could be opened.
        else:

            # The number of frames of the trajectory.
            nSteps = len(netcdfInputFile.variables['time'])

            # Case of a periodic universe.
            if netcdfInputFile.variables.has_key('cell_lengths'):
                # The box dimensions are stored in |box| and converted into nanometers.
                box = netcdfInputFile.variables['cell_lengths']

                # The box angles are stored in |angles| and converted into radians.
                angles = netcdfInputFile.variables['cell_angles']

                # The cell parameters a, b, c, alpha, beta and gamma are stored in |cellParams| list.
                cellParams = list(box[0]*Units.Ang) + list(angles[0]*Units.deg)

                # A MMTK periodic universe is created out of this parameters.
                universe = ParallelepipedicPeriodicUniverse(basisVectors(cellParams))

            # Case of a non periodic universe.
            else:
                # A MMTK infinite universe is created.
                universe = InfiniteUniverse()

        # The PDB input file is opened for reading.
        try:
            # The full contents of the PDB file.
            conf = PDBConfiguration(self.pdbFile)

        # The file could not be opened. Raises an error.
        except:
            raise Error('The PDB file %s could not be opened for reading.' % self.pdbFile)

        # The file could be opened.
        else:

            # |molecules| is a Collection of all objects contained in the PDB file.
            molecules = conf.createAll()

            # The objects are introduced into the universe.
            universe.addObject(molecules)

        # A MMTK trajectory is opened for writing.
        trajectory = Trajectory(universe, self.outputFile, 'w', 'Converted from AMBER NetCDF file.')

        # A frame generator is created.
        snapshot = SnapshotGenerator(universe, actions=[TrajectoryOutput(trajectory, ["all"], 0, None, 1)])

        # The total number of atoms of the universe.
        num = molecules.numberOfAtoms()

        # Loop over the trajectory frames.
        for i in range(nSteps):

            # The time corresponding to the frame.
            t = float(i+1)*self.timeStep

            # Case of a periodic universe. The cell parameters for the current frame are calculated.
            if universe.is_periodic:
                cellParams = list(box[i]*Units.Ang) + list(angles[i]*Units.deg)
                b = basisVectors(cellParams)
                cell = N.zeros((9,), typecode = N.Float)
                cell[0:3] = b[0]
                cell[3:6] = b[1]
                cell[6:9] = b[2]

            # Otherwise, the cell parameters are set to None.
            else:
                cell = None

            # The coordinates are converted in nanometers (they are stored in Angstrom in Amber).
            coordinates = netcdfInputFile.variables['coordinates'][i][:num]*Units.Ang

            # The universe configuration is updated with the configuration of the current frame.
            universe.setConfiguration(Configuration(universe, coordinates, cell))

            # A call to the snapshot generator produces the step corresponding to the current frame.
            snapshot(data = {'time': t})

        # The MMTK trajectory file is closed.
        trajectory.close()

        # The Amber NetCDF file is closed.
        netcdfInputFile.close()

class CASTEPConverter(object):
    """Converts a CASTEP Trajectory into a MMTK NetCDFFile. 
    """

    def __init__(self, castep_file, output_file):
        """The constructor.

        @param castep_file: the CASTEP file name of the trajectory to convert.
        @type castep_file: string

        @param output_file: the MMTK NetCDF output file name of the trajectory to convert.
        @type output_file: string

        @note: calling the constructor will actually perform the conversion.
        """

        # The arguments are copied to instance attributes.
        self.castep_file = castep_file
        self.output_file = output_file

        # Do the conversion.
        self.__convert()

    def __convert(self):
        """Performs the actual conversion.
        """
        
        # The Amber NetCDF input file is opened for reading.
        try:
            castep = open(self.castep_file, 'r')

        # Could not be opened. Raises an error.
        except:
            raise Error('The CASTEP file %s could not be opened for reading.' % self.castep_file)

        # The time step. Will be extracted from the CASTEP file.
        time_step   = None
        
        # A list of the a, b and c vectors. Will be extracted from the CASTEP file.
        abc_vectors = []
        
        # The number of atoms. Will be extracted from the CASTEP file.
        n_atoms     = None

       # Start to read the CASTEP file until the first frame in order to get the time step, 
       # the a, b, c vectors and the number of atoms.
        while True:

            # A line is read.
            line = castep.readline().lower()

            # If the line contains 'time step', this is the one where the time step is defined.
            if 'time step' in line:
                # The time step is extracted.
                time_step = float(re.findall('time step\s*:\s*(.*)\s*ps', line)[0])
                continue
        
            # If the line contains 'Real Lattice', the three next lines will contains the direct matrix.
            if 'real lattice' in line:
                # Resp. a, b and c vectors are extracted from the three next lines.
                for i in range(3):
                    line = castep.readline().split()[:3]
                    abc_vectors.append(Units.Ang*Vector([float(v) for v in line]))
                continue
            
            # If the line contains 'Total number of ions in cell', this is the one where the 
            # number of atoms is defined.
            if 'total number of ions in cell' in line:
                n_atoms = int(re.findall('total number of ions in cell\s*=\s*(.*)', line)[0])

            # If the line contains 'fractional coordinates of atoms', this is beginning of 
            # a frame. Stop reading the header.
            if 'fractional coordinates of atoms' in line:
                break
    
            # If the EOF was already reached here, the file is corrupted because no frame
            # could be found. 
            if not line:
                raise Error('The file %s does not contains any frame.' % self.castep_file)

        # If no time step could be found, use a default one of 1 fs.
        if time_step is None:
            time_step = 0.001
            
        # There must be a, b and c vectors in a CASTEP file.
        if not abc_vectors:
            raise Error('The file %s does not contains any PBC information.' % self.castep_file)

        # The number of atoms must be defined in a CASTEP file.
        if n_atoms is None:
            raise Error('The number of atoms is not defined in the CASTEP file %s.' % self.castep_file)

        # The universe that will store the system is created.
        universe = ParallelepipedicPeriodicUniverse(abc_vectors, None)

        # The two next lines following the frame declaration are not useful.
        line = castep.readline()
        line = castep.readline()

        # The n_atoms following lines will define the atomic contents of the system.
        # Note that this frame is not an actual MD frame. It is used only to declare the
        # atomic contents of the system.
        for i in range(n_atoms):
            
            line = castep.readline().split()
            
            # The symbol and the index of the atom are extracted from the line.
            symbol, idx = line[1:3]
            
            # An atom is created out of this symbol and added to the universe.
            universe.addObject(Atom(symbol, name = symbol + idx))

        # The output NetCDF trajectory is created.
        trajectory = Trajectory(universe, self.output_file, mode = 'w', comment = 'Converted from CASTEP file')
        
        # A snapshot creator is instanciated.
        snapshot = SnapshotGenerator(universe, actions = [TrajectoryOutput(trajectory,["all"], 0, None, 1)])

        comp = 0
        # Now parse the rest of the file and extract each frame.
        while True:
    
            # A new line is read.
            line = castep.readline().lower()

            # If the line contains 'fractional coordinates of atoms', this corresponds
            # to a new frame.
            if 'fractional coordinates of atoms' in line:
                
                try:
                    # The two next lines following the frame declaration are not useful.
                    line = castep.readline()
                    line = castep.readline()
                
                    # The n_atoms following lines will define the atom coordinates for that frame.
                    for i in range(n_atoms):

                        line = castep.readline().split()
                    
                        # The coordinates of atom i are extracted and stored into a Scientific Vector.
                        coord = Vector([float(v) for v in line[3:6]])
                    
                        # The coordinates are converted from box to real coordinates.
                        coord = universe.boxToRealCoordinates(coord)
                    
                        # The coordinates are registed in the universe current configuration.
                        universe.atomList()[i].setPosition(coord)

                    # The whole configuration is folded in to the simulation box.
                    universe.foldCoordinatesIntoBox()

                    # A snapshot is created out of the current configuration.
                    snapshot(data = {'time': comp*time_step})
        
                    comp += 1
                    
                except:
                    
                    # The CASTEP file is closed.
                    castep.close()

                    # The output trajectory is closed.
                    trajectory.close()

                    raise Error('Error when parsing frame %d. Will stop the conversion here.' % comp + 1)
    
            # This is the EOF condition. The conversion is stopped here.
            if not line:
                break
    
        # The CASTEP file is closed.
        castep.close()

        # The output trajectory is closed.
        trajectory.close()

class CHARMMConverter(object):
    """Converts a CHARMM Trajectory into a MMTK NetCDFFile. 

    @note: this code is based on the original converter written by Konrad Hinsen.
    """

    def __init__(self, pdbFile, dcdFile, outputFile):
        """The constructor.

        @param pdbFile: the CHARMM PDB file name of one frame.
        @type pdbFile: string

        @param dcdFile: the CHARMM DCD file name of the trajectory to convert.
        @type dcdFile: string

        @param outputFile: the MMTK NetCDF trajectory output file name of the trajectory to convert.
        @type outputFile: string

        @note: calling the constructor will actually perform the conversion.
        """

        # The arguments are copied to instance attributes.
        self.pdbFile = pdbFile
        self.dcdFile = dcdFile
        self.outputFile = outputFile

        # Do the conversion.
        self.__convert()

    def __convert(self):
        """Performs the actual conversion.
        """

        # Open the DCD trajectory file.
        dcd = DCDFile(self.dcdFile)

        # The starting step number.
        step = dcd.istart

        # The step increment.
        step_increment = dcd.nsavc

        # The MD time steps.
        dt = dcd.delta

        # The cell parameters a, b, c, alpha, beta and gamma (stored in |unit_cell|) 
        # and the x, y and z values of the first frame.
        unit_cell, x, y, z = dcd.readStep()

        # The DCD file does not contain any periodic boundary condition data.
        if not dcd.has_pbc_data:
            universe = InfiniteUniverse()

        # The DCD file contains periodic boundary condition datas. Create a parallelepipedic universe out of 
        # the first unit cell defined for the first step of the trajectory.
        else:
            universe = ParallelepipedicPeriodicUniverse(basisVectors(unit_cell))

        # Create all objects from the PDB file. The PDB file must match the
        # the DCD file (same atom order).
        conf = PDBConfiguration(self.pdbFile)

        # |molecules| is a Collection of all objects contained in the PDB file.
        molecules = conf.createAll()

        # The objects are introduced into the universe.
        universe.addObject(molecules)

        # A MMTK trajectory is opened for writing.
        trajectory = Trajectory(universe, self.outputFile, mode = 'w', comment = dcd.title)

        # A frame generator is created.        
        snapshot = SnapshotGenerator(universe, actions=[TrajectoryOutput(trajectory, ["all"], 0, None, 1)])

        # The array that will store the frame configuration.
        conf_array = universe.configuration().array

        # Iterates over the frame stored in the DCD file.
        while True:
            # The coordinates of the current frame are stored in |conf_array|.
            conf_array[:, 0] = x
            conf_array[:, 1] = y
            conf_array[:, 2] = z

            # Case of a periodic universe. Its shape is updated with the cell parameters of the current frame.
            if universe.is_periodic:
                universe.setShape(basisVectors(unit_cell))

            # The time corresponding to the frame.
            t = step*dt

            step_data = {'time': t, 'step': step}            

            # A call to the snapshot generator produces the step corresponding to the current frame.
            snapshot(data = step_data)

            step += step_increment

            # Tries to read the cell parameters and coordinates of the next frame.
            try:
                # They could be read, so store them in |unit_cell|, |x|, |y| and |z| variables.
                unit_cell, x, y, z = dcd.readStep()

            # The next step could not be read. This is the end of the file.
            except EndOfFile:
                break

        # Close the output trajectory
        trajectory.close()

class DL_POLYConverter(object):
    """Converts a DL_POLY Trajectory into a MMTK NetCDFFile. 
    """

    def __init__(self, fieldFile, historyFile, outputFile, specialAtoms = None):
        """The constructor.

        @param fieldFile: the DL_POLY FIELD file name of the trajectory to convert.
        @type fieldFile: string

        @param historyFile: the DL_POLY HISTORY file name of the trajectory to convert.
        @type historyFile: string

        @param outputFile: the MMTK NetCDF output file name.
        @type outputFile: string

        @param specialAtoms: dictionnary of the form {s1 : e1, s2 : e2 ...} where 's1', 's2' ... and 'e1', 'e2' ... 
            are respectively the DL_POLY name and the symbol of atoms 1, 2 ...
        @type specialAtoms: dict

        @note: calling the constructor will actually perform the conversion.
        """

        # The arguments are copied to instance attributes.
        self.fieldFile = fieldFile
        self.historyFile = historyFile
        self.outputFile = outputFile

        # The |specialAtoms| dictionnary is copied.
        if isinstance(specialAtoms, dict):
            self.specialAtoms = copy(specialAtoms)

        else:
            self.specialAtoms = {}

        # All the items of the special atoms dictionnary are lowerized. 
        for k,v in self.specialAtoms.items():
            if k.lower() != k:
                self.specialAtoms[k.lower()] = v.lower()
                del self.specialAtoms[k]
            else:
                self.specialAtoms[k] = v.lower()

        # Do the conversion.
        self.__convert()

    def __readFIELD(self):
        """Reads the DL_POLY FIELD file that contains informations about the system.

        @return: a list of tuples where each tuple stores the main information about each molecule declared in the FIELD file.
        @rtype: list
        """

        # The FIELD file is opened for reading, its contents stored into |lines| and then closed.
        try:
            fieldFile = open(self.fieldFile, 'r')
        except:
            raise Error('The FIELD file %s could not be opened for reading.' % self.fieldFile)

        lines = fieldFile.readlines()
        fieldFile.close()

        # Parses the contents of the FIELD file until the keywords 'molecules' or 'molecular types' have been found.
        while 1:

            # If the contents was parsed without finding any lines matching the molecules or molecular types 
            # keyword, raises an error.
            if not lines:
                raise Error('The keyword "molecules" or "molecular types" was not found in your DL_POLY/FIELD file.')

            # Check whether the keyword 'molecules' or ''molecular types has been matched.
            if ('MOLECULES' in lines[0].upper()) or ('MOLECULAR TYPES' in lines[0].upper()):
                try:
                    # |nmolecules| gives the number of different type of molecules in the system.
                    nmolecules = int(lines[0].split()[-1])
                except:
                    raise Error('The line containing the keywords "molecules" or "molecular types" must end with an integer indicating the number of molecules of the system.')
                else:
                    break

            # Otherwise, loops to the next line.
            else:            
                lines = lines[1:]

        lines = lines[1:]

        # |molecules| is a list that stores the main informations about the molecular system.
        molecules = []
        # Loop over the number of molecules declared in the FIELD file.
        for i in range(nmolecules):
            # |name| stores the name of the molecule.
            name = lines[0].strip()

            # |nummols| stores the number of such molecules.
            nummols = int(lines[1].split()[1])

            # |atoms| stores the number of atoms of that molecule.
            atoms = int(lines[2].split()[1])

            lines = lines[3:]

            # This list will stores tuples of the form (element,name) for each atom of the molecule appended in 
            # the order they are displayed in the FIELD file.
            atomicContents = []

            # Loop until the block defining atoms is not finished.
            while atoms > 0:
                # Each atom line is splitted in order to fetch the different atomic parameter
                # Take care the second and third values must be real, while the fourth and
                # fifth ones, although optional, must be integer.
                atomLine = lines[0].split()

                if len(atomLine) < 3:
                    raise Error('An atom declaration record must contain at least three element.')

                # |atom_name| stores the names of the atom.
                atom_name = atomLine[0].strip().lower()
                
                # The atomic mass of the current atom.
                atom_mass = string.strip("%6.1f" % float(atomLine[1]))

                # Determine what the element corresponding to that atom. Sets it first to the atom name.
                # Checks whether the atom name is a key of |self.specialAtoms| dictionnary.
                if atom_name in self.specialAtoms.keys():
                    # If so, the element is the corresponding value.
                    element = self.specialAtoms[atom_name]
                    
                # Otherwise, check whether its atomc mass is stored as a key of |M_to_Symbol| dictionnary.
                elif atom_mass in M_to_Symbol.keys():
                    element = M_to_Symbol[atom_mass]

                # Otherwise, tries to guess it from its name.
                else:

                    # By default the element is the name of the atom.
                    element = atom_name

                    # Tries to generate Atom out of the element string until one could be generated.
                    while len(element) > 0:

                        if not element:
                            raise Error('The element corresponding to atom %s could not be found in the MMTK Database.' % atom_name)

                        # Tries to generate an atom of type |element|.
                        try:
                            at = Atom(element)

                        # Failed, retries by removing the last character of |element|.
                        except:
                            element = element[:-1]

                        # Successful, the element corresponding to the atom is found.
                        else:
                            del at
                            break

                # The second entry of an atomic line must be a float indicating the atom weight.
                try:
                    atomweight = float(atomLine[1])

                # Otherwise, raises an error.
                except:
                    raise Error('The second value of the atomic record (atom weight) must be a number.')

                # The third entry of an atomic line must be a float indicating the atom charge.
                try:
                    atomchg = float(atomLine[2])

                # Otherwise, raises an error.
                except:
                    raise Error('The third value of the atomic record (atom charge) must be a number.')

                # The fourth entry of an atomic line must be an integer indicating the number of repeat of this atom.
                try:
                    nrepeat = max(int(atomLine[3]), 1)

                # Otherwise sets it to 1.
                except:
                    nrepeat = 1

                # The atomicContents list is updated with |nrepeat| tuples.
                atomicContents.extend([(element, atom_name)]*nrepeat)

                # The number of remaining atoms is updated.
                atoms = atoms - nrepeat

                # Remove the current line and treat the next one.
                lines = lines[1:]

            # |constraints| is a list of tuples storing the bond contraints found in the molecule.
            # Its element are of the form (atom 1 index, atom 2 index, bond length).
            constraints = []

            # Parses the contents of the FIELD file until the keywords 'CONSTRAINTS' have been found.
            # If the keyword 'FINISH' is found before, then loops to the next molecules.
            # This closes the definition of the current molecule.
            while (not re.match('FINISH', lines[0].upper())):

                # If the keyword 'finish' is not found, raises an error because it is compulsory to close a molecule definition block.
                if not lines:
                    raise Error('The keyword "finish" was not found in your DL_POLY/FIELD file. It is neccasry to close a molecular block.')

                if re.match('CONSTRAINTS',lines[0].upper()):
                    try:
                        # |nc| is the number of contrained bonds in the molecule
                        nc = int(lines[0].split()[1])
                    except:
                        raise Error('The keyword "constraints" must be followed by an integer indicating the number of constrained bonds.')
                    else:
                        for co in range(nc):
                            lines = lines[1:]
                            l = lines[0].split()
                            i1 = int(l[0])-1
                            i2 = int(l[1])-1
                            d = float(l[2])*Units.Ang
                            constraints.append((i1, i2, d))

                lines = lines[1:]

            # The |molecules| list updated by appending a tuple storing the main information about the molecule.
            # These informations are its name, its number, its atomic contents and its constraints.
            molecules.append([name, nummols, atomicContents, constraints])

            lines = lines[1:]

        return molecules

    def __readHISTORY(self):
        """Reads the HISTORY file and write the output trajectory.
        """

        try:
            # The HISTORY file is opened for reading.
            historyFile = open(self.historyFile, 'r')
        except:
            raise Error('The HISTORY file %s could not be opened for reading.' % self.historyFile)

        # The first line of the HISTORY file contains informations about the system.
        header = historyFile.readline().strip()

        # The second line of the HISTORY file gives some infos about the trajectory.
        # It must be made of of 3 integers.
        # The first one refers to the trajectory key (O = only coordinates, 1 = only velocities, 2 = both)
        # The second one refers to the periodic boundary key.
        # The third one gives the total number of atoms.
        keytrj, imcon, natms = FortranLine(historyFile.readline(), '3I10')
        # Checks that the type of periodic box is supported by MMTK. 
        if imcon >= 4:
            raise Error('Box shape not supported by MMTK')

        # |imcon| = 0 refers to non periodic boundary in DL_POLY. So creates an infinite universe.
        if imcon == 0:
            universe = InfiniteUniverse()

        # Otherwise, creates a periodic universe.
        else:
            universe = ParallelepipedicPeriodicUniverse()

        # Loops over the molecules declared in the FIELD file.
        for molname, nummols, atomicContents, constraints in self.molecules:
            # Loops over the number of molecules of of type |molname|.
            for comp in range(nummols):
                number = 0

                # This list will contains the MMTK instance of the atoms of the molecule.
                atom_objects = []

                # Loops over the atom of the molecule.
                for element, name in atomicContents:
                    # The atom is created.
                    a = Atom(element, name = name+str(number))
                    a.number = number
                    # And appended to |atom_objects| list.
                    atom_objects.append(a)
                    number = number + 1

                # If the length of ||atom_objects list is 1, then the molecule is actually a single atom.
                if len(atom_objects) == 1:
                    # Add the atom to the universe directly.
                    universe.addObject(atom_objects[0])

                # Otherwise, gather first the atoms of the molecules into an MMTK AtomCluster and append that
                # cluster to the universe.
                else:

                    # The AtomCluster is created.
                    ac = AtomCluster(atom_objects, name = molname)

                    # Loops over the constraints, if some were defined.
                    for i1, i2, d in constraints:
                        # Add the distances contraints to the definition of the AtomCluster.
                        ac.addDistanceConstraint(atom_objects[i1], atom_objects[i2], d)

                    # Add the AtomCluster to the universe.
                    universe.addObject(ac)

        # If |keytrj| > 0, then the HISTORY file contains also the velocities.
        if keytrj > 0:
            universe.initializeVelocitiesToTemperature(0.)

        # A MMTK trajectory is opened for writing.
        trajectory = Trajectory(universe, self.outputFile, mode = 'w', comment = header)

        # A frame generator is created.
        snapshot = SnapshotGenerator(universe, actions = [TrajectoryOutput(trajectory, ["all"], 0, None, 1)])

        # The current configuration.
        conf = universe.configuration()

        # The current velocities.
        vel = ParticleVector(universe)

        # |grad| is a ParticleVector MMTK object that associates a vector to each element of the universe.
        # That vector is the force associated to that element at a given time step.
        # This should not be of any use in nMoldyn but, as DLPOLY can provide it, it is done just in case ...
        grad = ParticleVector(universe)

        # Some predined format.
        history_timestep_line = FortranFormat('A8,4I10,F12.6')
        history_pbc_line = FortranFormat('3G12.4')

        # Loop until all the frames stored in the HISTORY file have been processed.
        while 1:

            # The HISTORY file was not closed so the reading continue from the third line.
            line = historyFile.readline()
            # An empty line means the EOF.
            if not line:
                break

            # The third line contains some time informations about the trajectory.
            data = line.split()

            # The current frame number.
            try:
                nstep = int(data[1])
            except:
                LogMessage('warning','The current frame number should be an integer.', ['gui','console'])
                break

            # The number of atoms in the configuration.
            try:
                natoms = int(data[2])
            except:
                LogMessage('warning','The number of atoms should be an integer.', ['gui','console'])
                break

            try:
                keytrj = int(data[3])
            except:
                LogMessage('warning','The trajectory key should be an integer.', ['gui','console'])
                break

            try:
                imcon = int(data[4])
            except:
                LogMessage('warning','The periodic boundary key should be an integer.', ['gui','console'])
                break

            # The trajectory time step.
            try:
                tstep = float(data[5])
            except:
                LogMessage('warning','The time step should be a float.', ['gui','console'])
                break

            step_data = {'time': nstep*tstep}

            # Case of a periodic universe. Sets the cell parameters.
            if imcon:

                    # This line refers to the basis vector 'a' in cartesian coordinates.
                line = historyFile.readline()
                try:
                    # The third line contains some time informations about the trajectory.
                    data = N.array([float(v) for v in line.split()], typecode = N.Float)
                    dirA = Vector(data)*Units.Ang
                except:
                    LogMessage('warning','The A vector should be declared using three floats.', ['gui','console'])
                    break

                # This line refers to the basis vector 'b' in cartesian coordinates.
                line = historyFile.readline()
                try:
                    # The third line contains some time informations about the trajectory.
                    data = N.array([float(v) for v in line.split()], typecode = N.Float)
                    dirB = Vector(data)*Units.Ang
                except:
                    LogMessage('warning','The B vector should be declared using three floats.', ['gui','console'])
                    break

                # This line refers to the basis vector 'c' in cartesian coordinates.
                line = historyFile.readline()
                try:
                    # The third line contains some time informations about the trajectory.
                    data = N.array([float(v) for v in line.split()], typecode = N.Float)
                    dirC = Vector(data)*Units.Ang
                except:
                    LogMessage('warning','The C vector should be declared using three floats.', ['gui','console'])
                    break

                cell = N.zeros((9,), typecode = N.Float)

                cell[0:3] = dirA
                cell[3:6] = dirB
                cell[6:9] = dirC

            # Case of non periodic universe. The cell parameters are set to None.
            else:
                cell = None

            # Loop over the total number of atoms.
            for i in range(natoms):

                historyFile.readline()

                # The configuration of atom |i| is updated.                
                conf.array[i] = [float(v) for v in historyFile.readline().split()]

                # Case where velocities are given in the HISTORY file.
                if keytrj > 0:
                    vel.array[i] = [float(v) for v in historyFile.readline().split()]
                    # Case where forces are also given in the HISTORY file.
                    if keytrj > 1:
                        grad.array[i] = [float(v) for v in historyFile.readline().split()]

            # The universe current configuration is updated with the corresponding frame of the HISTORY file.
            universe.setConfiguration(Configuration(universe, conf.array*Units.Ang, cell))

            # If |keytrj| > 0, the velocities are also updated.
            if keytrj > 0:
                universe.setVelocities(vel*Units.Ang/Units.ps)

            # If |keytrj| > 1, the are also introduced in the MMTK trajectory.
            if keytrj > 1:
                N.multiply(grad.array, -Units.amu*Units.Ang/Units.ps**2, grad.array)
                step_data['gradients'] = grad

            # A call to the snapshot generator produces the step corresponding to the current frame.
            snapshot(data = step_data)

        # The MMTK trajectory is closed.
        trajectory.close()

        # The HISTORY file is closed.
        historyFile.close()

    def __convert(self):
        """Performs the actual conversion.
        """

        # |self.molecules| is a list that contains the main informations about each molecular type found in the system.
        # Each element of that list corresponds to a molecular type found in the system and is the tuplet that
        # contains respectively the name of the molecular type, the number of molecules of that molecular type, the
        # list of atoms of one molecule, and the bond of the molecule with a constraints.
        self.molecules = self.__readFIELD()

        self.__readHISTORY()

class LAMMPSConverter(object):
    """Converts a LAMMPS Trajectory into a MMTK NetCDFFile. 
    """

    def __init__(self, lammps_cfg_file, lammps_file, output_file):
        """The constructor.

        @param lammps_file: the LAMMPS file name of the trajectory to convert.
        @type lammps_file: string

        @param output_file: the MMTK NetCDF output file name of the trajectory to convert.
        @type output_file: string

        @note: calling the constructor will actually perform the conversion.
        """

        # The arguments are copied to instance attributes.
        self.lammps_cfg_file = lammps_cfg_file
        self.lammps_file = lammps_file
        self.output_file = output_file

        # Do the conversion.
        self.__convert()

    def __convert(self):
        """Performs the actual conversion.
        """
        
        try:
            lammps_cfg = open(self.lammps_cfg_file, 'r')
        except:
            raise Error('The LAMMPS config file %s could not be opened for reading.' % self.lammps_cfg_file)

        nAtoms = None
        
        nAtomTypes = None
        
        idToElement = {}
                
        while True:
            
            line = lammps_cfg.readline()
            
            if not line:
                break
            
            line = line.upper().strip()
                        
            if nAtoms is None:
                match = re.findall("(\d+) ATOMS$",line)             
                if match:
                  nAtoms = int(match[0])
                  
            if nAtomTypes is None: 
                match = re.findall("(\d+) ATOM TYPES$",line)             
                if match:
                  nAtomTypes = int(match[0])

            if line == "MASSES":
                
                if nAtomTypes is None:
                    raise Error("Did not find the number of atom types.")
                                    
                line = lammps_cfg.readline()
                                
                for i in range(nAtomTypes):
                    id, mass = lammps_cfg.readline().strip().split()                    
                    idToElement[id] = M_to_Symbol[str(round(float(mass),1))]
                    
            if re.findall("^ATOMS$", line):
                
                if not idToElement:
                    raise Error("Did not find the type to mass table.")
            
                line = lammps_cfg.readline()
                                
                # The universe that will store the system is created.
                universe = ParallelepipedicPeriodicUniverse(None, None)
                ac = []                
                for i in range(nAtoms):
                    line = lammps_cfg.readline().split()
                    element = idToElement[line[2]]                                                            
                    at = Atom(element, name = "%s_%s" % (element, i+1))
                    ac.append(at)
                        
                clust = AtomCluster(ac, name = "clust")
                    
                # An atom is created out of this symbol and added to the universe.
                universe.addObject(clust)
                               
        lammps_cfg.close()
        
        # The output NetCDF trajectory is created.
        trajectory = Trajectory(universe, self.output_file, mode = 'w', comment = 'Converted from CASTEP file')
        
        # A snapshot creator is instanciated.
        snapshot = SnapshotGenerator(universe, actions = [TrajectoryOutput(trajectory,["all"], 0, None, 1)])
                            
        # The LAMMPS input file is opened for reading.
        try:
            lammps = open(self.lammps_file, 'r')

        # Could not be opened. Raises an error.
        except:
            raise Error('The LAMMPS file %s could not be opened for reading.' % self.lammps_file)

        # The time step. Will be extracted from the CASTEP file.
        timeStep   = None
        
        abcVectors = None
                
        while True:

            # A line is read.
            line = lammps.readline()
                                    
            if not line: 
                break
            
            line = line.upper()

            itemType = ""
            
            if line[:4] == "ITEM":
                itemType = "".join(line.split(":")[1:]).replace(" ","").strip()
                         
            if itemType == "BOXBOUNDS":
                abcVectors = N.zeros((9), typecode = N.Float)
                
                temp = [float(v) for v in lammps.readline().split()]
                if len(temp) == 2:
                    xlo, xhi = temp
                    xy = 0.0
                elif len(temp) == 3:
                    xlo, xhi, xy = temp
                else:
                    raise Error("Bad format for A vector components")

                temp = [float(v) for v in lammps.readline().split()]
                if len(temp) == 2:
                    ylo, yhi = temp
                    xz = 0.0
                elif len(temp) == 3:
                    ylo, yhi, xz = temp
                else:
                    raise Error("Bad format for B vector components")

                temp = [float(v) for v in lammps.readline().split()]
                if len(temp) == 2:
                    zlo, zhi = temp
                    yz = 0.0
                elif len(temp) == 3:
                    zlo, zhi, yz = temp
                else:
                    raise Error("Bad format for C vector components")
                              
                # The ax component.                                      
                abcVectors[0] = xhi - xlo
                
                # The bx and by components.                                      
                abcVectors[3] = xy
                abcVectors[4] = yhi - ylo
                
                # The cx, cy and cz components.                                      
                abcVectors[6] = xz
                abcVectors[7] = yz
                abcVectors[8] = zhi - zlo
                                                                
                abcVectors *= Units.Ang
                
            elif itemType == "TIMESTEP":
                temp = lammps.readline()
                timeStep = Units.fs*float(temp)
                
            elif itemType[0:5] == "ATOMS":
                                                                                
                universe.setCellParameters(abcVectors)

                conf = universe.configuration()
                                                            
                for i in range(nAtoms):
                    temp = lammps.readline().split()                    
                    id = int(temp[0])                    
                    coord = universe.boxToRealCoordinates(N.array([float(v) for v in temp[2:5]]) - 0.5)
                      
                    conf.array[id-1,:] = coord                    

                # The whole configuration is folded in to the simulation box.
                universe.foldCoordinatesIntoBox()

                # A snapshot is created out of the current configuration.
                snapshot(data = {'time': timeStep})
                                                                 
        # The LAMMPS trajectory file is closed.
        lammps.close()

        # The output trajectory is closed.
        trajectory.close()

class MaterialsStudioConverter(object):
    """Converts a MaterialsStudio Discover or Forcite Trajectory into a MMTK NetCDFFile. 
    """

    def __init__(self, module, xtdxsdFile, histrjFile, outputFile):
        """The constructor.

        @param module: the MaterialsStudio module the input trajectory is coming from.
        @type module: string being one of 'Discover' or 'Forcite'

        @param xtdxsdFile: the MaterialsStudio XTD or XSD file name of the trajectory to convert.
        @type xtdxsdFile: string

        @param histrjFile: the MaterialsStudio HIS (Discover) or TRJ (Forcite) file name of the trajectory to convert.
        @type histrjFile: string

        @param outputFile: the  MMTK NetCDF trajectory output file name.
        @type outputFile: string

        @note: calling the constructor will actually perform the conversion.
        """

        # The arguments are copied to instance attributes.
        self.module = module.lower()
        self.xtdxsdFile = xtdxsdFile
        self.histrjFile = histrjFile
        self.outputFile = outputFile        

        # Do the conversion.
        self.__convert()

    def createCluster(self, atDict, clust = None):
        """Retrieves the complete arborescence of the atoms linked to a given atom.

        @param atDict: a dictionnary storing some information about a given atom.
        @type atDict: dictionnary

        @param clust: a list storing the atom dictionnaries for the cluster under construction.
        @type clust: dictionnary

        @return: a list storing the atom dictionnaries for the cluster under construction.
        @rtype: list

        @note: this is a recursive function.
        """

        if clust is None:
            clust = []

        if atDict.has_key('cluster'):
            return

        atDict['cluster'] = True

        clust.append(atDict)

        for cAt in atDict['connections']:
            self.createCluster(self.atoms[cAt], clust)

        return clust

    def readXTDFile(self):
        """Reads the Materials Studio XTD or  XSD file and set up the universe from which the NetCDF 
        MMTK trajectory will be written. 

        @note: the XTD and XSD file are xml file.
        """        

        # Parse the file.
        xmlFile = xml.dom.minidom.parse(self.xtdxsdFile)

        # Does the system contain PBC ?
        symmetry = xmlFile.getElementsByTagName('SymmetrySystem')

        # The system has PBC
        if symmetry:
            # Search the Tag 'IdentityMapping' to characterize the system.
            systemDef = xmlFile.getElementsByTagName('IdentityMapping')[0]
            imageMapping = xmlFile.getElementsByTagName('ImageMapping')
            # Search for the Tag 'SpaceGroup'.
            spaceGroup = systemDef.getElementsByTagName('SpaceGroup')[0]

            # If found, extract the Tags 'AVector', 'BVector', 'CVector' and fills the 
            # |self.normalizedCell| attribute array with them.
            self.cell = N.zeros((9,), typecode = N.Float)
            self.cell[0:3] = Vector([float(v) for v in spaceGroup.getAttribute('AVector').split(',')])
            self.cell[3:6] = Vector([float(v) for v in spaceGroup.getAttribute('BVector').split(',')])
            self.cell[6:9] = Vector([float(v) for v in spaceGroup.getAttribute('CVector').split(',')])
            self.cell *= Units.Ang
            
            self.normalizedCell = N.zeros((9,), typecode = N.Float)
            self.normalizedCell[0:3] = Vector(self.cell[0:3]).normal()
            self.normalizedCell[3:6] = Vector(self.cell[3:6]).normal()
            self.normalizedCell[6:9] = Vector(self.cell[6:9]).normal()
            
            # Creates a periodic universe.
            self.universe = ParallelepipedicPeriodicUniverse()

        else:
            LogMessage('warning',\
                       'No SpaceGroup tag found in the xtd/xsd file.\nThe MMTK universe created will be an infinite universe.',\
                       ['gui','console'])
            # Search the Tag 'Molecule' to characterize the system.
            systemDef = xmlFile.getElementsByTagName('Molecule')[0]
            self.universe = InfiniteUniverse()
            imageMapping = None
            self.cell    = None
            
            
        atomsMapping = {}
        
        # Search the Tag 'Atom3d'.
        atomTag = systemDef.getElementsByTagName('Atom3d')        
        # Loop over the atoms stored int he xtd/xsd file.
        for a in atomTag:            
            # The atom index.
            atID = int(a.getAttribute('ID'))
            atomsMapping[atID] = atID
        
        if imageMapping:
            for img in imageMapping:
                # Search the Tag 'Atom3d'.
                atomsTag = img.getElementsByTagName('Atom3d')
                for a in atomsTag:                
                    # The atom index.
                    atID = int(a.getAttribute('ID'))
                    atomsMapping[atID] = int(a.getAttribute('ImageOf'))                    

        bonds = {}
        # Search for the Tag 'Bond'.
        bondTag = systemDef.getElementsByTagName('Bond')
        # Loop over the bonds.
        for b in bondTag:
            bondID = int(b.getAttribute('ID'))
            bonds[bondID] = [atomsMapping[int(v)] for v in b.getAttribute('Connects').split(',')]
                            
        if imageMapping:
            for img in imageMapping:
                bondTag = img.getElementsByTagName('Bond')
                for b in bondTag:
                    bondID = int(b.getAttribute('ID'))
                    bonds[bondID] = bonds[int(b.getAttribute('ImageOf'))]

        comp = 0
        
        self.atoms = {}
        
        # Search the Tag 'Atom3d'.
        atomTag = systemDef.getElementsByTagName('Atom3d')
       
        # Loop over the atoms stored int he xtd/xsd file.
        for a in atomTag:
            
            atomInfo = {}
            
            # The atom index.
            atomID = int(a.getAttribute('ID'))
            
            # The serial number.
            atomInfo['id'] = atomID
                        
            # The serial number.
            atomInfo['serial_number'] = comp
            
            # The corresponding element.
            atomInfo['element'] = str(a.getAttribute('Components')).split(',')[0].strip()
                        
            # The corresponding coordinates converted in nanometers.
            atomInfo['coord'] = N.array([float(v) for v in a.getAttribute('XYZ').split(',')], typecode = N.Float)
                
            # Try to retrieve the name of the atom by using the Tag 'Name'.
            name = str(a.getAttribute('Name')).strip()
            
            # Case where the atom has no tag 'Name'
            if len(name) == 0:
                # Try to build an atom name using the tag 'ForcefieldType'
                name = str(a.getAttribute('ForcefieldType')).strip()
                    
                # If the atom also do not have a tag 'ForcefieldType', use its element name.
                if len(name) == 0:
                    atomInfo['name'] = atomInfo['element'] + '_el'
                        
                else:
                    atomInfo['name'] = name + '_ff'                        
                    
            else:                
                atomInfo['name'] = name
            
            # This list will store the indexes of the atoms connected to the current atom.
            connections = []
            # Check whether the atom has some connections.
            if a.hasAttribute('Connections'):
                # List of the indexes of the bond in which the atom is involved.
                connectedBonds = [int(v) for v in a.getAttribute('Connections').split(',')]
                # Loop over these bonds.
                for b in connectedBonds:
                    # This ensure that the bonds is a covalent bond and not a H bond.
                    if bonds.has_key(b):
                        # Loop over the indexes of the two atoms involved in the bond.
                        for bondedAtom in bonds[b]:
                            # Adds the atoms connected to the current atom to the |connections| list.
                            if bondedAtom != atomID:
                                connections.append(bondedAtom)
                                break
            
            atomInfo['connections'] = connections
            
            self.atoms[atomID] = atomInfo
                        
            comp += 1
                                                    
        self.clusters = []
                
        # Loop over the dictionnaries of |self.atomContents| list in order to recreate the molecular arborescence.
        for a in self.atoms.values():
                        
            # If the atom dictionnary has already the key 'cluster', then it has already been clusterized. So, skip it.
            if a.has_key('cluster'):
                continue
            
            # Otherwise, retrieves recursively the complete arborescence of the atoms linked to that atom.
            else:
                
                clust = self.createCluster(a)

                # Sort the cluster by ascending first by atom name then by atom index.
                clust = sorted(clust, key = operator.itemgetter("element","id"))

                # Dictionnary whose keys are the elements found in the system and the values are the corresponding number.
                elements = {}

                for at in clust:
                                 
                    # The |atomNames| is updated.
                    if elements.has_key(at['element']):
                        elements[at['element']] += 1
                    else:
                        elements[at['element']] = 1

                # A name for the cluster is defined.
                cName = ''
                for k,v in sorted(elements.items()):
                    cName += '%s%d' % (k,v)
                                        
                cId = clust[0]["id"]
                                       
                # Appends the cluster to |self.clusters| list.
                self.clusters.append({"cname" : cName, "cid" : cId, "clust" : clust})

        # The XML file is closed.
        xmlFile.unlink()
        
        self.clusters = sorted(self.clusters, key = operator.itemgetter("cname","cid"))
                        
        # Loops over the number of clusters.
        for c in range(len(self.clusters)):
            
            atomNames = {}

            # List of the MMTK atoms forming the cluster.
            atomList = []

            # Loop over the atom dictionnaries of the cluster.
            for at in self.clusters[c]["clust"]:
                                                    
                # Handles the case where several atoms have the same name in the molecule.
                name = at['name']
                
                if atomNames.has_key(name):                    
                    atomNames[name] += 1 
                    
                else:
                    atomNames[name] = 1 
                    
                name += '_%d' % atomNames[name]
                    
                # The |atoms| list is updated with the MMTK Atom instance of the atom.
                atomList.append(Atom(at['element'], name = name))
                
            # An AtomCluster is made out of the |atoms| list and added to the universe.
            self.universe.addObject(AtomCluster(atomList, name = self.clusters[c]["cname"]))

    def readHISFile(self):
        """Reads a Materials Studio HIS file and builds the NetCDF trajectory file.
        """

        # This array will stores the coordinates of the selected atoms.
        selCoordinates = N.zeros((self.universe.numberOfAtoms(),3), typecode = N.Float)

        # This array will stores the velocities of the selected atoms.
        selVelocities = N.zeros((self.universe.numberOfAtoms(),3), typecode = N.Float)

        # The MMTK trajectory is opened for writing.
        trajectory = Trajectory(self.universe, self.outputFile, mode = 'w', comment = 'From Discover trajectory file')

        # A frame generator is created.
        snapshot = SnapshotGenerator(self.universe, actions = [TrajectoryOutput(trajectory,["all"], 0, None, 1)])

        # Read the HISTORY binary file.
        hisFile = open(self.histrjFile,'rb')

        # Record 1
        rec = '!4xi8x'
        recSize = calcsize(rec)
        hisFile.read(recSize)

        # Record 2
        rec = '!' + '4s' * 20 + 'd8x'
        recSize = calcsize(rec)        
        data = struct.unpack(rec, hisFile.read(recSize))
        Vershn = data[20]

        # Some intermediate information.
        LogMessage('info','Vershn -- > %s' % Vershn,['file'])

        # Record 3
        rec = '!' + '4s' * 20 + '8x'
        recSize = calcsize(rec)
        hisFile.read(recSize)

        # Record 4
        rec = '!' + '4s' * 20 + '8x'
        recSize = calcsize(rec)
        hisFile.read(recSize)

        # Record 5
        rec = '!i'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        NAtTyp = data[0]

        rec = '!' + '4s' * NAtTyp + '%sd8x' % NAtTyp
        recSize = calcsize(rec)
        hisFile.read(recSize)

        # Record 6        
        rec = '!i'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        NNmRes = data[0]

        rec = '!' + '4s' * NNmRes + '8x'
        recSize = calcsize(rec)
        hisFile.read(recSize)

        # Record 7
        rec = '!i'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        NAtoms = data[0]

        rec = '!%si' % NAtoms
        recSize = calcsize(rec)
        hisFile.read(recSize)

        if Vershn < 2.9:
            rec = '!' + '4s' * NAtoms

        else:
            rec = '!' + '5s' * NAtoms

        recSize = calcsize(rec)
        hisFile.read(recSize)

        hisFile.read(calcsize('!8x'))

        # Record 8
        rec = '!ii'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))

        NAtMov = data[1]

        if Vershn >= 2.6:
            rec = '!%si' % NAtMov
            recSize = calcsize(rec)
            data = struct.unpack(rec, hisFile.read(recSize))
            movableatoms = data
            LogMessage('info','movableatoms -- > %s' % str(movableatoms),['file'])

        hisFile.read(calcsize('!8x'))

        # Record 9
        rec = '!i'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        NMol = data[0]

        rec = '!%si%si' % (NMol,NMol)
        recSize = calcsize(rec)
        hisFile.read(recSize)

        hisFile.read(calcsize('!8x'))

        # Record 10
        rec = '!i'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        NRes = data[0]

        rec = '!%si%si' % (NRes*2,NRes)
        recSize = calcsize(rec)
        hisFile.read(recSize)

        hisFile.read(calcsize('!8x'))

        # Record 11
        rec = '!i'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        NBonds = data[0]

        if NBonds > 0 :
            rec = '!%si' % 2*NBonds
            recSize = calcsize(rec)
            hisFile.read(recSize)

        hisFile.read(calcsize('!8x'))

        # Record 12
        rec = '!4149di4di6d6i8x'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))

        # Record 13
        rec = '!idii8x'
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        NEner, timestep, frequency, startingstep = data

        LogMessage('info','NEner -- > %s' % NEner,['file'])
        LogMessage('info','timestep -- > %s fs' % timestep,['file'])
        LogMessage('info','frequency -- > %s step(s)' % frequency,['file'])
        LogMessage('info','startingstep -- > %s' % startingstep,['file'])

        # Record 14
        nrjBlockSize = 59+5*NMol+NEner+NMol*NEner
        rec = '!%sd8x' % nrjBlockSize
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))
        totalenergy = data[0]

        LogMessage('info','totalenergy -- > %s' % totalenergy,['file'])

        nAtoms3 = 3*NAtoms 

        # Record 15
        rec = '!%sf8x' % nAtoms3
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))        
        firstatomxyz = data[0:3]

        LogMessage('info','firstatomxyz -- > %s' % str(firstatomxyz),['file'])

        # Record 16
        rec = '!%sf8x' % nAtoms3
        recSize = calcsize(rec)
        data = struct.unpack(rec, hisFile.read(recSize))        
        firstatomvel = data[0:3]

        LogMessage('info','firstatomvel -- > %s' % str(firstatomvel),['file'])

        # Frame record
        frame = 0
        while True:

            try:
                # Record N
                rec = '!i8x'
                recSize = calcsize(rec)
                hisFile.read(recSize)

                # Record N + 1
                rec = '!%sd8x' % nrjBlockSize
                recSize = calcsize(rec)
                hisFile.read(recSize)

                # Record N + 2
                rec = '!6d9d8x'
                recSize = calcsize(rec)
                data = struct.unpack(rec, hisFile.read(recSize))

                self.cell = N.array(data[6:], typecode = N.Float)*Units.Ang                
                LogMessage('info','cell -- > %s' % str(self.cell),['file'])

                # Record N+3
                rec = '!%sf8x' % nAtoms3
                recSize = calcsize(rec)
                data = struct.unpack(rec, hisFile.read(recSize))
                xyz = N.array(data, typecode = N.Float).reshape((NAtoms,3))*Units.Ang

                # Record N+4
                rec = '!%sf' % nAtoms3
                recSize = calcsize(rec)
                data = struct.unpack(rec, hisFile.read(recSize))
                vel = N.array(data, typecode = N.Float).reshape((NAtoms,3))*Units.Ang

                hisFile.read(calcsize('!8x'))

            except:
                break

            else:
                frame += 1

                comp = 0
                # Loop over the cluster.
                for cluster in self.clusters:
                    # Loop over the atoms of the cluster.
                    for atom in cluster["clust"]:                        
                        selCoordinates[comp,:] = xyz[atom['serial_number'],:]
                        selVelocities[comp,:] = vel[atom['serial_number'],:]
                        comp += 1

                # The configuration is updated with the coordinates of the current frame.
                self.universe.setConfiguration(Configuration(self.universe, selCoordinates, self.cell))

                # The velocities are updated with the velocities of the current frame.
                self.universe.setVelocities(ParticleVector(self.universe, selVelocities))

                # A call to the snapshot generator produces the step corresponding to the current frame.
                snapshot(data = {'time': (startingstep + frequency*frame)*timestep*Units.fs})

        # The MMTK trajectory is closed.
        trajectory.close()

        # The HIS file is closed.
        hisFile.close()

    def readTRJFile(self):
        """Reads a Materials Studio HIS file and fills up the NetCDF trajectory file.
        """

        # This array will stores the coordinates of the selected atoms.
        selCoordinates = N.zeros((self.universe.numberOfAtoms(),3), typecode = N.Float)

        # This array will stores the velocities of the selected atoms.
        selVelocities = N.zeros((self.universe.numberOfAtoms(),3), typecode = N.Float)

        # The MMTK trajectory is opened for writing.
        trajectory = Trajectory(self.universe, self.outputFile, mode = 'w', comment = 'From Forcite trajectory file')

        # A frame generator is created.
        snapshot = SnapshotGenerator(self.universe, actions = [TrajectoryOutput(trajectory,["all"], 0, None, 1)])

        # Read the his file.
        trjFile = open(self.histrjFile,'rb')

        # Record 1
        rec = '!4x4s20i8x'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))

        HDR = data[0]
        LogMessage('info','HDR -- > %s' % HDR,['file'])

        ICNTRL = data[1:]
        LogMessage('info','ICNTRL -- > %s' % str(ICNTRL),['file'])

        if ICNTRL[0] == 2000:
            prec = 'f'

        elif ICNTRL[0] == 2010:
            prec = 'd'

        elif ICNTRL[0] == 3000:
            prec = 'd'

        else:
            LogMessage('warning','Unknown trj version number: %s.' % ICNTRL[0],['gui'])
            return

        # Diff with doc --> NTRJTI and TRJTIC not in doc
        rec = '!i'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))
        NTRJTI, = data

        rec = '!' + '80s'* NTRJTI
        recSize = struct.calcsize(rec)
        trjFile.read(recSize)

        trjFile.read(calcsize('!8x'))

        # Record 2
        rec = '!i'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))
        NEEXTI, = data

        rec = '!' + '80s'* NEEXTI
        recSize = struct.calcsize(rec)
        trjFile.read(recSize)

        trjFile.read(calcsize('!8x'))

        # Record 3
        rec = '!8i8x'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))

        PERTYPE, MOLXTL, LCANON, DEFCEL, PRTTHRM, LNOSE, LNPECAN, LTMPDAMP = data

        LogMessage('info','PERTYPE -- > %s' % PERTYPE,['file'])
        LogMessage('info','MOLXTL -- > %s' % MOLXTL,['file'])
        LogMessage('info','LCANON -- > %s' % LCANON,['file'])
        LogMessage('info','DEFCEL -- > %s' % DEFCEL,['file'])
        LogMessage('info','PRTTHRM -- > %s' % PRTTHRM,['file'])
        LogMessage('info','LNOSE -- > %s' % LNOSE,['file'])
        LogMessage('info','LNPECAN -- > %s' % LNPECAN,['file'])
        LogMessage('info','LTMPDAMP -- > %s' % LTMPDAMP,['file'])

        # Record 4
        rec = '!i'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))
        NFLUSD, = data

        # Diff with doc 8x at the end
        rec = '!%si%si' % (NFLUSD, NFLUSD) + '8s'* NFLUSD + '8x'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))

        MVATPF = data[0:NFLUSD]
        NATPFU = data[NFLUSD:2*NFLUSD]
        DECUSD = data[2*NFLUSD:3*NFLUSD]

        LogMessage('info','MVATPF -- > %s' % str(MVATPF),['file'])
        LogMessage('info','NATPFU -- > %s' % str(NATPFU),['file'])
        LogMessage('info','DECUSD -- > %s' % str(DECUSD),['file'])

        rec = '!i'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))

        TOTMOV = data[0]
        
        LogMessage('info','TOTMOV -- > %s' % str(TOTMOV),['file'])

        rec = '!%si8x' % TOTMOV
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))

        MVOFST = [d-1 for d in data]

        LogMessage('info','MVOFST -- > %s' % str(MVOFST),['file'])

        # Record 4a
        rec = '!i'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))

        LEEXTI, = data
        
        rec = '!' + '%ss'% LEEXTI
        recSize = struct.calcsize(rec)
        trjFile.read(recSize)

        trjFile.read(calcsize('!8x'))

        # Record 4b
        rec = '!i'
        recSize = struct.calcsize(rec)
        data = struct.unpack(rec, trjFile.read(recSize))

        LPARTI, = data
        
        rec = '!' + '%ss'% LPARTI
        recSize = struct.calcsize(rec)
        trjFile.read(recSize)

        trjFile.read(calcsize('!8x'))
        
        while True:

            try:            

                # Frame information
                # Record 1
                if ICNTRL[0] == 2000:
                    rec = '!%si33%s5i8x' % (prec, prec)
                    velInd = 37

                elif ICNTRL[0] == 2010:
                    rec = '!%si57%s6i8x' % (prec, prec)
                    velInd = 61

                elif ICNTRL[0] == 3000:
                    rec = '!%si58%s6i8x' % (prec, prec)
                    velInd = 62

                recSize = struct.calcsize(rec)
                data = struct.unpack(rec, trjFile.read(recSize))
                
                CurrentTime     = data[0]
                itstep          = data[1]
                Temp            = data[2]
                AvgTemp         = data[3]
                TimeStep        = data[4]
                InitialTemp     = data[5]
                FinalTemp       = data[6]
                TotalPE         = data[7]
                VelocityWritten = data[velInd]
                                
                if ICNTRL[0] == 2000:
                    ForcesWritten   = False

                elif ICNTRL[0] == 2010:
                    ForcesWritten   = data[62]
                    
                elif ICNTRL[0] == 3000:
                    ForcesWritten   = data[63]                                    

                # Record 2
                rec = '!12%s8x' % prec
                recSize = struct.calcsize(rec)
                data = struct.unpack(rec, trjFile.read(recSize))

                Press             = data[0]
                Volume            = data[1]
                Junk              = data[2:5]
                GyrationRadius    = data[5]
                AvgPress          = data[6]
                AvgVolume         = data[7]
                Junk              = data[8:11]
                AvgGyrationRadius = data[11]
                                
                # Record 3
                if LCANON:
                    rec = '!4%s8x' % prec
                    recSize = struct.calcsize(rec)
                    data = struct.unpack(rec, trjFile.read(recSize))

                    if LNOSE:
                        snose       = data[0]
                        snoseh      = data[1]
                        dsstot      = data[2]
                        dqcanonNose = data[3]

                    else:            
                        signose   = data[0]
                        zfrict    = data[1]
                        dzprfrict = data[2]
                        dqcanon   = data[3]
                        
                # Record 4
                if PERTYPE > 0:
                    rec = '!22%s8x' % prec
                    recSize = struct.calcsize(rec)
                    data = struct.unpack(rec, trjFile.read(recSize))

                    a = data[2]*Units.Ang
                    b = data[3]*Units.Ang
                    c = data[4]*Units.Ang

                    self.cell[0:3] = a*self.normalizedCell[0:3]
                    self.cell[3:6] = b*self.normalizedCell[3:6]
                    self.cell[6:9] = c*self.normalizedCell[6:9]
                    
                # Record 5
                if PERTYPE > 0:
                    rec = '!i14%s8x' % prec
                    recSize = struct.calcsize(rec)
                    data = struct.unpack(rec, trjFile.read(recSize))

                    NUnitCellAtoms = data[0]
                    Junk           = data[1:15]

                # Record 6
                if LNPECAN:
                    rec = '!3%s8x' % prec
                    recSize = struct.calcsize(rec)
                    trjFile.read(recSize)

                # Record 7
                if LTMPDAMP:
                    rec = '!%s8x' % prec
                    recSize = struct.calcsize(rec)
                    trjFile.read(recSize)

                # Record 8
                rec = '!%s%s' % (TOTMOV, prec)
                recSize = struct.calcsize(rec)
                data = struct.unpack(rec, trjFile.read(recSize))

                XCOORD = N.array(data[0:TOTMOV], typecode = N.Float)*Units.Ang
                trjFile.read(calcsize('!8x'))

                # Record 9
                rec = '!%s%s' % (TOTMOV, prec)
                recSize = struct.calcsize(rec)
                data = struct.unpack(rec, trjFile.read(recSize))

                YCOORD = N.array(data[0:TOTMOV], typecode = N.Float)*Units.Ang
                trjFile.read(calcsize('!8x'))

                # Record 10
                rec = '!%s%s' % (TOTMOV, prec)
                recSize = struct.calcsize(rec)
                data = struct.unpack(rec, trjFile.read(recSize))

                ZCOORD = N.array(data[0:TOTMOV], typecode = N.Float)*Units.Ang  
                trjFile.read(calcsize('!8x'))                

                if VelocityWritten:
                    # Record 11
                    rec = '!%s%s' % (TOTMOV, prec)
                    recSize = struct.calcsize(rec)
                    data = struct.unpack(rec, trjFile.read(recSize))

                    XVelocity = N.array(data[0:TOTMOV], typecode = N.Float)*Units.Ang
                    trjFile.read(calcsize('!8x'))                

                    # Record 12
                    rec = '!%s%s' % (TOTMOV, prec)
                    recSize = struct.calcsize(rec)
                    data = struct.unpack(rec, trjFile.read(recSize))

                    YVelocity = N.array(data[0:TOTMOV], typecode = N.Float)*Units.Ang
                    trjFile.read(calcsize('!8x'))                

                    # Record 13
                    rec = '!%s%s' % (TOTMOV, prec)
                    recSize = struct.calcsize(rec)
                    data = struct.unpack(rec, trjFile.read(recSize))

                    ZVelocity = N.array(data[0:TOTMOV], typecode = N.Float)*Units.Ang
                    trjFile.read(calcsize('!8x'))                

                if ForcesWritten:
                    # Record 14
                    rec = '!%s%s' % (TOTMOV, prec)
                    recSize = struct.calcsize(rec)
                    trjFile.read(recSize)

                    trjFile.read(calcsize('!8x'))                

                    # Record 15
                    rec = '!%s%s' % (TOTMOV, prec)
                    recSize = struct.calcsize(rec)
                    trjFile.read(recSize)

                    trjFile.read(calcsize('!8x'))                

                    # Record 16
                    rec = '!%s%s' % (TOTMOV, prec)
                    recSize = struct.calcsize(rec)
                    trjFile.read(recSize)

                    trjFile.read(calcsize('!8x'))                

            except:
                break

            else:

                self.universe.setCellParameters(self.cell)
                
                comp = 0
                # Loop over the clusters.
                for cluster in self.clusters:
                    
                    # Loop over the atom of the cluster.
                    for atom in cluster["clust"]:

                        # If the atom is not fixed then uses its coordinates and velocities in the current frame.
                        if atom['serial_number'] in MVOFST:
                            ind = MVOFST.index(atom['serial_number'])
                            selCoordinates[comp,0] = XCOORD[ind]
                            selCoordinates[comp,1] = YCOORD[ind]
                            selCoordinates[comp,2] = ZCOORD[ind]
                            
                            if VelocityWritten:
                                selVelocities[comp,0]  = XVelocity[ind]
                                selVelocities[comp,1]  = YVelocity[ind]
                                selVelocities[comp,2]  = ZVelocity[ind]

                        # Otherwise, uses its initial coordinates stored in box coordinates in the ctd file.
                        # The velocities for such atoms being obviously null.
                        else:
                            
                            if PERTYPE > 0:
                                selCoordinates[comp,:] = self.universe.boxToRealCoordinates(Vector(atom['coord']))
                                
                            else:
                                selCoordinates[comp,:] = atom['coord']*Units.Ang

                        comp += 1
                                                
                # The configuration is updated with the coordinates of the current frame.
                self.universe.setConfiguration(Configuration(self.universe, selCoordinates))
                
                # Sets the velocities if they were present in the file.
                if VelocityWritten:    
                    # The velocities are updated with the velocities of the current frame.
                    self.universe.setVelocities(ParticleVector(self.universe, selVelocities))

                # A call to the snapshot generator produces the step corresponding to the current frame.
                snapshot(data = {'time': CurrentTime})

        # The MMTK trajectory is closed.
        trajectory.close()

        # The TRJ file is closed.
        trjFile.close()

    def __convert(self):
        """Performs the actual conversion.
        """

        # Parses the XTD file and retrieves the arborescence of the different molecules of the system.
        self.readXTDFile()

        # Case of an input trajectory written with 'Discover'.
        if self.module == 'discover':
            self.readHISFile()

        # Case of an input trajectory written with 'Forcite'.
        elif self.module == 'forcite':
            self.readTRJFile()

class NAMDConverter(object):
    """Converts a NAMD Trajectory into a MMTK NetCDFFile. 

    @note: this code is based on the original converter written by Konrad Hinsen.
    """

    def __init__(self, pdbFile, dcdFile, xstFile, outputFile):
        """The constructor.

        @param pdbFile: the NAMD PDB file name of one frame of the trajectory to convert.
        @type pdbFile: string

        @param dcdFile: the NAMD DCD file name of the trajectory to convert.
        @type dcdFile: string

        @param xstFile: the NAMD XSTfile name of the trajectory to convert.
        @type xstFile: string

        @param outputFile: the MMTK NetCDF trajectory output file name.
        @type outputFile: string

        @note: calling the constructor will actually perform the conversion.
        """

        # The arguments are copied to instance attributes.
        self.pdbFile = pdbFile
        self.dcdFile = dcdFile
        self.xstFile = xstFile
        self.outputFile = outputFile

        # Do the conversion.
        self.__convert()

    def __convert(self):
        """Performs the actual conversion."""

        # Process the optional XST file if one was loaded.
        if self.xstFile is not None:

            # The file is opened for reading.
            xst = open(self.xstFile, 'r')

            # This list will contains the 9 cell parameters (ax,ay,az,bx,by,bz,cx,cy,cz) gathered for each 
            # time step of the MD.
            temp = []
            # Loop over the lines of the XST file that contains the cell parameters at each step of the MD.
            for pLines in xst.readlines()[3:]:
                # The 9 cell parameters are extended to the |temp| list.
                temp.extend(pLines.strip().split()[1:10])

            # The XST file is closed.
            xst.close()

            # An array is made out of |temp| list with an automatic float and nanometer conversion.
            a = N.array(temp, typecode = N.Float)*Units.Ang
            # The array is reshape as (number of time steps,3,3)
            a = N.reshape(a,(a.shape[0]/9,3,3))

            # This list will contains the nested lists of a,b and c vectors at time step of the MD.
            xstData = []            

            # Loop over the number of time steps.
            for comp in range(a.shape[0]):      
                # The Vector a, b and c are built and appended as nested list in |xstData| list.
                xstData.append([Vector(v) for v in a[comp,:,:]])

            # A periodic universe is created using the cell parameters of the first frame.                        
            universe = ParallelepipedicPeriodicUniverse(xstData[0])

        # Otherwise, create an infinite universe.
        else:
            universe = InfiniteUniverse()

        # Open the DCD trajectory file for reading.
        dcd = DCDFile(self.dcdFile)

        # The starting step number.
        step = dcd.istart

        # The step increment.
        step_increment = dcd.nsavc

        # The MD time steps.
        dt = dcd.delta

        # The cell parameters a, b, c, alpha, beta and gamma (stored in |unit_cell|) 
        # and the x, y and z values of the first frame.
        unit_cell, x, y, z = dcd.readStep()

        # Create all objects from the PDB file. The PDB file must match the
        # the DCD file (same atom order).
        conf = PDBConfiguration(self.pdbFile)

        # |molecules| is a Collection of all objects contained in the PDB file.
        molecules = conf.createAll()

        # The objects are introduced into the universe.
        universe.addObject(molecules)

        # A MMTK trajectory is opened for writing.
        trajectory = Trajectory(universe, self.outputFile, mode = 'w', comment = dcd.title)

        # A frame generator is created.        
        snapshot = SnapshotGenerator(universe, actions=[TrajectoryOutput(trajectory, ["all"], 0, None, 1)])

        # The array that will store the frame configuration.
        conf_array = universe.configuration().array

        # An internal counter to fetch the cell parameters cooresponding to the frame under process in the |xstData| list.
        comp = 0

        # Iterates over the frame stored in the DCD file.
        while True:
            # The coordinates of the current frame are stored in |conf_array|.
            conf_array[:, 0] = x
            conf_array[:, 1] = y
            conf_array[:, 2] = z

            # Case of a periodic universe. Its shape is updated with the cell parameters of the current frame by
            # taking the corresponding value in the |xstData| list. XST datas are more trustable than the values 
            # contained in the DCD file.
            if universe.is_periodic:
                universe.setShape(xstData[comp])

            # The time corresponding to the frame.
            t = step*dt

            step_data = {'time': t, 'step': step}

            # A call to the snapshot generator produces the step corresponding to the current frame.
            snapshot(data = step_data)

            step += step_increment

            # Tries to read the cell parameters and coordinates of the next frame.
            try:
                # They could be read, so store them in |unit_cell|, |x|, |y| and |z| variables.
                unit_cell, x, y, z = dcd.readStep()
                comp += 1

            # The next step could not be read. This is the end of the file.
            except EndOfFile:
                break

        # Close the output trajectory
        trajectory.close()

class VASPConverter(object):
    """Converts a VASP Trajectory into a MMTK NetCDFFile."""

    def __init__(self, xdatcarFile = None, contcarFile = None, atomSymbols = None, timeStep = None, outputFile = "vasp_to_mmtk.nc"):
        """The constructor.

        @param xdatcarFile: the VASP XDATCAR file name of the trajectory to convert.
        @type xdatcarFile: string

        @param contcarFile: the VASP CONTCAR or POSCAR file name of the trajectory to convert.
            Not used in VASP5.
        @type contcarFile: string
        
        @param atomSymbols: List of the element names (string) in the order they appear in the 
            trajectory. Not used for VASP 5.
        @type atomSymbols: list

        @param timeStep: the time step of the MD. Not used for VASP 4.
        @type timeStep: float

        @param outputFile: the name of MMTK NetCDF trajectory output file.
        @type outputFile: string        
        """
        
        # The arguments are copied to instance attribute.
        self.xdatcarFile = xdatcarFile
        self.contcarFile = contcarFile
        self.atomSymbols = atomSymbols
        self.timeStep = timeStep
        self.outputFile = outputFile
                    
        # Do the conversion.
        self.__convert()
        
    def read_contcar(self, contcarFile):
        """Read and extracts information from a CONTCAR file.
        """
        
        # The VASP/CONTCAR is opened for reading.
        try:
            contcar = open(contcarFile, 'r').readlines()
        except:
            raise Error('The CONTCAR file %s could not be opened for reading.' % self.contcarFile)
        
        # The first line is read. It may contains the name of the system.
        self.userInfo = contcar[0]

        # Pick up the scale factor that is the first element of the second line of CONTCAR file.
        self.scaleFactor = float(contcar[1].split()[0])

        # The direct basis vectors.
        self.vectA = Vector([float(v) for v in contcar[2].split()])*Units.Ang*self.scaleFactor
        self.vectB = Vector([float(v) for v in contcar[3].split()])*Units.Ang*self.scaleFactor
        self.vectC = Vector([float(v) for v in contcar[4].split()])*Units.Ang*self.scaleFactor
        
        self.atomNumbers = [int(v) for v in contcar[5].split()]
        
        self.nAtoms = sum(self.atomNumbers)
        
        # This dictionnary stores the number of atoms (values) for each atom element (keys).        
        self.composition = dict(zip(self.atomSymbols, self.atomNumbers))
        
    def __convert(self):
        """Performs the actual conversion."""
        
        # Open the XDATCAR file for reading.
        try:
            xdatcar = open(self.xdatcarFile, 'r').readlines()
        except:
            raise Error('The XDATCAR file %s could not be opened for reading.' % self.xdatcarFile)
        
        # Case of a VASP4 file.        
        if not xdatcar[5].strip():
            
            self.read_contcar(self.contcarFile)

            # The last entry of the second line of the XDATCAR file gives the timestep in seconds. Its is converted in ps.
            self.timeStep = float(xdatcar[1].split()[-1])*Units.s
            
            header = 5

        # Case of a VASP5 file.        
        else:
            
            self.userInfo = xdatcar[0].strip()
            
            # Pick up the scale factor that is the first element of the second line of CONTCAR file.
            self.scaleFactor = float(xdatcar[1].split()[0])

            # The direct basis vectors.
            self.vectA = Vector([float(v) for v in xdatcar[2].split()])*Units.Ang*self.scaleFactor
            self.vectB = Vector([float(v) for v in xdatcar[3].split()])*Units.Ang*self.scaleFactor
            self.vectC = Vector([float(v) for v in xdatcar[4].split()])*Units.Ang*self.scaleFactor
            
            self.atomSymbols = [v.strip() for v in xdatcar[5].split()]
            
            self.atomNumbers = [int(v) for v in xdatcar[6].split()]

            self.nAtoms = sum(self.atomNumbers)
            
            # This dictionnary stores the number of atoms (values) for each atom element (keys).        
            self.composition = dict(zip(self.atomSymbols, self.atomNumbers))
            
            header = 7
                                                    
        # In VASP, simulations are always done with periodic conditions. So creates directly a periodic universe.
        universe = ParallelepipedicPeriodicUniverse((self.vectA, self.vectB, self.vectC), None)
                
        # Loop over the number of element of the system.
        for a_name in self.atomSymbols:
            # Loop over the number of atoms per element.
            for i in range(self.composition[a_name]):
                universe.addObject(Atom(a_name, name = a_name + str(i+1)))

        # If that number differs from the one guessed CONTCAR file, raises an error. There might some unconsistencies 
        # between the two files.
        if self.nAtoms != universe.numberOfAtoms():
            raise Error('Wrong number of atoms.')
                                
        # The MMTK trajectory is opened for writing.
        trajectory = Trajectory(universe, self.outputFile, mode = 'w', comment = self.userInfo)

        # A frame generator is created.                
        snapshot = SnapshotGenerator(universe, actions = [TrajectoryOutput(trajectory,["all"], 0, None, 1)])
        
        bSize = self.nAtoms+1
        
        nFrames = (len(xdatcar) - header)/bSize
        
        for fInd in range(nFrames):
            
            start = header + fInd*bSize + 1
            end = header + (fInd+1)*bSize

            # The time corresponding to the frame.
            t = fInd*self.timeStep
            
            for lInd in range(start, end):

                # Conversion of the coordinates into a Scientific Vector object
                coord = Vector([float(v) for v in xdatcar[lInd].split()[0:3]])
                
                # In VASP the coordinates are given in box coordinates. They must be converted in real coordinates.
                coord = universe.boxToRealCoordinates(coord)

                # The universe is updated with the current position of the atom.
                universe.atomList()[lInd - start].setPosition(coord)

            # The real coordinates are foled then into the simulation box (-L/2,L/2). 
            universe.foldCoordinatesIntoBox()

            # A call to the snapshot generator produces the step corresponding to the current frame.
            snapshot(data = {'time': t})

        # The output trajectory is closed.
        trajectory.close()

class VASPBackConverter(object):
    """Converts a MMTK NetCDF file into a VASP Trajectory."""

    def __init__(self, mmtkFile, contcarFile, xdatcarFile):
        """The constructor.

        @param mmtkFile: the name of MMTK NetCDF trajectory to convert.
        @type mmtkFile: string

        @param contcarFile: the VASP CONTCAR output file name.
        @type contcarFile: string

        @param xdatcarFile: the VASP XDATCAR output file name.
        @type xdatcarFile: string
        """

        # The arguments are copied to instance attribute.
        self.mmtkFile = mmtkFile       
        self.contcarFile = contcarFile
        self.xdatcarFile = xdatcarFile

        # Do the conversion.
        self.__convert()

    def __convert(self):
        """Performs the actual conversion."""

        # The MMTK trajectory is opened.
        trajectory =  Trajectory(None, self.mmtkFile, 'r')

        # The universe
        universe = trajectory.universe

        # The number of atoms in the universe.
        nAtoms = universe.numberOfAtoms()

        # The number of steps of the trajectory.
        nStep = len(trajectory)

        # A nested-list where each nested list is of the form [atom symbol, atom index].
        temp = [(at.symbol.capitalize(),at.index) for at in universe.atomList()]
        # The list is sorted according the atomic symbols.
        temp.sort()

        # The atom symbols.
        atomSymbols = []
        # The atom indexes.
        atomIndexes = []

        # A list that stores the different atom symbols found in the universe.
        symbols = []
        for symbol, index in temp:
            atomSymbols.append(symbol)
            atomIndexes.append(index)
            if symbol not in symbols:
                symbols.append(symbol)

        # A nested-list where each nested-list is of the form [symbol, number of atoms found for that symbol]
        composition = []
        for s in symbols:
            composition.append([s,atomSymbols.count(s)])

        # The brute is formula is used as the system name.
        molName = ''.join([s+str(c) for s,c in composition])

        # The CONTCAR file is opened for writing.
        contcar = open(self.contcarFile, 'w')

        # The first line contains informations about the name of the system and its atomic composition.
        contcar.write('Molecular system %s:' % molName)
        [contcar.write(' %s' % s) for s in symbols]
        contcar.write('\n')

        # The scale factor is set to 1.
        contcar.write(' %22.16f\n' % 1)

        # The cell matrix is written.
        for v in universe.basisVectors():
            contcar.write(' %22.16f%22.16f%22.16f\n' % tuple(v))

        # The number of atoms of each type is written to the CONTCAR file.
        [contcar.write('%4d' % el[1]) for el in composition]        
        contcar.write('\nDirect\n')

        # The last frame is selected.
        universe.setFromTrajectory(trajectory, -1)        
        conf = universe.configuration()

        # then converted to box coordinates and written to the CONTCAR file.
        for ind in atomIndexes:
            bCoord = universe.realToBoxCoordinates(conf[ind])
            contcar.write(' %20.16f%20.16f%20.16f\n' % tuple(bCoord))

        # The CONTCAR file is closed.
        contcar.close()

        # The time step is computed.
        timeStep = trajectory.time[1] - trajectory.time[0]

        # And converted to seconds.
        timeStep /= Units.s

        # The XDATCAR file is opened for writing.
        xdatcar = open(self.xdatcarFile, 'w')

        # The first line of the XDATCAR.
        xdatcar.write('%4d%4d %5d\n' % (nAtoms, nAtoms, nStep))
        xdatcar.write(' %g %g %g %g %g\n' % (1,1,1,1,timeStep))
        xdatcar.write('%21.12f\n' % 1)

        xdatcar.write('  CAR\n')
        xdatcar.write(' unknown system\n\n')

        # All the frame converted to box coordinates are written to the XDATCAR file block by block.
        for i in range(nStep):

            universe.setFromTrajectory(trajectory, i)        
            conf = universe.configuration()

            for ind in atomIndexes:
                bCoord = universe.realToBoxCoordinates(conf[ind])
                xdatcar.write('  %11.8f %11.8f %11.8f\n' % tuple(bCoord))
            xdatcar.write('\n')

        # The XDATCAR file is closed.
        xdatcar.close()

        # The MMTK trajectory is closed.
        trajectory.close()

        return

class PDBConverter(object):
    """Converts a PDB Trajectory into a MMTK NetCDFFile. 

    """

    def __init__(self, pdbFile, outputFile, make_periodic = False):
        """The constructor.

        @param pdbFile: the PDB file name of the single frame to convert.
        @type pdbFile: string
        
        @param make_periodic: if True the trajectory will be made periodic even if no PBC was found in the PDB file.
        @type pdbFile: boolean

        @param outputFile: the MMTK NetCDF trajectory output file name.
        @type outputFile: string

        @note: calling the constructor will actually perform the conversion.
        """

        # The arguments are copied to instance attributes.
        self.pdbFile = pdbFile
        self.make_periodic = make_periodic
        self.outputFile = outputFile

        # Do the conversion.
        self.__convert()

    def __convert(self):
        """Performs the actual conversion."""
                
        try:
            pdb_config = PDBConfiguration(self.pdbFile)
        except:
            raise Error('The PDB file couls not be parsed properly.')

        # Create the universe corresponding to the PDB file. If no PBC is found in the PDB then 
        # the universe is an InfiniteUniverse.
        universe = pdb_config.createUnitCellUniverse()

        # The chemical contents of the PDB file.
        objects = pdb_config.createAll(None, 1)

        # The universe is filled with it.
        universe.addObject(objects)

        # Case where the universe has no PBC and the user decied to make it periodic.
        if self.make_periodic and not universe.is_periodic:
            # The coordinates of the low left corner.
            min_corner = Vector(universe.configuration().array.min(axis = 0))
            # The coordinates of the up right corner.
            max_corner = universe.configuration().array.max(axis = 0)            
            # The diagonal vector between those two corners.
            diff = max_corner - min_corner
            
            # The a, b and c vectors are created.
            a = Vector(diff[0],0,0)
            b = Vector(0,diff[1],0)
            c = Vector(0,0,diff[2])
            
            # The old infinite universe is removed.
            del universe
            
            # An parallelepipedic universe is created out of a, b and c vectors.
            universe = ParallelepipedicPeriodicUniverse([a,b,c])
            
            # And filled with the chemical contents of the PDB file.
            universe.addObject(objects, True)
            
            center = min_corner + Vector(a.length()/2,b.length()/2,c.length()/2)
            
            # Finally the whole universe is translated so as its origin to be in (0,0,0)
            universe.translateBy(-center)

        # A trajectory is opened for writing the single frame.
        trajectory = Trajectory(universe, self.outputFile, mode = "w", comment = "PDB single frame trajectory")

        # A snapshot is created.
        snapshot = SnapshotGenerator(universe, actions = [TrajectoryOutput(trajectory, ["all"], 0, None, 1)])

        # And done with a dummy time value of 1 ps created just for compatibility with 'normal' NetCDF trajectory.
        snapshot(data = {'time': 1.0})

        # The output trajectory is closed.
        trajectory.close()    