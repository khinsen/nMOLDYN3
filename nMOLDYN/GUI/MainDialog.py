"""This is where the main window of nMOLDYN is defined.

Classes:
    * MainDialog: The class that defines the nMOLDYN GUI main window and its associated actions.
"""

# The python distribution modules.
import glob
import inspect
import platform
import os
import re
import subprocess
import sys
import webbrowser

# The ScientificPython modules.
from Scientific import N
from Scientific.IO.NetCDF import NetCDFFile, _NetCDFFile

# The MMTK distribution modules.
from MMTK.Trajectory import Trajectory, trajectoryInfo, TrajectorySet

from Tkinter import *
from tkFileDialog import askopenfilename
from tkMessageBox import askyesno

# The nMOLDYN modules.
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Preferences import PREFERENCES, savePreferencesFile
from nMOLDYN.Analysis.Dynamics import *
from nMOLDYN.Analysis.NMR import *
from nMOLDYN.Analysis.Scattering import *
from nMOLDYN.Analysis.Structure import *
from nMOLDYN.Chemistry.Chemistry import hierarchizeUniverse
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import load_netcdf_file, load_trajectory_file
from nMOLDYN.Core.Logger import LogCreate, LogMessage
from nMOLDYN.Core.Monitoring import checkForNewVersion, downloadVersion, searchProgram

from nMOLDYN.GUI.AnalysisDialog import AnalysisDialog
from nMOLDYN.GUI.AnimationDialog import AnimationDialog
from nMOLDYN.GUI.StabilityTestsDialog import StabilityTestsDialog
from nMOLDYN.GUI.ASCIIToNetCDFConversionDialog import ASCIIToNetCDFConversionDialog
from nMOLDYN.GUI.PreferencesDialog import PreferencesDialog
from nMOLDYN.GUI.GeneralInformationsDialog import GeneralInformationsDialog
from nMOLDYN.GUI.NetCDFToASCIIConversionDialog import NetCDFToASCIIConversionDialog
from nMOLDYN.GUI.PDBSnapshotGeneratorDialog import PDBSnapshotGeneratorDialog
from nMOLDYN.GUI.TrajectoryConversionDialog import *
from nMOLDYN.GUI.ViewEffectiveModeDialog import ViewEffectiveModeDialog
from nMOLDYN.GUI.Widgets import ComboText

# The platform name.
PLATFORM = platform.system().upper()

class MainDialog(Tk):
    """
    This is the base class for the nMoldyn GUI. It launches the main window of nMoldyn from which
    different menus can be accessed.
    """

    def __init__(self, netcdf_filename=None):
        """
        The constructor. Displays the main window.

        @param netcdf_filename: a string specifying a NetCDF file, nMOLDYN should be started with.
        @rtype: string
        """
                
        Tk.__init__(self)

        self.option_add('*Dialog.msg.wrapLength', '30i')

        if PLATFORM == 'WINDOWS':
            self.wm_iconbitmap(default=os.path.join(GVAR['nmoldyn_path'], 'win32_files', 'nMOLDYN_Logo.ico'))

        self.title('nMOLDYN')

        self.netcdf_filename = netcdf_filename

        # Loads all the analysis setup from their respective ASCII files.
        self.loadAnalysisDatabase()

        self.initial_focus = self.body(self)

        self.grab_set()

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        # Tcl/Tk specificities for MacOS platform.
        if PLATFORM == 'DARWIN':
            try:
                self.tk.call('package', 'require', 'tile')
                self.tk.call('namespace', 'import', '-force', 'ttk::*')
                self.tk.call('tile::setTheme', 'aqua')
            except TclError:
                pass

        self.resizable(width=NO, height=NO)

        self.geometry("+20+20")

        self.focus_set()

        LogCreate(['console','file','gui'])

        # Search for the presence of external used (but not compulsory) by nMOLDYN.
        self.search_external_programs()

        # Wait 1s before checking for missing program.
        self.after(1000, self.define_valid_menus)
        
        if PREFERENCES['warning_check_for_new_version'] == 'yes':
            checkForNewVersion(fromGUI=True)
        
        self.wait_window(self)

    def body(self, master):

        mainMenu = Menu(master)

        if PLATFORM == 'DARWIN':
            accel = 'Cmd-'
            bindKey = 'Command-'
            menu = Menu(mainMenu, name='apple')
            mainMenu.add_cascade(label='nMOLDYN', menu=menu)
            menu.add_command(label='About nMOLDYN', command=self.aboutNMOLDYN)
            
        else:
            accel = 'Ctrl-'
            bindKey = 'Control-'

        # ####################
        # The 'File' Menu.
        # ####################
        self.fileMenu = Menu(mainMenu, tearoff=0)

        # Opens a netcdf file.
        self.fileMenu.add_command(label='Load NetCDF file', command=self.load_input_file, accelerator=accel + 'O')
        master.bind('<%s>' % (bindKey + 'o',), self.load_input_file)

        self.fileMenu.add_separator()

        conversionMenu = Menu(self.fileMenu, tearoff=0)
        self.fileMenu.add_cascade(label='Trajectory conversion', menu=conversionMenu)
        conversionMenu.add_command(label='Amber NetCDF to MMTK', command=lambda : AmberNetCDFConverterDialog(master))
        conversionMenu.add_command(label='CASTEP to MMTK', command=lambda : CASTEPConverterDialog(master))
        conversionMenu.add_command(label='CHARMM/X-PLOR to MMTK', command=lambda : CHARMMConverterDialog(master))
        conversionMenu.add_command(label='DL_POLY to MMTK', command=lambda : DL_POLYConverterDialog(master))
        conversionMenu.add_command(label='LAMMPS to MMTK', command=lambda : LAMMPSConverterDialog(master))

        materialsStudioMenu = Menu(conversionMenu, tearoff=0)
        conversionMenu.add_cascade(label='MaterialsStudio', menu=materialsStudioMenu)
        materialsStudioMenu.add_command(label='Discover to MMTK', command=lambda : MaterialsStudioConverterDialog(master, 'Discover'))
        materialsStudioMenu.add_command(label='Forcite to MMTK', command=lambda : MaterialsStudioConverterDialog(master, 'Forcite'))

        conversionMenu.add_command(label='NAMD to MMTK', command=lambda : NAMDConverterDialog(master))

        conversionMenu.add_command(label='VASP to MMTK', command=lambda : VASPConverterDialog(master))
        
        conversionMenu.add_separator()

        conversionMenu.add_command(label='MMTK to VASP4', command=lambda : VASPBackConverterDialog(master))
        
        conversionMenu.add_separator()
        
        conversionMenu.add_command(label='PDB to MMTK', command=lambda : PDBConverterDialog(master))

        # Opens a dialog to extract a selected frame from a given trajectory into a PDB file.
        self.fileMenu.add_command(label='Frame snapshot', command=self.extractTrajectoryFrame, accelerator=accel + 'R')
        master.bind('<%s>' % (bindKey + 'r',), self.extractTrajectoryFrame)

        self.fileMenu.add_separator()

        # Opens the NetCDF to ASCII file conversion dialog.
        self.fileMenu.add_command(label='Convert NetCDF to ASCII', command=self.convertNetCDFToASCII, accelerator=accel + 'U')
        master.bind('<%s>' % (bindKey + 'u',), self.convertNetCDFToASCII)

        # Opens the NetCDF to ASCII file conversion dialog.
        self.fileMenu.add_command(label='Convert ASCII to NetCDF', command=self.convertASCIIToNetCDF, accelerator=accel + 'I')
        master.bind('<%s>' % (bindKey + 'i',), self.convertASCIIToNetCDF)

        self.fileMenu.add_separator()

        # This open the configuration file editor used to set up some of the nMOLDYN environment variables.
        self.fileMenu.add_command(label='Preferences', command=self.setPreferences, accelerator=accel + 'E')
        master.bind('<%s>' % (bindKey + 'e',), self.setPreferences)

        self.fileMenu.add_separator()

        # Close the nMOLDYN apllication.
        self.fileMenu.add_command(label='Quit', command=self.cancel, accelerator=accel + 'Q')
        master.bind('<%s>' % (bindKey + 'q',), self.cancel)

        # This submenu is popped up when clicking on 'File' item of the main menu.
        mainMenu.add_cascade(label="File", menu=self.fileMenu)

        # #####################
        # The 'Analysis' Menu.
        # #####################
        self.analysisMenu = Menu(mainMenu, tearoff=0)

        self.buildAnalysisMenu(self.analysisMenu, self.menuHierarchy)

        # This submenu is popped up when clicking on |Analysis| item of the main menu.
        mainMenu.add_cascade(label="Analysis", menu=self.analysisMenu)

        for i in range(len(self.menuHierarchy)):
            self.analysisMenu.entryconfig(i, state=DISABLED)

        # #####################
        # The 'View' Menu.
        # #####################
        self.viewMenu = Menu(mainMenu, tearoff=0)

        # The 'Plot' entry. Launches the plot dialog for NetCDF file plotting.
        self.viewMenu.add_command(label='Plot', command=self.plotNetCDF, accelerator=accel + 'L')
        master.bind('<%s>' % (bindKey + 'l',), self.plotNetCDF)

        # The 'Animation' entry. Launches the Animation dialog.
        self.viewMenu.add_command(label='Animation', command=self.animateTrajectory, accelerator=accel + 'M')
        master.bind('<%s>' % (bindKey + 'm',), self.animateTrajectory)

        # The 'View Effective mode'. Launches the effective mode viewer dialog.
        self.viewMenu.add_command(label='Effective mode', command=self.viewEffectiveMode, accelerator=accel + 'K')
        master.bind('<%s>' % (bindKey + 'k',), self.viewEffectiveMode)

        # This submenu is popped up when clicking on |View| item of the main menu.
        mainMenu.add_cascade(label="View", menu=self.viewMenu)

        # #####################
        # The 'Help' Menu.
        # #####################
        self.helpMenu = Menu(mainMenu, tearoff=0)

        # This is the item 'nMOLDYN users guide'. Clicking on itopen the nMOLDYN users guide pdf file.
        self.helpMenu.add_command(label="Documentation", command=self.displayDocumentation, accelerator=accel + 'H')
        master.bind('<%s>' % (bindKey + 'h',), self.displayDocumentation)

        # This is the item 'nMOLDYN mailig list'.
        self.helpMenu.add_command(label="Mailing List", command=self.displayMailingList)

        # This is the item 'Technical notes'. Clicking on it opens the nMOLDYN API html files.
        self.helpMenu.add_command(label='API', command=self.displayAPI, accelerator=accel + 'T')
        master.bind('<%s>' % (bindKey + 't',), self.displayAPI)

        self.helpMenu.add_separator()

        # This is the item 'nMOLDYN benchmark'. Launches a toplevel window where the user can perform some analysis
        # benchmarks.
        self.helpMenu.add_command(label='Analysis benchmark', command=self.analysisBenchmark)

        self.helpMenu.add_separator()

        # This is the item 'nMOLDYN benchmark'. Launches a toplevel window where the user can perform some analysis
        # benchmarks.
        self.helpMenu.add_command(label='Check for updates', command=self.checkForUpdates)

        self.helpMenu.add_separator()

        self.helpMenu.add_command(label='Version history', command=self.displayVersionHistory, accelerator=accel + 'V')
        
        # This is the item 'About nMOLDYN'.
        # Clicking on it pops up a toplevel window containing some general informations about nMOLDYN.
        if PLATFORM != 'DARWIN':
            self.helpMenu.add_command(label='About nMOLDYN', command=self.aboutNMOLDYN, accelerator=accel + 'B')
        master.bind('<%s>' % (bindKey + 'b',), self.aboutNMOLDYN)

        # This submenu is popped up when clicking on |Help| item of the main menu.
        mainMenu.add_cascade(label="Help", menu=self.helpMenu)

        # Displays the menu
        master.config(menu=mainMenu)

        mainPanel = Frame(master)
        mainPanel.grid(row=0, column=0)

        # Display of an image of a neutron on the left side of the main window.
        self.neutronGIF = PhotoImage(file=os.path.join(GVAR['nmoldyn_path'], 'Images', 'Neutron.gif'))

        leftFrame = Frame(mainPanel)
        leftFrame.grid(row=0, column=0)
        Label(leftFrame, image=self.neutronGIF).grid(row=0)
        Label(leftFrame, image=self.neutronGIF).grid(row=1)
        Label(leftFrame, image=self.neutronGIF).grid(row=2)
        Label(leftFrame, image=self.neutronGIF).grid(row=3)
        Label(leftFrame, image=self.neutronGIF).grid(row=4)

        self.info = ComboText(mainPanel)
        self.info.config({'width' : 500, 'relief' : SUNKEN})
        self.info.text.config({'wrap' : WORD, 'width' : 80})
        self.info.grid(row=0, column=1, sticky="NSEW")

        self.info.text.tag_config('12', font=('Courier', 12))
        self.info.text.tag_config('12B', font=('Courier', 12, 'bold'))
        self.info.text.tag_config('r12B', foreground='red', font=('Courier', 12, 'bold'))
        self.info.text.tag_config('12I', font=('Courier', 12, 'normal', 'italic'))

        # Display of an image of a protein on the left side of the main window.
        self.proteinGIF = PhotoImage(file=os.path.join(GVAR['nmoldyn_path'], 'Images', 'Protein.gif'))
        rightFrame = Frame(mainPanel)
        rightFrame.grid(row=0, column=2)
        Label(rightFrame, image=self.proteinGIF).grid(row=0)
        Label(rightFrame, image=self.proteinGIF).grid(row=1)
        Label(rightFrame, image=self.proteinGIF).grid(row=2)
        Label(rightFrame, image=self.proteinGIF).grid(row=3)
        Label(rightFrame, image=self.proteinGIF).grid(row=4)

        mainPanel.grid_rowconfigure(0, weight=1)
        self.info.grid_rowconfigure(0, weight=1)
        self.info.text.grid_rowconfigure(0, weight=1)

        if self.netcdf_filename is not None:
            self.load_input_file(filename=self.netcdf_filename)
        
        return None

    def cancel(self, event=None):
        self.destroy()
        sys.exit(0)

    def load_input_file(self, event=None, filename=None):
        """
        This method is launched when the user clicks on the |Load NetCDF| of the |File| menu.
        It loads the NetCDF file and displays its main informations in the information window.
        """

        if filename is None:
            # A file dialog box is opened for browsing the 'nc' or 'ncs' files.
            # A 'ncs' file is a text file with the names of several 'nc' files.
            filename = askopenfilename(filetypes=[("NetCDF trajectory file", "*.nc"), ("NetCDF trajectory set", "*.ncs")], title="Open NetCDF",
                                       initialdir=PREFERENCES['trajfile_path'])

        # Do nothing if no trajectory has been selected.
        if not filename:
            return

        # The information window is cleaned up.
        self.info.cleanup()
                        
        traj = load_trajectory_file(filename)
                
        if isinstance(traj, Trajectory):

            GVAR['current_traj'] = filename

            # Detects the occurence of chemical fragments or atoms in the whole universe and indicates whether they are deuterable or not.
            self.chemicalObjectInfo, self.chemicalObjectNames = hierarchizeUniverse(traj.universe)

            try:

                firstStep = traj.time[0]
                lastStep = traj.time[-1]
                timeStep = traj.time[1] - traj.time[0]

            except:
                firstStep = 0.0
                lastStep = float(len(traj))
                timeStep = 1.0

            self.info.insert(contents='FILE TYPE: MMTK NETCDF TRAJECTORY FILE\n\n', tag='12B')

            # Some trajectory informations are displayed on the information window.
            self.info.insert(contents=trajectoryInfo(filename), tag='12B')
            self.info.insert(contents='\n\nMMTK objects found in the universe: ', tag='12')

            # Write into the info panel the details about the nmoldyn names found in the universe of the loaded trajectory.
            for k, v in sorted(self.chemicalObjectInfo.items()):
                # Do not consider the '*' nmoldyn name. Not interesting.
                if k == '*':
                    continue
                self.info.insert(contents='\n\t- %d %s (%s)' % (v['number'], k, v['objectclass']), tag='12')

            self.info.insert(contents='\n\nNumber of frames: ', tag='12')
            self.info.insert(contents=str(len(traj)), tag='12B')
            self.info.insert(contents='\nStarting at: ', tag='12')
            self.info.insert(contents=str(round(firstStep, 3)) + ' ps', tag='12B')
            self.info.insert(contents='\nEnding at: ', tag='12')
            self.info.insert(contents=str(round(lastStep, 3)) + ' ps', tag='12B')
            self.info.insert(contents='\nTime step: ', tag='12')
            self.info.insert(contents=str(round(timeStep, 3)) + ' ps', tag='12B')
            self.info.insert(contents='\n\nUniverse size: ', tag='12')
            self.info.insert(contents=str(len(traj.universe.atomList())), tag='12B')
            self.info.insert(contents='\n\nUniverse topology: ', tag='12')
            self.info.insert(contents=traj.universe.__class__.__name__, tag='12B')

            directBasis = traj.universe.basisVectors()
            if directBasis is not None:
                self.info.insert(contents='\n\nDirect basis:', tag='12')
                self.info.insert(contents='\n%8s\t%8s\t%8s' % tuple(['i', 'j', 'k']), tag='12I')
                for directVector in directBasis:
                    self.info.insert(contents='\n%8.3f\t%8.3f\t%8.3f' % tuple(directVector), tag='12I')

            reciprocalBasis = traj.universe.reciprocalBasisVectors()
            if reciprocalBasis is not None:
                self.info.insert(contents='\n\nReciprocal basis:', tag='12')
                self.info.insert(contents='\n%8s\t%8s\t%8s' % tuple(['i', 'j', 'k']), tag='12I')
                for reciprocalVector in reciprocalBasis:
                    self.info.insert(contents='\n%8.3f\t%8.3f\t%8.3f' % tuple(reciprocalVector), tag='12I')

            for i in range(len(self.menuHierarchy)):
                self.analysisMenu.entryconfig(i, state=NORMAL)
                
        elif isinstance(traj, _NetCDFFile):
            
            netcdf = load_netcdf_file(filename)
            
            if netcdf is None:
                raise Error("Error when loading netcdf file %s" % filename)

            for i in range(len(self.menuHierarchy)):
                self.analysisMenu.entryconfig(i, state=DISABLED)

            self.info.insert(contents='FILE TYPE: NETCDF FILE\n\n', tag='12B')

            # The list of the NetCDF string attributes that can be visualized to get informations about the loaded NetCDF file.
            # The list is initialized with the first one that should appear if present in the NetCDF file.
            netCDFStringAttributes = ['title', 'history', 'comment', 'jobinfo', 'jobinput']
            for attr in inspect.getmembers(netcdf):
                if isinstance(attr[1], str) & (attr[0] not in netCDFStringAttributes):
                    netCDFStringAttributes.append(attr[0])

            self.info.insert(contents='Information about NetCDF file:\n', tag='12B')
            self.info.insert(contents='%s.\n\n' % filename, tag='12')

            # Loop over the netcdf file string attributes that will be displayed for information purpose.
            for h in netCDFStringAttributes:
                # Checks that the NetCDF has this attributes.
                if hasattr(netcdf, h):
                    self.info.insert(contents='%s:\n' % h.capitalize(), tag='12B')
                    self.info.insert(contents=getattr(netcdf, h) + '\n\n', tag='12')

            if hasattr(netcdf, 'dimensions'):
                self.info.insert(contents='\n\nNetCDF dimensions found:', tag='12B')
                for d in netcdf.dimensions.items():
                    self.info.insert(contents='\n\n\t-%s ----> ' % d[0], tag='12')
                    self.info.insert(contents='%s' % d[1], tag='r12B')

            if hasattr(netcdf, 'variables'):
                self.info.insert(contents='\n\nNetCDF variables found:', tag='12B')
                for v in netcdf.variables.items():
                    self.info.insert(contents='\n\n\t-%s:' % v[0], tag='12B')
                    self.info.insert(contents='\n\t\t+ Dimension(s) ----> ', tag='12')
                    self.info.insert(contents=str(v[1].shape), tag='r12B')
                    if hasattr(v[1], 'typecode'):
                        self.info.insert(contents='\n\t\t+ Typecode ----> ', tag='12')
                        self.info.insert(contents=v[1].typecode(), tag='r12B')
                    else:
                        self.info.insert(contents='\n\t\t+ Typecode ----> ', tag='12')
                        self.info.insert(contents='Undefined', tag='r12B')

                    if hasattr(v[1], 'units'):
                        self.info.insert(contents='\n\t\t+ Units ----> ', tag='12')
                        self.info.insert(contents=v[1].units, tag='r12B')
                    else:
                        self.info.insert(contents='\n\t\t+ Units ----> ', tag='12')
                        self.info.insert(contents='Undefined', tag='r12B')

                    try:
                        if v[1].typecode() == N.Character:
                            continue

                        if len(v[1].getValue().shape) > 2:
                            raise

                        mini = min(v[1].getValue().ravel())
                        maxi = max(v[1].getValue().ravel())

                    except:
                        continue

                    else:
                        self.info.insert(contents='\n\t\t+ Range ----> ', tag='12')
                        self.info.insert(contents='%.5f to %.5f' % (mini, maxi), tag='r12B')

        else:
            raise Error("Error when opening the file %s" % filename)

    def extractTrajectoryFrame(self, event=None):
        """
        This method pops up a dialog from where the user can extract a PDB file from a NetCDF
        trajectory frame.
        """

        PDBSnapshotGeneratorDialog(self, 'PDB frame extraction')

    def convertNetCDFToASCII(self, event=None):
        """
        This method pops up a dialog where the user can proceed to a conversion from a file in NetCDF format
        to ASCII format.
        """

        NetCDFToASCIIConversionDialog(self, 'NetCDF to ASCII conversion')

    def convertASCIIToNetCDF(self, event=None):
        """
        This method pops up a dialog where the user can proceed to a conversion from a file in ASCII format
        to NetCDF format.
        """

        ASCIIToNetCDFConversionDialog(self, 'ASCII/CDL to NetCDF conversion')

    def search_external_programs(self):
        """This method will search the path for the executables of external programs used by nMOLDYN that are not defined in the
        Preferences.
        """

        # Sets the dictionnary storing the external programs to search with their associated executable name.
        if PLATFORM == 'WINDOWS':
            externalPrograms = {'vmd' : 'vmd.exe', 'ncdump' : 'ncdump.exe', 'ncgen' : 'ncgen.exe', 'acroread' : 'AcroRd32.exe'}

        elif PLATFORM == 'DARWIN':
            externalPrograms = {'vmd' : 'vmd', 'ncdump' : 'ncdump', 'ncgen' : 'ncgen'}

        else:
            externalPrograms = {'vmd' : 'vmd', 'ncdump' : 'ncdump', 'ncgen' : 'ncgen', 'acroread' : 'acroread'}

        # if the external programs path are not set, nMOLDYN tries to find them in classical locations.
        for prog, progName in externalPrograms.items():
            
            path = PREFERENCES[prog + '_path']

            if not os.path.isfile(path):

                LogMessage('info', 'Searching for %s external executable. Please wait ...' % prog, ['console'])

                progPath = searchProgram(progName)

                if progPath is not None:
                    PREFERENCES[prog + '_path'] = progPath
                    LogMessage('info', '%s found in %s. Your preferences will be modified for the running session.' % (prog, progPath), ['console'])

                else:
                    LogMessage('info', '%s not found.' % prog, ['console'])
                
            else:
                LogMessage('info', 'The path for %s is %s' % (prog,path), ['console'])

    def define_valid_menus(self):
        """This method will (un)validate some of the GUI menu depending on the value of some
        configuration and preferences variables.
        """
        
        missingPrograms = []

        try:
            import nMOLDYN.GUI.PlotNetCDFVariableDialog
        except:
            missingPrograms.append(['matplotlib library','Plot menu disabled'])
            self.viewMenu.entryconfig(0, state=DISABLED)
        else:
            self.viewMenu.entryconfig(0, state=NORMAL)
            
        # Check that the path to ncdump executable is set.    
        if os.path.exists(PREFERENCES['ncdump_path']):
            self.fileMenu.entryconfig(5, state=NORMAL)

        else:
            missingPrograms.append(['ncdump','NetCDF->ASCII conversion disabled'])
            self.fileMenu.entryconfig(5, state=DISABLED)

        # Check that the path to ncgen executable is set.    
        if os.path.exists(PREFERENCES['ncgen_path']):
            self.fileMenu.entryconfig(6, state=NORMAL)

        else:
            missingPrograms.append(['ncgen','ASCII->NetCDF conversion disabled'])
            self.fileMenu.entryconfig(6, state=DISABLED)

        # On Linux, check that the path to acrobat reader executable is set. 
        if (PLATFORM == 'LINUX') & (PREFERENCES['documentation_style'] == 'pdf'):
            if os.path.exists(PREFERENCES['acroread_path']):
                self.helpMenu.entryconfig(0, state=NORMAL)

            else:
                missingPrograms.append(['acrobat reader','Documentation not viewable in pdf format.'])
                self.helpMenu.entryconfig(0, state=DISABLED)
                
        # Check that the path to vmd executable is set.    
        if os.path.exists(PREFERENCES['vmd_path']):
            self.viewMenu.entryconfig(1, state=NORMAL)
            self.viewMenu.entryconfig(2, state=NORMAL)

        else:
            missingPrograms.append(['vmd', 'Trajectory viewer disabled.'])
            self.viewMenu.entryconfig(1, state=DISABLED)
            self.viewMenu.entryconfig(2, state=DISABLED)
                                                                
    def setPreferences(self, event=None):
        """
        This method pops up a dialog from where the user can edit the nMOLDYN configuration file.
        """

        PreferencesDialog(self, 'Preferences')
        
        self.define_valid_menus()

    def analysisDialog(self, analysis):
        """
        """
        
        AnalysisDialog(self, analysis, self.chemicalObjectInfo)

    def plotNetCDF(self, event=None):
        """
        This method pops up a dialog from where the user can display any numeric 2D or 3D NetCDF variables.
        """

        try:
            from nMOLDYN.GUI.PlotNetCDFVariableDialog import PlotNetCDFVariableDialog
        except:
            return
        else:
            PlotNetCDFVariableDialog(self, 'Plot NetCDF')

    def animateTrajectory(self, event=None):
        """
        This method pops up a dialog from where the user can animate a trajectory. If a trajectory has been loaded
        for analysis this will be the default one. Otherwise the user can still browse one from the dialog.
        The animation requires VMD.
        """

        AnimationDialog(self, 'Trajectory animation')

    def viewEffectiveMode(self, event=None):
        """
        This method pops up a dialog from where the user can animate an effective mode coming from a QHA
        analysis.
        The animation require VMD.
        """

        ViewEffectiveModeDialog(self, 'Effective mode viewer')

    def traceAnalysis(self, event=None):
        """
        This method pops up a dialog from where the user can check the march of the running jobs.
        The dialog can be updated dynamically by pressing its button 'refresh'.
        """

        CheckJobsStatusDialog(self, 'Job progress')

    def analysisBenchmark(self):
        """
        This method pops up a dialog from where the user can perform some analysis benchmark. The benchmark is
        done between the current version and a reference version that is the version 2.2.5 the last official release
        of nMOLDYN.
        """
        
        if not os.path.exists(os.path.join(GVAR['nmoldyn_path'], 'Tests')):
            LogMessage('warning', 'No tests available on your nMOLDYN installation.', ['gui'])
            return

#        AnalysisBenchmarkDialog(self, 'Analysis benchmark')
        StabilityTestsDialog(self, 'Analysis benchmark')

    def checkForUpdates(self, event=None):
        """This method will search for a new version of nMOLDYN and propose it for download if one is found.
        """

        newVersion = checkForNewVersion()
        
        if newVersion is not None:

            q = askyesno(title='Question', \
                         message='nMOLDYN %s is available.\n\n\
Do you want to download it ?' % newVersion)
            if q:                                
                downloadVersion(newVersion)
                
        else:
            LogMessage('info', 'No new version available', ['gui'])

    def displayVersionHistory(self, event=None):
        """This method displays the nMOLDYN CHANGELOG file on a Tkinter Text widget.
        """
        
        # The name of the CHANGELOG file.
        changeLogFile = os.path.join(GVAR['nmoldyn_path'], 'Doc', 'CHANGELOG')
        changeLog = ''.join(file(changeLogFile).read().split('\r'))
        
        # Create a new window with a text widget showing the file contents.
        new_window = Toplevel()
        new_window.title("Version history")
        text = Text(new_window, wrap=WORD, width=100)
        text.grid(column=0, row=0, sticky=W + N + S)
        scrollbar = Scrollbar(new_window)
        scrollbar.grid(column=1, row=0, sticky=E + N + S)
        scrollbar.config(command=text.yview)
        text.config(yscrollcommand=scrollbar.set)
        new_window.grid_columnconfigure(0, weight=1)
        new_window.grid_columnconfigure(1, weight=0)
        new_window.grid_rowconfigure(0, weight=1)
        text.config(state=NORMAL)
        text.insert(END, changeLog)
        text.config(state=DISABLED)

    def displayDocumentation(self, event=None):
        """
        This methode opens the nMOLDYN pdf users guide.
        The users guide was written by E. Pellegrini, V. Calandrini, P. Calligari, K. Hinsen and G.R. Kneller.
        """

        if not os.path.exists(os.path.join(GVAR['nmoldyn_path'], 'Doc', 'UsersGuide')):
            LogMessage('warning', 'No documentation found on your nMOLDYN installation.', ['gui'])
            return

        if PREFERENCES['documentation_style'] == 'html':
            webbrowser.open(os.path.join(GVAR['nmoldyn_path'], 'Doc', 'UsersGuide', 'HTML', 'index.html'))

        elif PREFERENCES['documentation_style'] == 'pdf':
            # The name of the users guide pdf file.
            pdfFile = os.path.join(GVAR['nmoldyn_path'], 'Doc', 'UsersGuide', 'PDF', 'nMOLDYN_ug.pdf')

            # Case of Win32
            if PLATFORM == 'WINDOWS':
                try:
                    # This will find automatically the appropriated application to open the PDF file. 
                    os.startfile(pdfFile)
                except:
                    raise Error('The file %s could not be opened. There might be no PDF reader installed.' % pdfFile)

            # Case of MacOS
            elif PLATFORM == 'DARWIN':
                try:
                    # This will find automatically the appropriated application to open the PDF file. 
                    subprocess.call(['open', pdfFile])
                except:
                    raise Error('The file %s could not be opened. There might be no PDF reader installed.' % pdfFile)
            else:

                try:
                    subprocess.call([PREFERENCES['acroread_path'], pdfFile])
                except:
                    raise Error('Unable to read the file %s\nwith the acrobat reader located in %s.' % (pdfFile, PREFERENCES['acroread_path']))
        else:
            raise Error('Unknow documentation format: %s.' % PREFERENCES['documentations_style'])

    def displayMailingList(self):
        """
        This methode opens the nMOLDYN mailing list.
        """

        webbrowser.open('http://sourcesup.cru.fr/forum/?group_id=194')

    def displayAPI(self, event=None):

        if not os.path.exists(os.path.join(GVAR['nmoldyn_path'], 'Doc', 'API')):
            LogMessage('warning', 'No API found on your nMOLDYN installation.', ['gui'])
            return

        if PREFERENCES['documentation_style'] == 'html':
            webbrowser.open(os.path.join(GVAR['nmoldyn_path'], 'Doc', 'API', 'HTML', 'index.html'))

        elif PREFERENCES['documentation_style'] == 'pdf':
            # The name of the users guide pdf file.
            pdfFile = os.path.join(GVAR['nmoldyn_path'], 'Doc', 'API', 'PDF', 'api.pdf')

            # Case of Win32
            if PLATFORM == 'WINDOWS':
                os.startfile(pdfFile)

            # Case of MacOS
            elif PLATFORM == 'DARWIN':
                subprocess.call(['open', pdfFile])

            else:

                try:
                    subprocess.call([PREFERENCES['acroread_path'], pdfFile])
                except:
                    raise Error('Unable to read the file %s\nwith the acrobat reader located in %s.' % (pdfFile, PREFERENCES['acroread_path']))

        else:
            raise Error('Unknow documentation format: %s.' % PREFERENCES['documentations_style'])

    def aboutNMOLDYN(self, event=None):
        """
        This method displays general informations about the program such as the developper, the main versions ...
        """

        GeneralInformationsDialog(self, 'About nMOLDYN')

    def loadAnalysisDatabase(self):
        """
        """

        path = os.path.join(GVAR['nmoldyn_path'], 'Database', 'Analysis', '*')
        db = glob.glob(path)

        self.menuHierarchy = {}

        for f in db:
            p = {}
            exec file(f) in None, p

            if not p.has_key('db_menu'):
                continue

            temp = self.menuHierarchy
            for m in p['db_menu']:
                if not temp.has_key(m):
                    temp[m] = {}
                temp = temp[m]

            else:
                temp[p['db_longname']] = f

    def buildAnalysisMenu(self, menu, hierarchy):
        """
        """

        for key in sorted(hierarchy.keys()):

            value = hierarchy[key]

            if isinstance(value, dict):
                submenu = Menu(menu, tearoff=0)
                menu.add_cascade(label=key, menu=submenu)
                self.buildAnalysisMenu(submenu, value)

            else:
                menu.add_command(label=key, command=lambda a=value: self.analysisDialog(a))


