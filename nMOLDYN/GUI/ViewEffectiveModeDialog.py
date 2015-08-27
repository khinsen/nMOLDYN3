"""This modules implements I{View --> Effective Mode} dialog.

Classes:
    * ViewEffectiveModeDialog: creates I{View --> Effective Mode} dialog used to 
      view an animation of the effective modes resulting from a Quasi Harmonic Analysis.
"""

# The python distribution modules
import copy
import os
import sys

# The Tcl/Tk modules
from tkFileDialog import askopenfilename
from Tkinter import *

# The ScientificPython modules
from Scientific import N as Num
from Scientific.IO.NetCDF import NetCDFFile

# The MMTK distribution modules
from MMTK import Configuration
from MMTK import Skeleton
from MMTK.ParticleProperties import ParticleVector
from MMTK.Visualization import definePDBViewer, viewSequenceVMD

# The nMOLDYN modules
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import ComboIntegerEntry, ComboFloatEntry, ComboFileBrowser, ComboListbox, ComboRadiobutton, PortableToplevel

class ViewEffectiveModeDialog(PortableToplevel):
    """Sets up a dialog used to visualize the effective modes resulting from a QHA analysis.
    """

    def __init__(self, parent, title = None):
        """The constructor.
        
        @param parent: the parent widget.
        
        @param title: a string specifying the title of the dialog.
        @type title: string
        """

        PortableToplevel.__init__(self, parent)        
        self.transient(parent)
        
        if title:
            self.title(title)

        self.parent = parent       
                        
        body = Frame(self)
        self.initial_focus = self.body(body)
        body.grid(row = 0, column = 0, sticky = EW)  

        self.buttonbox()        
        
        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.resizable(width = NO, height = NO)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50, parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)            

    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """

        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget for the file browser.
        self.fileBrowser = ComboFileBrowser(settingsFrame,\
                                            frameLabel = "QHA input file",\
                                            tagName = 'view_effective_modes_qha_input_file',\
                                            contents = '',\
                                            save = False,\
                                            command = self.openNetCDFFile,\
                                            filetypes = [("NetCDF file", ".nc"),])
        self.fileBrowser.grid(row = 0, column = 0, sticky = EW, padx = 2, pady = 2)
        self.fileBrowser.grid_columnconfigure(0, weight = 1)
        self.fileBrowser.entry.bind('<Return>', self.openNetCDFFile)

        # The combo listbox that will contain the X variables.
        self.selectedModeLb = ComboListbox(settingsFrame,\
                                           frameLabel = 'Quasi-Harmonic mode',\
                                           tagName = 'quasi_harmonic_mode',\
                                           contents = [])
        self.selectedModeLb.lb.config({'exportselection' : 0, 'width' : 22, 'height' : 8, 'selectmode' : MULTIPLE})
        self.selectedModeLb.grid(row = 1, column = 0, sticky = EW, padx = 2, pady = 2)
        self.selectedModeLb.grid_columnconfigure(0, weight = 1)
        
        # The combo widget to set the number of frames for the animation.
        self.nFramesEntry = ComboIntegerEntry(settingsFrame,\
                                              frameLabel = 'Number of frames',\
                                              tagName = 'view_effective_modes_number_of_frames')
        self.nFramesEntry.grid(row = 2, column = 0, sticky = EW, padx = 2, pady = 2)
        self.nFramesEntry.grid_columnconfigure(0, weight = 1)

        # The combo widget to set the amplitude of the effective mode to view.
        self.amplitudeEntry = ComboFloatEntry(settingsFrame,\
                                              frameLabel = 'Amplitude (in nm)',\
                                              tagName = 'view_effective_modes_amplitude')
        self.amplitudeEntry.grid(row = 3, column = 0, sticky = EW, padx = 2, pady = 2)
        self.amplitudeEntry.grid_columnconfigure(0, weight = 1)
        
        return None        

    def buttonbox(self):
        """
        Add standard button box.
        """

        # The frame that contains the 'Cancel' and 'OK' buttons.
        box = LabelFrame(self, text = 'Actions', bd = 2, relief = GROOVE)
        box.grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3)
        box.grid_columnconfigure(0, weight = 1)

        w = Button(box, text = "Cancel", width=10, command = self.cancel)
        w.grid(row = 0, column = 0, sticky = E)
        w = Button(box, text = "OK", width=10, command = self.ok, default=ACTIVE)
        w.grid(row = 0, column = 1, sticky = E)
        
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

    # Standard button semantics.
    def ok(self, event = None):

        if not self.validate():
            self.initial_focus.focus_set()
            return

        self.update_idletasks()

        self.apply()
        
    def cancel(self, event=None):

        # Put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    # Command hooks
    def validate(self):
        
        try:
                        
            if not self.selectedModeLb.lb.curselection():
                LogMessage('warning','Please select a vibration mode.',['gui'])
                raise

            self.selectedMode = [int(v) - 1 for v in self.selectedModeLb.lb.curselection()]

            self.amplitude = self.amplitudeEntry.getValue()
            self.nFrames = self.nFramesEntry.getValue()
            
            if self.amplitude <= 0.0:
                raise
            
            if self.nFrames < 0:
                raise
                        
        except:
            LogMessage('warning','Bad input. Please try again.',['gui'])
            return False
            
        return True
    
    def apply(self):

        try:
            if os.path.exists(PREFERENCES['vmd_path']):
                definePDBViewer('vmd', PREFERENCES['vmd_path'])
            
            else:
                raise
            
        except:
            raise Error('Error when defining the PDB viewer from %s path.' % PREFERENCES['vmd_path'])
            
        try:
        
            local = {}
            skeleton = eval(self.description, vars(Skeleton), local)
            universe = skeleton.make({}, self.avgStruct)
            universe.setCellParameters(self.cell)

            avg = Configuration(universe, self.avgStruct)
            pseudoTraj = [avg]
            
            for comp in range(self.nFrames):
                vibr = copy.copy(avg)
                for selMode in self.selectedMode:
                    dx = copy.copy(self.dx[selMode])
                    dx.shape = (universe.numberOfAtoms(), 3)
                    d = ParticleVector(universe, dx)
                    vibr += self.amplitude*Num.sin(2.0*Num.pi*float(comp)/self.nFrames)*d
                    
                pseudoTraj.append(vibr)

            viewSequenceVMD(universe, pseudoTraj, periodic = 1)
            
        except:
            raise Error('Error when animating the selected mode(s).')
        
    def openNetCDFFile(self, event = None):
        """
        This method open the NetCDF that contains the effective modes.
        Arguments:
            - event: Tkinter event.
        """

        # Case where the user enters a file name directly in the entry widget without using the browser.
        if event is not None:
            if event.widget == self.fileBrowser.entry:
                filename = self.fileBrowser.getValue()
            else:
                return
            
        else:
            # The name of the NetCDF file to load.
            filename = askopenfilename(parent = self,\
                                       filetypes = [('NetCDF file','*.nc')],\
                                       initialdir = PREFERENCES['trajfile_path'])

        # The file must exist.
        if filename:
            try:
                self.netcdf = NetCDFFile(filename, 'r')
                
            except IOError:
                LogMessage('warning','Problem when reading the NetCDF file.',['gui'])
                self.fileBrowser.setValue('')
                self.selectedModeEntry.setValue('')
                self.nFramesEntry.setValue('')
                self.amplitudeEntry.setValue('')
                
            try:
                self.description = self.netcdf.variables['description'][:].tostring()

                # The frequencies values.
                self.omega = self.netcdf.variables['omega'].getValue()
                
                # The displacements.
                self.dx = self.netcdf.variables['dx'].getValue()
                
                # The average structure.
                self.avgStruct = self.netcdf.variables['avgstruct'].getValue()
                                
            except KeyError:
                LogMessage('warning','The NetCDF file %s miss some QHA analysis keywords.' % filename,['gui'])
                
            self.fileBrowser.setValue(filename)

            self.selectedModeLb.lb.delete(0, END)

            for i in range(len(self.omega)):
                ome = self.omega[i]
                self.selectedModeLb.lb.insert(END, 'Mode %s (%s cm-1)' % (i+1, ome))

            self.nFramesEntry.setValue(10)
            self.amplitudeEntry.setValue(0.1)

            try:
                self.cell = self.netcdf.variables['box_size'][:]
            except KeyError:
                self.cell = None
                        
            self.netcdf.close()

        return 'break'
