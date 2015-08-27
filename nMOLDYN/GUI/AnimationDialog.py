"""This modules implements I{View-->Animation} dialog.

Classes:
    * AnimationDialog: creates I{View-->Animation} dialog used to view an animation of a trajectory.
"""

# The python distribution modules
import os
import sys

# The Tcl/Tk modules
from tkFileDialog import askopenfilename
from Tkinter import *

# The ScientificPython modules

# The MMTK distribution modules
from MMTK.Trajectory import Trajectory
from MMTK.Visualization import definePDBViewer, viewTrajectory

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import ComboFileBrowser, ComboIntegerEntry, ComboRadiobutton, PortableToplevel

class AnimationDialog(PortableToplevel):
    """Sets up a dialog used to visualize a loaded trajectory. 
    
    @note: if a trajectory has been previously loaded in nMOLDYN for simulation purposes this will be the one 
    proposed for visualization by default. Otherwise, the user can still choose a trajectory to visualize from
    the dialog.
    """

    def __init__(self, parent, title = None):
        """The constructor.
        
        @param parent: the parent widget.
        
        @param title: a string specifying the title of the dialog.
        @type title: string

        @param trajectory: an opened MMTK trajectory.
        @type trajectory: an instance of MMTK.Trajectory.Trajectory class.
        """
        
        PortableToplevel.__init__(self, parent)        
        self.transient(parent)
        
        if title:
            self.title(title)

        self.parent = parent       

        self.trajectory = GVAR['current_traj']

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
                                            frameLabel = "MMTK Trajectory file",\
                                            tagName = 'animation_mmtk_trajectory_file',\
                                            contents = '',\
                                            save = False,\
                                            command = self.openTrajectory)
        self.fileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.fileBrowser.grid_columnconfigure(0, weight = 1)
        self.fileBrowser.entry.bind('<Return>', self.openTrajectory)

        # The combo widget for the first step setting.
        self.firstStepEntry = ComboIntegerEntry(settingsFrame,\
                                                frameLabel = 'First step',\
                                                tagName = 'animation_first_step')
        self.firstStepEntry.grid(row = 1, column = 0, sticky = EW, padx = 2, pady = 2)
        self.firstStepEntry.grid_columnconfigure(0, weight = 1)

        # The combo widget for the last step setting.
        self.lastStepEntry = ComboIntegerEntry(settingsFrame,\
                                               frameLabel = 'Last step',\
                                               tagName = 'animation_last_step')
        self.lastStepEntry.grid(row = 2, column = 0, sticky = EW, padx = 2, pady = 2)
        self.lastStepEntry.grid_columnconfigure(0, weight = 1)

        # The combo widget for the skip step setting.
        self.skipStepEntry = ComboIntegerEntry(settingsFrame,\
                                               frameLabel = 'Skip step',\
                                               tagName = 'animation_skip_test')
        self.skipStepEntry.grid(row = 3, column = 0, sticky = EW, padx = 2, pady = 2)
        self.skipStepEntry.grid_columnconfigure(0, weight = 1)

        if self.trajectory is not None:
            # The control variables are updated with the informations coming from the loaded trajectory.
            self.fileBrowser.setValue(self.trajectory)
            self.firstStepEntry.setValue(1)
            self.lastStepEntry.setValue(1)
            self.skipStepEntry.setValue(1)

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

        self.withdraw()
        self.update_idletasks()

        self.apply()
        
        self.cancel()

    def cancel(self, event=None):

        # Put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    # Command hooks
    def validate(self):
                    
        if not self.trajectory:
            raise Error('No MMTK trajectory file loaded for animation.')
            
        if self.trajectory != self.fileBrowser.getValue():
            raise Error('Mismatch between the loaded trajectory and the displayed MMTK trajectory file name.\
You should validate that file by pressing Return on its corresponding entry.')
            
        self.viewableAtoms = self.trajectory.universe.atomsWithDefinedPositions()
            
        if not self.firstStepEntry.getValue():
            self.first = 0
            LogMessage('warning', 'First frame undefined. Frame 1 will be taken by default.', ['gui'])
        else:
            self.first = self.firstStepEntry.getValue() - 1

        if self.first < 0:
            raise Error('First frame not within [1,%d].' % len(self.trajectory))

        if not self.lastStepEntry.getValue():
            self.last = 1
            LogMessage('warning', 'Last frame undefined. Frame 1 will be taken by default.', ['gui'])
        else:
            self.last = self.lastStepEntry.getValue()

        if self.last > len(self.trajectory):
            raise Error('Last frame not within [1,%d].' % len(self.trajectory))

        if not self.skipStepEntry.getValue():
            self.skip = 1
            LogMessage('warning', 'Frame skip undefined. 1 will be taken by default.', ['gui'])
        else:
            self.skip = self.skipStepEntry.getValue()

        if (self.skip < 1) or (self.skip >= len(self.trajectory)):
            raise Error('Last frame not within [1,%d[.' % len(self.trajectory))
            
        return True
    
    def apply(self):

        try:
            if os.path.exists(PREFERENCES['vmd_path']):
                definePDBViewer('vmd', PREFERENCES['vmd_path'])            
            else:
                raise            
        except:
            raise Error('Error when defining VMD molecular viewer from %s path.' % PREFERENCES['vmd_path'])

        try:
            # This MMTK function calls VMD to visualize the trajectory. 
            viewTrajectory(trajectory = self.trajectory, first = self.first, last = self.last, skip = self.skip,\
                           subset = self.viewableAtoms)
        except:
            raise Error('Error when animating the trajectory.')
    
    def openTrajectory(self, event = None):
        """
        The method is called when the user clicks on the 'Browse' button of the trajectory visualization dialog.
        It opens a file browser. After the file selection some of the dialog widgets are updated with the informations
        coming from the loaded trajectory.
        
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
                # The trajectory is loaded.                
                self.trajectory = Trajectory(None, filename, 'r')
                
            except IOError:
                LogMessage('warning', 'Problem when reading the trajectory.',['gui'])
                self.fileBrowser.setValue('')
                self.firstStepEntry.setValue('')
                self.lastStepEntry.setValue('')
                self.skipStepEntry.setValue('')
                
            else:
                # The control variables are updated with the informations about the loaded trajectory.
                self.fileBrowser.setValue(filename)
                self.firstStepEntry.setValue(1)
                self.lastStepEntry.setValue(1)
                self.skipStepEntry.setValue(1)
                
        return 'break'
