"""This modules implements I{View-->Animation} dialog.

Classes:
    * PyroServerDialog: creates the dialog used to set up the Pyro server.
"""

# The python distribution modules
import os
import subprocess
import sys
import traceback

# The Tcl/Tk modules
from tkFileDialog import askopenfilename
from Tkinter import *

# The nMOLDYN modules
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Core.Monitoring import getCPUInfo
from nMOLDYN.GUI.Widgets import ComboFileBrowser, ComboLabel, ComboRadiobutton, ComboSpinbox, ComboText, PortableToplevel

class PyroServerDialog(PortableToplevel):
    """Sets up a dialog used to set up the Pyro server. 
    """

    def __init__(self, parent, widget):
        """The constructor.
        
        @param parent: the parent widget.
        """

        PortableToplevel.__init__(self, parent)        
        self.transient(parent)
        
        self.title('Pyro server')

        self.parent = parent
        
        self.selectedParallelMode = self.defaultSelectedParallelMode = widget.getValue()
                        
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
        
        self.settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        self.settingsFrame.grid(row = 0, column = 0, sticky = 'NEW', padx = 3, pady = 3)
        self.settingsFrame.grid_columnconfigure(0, weight = 1)

        self.parallelMode = ComboRadiobutton(self.settingsFrame,\
                                                 frameLabel = "Running mode",\
                                                 tagName = 'pyro_server',\
                                                 contents = ["monoprocessor", "multiprocessor", "cluster"],\
                                                 layout = (3,1))
        
        for rb in self.parallelMode.radio:
            rb.config(command = self.refreshParallelModeInfo)
            
        self.parallelMode.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = 'NSEW')
        self.parallelMode.grid_rowconfigure(10, weight = 1)
                
        self.parallelModeInfo = ComboText(self.settingsFrame,\
                                          frameLabel = 'Informations')
        self.parallelModeInfo.text.config({'width' : 60, 'height' : 15, 'wrap' : WORD})
        self.parallelModeInfo.grid(row = 0, column = 1, sticky = 'NSEW', padx = 2, pady = 2)
        self.parallelModeInfo.grid_columnconfigure(0, weight = 1)
        
        self.refreshParallelModeInfo()
        
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
    def ok(self):
                
        if not self.validate():
            self.initial_focus.focus_set()
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        # Put focus back to the parent window
        self.parent.focus_set()
        
        # The parent takes all the forthcoming events.
        self.parent.grab_set()
        
        self.destroy()
        
    def cancel(self, event=None):

        self.selectedParallelMode = self.defaultSelectedParallelMode

        # Put focus back to the parent window
        self.parent.focus_set()
        
        # The parent takes all the forthcoming events.
        self.parent.grab_set()
        
        self.destroy()

    # Command hooks
    def validate(self):        
        return True
    
    def apply(self):
        """Builds the pyro server configuration string. This can be:
            -'monoprocessor' for monoprocessor running mode
            -'multiprocessor::hostname:number of allocated cpus'
            -'cluster' for cluster running mode.
        """

        self.selectedParallelMode = self.parallelMode.getValue()

        if self.parallelMode.getValue() == 'multiprocessor':
            self.selectedParallelMode += ':%s' % str(self.numberOfSelectedProcs.get())
        
    def getValue(self):
        """Returns the pyro server specification string.
        """

        return self.selectedParallelMode        
                                                                                
    def refreshParallelModeInfo(self):
        """Updates the text widget that contains the informations about the pyro server.
        """
        
        self.parallelModeInfo.cleanup()
                    
        if self.parallelMode.getValue() == 'monoprocessor':            
            self.parallelModeInfo.insert(END, 'The analysis will be run in monoprocessor mode without the use of \
Pyro server.')
                    
        elif self.parallelMode.getValue() == 'multiprocessor':

            self.parallelModeInfo.insert(END, 'The analysis will be run in multiprocessor mode on your local machine.\n\n\
The Pyro settings will be automatically done once the number of processors you would like to allocate to run the analysis \
will be selected.\n\n')

            nTotal, nLoaded, nFree = getCPUInfo()
                        
            self.parallelModeInfo.text.window_create(END,\
                                                     window = Label(self.parallelModeInfo.text,
                                                                    text = '# cpus',\
                                                                    width = 15,\
                                                                    bg = 'white',\
                                                                    anchor = CENTER))
            self.parallelModeInfo.text.window_create(END,\
                                                     window = Label(self.parallelModeInfo.text,
                                                                    text = '# free cpus',\
                                                                    width = 15,\
                                                                    bg = 'white',\
                                                                    anchor = CENTER))
            
            self.parallelModeInfo.text.window_create(END,\
                                                     window = Label(self.parallelModeInfo.text,
                                                                    text = '# allocated cpus',\
                                                                    width = 15,\
                                                                    bg = 'white',\
                                                                    anchor = CENTER))
            
            self.parallelModeInfo.insert(END, '\n')
            
            self.parallelModeInfo.text.window_create(END,\
                                                     window = Label(self.parallelModeInfo.text,
                                                                    text = str(nTotal),\
                                                                    width = 15,\
                                                                    relief = GROOVE,\
                                                                    anchor = CENTER))

            self.parallelModeInfo.text.window_create(END,\
                                                     window = Label(self.parallelModeInfo.text,
                                                                    text = str(nFree),\
                                                                    width = 15,\
                                                                    relief = GROOVE,\
                                                                    anchor = CENTER))

            selProcs = [str(v) for v in range(1, nTotal + 1)]
            
            self.numberOfSelectedProcs = Spinbox(self.parallelModeInfo.text,\
                                                 values = tuple(selProcs),\
                                                 wrap = True,\
                                                 bg = 'white',\
                                                 width = 15)
            
            self.parallelModeInfo.text.window_create(END, window = self.numberOfSelectedProcs)

            [self.numberOfSelectedProcs.invoke('buttonup') for v in range(nFree - 1)]

            self.parallelModeInfo.insert(END, '\n')

        elif self.parallelMode.getValue() == 'cluster':
            self.parallelModeInfo.insert(END,
"""The analysis will be run in cluster mode.

You will have to:

* save this analysis as a <scriptname>

* run the command

\tpython <scriptname> -t <taskname>

  on one node of your cluster

* run the command:

\ttask_manager -w slave <taskname>

  once per compute process on some node of your cluster.

In the above command lines, <scriptname> is the name you give to
the analysis script and <taskname> is a unique name for the analysis
that must be different from the names of all other tasks running
at the same time.""")
