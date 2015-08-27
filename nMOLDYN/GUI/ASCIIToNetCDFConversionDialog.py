"""This modules implements I{File-->Convert ASCII to NetCDF} dialog.

Classes:
    * ASCIIToNetCDFConversionDialog: creates I{File-->Convert ASCII to NetCDF} dialog used to 
      convert a file in ASCII format to a file in NetCDF format.
"""

# The python distribution modules
import copy
import os
import sys

# The Tcl/Tk modules
from tkFileDialog import askopenfilename
from Tkinter import *

# The ScientificPython modules
from Scientific.IO.NetCDF import NetCDFFile

# The nMOLDYN modules
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import convertASCIIToNetCDF
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import *

class ASCIIToNetCDFConversionDialog(PortableToplevel):
    """Sets up a dialog from where the user can convert a file with numeric data in ASCII or CDL format to NetCDF format.
    
    The ASCII file may contain some comments introduced with the # character. These comments will also be written 
    in the NetCDF output file (|comment| attribute). The numeric datas have to be organized by column. The only 
    restriction is that all the columns should have the same length.
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
        self.inputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                 frameLabel = "ASCII input file",\
                                                 tagName = 'convert_ascii_to_netcdf_ascii_input_file',\
                                                 contents = '',\
                                                 save = False,\
                                                 command = self.openASCIIFile)
        self.inputFileBrowser.grid(row = 0, column = 0, sticky = EW, padx = 2, pady = 2)
        self.inputFileBrowser.grid_columnconfigure(0, weight = 1)
        self.inputFileBrowser.entry.bind('<Return>', self.openASCIIFile)

        # The combo widget for the file browser.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "NetCDF output file",\
                                                  tagName = 'convert_ascii_to_netcdf_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True,\
                                                  filetypes = [("NetCDF file", ".nc")])
        self.outputFileBrowser.grid(row = 3, column = 0, sticky = EW, padx = 2, pady = 2)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
        
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
        # An intput file must have been set.
        self.inputFile = self.inputFileBrowser.getValue()
        if not self.inputFile:
            raise Error('Please enter an ASCII input file.')

        self.outputFile = self.outputFileBrowser.getValue()        
        # An output file must have been set.
        if not self.outputFile:
            raise Error('Please enter a NetCDF output file.')
                        
        return True

    # Command hooks
    def apply(self):
        """
        This method is called when the user clicks on the OK button of the conversion dialog. It performs the 
        conversion from the loaded NetCDF file to the selected ASCII file.
        """

        convertASCIIToNetCDF(self.inputFile, self.outputFile)
        LogMessage('info', 'Conversion successful', ['gui'])

    def openASCIIFile(self, event = None):
        """
        This method/callback is called when the user press Return on the entry of the input file browser 
        or browse directlry from the file browser. It will set the filebrowser entry to the name of the browsed
        file and propose and set a name for the output file based on the basename of the browsed file.
        """

        # Case where the user enters a file name directly in the entry widget without using the browser.
        if event is not None:
            if event.widget == self.inputFileBrowser.entry:
                filename = self.inputFileBrowser.getValue()
            else:
                return
            
        else:
            # The name of the ASCII file to load.
            filename = askopenfilename(parent = self,\
                                       filetypes = [('CDL file','.cdl'), ('All files', '.*')],\
                                       initialdir = PREFERENCES['outputfile_path'])

        if filename:            
            # The filebrowser entry is updated with the loaded filename.
            self.inputFileBrowser.setValue(filename)

            # A default filename for the output NetCDF file is built.
            self.outputFileBrowser.setValue(os.path.splitext(filename)[0] + '.nc')
            
        return 'break'
