"""This modules implements I{File-->Convert NetCDF to ASCII} dialog.

Classes:
    * NetCDFToASCIIConversionDialog: creates I{File-->Convert NetCDF to ASCII} dialog used to 
      convert a file in NetCDF format to a file in ASCII format.
"""

# The python distribution modules
import copy
import inspect
import os
import re
import sys

# The Tcl/Tk modules
from tkFileDialog import askopenfilename
from Tkinter import *

# The ScientificPython modules
from Scientific import N as Num
from Scientific.IO.NetCDF import _NetCDFFile, NetCDFFile

from MMTK.Trajectory import Trajectory

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import convertNetCDFToASCII
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import *

class NetCDFToASCIIConversionDialog(PortableToplevel):
    """Sets up a dialog that allows the conversion of any numeric variables present in a 
    NetCDF file into an ASCII file.
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

        self.netcdf = None

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
                                                 frameLabel = "NetCDF input file",\
                                                 tagName = 'convert_netcdf_to_ascii_netcdf_input_file',\
                                                 contents = '',\
                                                 save = False,\
                                                 command = self.openNetCDFFile,\
                                                 filetypes = [("NetCDF file", ".nc"),])
        self.inputFileBrowser.entry.bind('<Return>', self.openNetCDFFile)
        self.inputFileBrowser.grid(row = 0, column = 0, sticky = EW, padx = 2, pady = 2)
        self.inputFileBrowser.grid_columnconfigure(0, weight = 1)
        
        self.variablesInfo = ComboText(settingsFrame,\
                                       frameLabel = 'Available variables',\
                                       tagName = 'convert_netcdf_to_ascii_available_variables')
        self.variablesInfo.grid(row = 1, column = 0, columnspan = 2, sticky = EW, padx = 2, pady = 2)
        self.variablesInfo.grid_columnconfigure(0, weight = 1)
        self.variablesInfo.text.config({'height' : 10, 'width' : 68})

        self.variablesInfo.text.window_create('1.0',\
                                              window = Label(self.variablesInfo.text,
                                                             text = 'Select',\
                                                             width = 13,\
                                                             bg = 'Grey'))

        self.variablesInfo.text.window_create('1.1',\
                                              window = Label(self.variablesInfo.text,
                                                             text = 'NetCDF Variable',\
                                                             width = 50,\
                                                             anchor = W,\
                                                             bg = 'white'))
        
        # The combo widget for the float precision.
        self.floatPrecisionEntry = ComboIntegerEntry(settingsFrame,\
                                                     frameLabel = 'Float precision',\
                                                     tagName = 'convert_netcdf_to_ascii_float_precision',\
                                                     contents = 9)
        self.floatPrecisionEntry.grid(row = 2, column = 0, sticky = EW, padx = 2, pady = 2)
        self.floatPrecisionEntry.grid_columnconfigure(0, weight = 1)
        
        # The combo widget for the double precision.
        self.doublePrecisionEntry = ComboIntegerEntry(settingsFrame,\
                                                      frameLabel = 'Double precision',\
                                                      tagName = 'convert_netcdf_to_ascii_double_precision',\
                                                      contents = 17)
        self.doublePrecisionEntry.grid(row = 3, column = 0, sticky = EW, padx = 2, pady = 2)
        self.doublePrecisionEntry.grid_columnconfigure(0, weight = 1)
        
        # The combo widget for the file browser.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "ASCII/CDL output file",\
                                                  tagName = 'convert_netcdf_to_ascii_ascii_output_file',\
                                                  contents = '',\
                                                  filetypes = [("CDL file", ".cdl"),],\
                                                  save = True)
        self.outputFileBrowser.grid(row = 4, column = 0, sticky = EW, padx = 2, pady = 2)
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

        self.update_idletasks()

        self.apply()

    def cancel(self, event=None):

        # Put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    # Command hooks
    def validate(self):
            
        # An input file must have been set.
        self.inputFile = self.inputFileBrowser.getValue()
        if not self.inputFile:
            raise Error('Please enter a NetCDF input file.')
                    
        # The float precision entry is checked and stored.
        self.floatPrecision = self.floatPrecisionEntry.getValue()
        if self.floatPrecision <= 0:
            raise Error('The float precision must > 0.')
        
        # The double precision entry is checked and stored.
        self.doublePrecision = self.doublePrecisionEntry.getValue()
        if self.doublePrecision <= 0:
            raise Error('The Double precision must > 0.')

        # An output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter an ASCII output file.')

        # The extension '.cdl' is appended to the output file name if necessary.
        if self.outputFile[-4:] != '.cdl':
            self.outputFile += '.cdl'
            
        # The selected variables are stored.
        self.selectedVariables = []
        for k, v in self.variables.items():
            if v.get():
                self.selectedVariables.append(k)
                
        if not self.selectedVariables:
            raise Error('No NetCDF variables selected for conversion.')
            
        return True
        
    # Command hooks
    def apply(self):
        """
        This method is called when the user clicks on the OK button of the conversion dialog. It performs the 
        conversion from the loaded NetCDF file to the selected ASCII/CDL file.
        """
        
        convertNetCDFToASCII(inputFile = self.inputFile,\
                             outputFile = self.outputFile,\
                             variables = self.selectedVariables,\
                             floatPrecision = self.floatPrecision,\
                             doublePrecision = self.doublePrecision)
            
        LogMessage('info', 'Conversion successful', ['gui'])
        
    
    def openNetCDFFile(self, event = None):
        """
        This method opens a NetCDF file and updates the dialog with the data read from that file.
        Arguments:
            -event: Tkinter event.
        """

        self.variablesInfo.text.config(state = NORMAL)
        self.variablesInfo.text.delete('2.0', END)
        self.variablesInfo.text.config(state = DISABLED)

        # Case where the user enters a file name directly in the entry widget without using the browser.
        if event is not None:
            if event.widget == self.inputFileBrowser.entry:
                filename = self.inputFileBrowser.getValue()
            else:
                return
            
        else:
        
            # The name of the NetCDF file to load.
            filename = askopenfilename(parent = self,\
                                       filetypes = [('NetCDF file','*.nc')],\
                                       initialdir = PREFERENCES['outputfile_path'])
            
        # The file must exist otherwise do nothing.
        if filename:
            self.netcdf = NetCDFFile(filename, 'r')
                    
            self.displayNetCDFContents()

            # The filebrowser entry is updated with the loaded filename.
            self.inputFileBrowser.setValue(filename)

            # A default filename for the output ASCII file is built.
            self.outputFileBrowser.setValue(os.path.splitext(filename)[0] + '.cdl')                
                
        return 'break'
    
    def selectVariable(self):
        
        for k, v in self.variables.items():
            if v.get():
                self.textVariables[k].set('selected')
                
            else:
                self.textVariables[k].set('unselected')
        

    def displayNetCDFContents(self):
        """
        This method display the variables found in the NetCDF file.
        """

        self.variables = {}                
        self.textVariables = {}

        self.variablesInfo.text.config(state = NORMAL)
        self.variablesInfo.text.delete('2.0', END)
        self.variablesInfo.text.config(state = DISABLED)

        comp = 0
        for varName in sorted(self.netcdf.variables.keys()):
            self.variables[varName] = BooleanVar()
            self.variables[varName].set(True)
            
            self.textVariables[varName] = StringVar()
            self.textVariables[varName].set('selected')
            
            self.variablesInfo.insert(END, '\n')

            self.variablesInfo.text.window_create(str(comp + 2)+'.0',\
                                                  window = Checkbutton(self.variablesInfo.text,\
                                                                       variable = self.variables[varName],\
                                                                       textvariable = self.textVariables[varName],\
                                                                       command = self.selectVariable,\
                                                                       anchor = W,\
                                                                       justify = LEFT,\
                                                                       bd = 1,\
                                                                       width = 10,\
                                                                       bg = 'Grey'))

            self.variablesInfo.text.window_create(str(comp + 2)+'.1',\
                                                  window = Label(self.variablesInfo.text,\
                                                                 text = varName,\
                                                                 width = 50,\
                                                                 anchor = W,\
                                                                 bg = 'white'))
            comp += 1
