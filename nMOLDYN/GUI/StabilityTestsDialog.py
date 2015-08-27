"""This modules implements I{Help-->nMOLDYN benchmark} dialog.

Classes:
    * StabilityTestsDialog: creates I{Help-->nMOLDYN benchmark} dialog used to run the stability tests 
      between nMOLDYN v2.2.5 and the current version of nMOLDYN.
"""

# The python distribution modules
import copy
import glob
import os
import re

# The Tcl/Tk modules
from Tkinter import *
from tkMessageBox import showinfo

# The nMOLDYN modules
from nMOLDYN.__pkginfo__ import __version__ as NMOLDYN_VERSION
from nMOLDYN.Analysis.Templates import *
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Core.Error import Error
from nMOLDYN.GUI.Widgets import PortableToplevel

ALL_TESTS = ('ARA','ARDCSF','ARDISF','AVACF','ADOS','CDOS','CVACF','DCSF','DISF','DISFGA','EISF','MSD','RCF')

class StabilityTestsDialog(PortableToplevel):
    """Builds the dialog for nMOLDYN analysis benchmarks.
    
    These benchmarks will guarantee the stability of the current version with a reference one that is the 
    version 2.2.5 the last official release of nMOLDYN.
    
    @note: these benchmarks should be performed every time the code source of nMOLDYN has been touched for 
        any reason.
    """

    def __init__(self, parent, title = None):
        """The constructor.
        
        @param parent: the parent widget.
        
        @param title: a string specifying the title of the dialog.
        @type title: string
        """
        
        # The base class is initialized.
        PortableToplevel.__init__(self, parent)

        # Will store all the sensible widgets ids.
        self.widgets = {}
        
        # Will store the name of the currently selected test.
        self.selectedTest = None

        # The dialog will be transiented to the parent dialog.                              
        self.transient(parent)
                
        # A title for the dialog is set.
        if title:
            self.title(title)
                    
        # This will make the dialog resizable.
        self.rowconfigure(0, weight = 1)                    
        self.columnconfigure(0, weight = 1)                    
                    
        # A copy of the parent dialog id.
        self.parent = parent
        
        # The body of the dialog is built.          
        self.body()       

        # The buttonbar of the dialog is built.
        self.buttonbox()        
        
        # The dialog grabs all incoming events.
        self.grab_set()

        # The cancel method is bound to the dialog closure event.
        self.protocol("WM_DELETE_WINDOW", self.cancel)

        # The geometry of the dialog.
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50, parent.winfo_rooty()+50))

        # The dialog takes the focus.
        self.focus_set()
        
        # The resizeDialog method is bound to any resizing/moving event.
        self.bind("<Configure>", self.resizeDialog)
        
        # Th dialog is updated.
        self.update_idletasks()

        # The minimum sizes for the dialog.
        self.minWidth  = 500
        self.minHeight = 500
        
        # The dialog event loop is started.
        self.wait_window(self)            

    def body(self):
        """Creates the dialog body.
        """
        
        self.widgets["main"] = Frame(self)
        self.widgets["main"].grid(row = 0, column = 0, sticky = 'nswe')  
        self.widgets["main"].rowconfigure(1, weight = 1000)  
        self.widgets["main"].columnconfigure(0, weight = 1)  
        
        testsFrame = LabelFrame(self.widgets["main"], text = "Tests Settings", bd = 2, relief = GROOVE)
        testsFrame.grid(row = 0, column = 0, sticky = 'nswe', padx = 3, pady = 3)
        testsFrame.columnconfigure(4, weight = 1000)
        
        self.widgets["testnames"] = Listbox(testsFrame, width = 15, exportselection = 0)
        self.widgets["testnames"].grid(row = 0, column = 0, sticky = "nswe", padx = 3, pady = 3)
        sb = Scrollbar(testsFrame, orient = VERTICAL)
        sb.grid(row = 0, column = 1, sticky = "nswe")
        self.widgets["testnames"].configure(yscrollcommand = sb.set)
        [self.widgets["testnames"].insert(END, t) for t in ALL_TESTS]
        self.widgets["testnames"].bind('<ButtonRelease-1>', self.selectTest)
        sb['command'] = self.widgets["testnames"].yview

        self.widgets["testids"] = Listbox(testsFrame, width = 10, selectmode = MULTIPLE)
        self.widgets["testids"].grid(row = 0, column = 2, sticky = "nswe", padx = 3, pady = 3)
        sb = Scrollbar(testsFrame, orient = VERTICAL)
        sb.grid(row = 0, column = 3, sticky = "nswe")
        self.widgets["testids"].configure(yscrollcommand = sb.set)
        sb['command'] = self.widgets["testids"].yview
        self.widgets["testids"].bind('<ButtonRelease-1>', self.selectTestIds)
        self.widgets["testids"].bind('<Control-a>', self.selectTestIds)

        self.widgets["testcontents"] = Text(testsFrame, height = 10, yscrollcommand = sb.set, state = DISABLED)
        self.widgets["testcontents"].grid(row = 0, column = 4, sticky = "nswe", padx = 3, pady = 3)
        sb = Scrollbar(testsFrame, orient = VERTICAL)
        sb.grid(row = 0, column = 5, sticky = "nswe")
        self.widgets["testcontents"].configure(yscrollcommand = sb.set)
        sb['command'] = self.widgets["testcontents"].yview
        sb = Scrollbar(testsFrame, orient = HORIZONTAL)
        sb.grid(row = 1, column = 4, sticky = "nswe")
        self.widgets["testcontents"].configure(xscrollcommand = sb.set)
        sb['command'] = self.widgets["testcontents"].xview

        resultsFrame = LabelFrame(self.widgets["main"], text = 'Results', bd = 2, relief = GROOVE)
        resultsFrame.grid(row = 1, column = 0, sticky = 'nswe', padx = 3, pady = 3)
        resultsFrame.rowconfigure(0, weight = 1)  
        resultsFrame.columnconfigure(0, weight = 1)  
        
        self.widgets["testresults"] = Text(resultsFrame, state = DISABLED)
        self.widgets["testresults"].grid(row = 0, column = 0, sticky = "nswe", padx = 3, pady = 3)
        sb = Scrollbar(resultsFrame, orient = VERTICAL)
        sb.grid(row = 0, column = 1, sticky = "nswe")
        self.widgets["testresults"].configure(yscrollcommand = sb.set)
        sb['command'] = self.widgets["testresults"].yview
        sb = Scrollbar(resultsFrame, orient = HORIZONTAL)
        sb.grid(row = 1, column = 0, sticky = "nswe")
        self.widgets["testresults"].configure(yscrollcommand = sb.set)
        sb['command'] = self.widgets["testresults"].xview
                        
    def buttonbox(self):
        """Creates the dialog button bar."""

        # The frame that contains the 'Cancel' and 'OK' buttons.
        self.widgets["buttons"] = LabelFrame(self, text = 'Actions', bd = 2, relief = GROOVE)
        self.widgets["buttons"].grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3)
        self.widgets["buttons"].columnconfigure(0, weight = 1000)
        self.widgets["buttons"].columnconfigure(1, weight = 1)

        w = Button(self.widgets["buttons"], text = "Cancel", width=10, command = self.cancel)
        w.grid(row = 0, column = 0, sticky = E)
                
        w = Button(self.widgets["buttons"], text = "OK", width=10, command = self.ok, default=ACTIVE)
        w.grid(row = 0, column = 1, sticky = E)
        
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

    def ok(self, event = None):
        """Called when OK button is pressed. Will start the tests.
        """
                
        # The test contents text widget is enabled for writing.                
        self.widgets["testresults"].config(state = NORMAL)
        self.widgets["testresults"].delete('0.0', END)
        self.widgets["testresults"].insert(END, "Running tests. Please wait ...")
        self.widgets["testresults"].config(state = DISABLED)
        
        # Updates the dialogs.        
        self.update_idletasks()
        
        # The listbox index of the selected test name.
        selectedTest = self.widgets["testnames"].curselection()
        
        # If no tests was selected, just return.
        if not selectedTest:
            return

        # The index is converted to the actual selected test name.
        selectedTest = self.widgets["testnames"].get(int(selectedTest[0]))

        # The corresponding analysis file is read.
        info = {}
        exec file(os.path.join(GVAR['nmoldyn_analysis'], selectedTest)) in None, info

        # The listbox indexes of the tests ids.
        selectedIds = self.widgets["testids"].curselection()

        # If no tests ids was selected, just return.
        if not selectedIds:
            return
                
        # The index is converted to integer.
        selectedIds = [int(self.widgets["testids"].get(i)) for i in selectedIds]
        
        # The analysis corresponding to the selected test is initialized.
        analysis = eval('%s_serial()' % info["db_internalname"])
                
        # And run in test mode.
        results = analysis.test(selectedIds)
        
        summary = results.summary()                        

        # The test contents text widget is enabled for writing.                
        self.widgets["testresults"].config(state = NORMAL)                
        self.widgets["testresults"].delete('0.0', END)
        self.widgets["testresults"].insert(END, summary)                
        self.widgets["testresults"].config(state = DISABLED)
        
    def cancel(self, event = None):
        """Called when the dialog is closed.
        """

        # Gives the focus back to the parent window
        self.parent.focus_set()
        
        # Close the dialog.
        self.destroy()

    def selectTest(self, event):
        """Called when the user selects a test name. It will display the available test ids for the 
        selected test.
        """
           
        # The listbox index of the selected test name.
        selectedTest = self.widgets["testnames"].curselection()
        
        # If no tests was selected, just return.
        if not selectedTest:
            return
                        
        # The index is converted to the actual selected test name.
        self.selectedTest = self.widgets["testnames"].get(int(selectedTest[0]))
                        
        # The corresponding analysis file is read.
        info = {}
        exec file(os.path.join(GVAR['nmoldyn_analysis'], self.selectedTest)) in None, info

        # Display the tests ids available for the selected test.
        self.widgets["testids"].delete(0, END)
        nTests = len(glob.glob(os.path.join(GVAR['nmoldyn_tests'], info["db_shortname"],"*_Reference.py")))
        [self.widgets["testids"].insert(END, i) for i in range(1, nTests + 1)]
        
    def selectTestIds(self, event):
        """Called when the user selects a set of test ids.
        """
                   
        # If there is no selected test name, do nothing.     
        if not self.selectedTest:
            return
                        
        if event.keysym == "a":
            self.widgets["testids"].selection_set(0, END)

        # The listbox indexes of the tests ids.
        selectedIds = self.widgets["testids"].curselection()

        # The test contents text widget is enabled for writing.                
        self.widgets["testcontents"].config(state = NORMAL)
        self.widgets["testcontents"].delete('0.0', END)
        
        nTests = len(selectedIds)
        
        # If only one id was selected, display in the test contents text widget the input file for the selected test.
        if nTests == 1:
            selectedIds = [int(self.widgets["testids"].get(i)) for i in selectedIds]
            testContents = os.path.join(GVAR['nmoldyn_tests'], self.selectedTest, "%s%d_Current.py" % (self.selectedTest, selectedIds[-1]))                
            self.widgets["testcontents"].insert(END, open(testContents, "r").read())
        
        # Otherwise, just insert an informative message.    
        else:
            self.widgets["testcontents"].insert(END, "Click on a single test id to see its contents.")

        self.widgets["testcontents"].config(state = DISABLED)
                                                
    def resizeDialog(self, event):
        """Called when the user resizes the dialog. It will fit the widgets of the dialog to the new size.
        """
        
        # Check that the resizing event comes from the main dialog.
        if event.widget == self:
                                    
            # The new dimensions can not be smaller than the minimum dimensions. 
            if event.width < self.minWidth or event.height < self.minHeight:
                self.geometry(newGeometry = "%dx%d" % (self.minWidth, self.minHeight))
                
                        