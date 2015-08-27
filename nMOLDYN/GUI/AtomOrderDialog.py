"""This modules implements the dialog used to select the atoms in a specific order.
"""

# The python distribution modules
import os
import re
import sys

# The Tcl/Tk modules
from Tkinter import *
from tkFileDialog import askopenfilename

# The nMOLDYN modules
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import ComboFileBrowser, ComboFrame, ComboListbox, ComboRadiobutton, ComboText, PortableToplevel
            
class AtomOrderDialog(PortableToplevel):
    """Sets up a dialog to perform the select the atom in a specific order.
            
    @note: inherits from Tkinter.PortableToplevel class.
    """

    def __init__(self, parent, atomNames):
        """The constructor.
        
        @param parent: the parent dialog.
        @type parent: a Tkinter PortableToplevel or Tk instance.
        
        @param atomNames: list of the MMTK atom name of the atoms of the universe.
        @type atomNames: list of strings
        """

        # The PortableToplevel is window is set up.
        PortableToplevel.__init__(self, parent)
        
        # And transiented to the parent dialog.
        self.transient(parent)
        
        # The title of the selection dialog will depend on the selection type.
        self.title('Atom order dialog')

        # The arguments are copied to instance attribute.
        self.parent = parent
        self.atomNames = atomNames
                    
        # The main frame of the dialog.
        body = Frame(self)
        
        # The body of the dialog is built.
        self.body(body)
        
        # and displayed.
        body.grid(row = 0, column = 0, sticky = EW)  

        # The button box at the bottom of the dialog is built.
        self.buttonbox()        
        
        # The dialog grab all events.
        self.grab_set()

        # Clicking on the up-right corner destroys the dialog and reset the |self.selectionString| to its default value.
        self.protocol("WM_DELETE_WINDOW", self.cancel)

        # This dialog can not be resized.
        self.resizable(width = NO, height = NO)

        # Its position is set from the one of the parent dialog and its dimension are (50,50).
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50, parent.winfo_rooty()+50))

        # The dialog takes the focus.
        self.focus_set()

        # And the event waiting loop is started.
        self.wait_window(self)            

    def body(self, master):
        """Builds up the dialog.
        
        @param master: the main frame of the dialog.
        @type master: a Tkinter Frame instance
        """

        # The frame that will host all the dialog settings widgets.
        self.settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        self.settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        self.settingsFrame.grid_columnconfigure(0, weight = 1)
                                            
        # The first listbox is for the object name selection.
        self.atomNamesListbox = ComboListbox(self.settingsFrame,\
                                             frameLabel = 'Atom names',\
                                             contents = self.atomNames)
        self.atomNamesListbox.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.atomNamesListbox.grid_columnconfigure(0, weight = 1)
        self.atomNamesListbox.lb.config(width = 20, exportselection = 0, selectmode = SINGLE)
        self.atomNamesListbox.lb.bind('<ButtonRelease-1>', self.selectAtomName)

        # The first listbox is for the object name selection.
        self.selectedAtomNamesListbox = ComboListbox(self.settingsFrame,\
                                            frameLabel = 'Selected atom names')
        self.selectedAtomNamesListbox.grid(row = 0, column = 1, padx = 2, pady = 2, sticky = EW)
        self.selectedAtomNamesListbox.grid_columnconfigure(0, weight = 1)
        self.selectedAtomNamesListbox.lb.config(width = 20, exportselection = 0, selectmode = SINGLE)
        self.selectedAtomNamesListbox.lb.bind('<KeyPress-Delete>', self.deleteAtomName)
                
    def buttonbox(self):
        """Adds the 'Cancel' and 'OK' buttons.
        """

        # The frame that will contain the 'Cancel' and 'OK' buttons.
        box = LabelFrame(self, text = 'Actions', bd = 2, relief = GROOVE)
        box.grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3)
        box.grid_columnconfigure(0, weight = 1)

        # The 'Cancel' button will close the dialog resting the selection string to its default value.
        w = Button(box, text = "Cancel", width=10, command = self.cancel)
        w.grid(row = 0, column = 0, sticky = E)
        
        # The 'OK' button will close the dialog keeping the selection string in its current state.
        w = Button(box, text = "OK", width=10, command = self.ok, default=ACTIVE)
        w.grid(row = 0, column = 1, sticky = E)
        
        # Bind 'Return' keyboard key to |self.ok| method. Emulates the 'OK' button.
        self.bind("<Return>", self.ok)

        # Bind 'Esc' keyboard key to |self.cancel| method. Emulates the 'Cancel' button.
        self.bind("<Escape>", self.cancel)

    def ok(self, event = None):
        """Updates the atom selection and closes the dialog.
        
        @param event: pressing the 'Return' keyboard key event if not None.
        @type event: a Tkinter event or None
        """
        
        # The selection string is set with the values stored in the 'Selected atoms' listbox.
        self.selectedAtomNames = ','.join(self.selectedAtomNamesListbox.lb.get(0, END))
                
        # Withdraws the dialog before closing it.
        self.withdraw()
        
        # Updates the display of the widget.
        self.update_idletasks()
        
        # Put the focus back on the parent dialog.
        self.parent.focus_set()
        
        # Close the dialog.
        self.destroy()

    def cancel(self, event = None):
        """Reset the atom selection and closes the dialog.
        """

        # The string that will store the selected atom names is set to an empty string.
        self.selectedAtomNames = ''
                
        # Put the focus back on the parent dialog.
        self.parent.focus_set()

        # Close the dialog.
        self.destroy()
                                
    def getValue(self):
        """Returns the string storing the selected atom names.

        Thanks to this method, this dialog can be used like any other ComboWidget for which the getValue
        allows to fetch their contents.
        """
        
        return self.selectedAtomNames
                                            
    def selectAtomName(self, event):
        """Displays into the 'Selected atom' listbox the selected atoms.
        
        @param event: left-clicking on one item of the 'Object name' listbox.
        @type event: a Tkinter event or None
        """
        
        # If no atom was selected in the 'Atom name' listbox, displays a warning and return.
        if not self.atomNamesListbox.lb.curselection():
            LogMessage('warning','Please select an atom name.',['gui'])
            return
        
        # The selected atom name.
        atomName = self.atomNamesListbox.lb.get(int(self.atomNamesListbox.lb.curselection()[0]))
        
        # If the selected atom is not in the 'Selected atoms' listbox, then inserts it.
        if atomName not in self.selectedAtomNamesListbox.lb.get(0, END):
            self.selectedAtomNamesListbox.lb.insert(END, atomName)

    def deleteAtomName(self, event):
        """Delete a selected atom name from the 'Selected atoms' listbox.

        @param event: pressing the 'Delete' key after have selected an atom name.
        @type event: a Tkinter event or None
        """

        # If no atom was selected in the 'Atom names' listbox, displays a warning and return.
        if not self.selectedAtomNamesListbox.lb.curselection():
            LogMessage('warning','Please select an atom name.',['gui'])
            return

        # Delete the selected atom name from the listbox.
        self.selectedAtomNamesListbox.lb.delete(int(self.selectedAtomNamesListbox.lb.curselection()[0]))
