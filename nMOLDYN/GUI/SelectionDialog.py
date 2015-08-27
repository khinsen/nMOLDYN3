"""This modules implements the selection dialog used to performs different types of selections.

The atom selection can be performed for various purposes such as selection of atoms for an analysis, selection
of hydrogen atoms to deuterate or selection of several group of atoms on which an analysis will be performed 
collectively.

Classes:
    * SelectionDialog: sets up a dialog to perform the selection of atoms or group of atoms prior to an analysis.
"""

# The python distribution modules
import copy
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
            
class SelectionDialog(PortableToplevel):
    """Sets up a dialog to perform the selection of atoms or group of atoms prior to an analysis.
    
    There are three kind of selections that can be performed in nMOLDYN:
        - subset selection: select a subset of atoms on which the analysis will be performed.
        - deuteration selection: select a subset of hydrogens atoms that will considered as deuterium 
        atoms during an analysis. Its used the same mechanism thatn subset selection.
        - group selection: select groups of atoms on which an analysis will be performed collectively.
        
    For each of these types of selection, one can perform the selection from three kind of medias:
        - from a file whose opened unit can be executed by Python. That file must assign:
            * 'pdb' variable to the path of a PDB file of the system.
            * 'subset' variable to a list of the indexes (starting from 1) of the atoms to select in the PDB file 
            in case of a subset selection.
            * 'deuteration' variable to a list of the indexes (starting from 1) of the atoms to select in the PDB file 
            in case of a deuteration selection.
            * 'group' variable to nested lists where each nested list contains the indexes (starting from 1) of the atoms 
            to select from the PDB file that will be included in a given group in case of a group selection.
            
        - from a valid python expression assigning the variable 'selection' to the selection.
        - from a selection string built directly from the dialog.
        
    @note: inherits from Tkinter.PortableToplevel class.
    """

    def __init__(self, parent, selectionType, univContents, widget):
        """The constructor.
        
        @param parent: the parent dialog.
        @type parent: a Tkinter PortableToplevel or Tk instance.
        
        @param selectionType: the kind of selection to perform.
        @type selectionType: string being one of 'subset', 'deuteration' or 'group'

        @param univContents: a dictionnary that hierarchizes the universe contents. Its keys are the nMOLDYN name 
        of the chemical objects of the universe and its values are subdictionnaries used to perform the selection. 
        Their keys/values represents respectively the selection keywords and their associated values:
            - 'name'/a python set of the different MMTK names under which the object can be found plus the '*' wildcard.
            - 'number'/the number of such objects in the universe.
        and some others specific to the class of the MMTK object:
            - object of Atom class:
                * 'objectclass'/'Atom'
                * 'groupinglevel'/['atom'] the level at which the atoms of the object can be grouped.                
            - object of AtomCluster class:
                * 'objectclass'/'AtomCluster'
                * 'groupinglevel'/['atom','cluster'] the level at which the atoms of the object can be grouped.
                * 'atomname'/a python set of the MMTK full names of the atoms of the object plus the '*' wildcard.
                * 'atomelement'/a python set of the chemical elements found in the object plus the '*' wildcard.
                
            - object of Molecule class:
                * 'objectclass'/'Molecule'
                * 'groupinglevel'/['atom','amine', 'hydroxy', 'methyl', 'thiol', 'molecule'] the level at which the atoms of the object can be grouped.
                * 'atomname'/a python set of the MMTK full names of the atoms of the object plus the '*' wildcard.
                * 'atomelement'/a python set of the chemical elements found in the object (e.g. hydrogen, carbon ...) plus the '*' wildcard.
                * 'chemfragment'/['amine', 'hydroxy', 'methyl', 'thiol'] the chemical fragments handled by nMOLDYN.
                
            - object of NucleotideChain class:
                * 'objectclass'/'NucleotideChain'
                * 'groupinglevel'/['atom', 'amine', 'residue', 'nucleicacid'] the level at which the atoms of the object can be grouped.
                * 'atomname'/a python set of the MMTK full names of the atoms of the object plus the '*' wildcard.
                * 'atomelement'/a python set of the chemical elements found in the object plus the '*' wildcard.
                * 'atomtype'/a python set of the MMTK atom type found in the object plus the '*' wildcard.
                * 'nuclname'/a python set of the nucleotide names found in the object plus the '*' wilcard.
                * 'nucltype'/a python set of the nucleotide types found in the object plus the '*' wilcard.
                * 'misc'/['backbone', 'bases'] some miscellaneous selectable parts of the object.
                
            - object of PeptideChain class:
                * 'objectclass'/'PeptideChain'
                * 'groupinglevel'/['atom', 'amine', 'hydroxy', 'methyl', 'thiol', 'residue','chain'] the level at which the atoms of the object can be grouped.
                * 'atomname'/a python set of the MMTK full names of the atoms of the object plus the '*' wildcard.
                * 'atomelement'/a python set of the chemical elements found in the object plus the '*' wildcard.
                * 'atomtype'/a python set of the MMTK atom type found in the object plus the '*' wildcard.
                * 'resname'/a python set of the residue names found in the object plus the '*' wilcard.
                * 'restype'/a python set of the residue types found in the object plus the '*' wilcard.
                * 'resclass'/a python set of the residue classes found in the object plus the '*' wilcard.
                * 'chemfragment'/['amine', 'c_alphas', 'hydroxy', 'methyl', 'thiol'] the chemical fragments handled by nMOLDYN.
                * 'misc'/['backbone', 'sidechains']  some miscellaneous selectable parts of the object.

            - object of Protein class:
                * 'objectclass'/'Protein'
                * 'groupinglevel'/['atom', 'amine', 'hydroxy', 'methyl', 'thiol', 'residue','chain', 'protein'] the level at which the atoms of the object can be grouped.
                * 'atomname'/a python set of the MMTK full names of the atoms of the object plus the '*' wildcard.
                * 'atomelement'/a python set of the chemical elements found in the object plus the '*' wildcard.
                * 'atomtype'/a python set of the MMTK atom type found in the object plus the '*' wildcard.
                * 'resname'/a python set of the residue names found in the object plus the '*' wilcard.
                * 'restype'/a python set of the residue types found in the object plus the '*' wilcard.
                * 'resclass'/a python set of the residue classes found in the object plus the '*' wilcard.
                * 'chemfragment'/['amine', 'c_alphas', 'hydroxy', 'methyl', 'thiol'] the chemical fragments handled by nMOLDYN.
                * 'misc'/['backbone', 'sidechains']  some miscellaneous selectable parts of the object.
                * 'chainname'/a python set of the chain names found in the object plus the '*' wilcard.                            
        @type univContents: dict
        
        @param widget: The ComboButton that triggers the opening of the dialog.
        @type widget: an instance of ComboButton class
        """

        # The PortableToplevel is window is set up.
        PortableToplevel.__init__(self, parent)
        
        # And transiented to the parent dialog.
        self.transient(parent)
        
        # The title of the selection dialog will depend on the selection type.
        self.title('%s selection dialog' % selectionType)

        # The arguments are copied to instance attribute.
        self.parent = parent
        self.selectionType = selectionType        
        self.univContents = copy.deepcopy(univContents)

        # |self.keywords| attribute is a list of all selection keywords arguments that can be found in the universe. 
        self.keywords = ['atomname','atomelement','atomtype','resname','restype','resclass',\
                         'nuclname','nucltype','chainname','chemfragment','misc']

        # The selection dictionnary on which the selection string will be built. It will be used to store the selection under process on a easily 
        # mutatable form. It will contain just one key/value pair being one of those:
        #    - 'filename'/name of the selection file for a file-based selection.
        #    - 'expression'/valid python expression for an expression-based selection
        #    - 'object'/list of subdictionnaries that stores the selection processed 
        # for the different objects. Its contains the following keys/values:
        #        * 'objectname'/the name of the object.
        #        * 'selection'/list of selection keywords/selection values dictionnaries separated by linker-strings.
        #        * 'groupinglevel'/the level at which will be grouped the atoms of the object. Only set in the case
        # of a group selection.
        self.selection = {}
                                
        # Initializes the |self.selectionString| attribute using the value stored in the entry associated with
        # the ComboButton from which the dialog was launched. This will also be the default value.
        self.selectionString = self.defaultSelectionString = widget.getValue()
        
        # Case of a group selection.
        if self.selectionType == 'group':
            # This dictionnary stores the default grouping levels corresponding to each class of chemical object handled
            # by MMTK.
            self.defaultGroupingLevels = {'Atom' : 'atom',\
                                          'AtomCluster' : 'cluster',\
                                          'Molecule' : 'molecule',\
                                          'NucleotideChain' : 'residue',\
                                          'PeptideChain' : 'residue',\
                                          'Protein' : 'residue',\
                                          'AllClass' : 'default'}
            
        # This attribute stores the selected selection keyword.
        self.selectedKeyword = None

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

        # A subframe containing the radio buttons used to set the selection media.
        self.selectionMediasRadio = ComboRadiobutton(self.settingsFrame,\
                                                     frameLabel = "Selection media",\
                                                     contents = ["from a selection file",\
                                                                 "from the loaded trajectory",\
                                                                 "from an expression string"],\
                                                    layout = (1,3),\
                                                    default = 1)
        # Each of the radiobutton is bound to the |self.changeSelectionMedia| method that will update the dialog
        # when the user changes of selection media.
        [r.config(width = 22, command = self.changeSelectionMedia) for r in self.selectionMediasRadio.radio]
        self.selectionMediasRadio.grid(row = 0, column = 0, padx = 2, pady = 5, sticky = EW)

        # A subframe whose contents depends on the selected selection media.
        self.selectionFrame = Frame(self.settingsFrame, bd = 2, relief = GROOVE, width = 1000, height = 300)
        self.selectionFrame.grid(row = 1, column = 0, padx = 3, pady = 3, sticky = 'NEW')
        self.selectionFrame.grid_propagate(0)
        self.selectionFrame.grid_columnconfigure(0, weight = 1)
        self.selectionFrame.grid_rowconfigure(10, weight = 1)

        # For a file-based selection media, this is a file browser.
        self.fileBrowser = ComboFileBrowser(self.selectionFrame,\
                                            frameLabel = 'Selection from a selection file',\
                                            command = self.selectFromFile)
        self.fileBrowser.entry.config(width = 100)
        self.fileBrowser.entry.bind('<Return>', self.selectFromFile)
        
        # For an object-based selection media, this is a set of listboxes from where one can build
        # the object-based selection string. This is that contents that is displayed by default when opening the dialog.
        self.objectBrowser = ComboFrame(self.selectionFrame, frameLabel = 'Selection from the loaded trajectory')
        self.objectBrowser.grid(row = 0, column = 0, padx = 3, pady = 3, sticky = EW)
        self.objectBrowser.grid_columnconfigure(0, weight = 1)
        self.objectBrowser.grid_columnconfigure(1, weight = 1)
        self.objectBrowser.grid_columnconfigure(2, weight = 1)
                                    
        # The first listbox is for the object name selection.
        self.objectName = ComboListbox(self.objectBrowser,\
                                       frameLabel = 'Object name',\
                                       contents = sorted(self.univContents.keys()))
        self.objectName.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.objectName.grid_columnconfigure(0, weight = 1)
        self.objectName.lb.config(width = 20, exportselection = 0, selectmode = SINGLE)        
        self.objectName.lb.bind('<ButtonRelease-1>', self.selectObjectName)

        # The second listbox is for the selection keyword selection.
        self.selectionKeyword = ComboListbox(self.objectBrowser, 'Selection keyword')
        self.selectionKeyword.grid(row = 0, column = 1, padx = 2, pady = 2, sticky = EW)
        self.selectionKeyword.grid_columnconfigure(0, weight = 1)
        self.selectionKeyword.lb.config(width = 20, exportselection = 0, selectmode = SINGLE)
        self.selectionKeyword.lb.bind('<ButtonRelease-1>', self.selectKeyword)
        
        # The third listbox is for the selection value selection.
        self.selectionValue = ComboListbox(self.objectBrowser, 'Selection value')
        self.selectionValue.grid(row = 0, column = 2, padx = 2, pady = 2, sticky = EW)
        self.selectionValue.grid_columnconfigure(0, weight = 1)
        self.selectionValue.lb.config(width = 40, exportselection = 0, selectmode = MULTIPLE)
        self.selectionValue.lb.bind('<ButtonRelease-1>', self.selectValue)
        self.selectionValue.lb.bind('<Control-a>', self.selectValue)

        # The fourth listbox is specific to group selection and is for the selection of the level at which the atoms
        # of the selected object will be grouped.
        if self.selectionType == 'group':
            self.groupingLevel = ComboListbox(self.objectBrowser,'Grouping level')
            self.groupingLevel.grid(row = 0, column = 3, padx = 2, pady = 2, sticky = EW)
            self.groupingLevel.grid_columnconfigure(0, weight = 1)
            self.groupingLevel.lb.config(width = 20, exportselection = 0, selectmode = SINGLE)
            self.groupingLevel.lb.bind('<ButtonRelease-1>', self.selectGroupingLevel)

        # This subframe will contain the button used to link the object-based selection string with logical operators plus
        # an undo and clear buttons.
        selectionButtons = ComboFrame(self.objectBrowser, frameLabel = 'Linkers')
        selectionButtons.grid(row = 1, column = 0, columnspan = 5, padx = 2, pady = 2, sticky = EW)

        # The 'Clear' button will reset the object-based selection string to its default value.
        Button(selectionButtons, text = 'Clear', width = 4, command = lambda : self.undo(True)).grid(row = 0, column = 0)

        # The 'Undo' button will set the object-based selection string to a previous state.
        Button(selectionButtons, text = 'Undo' , width = 4, command = lambda : self.undo(False)).grid(row = 0, column = 1)

        # The '(' button will introduce a ( in the object-based selection string.
        Button(selectionButtons, text = '('    , width = 4, command = lambda : self.appendLinker('(')).grid(row = 0, column = 2)

        # The ')' will introduce a ) in the object-based selection string.
        Button(selectionButtons, text = ')'    , width = 4, command = lambda : self.appendLinker(')')).grid(row = 0, column = 3)
        
        # The 'AND' will introduce a logical AND operator in the object-based selection string.
        Button(selectionButtons, text = 'AND'  , width = 4, command = lambda : self.appendLinker('AND')).grid(row = 0, column = 4)

        # The 'OR' will introduce a logical OR operator in the object-based selection string.
        Button(selectionButtons, text = 'OR'   , width = 4, command = lambda : self.appendLinker('OR')).grid(row = 0, column = 5)

        # For a expression-based selection media, this is a text widget from where the user can enter a valid python
        # expression.
        self.expressionBrowser = ComboText(self.selectionFrame, 'Selection from an expression string')
        self.expressionBrowser.text.config(height = 10, state = NORMAL)
        self.expressionBrowser.text.bind('<Return>', self.selectFromExpression)
        
        # This frame will contain a text widget that will display a preview of the selection string under construction.
        self.selectionPreview = ComboText(self.settingsFrame, 'Selection preview')
        self.selectionPreview.grid(row = 2, column = 0, padx = 3, pady = 3, sticky = EW)
        self.selectionPreview.grid_columnconfigure(0, weight = 1)
        self.selectionPreview.text.config(height = 10, font = ('Courier', '12'))
        
        # Some tags are defined for that text widget. They will allow some specific words to be highlighted.
        self.selectionPreview.text.tag_config('media'        , foreground = 'blue'  , font = ('Courier','12','bold'))
        self.selectionPreview.text.tag_config('italic'       , foreground = 'black' , font = ('Courier','12','normal', 'italic'))
        self.selectionPreview.text.tag_config('bold'         , foreground = 'black' , font = ('Courier','12','bold'))
        self.selectionPreview.text.tag_config('objectname'   , foreground = 'red'   , font = ('Courier','12','bold'))
        self.selectionPreview.text.tag_config('brace'        , foreground = 'black' , font = ('Courier','12','bold'))
        self.selectionPreview.text.tag_config('logical'      , foreground = 'purple', font = ('Courier','12','bold'))

        # The default selection string is displayed.
        self.displaySelectionString()        
        
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
        """Closes the selection editor dialog and do a final update of the selection string.
        
        @param event: pressing the 'Return' keyboard key event if not None.
        @type event: a Tkinter event or None
        """

        # Performs a last checking of the selection string.
        # Case, where something is wrong with the selection string. Displaya a warning and resets everything.
        if not self.checkSelectionString():
            LogMessage('warning', 'This is not a valid selection string.', ['gui'])
            self.undo(True)
            return

        # Withdraws the dialog before closing it.
        self.withdraw()
        
        # Updates the display of the widget.
        self.update_idletasks()
        
        # Put the focus back on the parent dialog.
        self.parent.focus_set()
        
        # The parent takes all the forthcoming events.
        self.parent.grab_set()
        
        # Close the dialog.
        self.destroy()

    def cancel(self, event = None):
        """Closes the dialog setting the selection string back to its default value.
        """

        # The selection string is set to its default value.
        self.selectionString = self.defaultSelectionString
        
        # Put the focus back on the parent dialog.
        self.parent.focus_set()

        # The parent takes all the forthcoming events.
        self.parent.grab_set()

        # Close the dialog.
        self.destroy()

    def checkSelectionString(self):
        """Checks the selection string.
        """

        # Case of a file-based selection. Nothing to chech there.
        if self.selection.has_key('filename'):
            return True
        
        # Case of an expression-based selection. Nothing to chech there.
        elif self.selection.has_key('expression'):
            return True
        
        # Case of an object-based selection. Some basic construction rules will be checked.
        else:

            # Reverse loop over the number of subdictionnaries of the |self.selection['object']| list.
            for comp in range(len(self.selection['object'])-1,-1,-1):

                # The current subdictionnary.
                subDictionnary = self.selection['object'][comp]

                # Check that its 'selection'-key value is not an empty list.
                if subDictionnary['selection']:
                
                    # If not, check that the last element of the list is not an 'AND' or 'OR' logical linker. 
                    # This would have no sense.
                    if subDictionnary['selection'][-1] in ['AND', 'OR']:
                        
                        # If so, delete it from the list.
                        del subDictionnary['selection'][-1]
                    
                        # If after this the list is empty, delete the current selection subdictionnary.
                        if not subDictionnary['selection']:
                            del self.selection['object'][comp]
                        
                # If it is empty, delete the current subdictionnary.
                else:
                    del self.selection['object'][comp]

            # If after all of this, the |self.selection['object']| is empty, raises an error. Something was certainly wrong 
            # with the processed selection.
            if not self.selection['object']:
                return False
                
            # Updates the selection string by rebuilding it from the refreshed |self.selection['object']| list.
            self.buildSelectionString()
            
            return True
                                
    def getValue(self):
        """Returns the constructed selection string.

        Thanks to this method, the selection dialog can be used like any other ComboWidget for which the getValue
        allows to fetch their contents.
        """
        return self.selectionString
        
    
    def changeSelectionMedia(self):
        """Changes the selection media from which the selection will be performed. 
        
        It can be either from a selection file, either from the loaded trajectory or from an expression string. 
        When changing selection media, the previous selection is cleared.
        """
        
        # The selected selection media.
        selectionMedia = self.selectionMediasRadio.getValue()
        
        # Suppress the widgets corresponding to the previously displayed selection media.
        for w in self.selectionFrame.winfo_children():
            w.grid_forget()

        # And displays the ones corresponding to newly selected selection media.        
        # Case of a file-based selection media.
        if selectionMedia == 'from a selection file':
            self.fileBrowser.grid(row = 0, column = 0, padx = 3, pady = 3, sticky = EW)
            self.fileBrowser.grid_columnconfigure(0, weight = 1)
            
        # Case of an object-based selection media.
        elif selectionMedia == 'from the loaded trajectory':
            self.objectBrowser.grid(row = 0, column = 0, padx = 3, pady = 3, sticky = EW)
            self.objectBrowser.grid_columnconfigure(0, weight = 1)
            self.objectBrowser.grid_columnconfigure(1, weight = 1)
            self.objectBrowser.grid_columnconfigure(2, weight = 1)

        # Case of an expression-based selection media.
        elif selectionMedia == 'from an expression string':
            self.expressionBrowser.grid(row = 0, column = 0, padx = 3, pady = 3, sticky = EW)
            self.expressionBrowser.grid_columnconfigure(0, weight = 1)
            
    def undo(self, clearAll):
        """This methods sets the selection to its previous state.
        
        @param clearAll: if True, reset the selection to its default state.
        @type clearAll: bool
        """
        
        # Case where |clearAll| argument is True. Deletes the |object| key from the self.selection dictionnary.
        if clearAll:                
            del self.selection['object']
            
        # Otherwise, removes the last subdictionnary of the |self.selection['object']| list.
        else:
            
            del self.selection['object'][-1]

            # If after this the list is empty, deletes the 'object' key from the |self.selection| dictionnary.
            if not self.selection['object']:
                del self.selection['object']
            
        # The 'Object name' listbox is cleared.
        self.objectName.lb.selection_clear(0, END)
        
        # The 'Selection keyword' listbox is cleared.
        self.selectionKeyword.cleanup()

        # And the current |self.selectedKeyword| reset to None.
        self.selectedKeyword = None

        # The 'Selection value' listbox is cleared.
        self.selectionValue.cleanup()

        # The 'Grouping level' listbox is cleared.
        if self.selectionType == 'group':
            self.groupingLevel.lb.selection_clear(0, END)
        
        # Builds the selection string out of the updated |self.selection| dictionnary.
        self.buildSelectionString()
       
        # Displays the selection string on the 'Selection preview' text widget.
        self.displaySelectionString()

    def buildSelectionString(self):
        """Builds the selection string.
        """

        # Case of a file-based selection, the selection string just contains the keyword 'filename' 
        # followed by the name of the selected file.
        if self.selection.has_key('filename'):
            self.selectionString = 'filename %s' % self.selection['filename']

        # Case of an expression-based selection, the selection string just contains the keyword 'expression' 
        # followed by the entered python expression.
        elif self.selection.has_key('expression'):
            self.selectionString = 'expression %s' % self.selection['expression']
        
        # Case of an object-based selection, the selection string is built from the subdictionnaries 
        # stored in |self.selection['object']| list.
        elif self.selection.has_key('object'):
        
            # The selection string will start with 'object '.
            self.selectionString = 'object '
            
            # The number of subdictionnaries.
            nSubdictionnaries = len(self.selection['object'])
            
            # Loop over the number of subdictionnaries.
            for comp in range(nSubdictionnaries):
                
                # The current subdictionnary.
                subdictionnary = self.selection['object'][comp]
                        
                # The object name corresponding to that subdictionnary.
                objectName = subdictionnary['objectname']
                
                # The keywords/values list corresponding to that subdictionary.
                selection = subdictionnary['selection']
                
                # Check that this list is not empty.
                if selection:
                                        
                    # Append the name of the object to the selection string.
                    self.selectionString += 'objectname %s ' % objectName
                                                    
                    # Loop over the selection list element.
                    for sel in selection:
                        # Case where it is a selection keywords/selection values dictionnary
                        if isinstance(sel,dict):
                            # The keys (i.e. selection keywords) are sorted.
                            sortedKeywords = sorted(sel.keys())
                            
                            # Loop over the number of selection keywords of the dictionnary.
                            for index in range(len(sortedKeywords)):
                                
                                # The current keyword.
                                keyword = sortedKeywords[index]
                                
                                # The current value. A list of selected values for the corresponding keyword.
                                value = sel[keyword]

                                # The selection string is updated.
                                self.selectionString += '%s %s ' % (keyword, ','.join(value))

                                # A logical 'OR' is always appended between two consecutive keywords.
                                if index < len(sortedKeywords) - 1:
                                    self.selectionString += 'OR '

                        # Case where it is a linker-string.
                        else:
                            # The linker is appended to the selection string.
                            self.selectionString += sel + ' '
                            
                # Case of a group selection, append the selected grouping level for the processed subdictionnary.
                if self.selectionType == 'group':
                    self.selectionString += 'groupinglevel %s ' % self.selection['object'][comp]['groupinglevel']
                    
                # A logical 'OR' is always appended between two consecutive subdictionnaries.
                if comp < nSubdictionnaries - 1:
                    self.selectionString += 'OR '
                                
            # The string is stripped.
            self.selectionString = self.selectionString.strip()
            
            # If the subdictionnaries did not contain any actual selection, reset the selection string to 
            # its default value.
            if self.selectionString == 'object':
                self.selectionString = self.defaultSelectionString
                
        # Otherwise, sets the selection string to its default value.
        else:
            self.selectionString = self.defaultSelectionString
                                
    def displaySelectionString(self):
        """Displays the selection string under process in the selection preview text widget.
        """
                
        # If the selection string is empty, do nothing.
        if not self.selectionString:
            return

        # The 'Selection preview' text widget is cleared.
        self.selectionPreview.cleanup()        
        
        # The selection string is inserted in the 'Selection preview' text widget.
        self.selectionPreview.insert(contents = self.selectionString)
        
        # The word 'objectname' is highlighed.
        self.selectionPreview.tag('objectname','objectname')

        # The word 'groupinglevel' is highlighted.
        self.selectionPreview.tag('groupinglevel','bold')

        # The selection media is highlighted.
        for k in self.selection.keys():
            self.selectionPreview.tag(k,'media')

        # The selection keywords are highlighted.
        for k in self.keywords:
            self.selectionPreview.tag(k,'italic')

        # The linker are highlighted.
        self.selectionPreview.tag('(','brace')
        self.selectionPreview.tag(')','brace')
        self.selectionPreview.tag('AND','logical')
        self.selectionPreview.tag('OR','logical')
        
    def selectFromExpression(self, event):
        """Checks and registers the entered expression-based selection string. 
        
        The expression must be a set of valid ';'-separated python instructions the last one being 'selection = ...'        
        as the selection string parser will search for the selection variables when executing the expression string.

        To refer to the universe just use the variable 'self.universe'.
        
        @param event: pressing the 'Return' keyboard key on the expression-based selection string editor.
        @type event: a Tkinter event
        """

        # The expression-based selection string.        
        expression = self.expressionBrowser.text.get('0.0', END).strip()
        
        # If it empty display a warning and does nothing.
        if not expression:
            LogMessage('warning','Please enter a valid expression.',['gui'])
            return

        # If it does not contain a 'selection =', displays a warning and resets the expression to the empty string.
        if not re.findall('selection\s*=', expression):
            LogMessage('warning','An expression-based selection string must contain an assignement to "selection" variable.',['gui'])
            self.expressionBrowser.setValue('')
            return
        
        # The expression seems OK, updates the |self.selection| dictionnary.
        self.selection = {'expression' : expression}
        
        # Builds the selection string out of the |self.selection| dictionnary.
        self.buildSelectionString()

        # Displays the selection string on the 'Selection preview' text widget.
        self.displaySelectionString()
        
        # This prevents the next binding 'Return' key pressing event to occur namely the call to |self.ok| method.
        return 'break'
            
    def selectFromFile(self, event = None):
        """Performs a selection from a file.
        
        @param event: pressing the 'Return' keyboard key event if not None.
        @type event: a Tkinter event or None
        """            

        # Case where the user entered directly the file name in the entry associated to the file browser and pressed 
        # 'Return' key.
        if event:
            # Gets the file name entered in the entry.
            filename = self.fileBrowser.getValue()
        
        # Case where the user selects the file from the browser.
        else:
            # The name of the NetCDF file to load.
            filename = askopenfilename(parent = self, filetypes = [('NMS file','.nms'), ('Python file', '.py'), ('All files', '.*')])
                
        # Checks that the file name is not an empty string.
        if not filename:
            # If so, displays a warning, resets the entry to '' and return.
            LogMessage('warning','Please enter a selection file name.',['gui'])
            self.fileBrowser.setValue('')
            return
        
        # Checks that the file exists.
        if not os.path.exists(filename):
            # If not, displays a warning, resets the entry and return.
            LogMessage('warning','Please enter a valid selection file name.',['gui'])
            self.fileBrowser.setValue('')
            return
                    
        # Update the entry associated to the file browser.
        self.fileBrowser.setValue(filename)
        
        # Update the |self.selection| dictionnary.
        self.selection = {'filename' : filename}
        
        # Builds the selection string out of the |self.selection| dictionnary.
        self.buildSelectionString()

        # Displays the selection string on the 'Selection preview' text widget.
        self.displaySelectionString()

        # This prevents the next binding 'Return' key pressing event to occur namely the call to |self.ok| method.
        return 'break'
        
    def selectObjectName(self, event):
        """Displays into the 'Selection keywords' listbox all the selection keywords corresponding to the selected 
        object name.
        
        @param event: left-clicking on one item of the 'Object name' listbox.
        @type event: a Tkinter event or None
        """
        
        # If no object was selected in the 'Object name' listbox, displays a warning and return.
        if not self.objectName.lb.curselection():
            LogMessage('warning','Please select an object name.',['gui'])
            return
        
        # The selected object name.
        objectName = self.objectName.lb.get(int(self.objectName.lb.curselection()[0]))
        
        # Case where the |self.selection| dictionnary does not have the 'object' key.
        if not self.selection.has_key('object'):
            
            # Sets it with to a list with a first subdictionnary to fill further.
            self.selection['object'] = [{'objectname' : objectName , 'selection' : []}]
            
        # Otherwise,
        else:
            # If the last subdictionnary of the |self.selection['object']| list contains no selection, replace it
            # by a brand new one.
            if self.selection['object'][-1]['selection']:
                self.selection['object'].append({'objectname' : objectName , 'selection' : []})
                
            # Otherwise, append to the list a new subdirectionnary to fill further.
            else:
                self.selection['object'][-1] = {'objectname' : objectName , 'selection' : []}

        # Case of a group-based selection.
        if self.selectionType == 'group':
            
            # The 'Grouping level' listbox is cleared.
            self.groupingLevel.cleanup()
            
            # And updated with the grouping levels available for the selected object name.
            for gl in self.univContents[objectName]['groupinglevel']:
                self.groupingLevel.lb.insert(END, gl)
            
            # The MMTK class of the selected object.
            objectClass = self.univContents[objectName]['objectclass']
            
            # The default grouping level for the corresponding MMTK class.
            defaultGroupingLevel = self.defaultGroupingLevels[objectClass]
            
            # The default grouping level is activated in the 'Grouping levels' listbox.
            index = self.univContents[objectName]['groupinglevel'].index(defaultGroupingLevel)            
            self.groupingLevel.lb.selection_set(index)
            
            # And sets as athe value of 'groupinglevel' key of the new subdictionnary.
            self.selection['object'][-1]['groupinglevel'] = self.groupingLevel.lb.get(index)

        # The 'Selection keyword' listbox is cleared.
        self.selectionKeyword.cleanup()
        
        # And updated with the selection keywords available for the selected object name.
        for k in sorted(self.univContents[objectName].keys()):
            
            # This keywords are not used for processing an object.
            if k in ['number','objectclass', 'groupinglevel']:
                continue
                
            self.selectionKeyword.lb.insert(END, k)

        # The 'Selection values' listbox is cleared.
        self.selectionValue.cleanup()
                        
    def selectGroupingLevel(self, event):
        """Selects a grouping level for the atoms of the object under process in case of a group selection.
        
        @param event: left-clicking on one item of the 'Grouping levels' listbox.
        @type event: a Tkinter event or None
        """
        
        # If no grouping level is selected, displays a warning and return.
        if not self.groupingLevel.lb.curselection():
            LogMessage('warning','Please select a grouping level.',['gui'])            
            return
            
        # The selected grouping level.
        groupingLevel = self.groupingLevel.lb.get(int(self.groupingLevel.lb.curselection()[0]))
        
        # Its value is set to the 'groupinglevel' key of the current subdictionnary.
        self.selection['object'][-1]['groupinglevel'] = groupingLevel
                
        # Builds the selection string out of the |self.selection| dictionnary.
        self.buildSelectionString()
       
        # Displays the selection string on the 'Selection preview' text widget.
        self.displaySelectionString()
            
    def selectKeyword(self, event):
        """Display in the 'Selection value' listbox the selection values available for the 
        selected selection keyword.

        @param event: left-clicking on one item of the 'Selection keyword' listbox.
        @type event: a Tkinter event or None
        """
        
        # If no keyword was selected in the 'Selection keyword' listbox, displays a warning and return.
        if not self.selectionKeyword.lb.curselection():
            LogMessage('warning','Please select a selection keyword.',['gui'])
            return

        # The selected selection keyword.
        self.selectedKeyword = self.selectionKeyword.lb.get(int(self.selectionKeyword.lb.curselection()[0]))
       
        # The 'Selection value' listbox is cleared.
        self.selectionValue.cleanup()
        
        # The current object name is extracted from the current subdictionnary.
        objectName = self.selection['object'][-1]['objectname']
        
        # List of the selection values available for the the selected selection keyword.
        values = self.univContents[objectName][self.selectedKeyword]
        
        # These selection values are displayed in the 'Selection value' listbox.
        self.selectionValue.insert(values)

        if self.selection['object'][-1]['selection']:
            
            # Loop over the number of available selection values.
            for index in range(len(values)):
            
                # The current selection value.
                value = values[index]

                # Check whether there is an ongoing keywords/selection values selection under process.
                if isinstance(self.selection['object'][-1]['selection'][-1],dict):
                    
                    # If so, checks whether the seleced keyword was previoulsy selected during that selection.
                    if self.selection['object'][-1]['selection'][-1].has_key(self.selectedKeyword):
                        
                        # If so, selects by default the value corresponding to that keyword in the 
                        # 'Selection values' listbox.
                        if value in self.selection['object'][-1]['selection'][-1][self.selectedKeyword]:
                            self.selectionValue.lb.selection_set(index)
                    
    def selectValue(self, event):
        """Updates the selection string.

        @param event: left-clicking on one item of the 'Selection values' listbox or pressing 'ctrl-a' keys.
        @type event: a Tkinter event or None
        """

        # If the event is Ctrl-a, selects everything.
        if event.keysym == 'a':        
            self.selectionValue.lb.selection_set(0, END)

        # If there is no ongoing selection keyword/selection values selection under process append 
        # an empty dictionnary that will store a new one to |self.selection['object'][-1]['selection']| list.
        if not self.selection['object'][-1]['selection']:
            self.selection['object'][-1]['selection'].append({})
            
        # Otherwise, if there is an ongoing selection,
        else:
            # BUT the last element of the |self.selection['object'][-1]['selection']| list is a linker then 
            # append an empty dictionnary that will store a new one to the list.
            if not isinstance(self.selection['object'][-1]['selection'][-1],dict):
                self.selection['object'][-1]['selection'].append({})            
                    
        # Case where some values have been selected.
        if self.selectionValue.lb.curselection():
            
            # List of the selected values.
            selectedValues = [self.selectionValue.lb.get(int(v)) for v in self.selectionValue.lb.curselection()]

            # The dictionnary corresponding to the ongoing selection keyword/selection values selection is updated
            # for the selected keyword with the corresponding selected values.
            self.selection['object'][-1]['selection'][-1][self.selectedKeyword] = selectedValues
           
        # Case where no value is selected. This can happen the first time a selection keyword is selected
        # or when a selection keyword has been selected but the selected values have been unselected.
        else:
            # The lastter case.
            if self.selection['object'][-1]['selection'][-1].has_key(self.selectedKeyword):
                # Delete the key corresponding to the selected keyword from the ongoing 
                # selection keywords/selection values dictionnary.
                del self.selection['object'][-1]['selection'][-1][self.selectedKeyword]
        
        # If the ongoing selection keywords/selection values is still an empty dictionnary delete it.
        if not self.selection['object'][-1]['selection'][-1]:
            del self.selection['object'][-1]['selection'][-1]
                    
        # Build the selection string out of the self.selection dictionnary.
        self.buildSelectionString()
       
        # Displays the selection string on the 'Selection preview' text widget.
        self.displaySelectionString()
        
    def appendLinker(self, linker):
        """Appends a linker to the selection string after having checked that it could be appended.
        
        @param linker: the linker to append.
        @type linker: string
        """
                       
        # Some basic checkings to prevent user from badly formatting a selection pattern string.
        # May not be sufficient but should throw most of the errors.

        # The linker must be one of '(', ')', 'AND', 'OR'.
        if linker not in ['(',')','AND','OR']:            
            return

        # The selection list of the last subdictionnary.
        selList = self.selection['object'][-1]['selection']
                    
        # If the selection pattern is empty, then the only possible linker to append is '('.
        if not selList:
            if linker != '(':
                return
        else:
            # Case of '('. Must follow '(', 'AND' or 'OR'.
            if linker == '(':
                if selList[-1] not in ['AND','OR','(']:
                    return
                
            # Case of ')'.
            #   -Should always be less ')' than '('.
            #   -Must not follow '(', 'AND' or 'OR'.
            elif linker == ')':
                if selList.count('(') <= selList.count(')'):
                    return
                if selList[-1] in ['AND', 'OR', '(']:
                    return
                
            # Case of 'AND'. Must not follow '(', 'AND' or 'OR'.
            elif linker == 'AND':
                if selList[-1] in ['AND','OR','(']:
                    return
                
            # Case of 'OR'. Must not follow '(', 'AND' or 'OR'.
            elif linker == 'OR':
                if selList[-1] in ['AND','OR','(']:
                    return

        # Everything seems OK, the linker can be appended.
        selList.append(linker)
        
        # The 'Selection value' listbox is cleared.
        self.selectionValue.lb.selection_clear(0, END)

        # Builds the selection string out of the |self.selection| dictionnary.
        self.buildSelectionString()
       
        # Displays the selection string on the 'Selection preview' text widget.
        self.displaySelectionString()
