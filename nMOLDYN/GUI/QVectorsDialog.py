"""The Q vectors dialog.

Classes:
    * QVectorsDialog: sets up a dialog to set the Q vectors generation.
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
from nMOLDYN.GUI.Widgets import ComboStringEntry, ComboIntegerEntry, ComboFloatEntry, ComboFileBrowser, ComboFrame, ComboListbox, ComboRadiobutton, ComboText, PortableToplevel
            
class QVectorsDialog(PortableToplevel):
    """Sets up a dialog to set the Q vectors generation..
            
    @note: inherits from Tkinter.PortableToplevel class.
    """

    def __init__(self, parent, widget):
        """The constructor.
        
        @param parent: the parent dialog.
        @type parent: a Tkinter PortableToplevel or Tk instance.
                
        @param widget: The ComboButton that triggers the opening of the dialog.
        @type widget: an instance of ComboButton class
        """
                                
        self.defaultSettings = {'qgeometry' : 'spatial', 'qshellwidth' : 1.0, 'qshellvalues' : '0.0:10.0:1.0', \
                                'qvectorspershell' : 50, 'hkls' : None,}
        
        self.widgets = {}
        self.widgets["qgeometry"] = StringVar()

        # The PortableToplevel is window is set up.
        PortableToplevel.__init__(self, parent)
        
        # And transiented to the parent dialog.
        self.transient(parent)
        
        # The title of the selection dialog will depend on the selection type.
        self.title('Q vectors dialog')

        # The arguments are copied to instance attribute.
        self.parent = parent
        
        self.widget = widget
                                                                    
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
                                        
        self.updateDialog(self.widget.getValue())
                
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
        
        self.qShellsFrame = LabelFrame(self.settingsFrame, text = 'Q shells', bd = 2, relief = GROOVE)
        self.qShellsFrame.grid(row = 0, column = 0, padx = 20, pady = 5, sticky = EW)
        
        self.widgets["qshellvalues"] = ComboStringEntry(self.qShellsFrame, frameLabel = "Q values (in nm-1)", tagName = 'qshell_values')
        self.widgets["qshellvalues"].grid(row = 0, column = 0, padx = 2, pady = 5, sticky = EW)

        self.widgets["qshellwidth"] = ComboFloatEntry(self.qShellsFrame, frameLabel = "Q shell width (in nm-1)", tagName = 'qshell_width')
        self.widgets["qshellwidth"].grid(row = 1, column = 0, padx = 2, pady = 5, sticky = EW)

        self.widgets["qvectorspershell"] = ComboIntegerEntry(self.qShellsFrame, frameLabel = "Q vectors per shell", tagName= 'qvectors_per_shell')
        self.widgets["qvectorspershell"].grid(row = 2, column = 0, padx = 2, pady = 5, sticky = EW)

        f = LabelFrame(self.settingsFrame, text = 'Q vectors type', bd = 2, relief = GROOVE)
        f.grid(row = 1, column = 0, padx = 20, pady = 5, sticky = EW)
        
        r = Radiobutton(f, text = "Spatial", variable = self.widgets["qgeometry"], value = "spatial", command = self.changeQGeometry)
        r.grid(row = 0, column = 0, padx = 2, pady = 5, sticky = W)
        
        r = Radiobutton(f, text = "Planar", variable = self.widgets["qgeometry"], value = "planar", command = self.changeQGeometry)
        r.grid(row = 1, column = 0, padx = 2, pady = 5, sticky = NW)
        
        self.planarFrame = Frame(f)
        self.planarFrame.grid(row = 1, column = 1, padx = 2, sticky = E)

        self.widgets["plan1"] = ComboStringEntry(self.planarFrame, frameLabel = "Direction 1 (h,k,l)", tagName = 'hkl1')
        self.widgets["plan1"].grid(row = 0, column = 0, padx = 2, pady = 5, sticky = EW)
        
        self.widgets["plan2"] = ComboStringEntry(self.planarFrame, frameLabel = "Direction 2 (h,k,l)", tagName = 'hkl2')
        self.widgets["plan2"].grid(row = 1, column = 0, padx = 2, pady = 5, sticky = EW)
        
        r = Radiobutton(f, text = "Axial", variable = self.widgets["qgeometry"], value = "axial", command = self.changeQGeometry)        
        r.grid(row = 2, column = 0, padx = 2, pady = 5, sticky = W)

        self.axialFrame = Frame(f)
        self.axialFrame.grid(row = 2, column = 1, padx = 2, pady = 5, sticky = E)

        self.widgets["axis"] = ComboStringEntry(self.axialFrame, frameLabel = "Direction 1 (h,k,l)", tagName = 'hkl1')
        self.widgets["axis"].grid(row = 0, column = 0, padx = 2, pady = 5, sticky = EW)
               
        r = Radiobutton(f, text = "Userdefined", variable = self.widgets["qgeometry"], value = "userdefined", command = self.changeQGeometry)        
        r.grid(row = 3, column = 0, padx = 2, pady = 5, sticky = W)

        self.userDefinedFrame = Frame(f)
        self.userDefinedFrame.grid(row = 3, column = 1, padx = 2, pady = 5, sticky = E)
                               
        self.widgets["hkls"] = ComboStringEntry(self.userDefinedFrame, frameLabel = "hkl(s)", tagName = 'hkls')
        self.widgets["hkls"].grid(row = 0, column = 0, padx = 2, pady = 5, sticky = EW)                
        
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
        
        self.buildQGenerator()

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
        """Closes the dialog.
        """
                
        self.qVectorGenerator = copy.copy(self.defaultSettings)
        
        # Put the focus back on the parent dialog.
        self.parent.focus_set()

        # The parent takes all the forthcoming events.
        self.parent.grab_set()

        # Close the dialog.
        self.destroy()
        
    def buildQGenerator(self): 

        self.currentSettings = {}
        
        self.currentSettings["qgeometry"] = self.widgets["qgeometry"].get()
        
        if self.currentSettings["qgeometry"] in ["spatial", "planar", "axial"]:
            
            self.currentSettings["qshellvalues"] = self.widgets["qshellvalues"].getValue()
            self.currentSettings["qshellwidth"] = self.widgets["qshellwidth"].getValue()
            self.currentSettings["qvectorspershell"] = self.widgets["qvectorspershell"].getValue()

            if self.currentSettings["qgeometry"] == "planar":
                
                self.currentSettings["hkls"] = "%s;%s" % (self.widgets["plan1"].getValue(),self.widgets["plan2"].getValue())
                
            elif self.currentSettings["qgeometry"] == "axial":
                
                self.currentSettings["hkls"] = self.widgets["axis"].getValue()
            
        else:
            
            self.currentSettings["hkls"] = self.widgets["hkls"].getValue()
                                                
    def getValue(self):
        """Returns the q vector generator.

        Thanks to this method, the dialog can be used like any other ComboWidget for which the getValue
        allows to fetch their contents.
        """
                
        return self.currentSettings 
    
    def changeQGeometry(self):
        """Updates the dialog when the user changes Q vector type (spatial, planar or axial).
        """
                
        if self.widgets["qgeometry"].get() == "spatial":
            self.setState(self.qShellsFrame, state = "normal")
            self.setState(self.planarFrame, state = "disabled")
            self.setState(self.axialFrame, state = "disabled")
            self.setState(self.userDefinedFrame, state = "disabled")

        elif self.widgets["qgeometry"].get() == "planar":
            self.setState(self.qShellsFrame, state = "normal")
            self.setState(self.planarFrame, state = "normal")
            self.setState(self.axialFrame, state = "disabled")
            self.setState(self.userDefinedFrame, state = "disabled")

        elif self.widgets["qgeometry"].get() == "axial":            
            self.setState(self.qShellsFrame, state = "normal")
            self.setState(self.planarFrame, state = "disabled")
            self.setState(self.axialFrame, state = "normal")
            self.setState(self.userDefinedFrame, state = "disabled")

        elif self.widgets["qgeometry"].get() == "userdefined":            
            self.setState(self.qShellsFrame, state = "disabled")
            self.setState(self.planarFrame, state = "disabled")
            self.setState(self.axialFrame, state = "disabled")
            self.setState(self.userDefinedFrame, state = "normal")
        
    def setState(self, widget, state = "disabled"):
        """Recursively set the state of all the children of a given widget.
        """
                                
        try:
            widget.config(state = state)
            
        except:
            pass
        
        for child in widget.winfo_children():            
            self.setState(child, state = state)
                        
    def updateDialog(self, settings):
        """
        """
                                
        self.currentSettings = eval(settings)

        if not isinstance(self.currentSettings, dict):
            self.currentSettings = copy.copy(self.defaultSettings)
        
        try:
            
            self.widgets["qgeometry"].set(self.currentSettings["qgeometry"])
        
            if self.currentSettings["qgeometry"] in ["spatial", "planar", "axial"]:

                self.widgets["qshellvalues"].setValue(self.currentSettings["qshellvalues"])
                self.widgets["qshellwidth"].setValue(self.currentSettings["qshellwidth"])
                self.widgets["qvectorspershell"].setValue(self.currentSettings["qvectorspershell"])

                if self.currentSettings["qgeometry"] == "planar":
                    plan1, plan2 = [v.strip() for v in self.currentSettings["hkls"].split(";")]
                    self.widgets["plan1"].setValue(plan1)
                    self.widgets["plan2"].setValue(plan2)
                
                elif self.currentSettings["qgeometry"] == "axial":
                    self.widgets["axis"].setValue(self.currentSettings["hkls"])
                        
            else:                                        
            
                self.widgets["hkls"].setValue(self.currentSettings["hkls"])
            
        except:
            self.updateDialog(self.defaultSettings)                                        
                                 
        self.changeQGeometry()
                                
