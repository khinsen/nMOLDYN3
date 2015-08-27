# -*- coding: iso-8859-1 -*-
"""This modules implements I{Help-->About nMOLDYN} dialog.

Classes:
    * CheckJobsStatusDialog: creates I{Help-->About nMOLDYN} dialog used to 
    display some general informations about nMOLDYN.
"""

import os
import webbrowser

# The Tcl/Tk modules
from Tkinter import *

from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.__pkginfo__ import __version__, __versiondate__
from nMOLDYN.GUI.Widgets import PortableToplevel

NMOLDYN_VERSION = __version__
NMOLDYN_VERSION_DATE = __versiondate__

# Some global variables used for the 'about nMOLDYN' generation.
aboutAuthors = """
nMOLDYN\n
A program for analyzing Molecular Dynamics Trajectories\n\n
Version: %s (%s)\n\n
Contacts:\n
\t- Eric Pellegrini (email: pellegrini@ill.fr)\n
\t- Konrad Hinsen (email: hinsen@cnrs-orleans.fr)\n
\t- Gerald Kneller (email: kneller@cnrs-orleans.fr)
""" % (NMOLDYN_VERSION, NMOLDYN_VERSION_DATE)

class GeneralInformationsDialog(PortableToplevel):
    """Sets up a dialog used to visualize some general informations about NMOLDYN.
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
    
        # The logo 'nMOLDYN'.
        self.nMOLDYNLogo = PhotoImage(file = os.path.join(GVAR['nmoldyn_path'], 'Images', 'nMOLDYN.gif'))                
        
        leftFrame = Frame(master)
        leftFrame.grid(row = 0, column = 0, padx = 10, sticky = NS)
        Label(leftFrame, image = self.nMOLDYNLogo).grid(row = 0, column = 0, pady = 10, sticky = NS)

        rightFrame = Frame(master, relief = GROOVE, bd = 2)
        rightFrame.grid(row = 0, column = 1)
        Label(rightFrame, text = aboutAuthors, justify = 'left').grid(row = 0, column = 0)
        Button(rightFrame, text = 'Visit nMOLDYN Homepage', width = 25, command = lambda : self.homepage('nMOLDYN')).grid(row = 1, column = 0)
        Button(rightFrame, text = 'Visit MMTK Homepage', width = 25, command = lambda : self.homepage('MMTK')).grid(row = 2, column = 0)        

        Label(rightFrame, text = "").grid(row = 3, column = 0)
        
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
        return True
    
    def apply(self):
        pass
    
    def homepage(self, website):
        
        if website == 'nMOLDYN':
            webbrowser.open('http://dirac.cnrs-orleans.fr/plone/software/nmoldyn/')
            
        elif website == 'MMTK':
            webbrowser.open('http://dirac.cnrs-orleans.fr/MMTK')
            
        else:
            return
    
