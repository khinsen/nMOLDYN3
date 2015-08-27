"""This modules implements I{File-->Preferences} dialog.

Classes:
    * PreferencesDialog: creates I{File-->Preferences} dialog used to to set up nMOLDYN PREFERENCES variables.
"""


# The python distribution modules
import os
import platform
import sys

# The Tcl/Tk modules
from tkFileDialog import askopenfilename, asksaveasfilename, askdirectory
from Tkinter import *

# The nMOLDYN modules
from nMOLDYN.Preferences import PREFERENCES, loadPreferencesFile, savePreferencesFile, setPreferences
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import *

# The platform name.
PLATFORM = platform.system().upper()

class PreferencesDialog(PortableToplevel):
    """Sets up a dialog used to change some PREFERENCES.
    """
    
    first = False

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
        settingsFrame.grid(row = 0, column = 0, sticky = 'NEW', padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)
        settingsFrame.grid_rowconfigure(0, weight = 1)

        self.preferences = {}        
        
        pagesFrame = Frame(settingsFrame)
        pagesFrame.grid(row = 0, column = 0, sticky = 'NSEW', padx = 3, pady = 3)
            
        pagesFrame.grid_rowconfigure(0, weight = 1)
        pagesFrame.grid_columnconfigure(0, weight = 1)
        pagesFrame.grid_columnconfigure(1, weight = 1)

        self.pages = ComboRadiobutton(pagesFrame,\
                                      frameLabel = "Preferences sections",\
                                      contents = ["File handling", "External programs", "Miscellaneous"],\
                                      layout = (3,1))
        [r.config(width = 16, command = self.changePage) for r in self.pages.radio]
        self.pages.grid(row = 0, column = 0, padx = 2, pady = 5, sticky = 'NSEW')
        self.pages.grid_rowconfigure(10, weight = 1)
        self.pages.grid_columnconfigure(0, weight = 1)
            
        self.preferencesFrame = Frame(settingsFrame, width = 500, height = 350)
        self.preferencesFrame.grid(row = 0, column = 1, padx = 3, pady = 3, sticky = 'NEW')
        self.preferencesFrame.grid_propagate(0)
        self.preferencesFrame.grid_columnconfigure(0, weight = 1)
        self.preferencesFrame.grid_rowconfigure(10, weight = 1)

        self.preferences['trajfile_path'] = ComboFileBrowser(self.preferencesFrame,\
                                                             frameLabel = "Trajectory file path",\
                                                             tagName = 'preferences_trajfile_path',\
                                                             contents = PREFERENCES['trajfile_path'],\
                                                             save = True,\
                                                             command = lambda k = 'trajfile_path': self.openDirectory(k))
        self.preferences['trajfile_path'].grid(row = 1, column = 0, padx = 2, pady = 5, sticky = 'NEW')
        self.preferences['trajfile_path'].grid_columnconfigure(0, weight = 1)

        self.preferences['outputfile_path'] = ComboFileBrowser(self.preferencesFrame,\
                                                               frameLabel = "Output file path",\
                                                               tagName = 'preferences_outputfile_path',\
                                                               contents = PREFERENCES['outputfile_path'],\
                                                               save = True,\
                                                               command = lambda k = 'outputfile_path': self.openDirectory(k))
        self.preferences['outputfile_path'].grid(row = 2, column = 0, padx = 2, pady = 5, sticky = 'NEW')
        self.preferences['outputfile_path'].grid_columnconfigure(0, weight = 1)

        if PLATFORM != 'DARWIN':
            self.preferences['acroread_path'] = ComboFileBrowser(self.preferencesFrame,\
                                                                 frameLabel = "Acrobat reader path",\
                                                                 tagName = 'preferences_acroread_path',\
                                                                 contents = PREFERENCES['acroread_path'],\
                                                                 save = True,\
                                                                 command = lambda k = 'acroread_path': self.openFilename(k))
            
        self.preferences['vmd_path'] = ComboFileBrowser(self.preferencesFrame,\
                                                        frameLabel = "VMD molecular viewer path",\
                                                        tagName = 'preferences_vmd_path',\
                                                        contents = PREFERENCES['vmd_path'],\
                                                        save = True,\
                                                        command = lambda k = 'vmd_path': self.openFilename(k))

        self.preferences['ncdump_path'] = ComboFileBrowser(self.preferencesFrame,\
                                                           frameLabel = "ncdump path",\
                                                           tagName = 'preferences_ncdump_path',\
                                                           contents = PREFERENCES['ncdump_path'],\
                                                           save = True,\
                                                           command = lambda k = 'ncdump_path': self.openFilename(k))

        self.preferences['ncgen_path'] = ComboFileBrowser(self.preferencesFrame,\
                                                          frameLabel = "ncgen path",\
                                                          tagName = 'preferences_ncgen_path',\
                                                          contents = PREFERENCES['ncgen_path'],\
                                                          save = True,\
                                                          command = lambda k = 'ncgen_path': self.openFilename(k))
        
        self.preferences['task_manager_path'] = ComboFileBrowser(self.preferencesFrame,\
                                                          frameLabel = "task_manager path",\
                                                          tagName = 'preferences_task_manager_path',\
                                                          contents = PREFERENCES['task_manager_path'],\
                                                          save = True,\
                                                          command = lambda k = 'task_manager_path': self.openFilename(k))
        
        PreferencesDialog.first = True
        
        self.preferences['progress_rate'] = ComboScale(self.preferencesFrame,\
                                                       frameLabel = "Progress rate",\
                                                       tagName = 'preferences_progress_rate',\
                                                       contents = [5,95,5])
        self.preferences['progress_rate'].setValue(PREFERENCES['progress_rate'])
            
        self.preferences['documentation_style'] = ComboRadiobutton(self.preferencesFrame,\
                                                                   frameLabel = "Documentation style",\
                                                                   tagName = 'preferences_documentation_style',\
                                                                   contents = ["html", "pdf"],\
                                                                   layout = (1,2))
        self.preferences['documentation_style'].setValue(PREFERENCES['documentation_style'])
                                        
        self.preferences['warning_check_for_new_version'] = ComboRadiobutton(self.preferencesFrame,\
                                                                              frameLabel = "check for new version at start up",\
                                                                              tagName = 'preferences_warning_check_for_new_version',\
                                                                              contents = ['yes', 'no'],\
                                                                              layout = (1,2))
        self.preferences['warning_check_for_new_version'].setValue(PREFERENCES['warning_check_for_new_version'].lower())

        self.preferences['database_version'] = ComboRadiobutton(self.preferencesFrame,\
                                                                frameLabel = "Atoms database version",\
                                                                tagName = 'preferences_database_version',\
                                                                contents = ['nmoldyn', 'mmtk'],\
                                                                layout = (1,2))
        self.preferences['database_version'].setValue(PREFERENCES['database_version'])
                        
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
        w = Button(box, text = "Save", width=10, command = self.savePreferences, default=ACTIVE)
        w.grid(row = 0, column = 1, sticky = E)
        w = Button(box, text = "Load", width=10, command = self.loadPreferences, default=ACTIVE)
        w.grid(row = 0, column = 2, sticky = E)
        w = Button(box, text = "OK", width=10, command = self.ok, default=ACTIVE)
        w.grid(row = 0, column = 3, sticky = E)
        
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
    
    def validate(self):
        return True
        
    def apply(self):
        """
        This method will set the nMOLDYN configuration.
        """
      
        prefs = {}
                     
        # The nMOLDYN configuration variables are set to the value defined in the dialog.
        for k in self.preferences.keys():
            if PREFERENCES.has_key(k):
                val = self.preferences[k].getValue() 
                if val is None:
                    prefs[k] = ''
                else:
                    prefs[k] = val
                    
        setPreferences(prefs)

    def changePage(self):
        
        page = self.pages.getValue()
        
        for w in self.preferencesFrame.winfo_children():
            w.grid_forget()
            
        if page == 'File handling':
            self.preferences['trajfile_path'].grid(row = 1, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['trajfile_path'].grid_columnconfigure(0, weight = 1)

            self.preferences['outputfile_path'].grid(row = 2, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['outputfile_path'].grid_columnconfigure(0, weight = 1)
                        
        elif page == 'External programs':
            if PLATFORM != 'DARWIN':
                self.preferences['acroread_path'].grid(row = 0, column = 0, padx = 2, pady = 5, sticky = 'NEW')
                self.preferences['acroread_path'].grid_columnconfigure(0, weight = 1)
            
            self.preferences['vmd_path'].grid(row = 1, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['vmd_path'].grid_columnconfigure(0, weight = 1)

            self.preferences['ncdump_path'].grid(row = 2, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['ncdump_path'].grid_columnconfigure(0, weight = 1)

            self.preferences['ncgen_path'].grid(row = 3, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['ncgen_path'].grid_columnconfigure(0, weight = 1)

            self.preferences['task_manager_path'].grid(row = 4, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['task_manager_path'].grid_columnconfigure(0, weight = 1)
            
        elif page == 'Miscellaneous':
            self.preferences['progress_rate'].grid(row = 0, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['progress_rate'].grid_columnconfigure(0, weight = 1)
            
            self.preferences['documentation_style'].grid(row = 1, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['documentation_style'].grid_columnconfigure(0, weight = 1)
        
            self.preferences['warning_check_for_new_version'].grid(row = 2, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['warning_check_for_new_version'].grid_columnconfigure(0, weight = 1)

            self.preferences['database_version'].grid(row = 3, column = 0, padx = 2, pady = 5, sticky = 'NEW')
            self.preferences['database_version'].grid_columnconfigure(0, weight = 1)
                                
    def openDirectory(self, prefKey):
        
        directory = askdirectory(parent = self, title = 'Select directory')
        
        if not directory:
            return

        self.preferences[prefKey].setValue(directory)
        
    def openFilename(self, prefKey):
        
        filename = askopenfilename(parent = self, title = 'Select filename')
        
        if not filename:
            return

        self.preferences[prefKey].setValue(filename)
                
    def savePreferences(self):
        """
        This method will save the preferences defined in the configuration dialog into a file.
        """

        try:
            # The Win32 config file directory.
            if PLATFORM == 'WINDOWS':
                cfgFileDir = os.path.join(os.environ['USERPROFILE'],'Application Data','nMOLDYN')
                initialName = 'nMOLDYN.ini'
            elif PLATFORM == 'DARWIN':
                cfgFileDir = os.path.join(os.environ['HOME'],'Library','Preferences','nMOLDYN')
                initialName = 'nMOLDYN.pref'
            else:
                cfgFileDir = os.environ['HOME']
                initialName = '.nMOLDYN'
        except KeyError:
            cfgFileDir = None
        else:
            if not os.path.exists(cfgFileDir):
                os.makedirs(cfgFileDir)

        filename = asksaveasfilename(parent = self,\
                                     title = 'Save nMOLDYN configuration file',\
                                     initialdir = cfgFileDir,\
                                     initialfile = initialName)
        if not filename:
            return

        pref = {}
                
        pref['trajfile_path'] = self.preferences['trajfile_path'].getValue()
        pref['outputfile_path'] = self.preferences['outputfile_path'].getValue()        
        pref['vmd_path'] = self.preferences['vmd_path'].getValue()
        pref['ncdump_path'] = self.preferences['ncdump_path'].getValue()
        pref['ncgen_path'] = self.preferences['ncgen_path'].getValue()
        pref['task_manager_path'] = self.preferences['task_manager_path'].getValue()
        pref['progress_rate'] = self.preferences['progress_rate'].getValue()
        pref['documentation_style'] = self.preferences['documentation_style'].getValue()
        pref['warning_check_for_new_version'] = self.preferences['warning_check_for_new_version'].getValue()
        pref['database_version'] = self.preferences['database_version'].getValue()

        if PLATFORM != 'DARWIN':
            pref['acroread_path'] = self.preferences['acroread_path'].getValue()
                    
        savePreferencesFile(filename, pref)
        
    def loadPreferences(self):
        """
        This method will load a preferences file and display its contents into the configuration dialog..
        """

        try:
            # The Win32 config file directory.
            if PLATFORM == 'WINDOWS':
                cfgFileDir = os.path.join(os.environ['USERPROFILE'],'Application Data','nMOLDYN')
                initialName = 'nMOLDYN.ini'
            elif PLATFORM == 'DARWIN':
                cfgFileDir = os.path.join('Library','Preferences','nMOLDYN')
                initialName = 'nMOLDYN.pref'
            # The POSIX config file directory.
            else:
                cfgFileDir = os.environ['HOME']
                initialName = '.nMOLDYN'
        except KeyError:
            cfgFileDir = None
        
        filename = askopenfilename(parent = self,\
                                   title = 'Open nMOLDYN configuration file',\
                                   initialdir = cfgFileDir,\
                                   initialfile = initialName)

        if not filename:
            return
        
        cfg = loadPreferencesFile(filename = filename)
        
        self.preferences.update(cfg)
