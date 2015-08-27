"""This modules implements I{File --> Trajectory conversion --> converter} dialog.

Classes:
    * ConverterDialog: a generic abstract class for all the converters.
    * AmberNetCDFConverterDialog: sets up I{File-->Trajectory conversion --> Amber NetCDF to MMTK} dialog used to 
      convert an AMBER trajectory to a MMTK NetCDF trajectory.
    * CASTEPConverterDialog: sets up I{File-->Trajectory conversion --> CASTEP to MMTK} dialog used to 
      convert a CASTEP trajectory to a MMTK NetCDF trajectory.      
    * CHARMMConverterDialog: sets up I{File-->Trajectory conversion --> CHARMM/X-PLOR to MMTK} dialog used to 
      convert a CHARMM or X-PLOR trajectory to a MMTK NetCDF trajectory.
    * DL_POLYConverterDialog: sets up I{File-->Trajectory conversion --> DL_POLY to MMTK} dialog used to 
      convert a DL_POLY trajectory to a MMTK NetCDF trajectory.      
    * LAMMPSConverterDialog: sets up I{File-->Trajectory conversion --> LAMMPS to MMTK} dialog used to 
      convert a LAMMPS trajectory to a MMTK NetCDF trajectory.      
    * MaterialsStudioConverterDialog: sets up I{File-->Trajectory conversion --> MaterialsStudio --> Discover to MMTK} or 
      I{File-->Trajectory conversion --> MaterialsStudio --> Forcite to MMTK} dialogs used to convert a Discover or a 
      Forcite trajectory to a MMTK NetCDF trajectory.
    * NAMDConverterDialog: sets up I{File-->Trajectory conversion --> NAMD to MMTK} dialog used to 
      convert a NAMD trajectory to a MMTK NetCDF trajectory.
    * PDBConverterDialog: sets up I{File-->Trajectory conversion --> PDB to MMTK} dialog used to 
      convert a single PDB frame to a MMTK NetCDF trajectory.
    * VASPBackConverterDialog: sets up I{File-->Trajectory conversion --> MMTK to VASP} dialog used to 
      convert a MMTK NetCDF trajectory back to a VASP trajectory.
    * VASPConverterDialog: sets up I{File-->Trajectory conversion --> VASP to MMTK} dialog used to 
      convert a VASP trajectory to a MMTK NetCDF trajectory.
"""

import re

# The Tcl/Tk modules
from tkFileDialog import askopenfilename
from Tkinter import *

# The nMOLDYN modules
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import AmberNetCDFConverter, CASTEPConverter, CHARMMConverter,\
                            DL_POLYConverter, LAMMPSConverter, MaterialsStudioConverter, \
                            NAMDConverter, PDBConverter, VASPBackConverter, VASPConverter
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import *

class ConverterDialog(PortableToplevel):
    """This is a generic class for all the trajectory converter dialogs.
    Its creates their top level window and define the Cancel and OK buttons.     
    """
    
    def __init__(self, parent, title):

        PortableToplevel.__init__(self, parent)        
        self.transient(parent)
        
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
        
        self.parent.grab_set()
        
        self.destroy()

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

class AmberNetCDFConverterDialog(ConverterDialog):
    """Sets up a dialog for the conversion from an Amber NetCDF trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires a AMber NetCDF and PDB file and the time step in ps between two frames.
    """

    def __init__(self, parent):
        """The constructor.
        
        @param parent: the parent widget.
        """

        ConverterDialog.__init__(self, parent, 'Amber --> MMTK NetCDF converter')
        
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
                    
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the PDB input file.
        self.pdbFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "PDB input file",\
                                               tagName = 'amber_converter_pdb_input_file',\
                                               contents = '',\
                                               filetypes = [('PDB', '*.pdb'), ('All files', '*')],\
                                               save = False)
        self.pdbFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = W)

        # The combo widget to browse the Amber NetCDF input file.
        self.amberNetCDFFileBrowser = ComboFileBrowser(settingsFrame,\
                                                       frameLabel = "Amber NetCDF input file",\
                                                       tagName = 'amber_converter_netcdf_input_file',\
                                                       contents = '',\
                                                       filetypes = [('Amber NetCDF', '*.nc'), ('All files', '*')],\
                                                       save = False)
        self.amberNetCDFFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.amberNetCDFFileBrowser.grid_columnconfigure(0, weight = 1)
        self.amberNetCDFFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        # The combo widget to select the time step in ps between two consecutive frames..
        self.timeStepEntry = ComboFloatEntry(settingsFrame,\
                                             frameLabel = "Time step (in ps)",\
                                             tagName = 'amber_converter_time_step',\
                                             contents = 1.0)
        self.timeStepEntry.grid(row = 2, column = 0, padx = 2, pady = 2, sticky = EW)
        self.timeStepEntry.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = 'amber_converter_mmtk_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)

        return None
        
    # Command hooks
    def validate(self):
                    
        self.pdbFile = self.pdbFileBrowser.getValue()
        # A PDB input file must have been set.
        if not self.pdbFile:
            raise Error('Please enter a PDB input file.')
        
        # A NetCDF input file must have been set.
        self.amberNetCDFFile = self.amberNetCDFFileBrowser.getValue()
        if not self.amberNetCDFFile:
            raise Error('Please enter an Amber NetCDF input file.')
                    
        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'

        self.timeStep = self.timeStepEntry.getValue()
        if self.timeStep <= 0.0:
            raise Error('Time step must be > 0.')

        return True
        
    def apply(self):
        """
        This method starts the conversion.
        """
        
        converter = AmberNetCDFConverter(self.pdbFile, self.amberNetCDFFile, self.outputFile, self.timeStep)
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
                    
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.amberNetCDFFileBrowser.getValue():
            self.outputFileBrowser.setValue(os.path.splitext(self.amberNetCDFFileBrowser.getValue())[0] + '_mmtk.nc')
                
class CASTEPConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a CASTEP trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires only the CASTEP file.
    """

    def __init__(self, parent):
        """
        The constructor.
        """

        ConverterDialog.__init__(self, parent, 'CASTEP --> MMTK NetCDF converter')
                
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
        
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the CASTEP input file.
        self.castepFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "CASTEP input file",\
                                                  tagName = 'castep_converter_castep_input_file',\
                                                  contents = '',\
                                                  filetypes = [('CASTEP input file', '*.castep'), ('All files', '*')],\
                                                  save = False)
        self.castepFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.castepFileBrowser.grid_columnconfigure(0, weight = 1)
        self.castepFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = 'castep_converter_mmtk_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
                
        return None        
                
    # Command hooks
    def validate(self):
                            
        # A DCD input file must have been set.
        self.castepFile = self.castepFileBrowser.getValue()
        if not self.castepFile:
            raise Error('Please enter a CASTEP input file.')
                    
        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'

        return True
        
    def apply(self):
        """
        This method starts the conversion.
        """

        converter = CASTEPConverter(self.castepFile, self.outputFile)
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
    
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.castepFileBrowser.getValue():
            self.outputFileBrowser.setValue(os.path.splitext(self.castepFileBrowser.getValue())[0] + '.nc')
            
class CHARMMConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a CHARMM/X-PLOR/NAMD trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires CHARMM/X-PLOR DCD and PDB files.
    """

    def __init__(self, parent):
        """The constructor.
        
        @param parent: the parent widget.
        """

        ConverterDialog.__init__(self, parent, 'CHARMM/X-PLOR --> MMTK NetCDF converter')

    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
        
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the PDB input file.
        self.pdbFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "PDB input file",\
                                               tagName = "charmm_converter_pdb_input_file",\
                                               contents = '',\
                                               filetypes = [('PDB', '*.pdb'), ('All files', '*')],\
                                               save = False)
        self.pdbFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.pdbFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the DCD input file.
        self.dcdFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "DCD input file",\
                                               tagName = "charmm_converter_dcd_input_file",\
                                               contents = '',\
                                               filetypes = [('DCD input file', '*.dcd'), ('All files', '*')],\
                                               save = False)
        self.dcdFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.dcdFileBrowser.grid_columnconfigure(0, weight = 1)
        self.dcdFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = "charmm_converter_mmtk_netcdf_output_file",\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
        
        return None
        
    # Command hooks
    def validate(self):
                    
        self.pdbFile = self.pdbFileBrowser.getValue()
        # A PDB input file must have been set.
        if not self.pdbFile:
            raise Error('Please enter a PDB input file.')
        
        # A DCD input file must have been set.
        self.dcdFile = self.dcdFileBrowser.getValue()
        if not self.dcdFile:
            raise Error('Please enter a DCD input file.')
                    
        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'

        return True
        
    def apply(self):
        """
        This method starts the conversion.
        """

        converter = CHARMMConverter(self.pdbFile, self.dcdFile, self.outputFile)
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
            
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.dcdFileBrowser.getValue():
            self.outputFileBrowser.setValue(os.path.splitext(self.dcdFileBrowser.getValue())[0] + '.nc')
        
class DL_POLYConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a DL_POLY trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires the DL_POLY FIELD and HISTORY files and additional information to specify 
    the atoms whose names is not sufficient to determine which element they are.
    """

    def __init__(self, parent):
        """The constructor.
        
        @param parent: the parent widget.
        """

        ConverterDialog.__init__(self, parent, 'DL_POLY --> MMTK NetCDF converter')
        
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
        
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the DL_POLY FIELD file.
        self.fieldFileBrowser = ComboFileBrowser(settingsFrame,\
                                                 frameLabel = "FIELD input file",\
                                                 tagName = "dlpoly_converter_field_input_file",\
                                                 contents = '',\
                                                 filetypes = [('FIELD files', 'FIELD*'),('All files', '*')],\
                                                 save = False)
        self.fieldFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.fieldFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the DL_POLY HISTORY file.
        self.historyFileBrowser = ComboFileBrowser(settingsFrame,\
                                                   frameLabel = "HISTORY input file",\
                                                   tagName = "dlpoly_converter_history_input_file",\
                                                   contents = '',\
                                                   filetypes = [('HISTORY files', 'HISTORY*'),('All files', '*')],\
                                                   save = False)
        self.historyFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.historyFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to specify the atoms whose names is not sufficient to determine which element they are.
        self.specialAtomsEntry = ComboStringEntry(settingsFrame,\
                                                  frameLabel = "Special atoms",\
                                                  tagName = 'dlpoly_converter_special_atoms',\
                                                  contents = '')
        self.specialAtomsEntry.grid(row = 2, column = 0, padx = 2, pady = 2, sticky = EW)
        self.specialAtomsEntry.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = 'dlpoly_converter_mmtk_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)

        return None
        
    # Command hooks
    def validate(self):
                    
        # A FIELD input file must have been set.
        self.fieldFile = self.fieldFileBrowser.getValue()
        if not self.fieldFile:
            raise Error('Please enter a FIELD input file.')
        
        # A HISTORY input file must have been set.
        self.historyFile = self.historyFileBrowser.getValue()
        if not self.historyFile:
            raise Error('Please enter a HISTORY input file.')

        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'

        self.specialAtoms = self.specialAtomsEntry.getValue()
        if not self.specialAtoms:
            self.specialAtoms = {}
                
        else:
            try:
                self.specialAtoms = dict(re.findall('(\w+):(\w+)',self.specialAtoms))
            except:
                raise Error('Wrong format for special atoms declaration. \
Must be of the form special atom name1:corresponding element,...')
            
        return True
    
    def apply(self):        

        converter = DL_POLYConverter(self.fieldFile, self.historyFile, self.outputFile, self.specialAtoms)            
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
        
class LAMMPSConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a LAMMPS trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires only the LAMMPS file.
    """

    def __init__(self, parent):
        """
        The constructor.
        """

        ConverterDialog.__init__(self, parent, 'LAMMPS --> MMTK NetCDF converter')
                
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
        
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the LAMMPS config file.
        self.lammpsCfgFileBrowser = ComboFileBrowser(settingsFrame,\
                                                     frameLabel = "LAMMPS config file",\
                                                     tagName = 'lammps_converter_lammps_config_file',\
                                                     contents = '',\
                                                     filetypes = [('LAMMPS config file', '*.config'), ('All files', '*')],\
                                                     save = False)
        self.lammpsCfgFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.lammpsCfgFileBrowser.grid_columnconfigure(0, weight = 1)
        self.lammpsCfgFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        # The combo widget to browse the LAMMPS input file.
        self.lammpsFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "LAMMPS input file",\
                                                  tagName = 'lammps_converter_lammps_input_file',\
                                                  contents = '',\
                                                  filetypes = [('LAMMPS input file', '*.lammps'), ('All files', '*')],\
                                                  save = False)
        self.lammpsFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.lammpsFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = 'lammps_converter_mmtk_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 2, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
                
        return None        
                
    # Command hooks
    def validate(self):

        self.lammpsCfgFile = self.lammpsCfgFileBrowser.getValue()
        if not self.lammpsCfgFile:
            raise Error('Please enter a LAMMPS config file.')
                            

        self.lammpsFile = self.lammpsFileBrowser.getValue()
        if not self.lammpsFile:
            raise Error('Please enter a LAMMPS input file.')
                    
        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'

        return True
        
    def apply(self):
        """
        This method starts the conversion.
        """

        converter = LAMMPSConverter(self.lammpsCfgFile, self.lammpsFile, self.outputFile)
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
    
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.lammpsCfgFileBrowser.getValue():
            self.outputFileBrowser.setValue(os.path.splitext(self.lammpsCfgFileBrowser.getValue())[0] + '.nc')

class MaterialsStudioConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a MaterialsStudio trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires MaterialsStudio XTD or XSD file and HIS or TRJ files according to the MaterialsStudio modules used to 
    produce the trajectory (HIS for Discover and TRJ for Forcite).
    """

    def __init__(self, parent, module):
        """The constructor.
        
        @param parent: the parent widget.
        """
        
        self.module = module.lower()

        ConverterDialog.__init__(self, parent, 'MaterialsStudio-%s --> MMTK NetCDF converter' % module)
        
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
                
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)            
        
        # The combo widget to browse the XTD or the XSD input file.
        self.xtdxsdFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "XTD/XSD input file",\
                                                  tagName = '%s_converter_xtd_xsd_input_file' % self.module,\
                                                  contents = '',\
                                                  filetypes = [('XTD input file', '*.xtd'), ('XSD input file', '*.xsd'), ('All files', '*')],\
                                                  save = False)
        self.xtdxsdFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.xtdxsdFileBrowser.grid_columnconfigure(0, weight = 1)

        if self.module == 'discover':
            # The combo widget to browse the HIS input file.
            self.histrjFileBrowser = ComboFileBrowser(settingsFrame,\
                                                      frameLabel = "HIS input file",\
                                                      tagName = 'discover_converter_his_input_file',\
                                                      contents = '',\
                                                      filetypes = [('HIS input file', '*.his'), ('All files', '*')],\
                                                      save = False)
        elif self.module == 'forcite':
            # The combo widget to browse the TRJ input file.
            self.histrjFileBrowser = ComboFileBrowser(settingsFrame,\
                                                      frameLabel = "TRJ input file",\
                                                      tagName = 'forcite_converter_trj_input_file',\
                                                      contents = '',\
                                                      filetypes = [('TRJ input file', '*.trj'), ('All files', '*')],\
                                                      save = False)
        
        self.histrjFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.histrjFileBrowser.grid_columnconfigure(0, weight = 1)
        self.histrjFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = '%s_converter_mmtk_netcdf_output_file' % self.module,\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
        
        return None
        
    # Command hooks
    def validate(self):
        
        # A XTD or XSD input file must have been set.
        self.xtdxsdFile = self.xtdxsdFileBrowser.getValue()
        if not self.xtdxsdFile:
            raise Error('Please enter a XTD or XSD input file.')

        # A HIS or TRJ input file must have been set.
        self.histrjFile = self.histrjFileBrowser.getValue()
        if not self.histrjFile:
            raise Error('Please enter a HIS or TRJ input file.')
                    
        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'
                        
        return True
    
    def apply(self):

        converter = MaterialsStudioConverter(self.module, self.xtdxsdFile, self.histrjFile, self.outputFile) 
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
        
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.histrjFileBrowser.getValue():
            self.outputFileBrowser.setValue(os.path.splitext(self.histrjFileBrowser.getValue())[0] + '.nc')
                                        
class NAMDConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a CHARMM/X-PLOR/NAMD trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires CHARMM/X-PLOR/NAMD DCD and PDB files and optionnaly the simulation box dimensions.
    """

    def __init__(self, parent):
        """
        The constructor.
        """

        ConverterDialog.__init__(self, parent, 'NAMD --> MMTK NetCDF converter')

    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
        
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the PDB input file.
        self.pdbFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "PDB input file",\
                                               tagName = 'namd_converter_pdb_input_file',\
                                               contents = '',\
                                               filetypes = [('PDB', '*.pdb'), ('All files', '*')],\
                                               save = False)
        self.pdbFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.pdbFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the DCD input file.
        self.dcdFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "DCD input file",\
                                               tagName = 'namd_converter_dcd_input_file',\
                                               contents = '',\
                                               filetypes = [('DCD input file', '*.dcd'), ('All files', '*')],\
                                               save = False)
        self.dcdFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.dcdFileBrowser.grid_columnconfigure(0, weight = 1)
        self.dcdFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        # The combo widget to browse the XST input file.
        self.xstFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "XST input file",\
                                               tagName = 'namd_converter_xst_input_file',\
                                               contents = '',\
                                               filetypes = [('XST input file', '*.xst'), ('All files', '*')],\
                                               save = False)
        self.xstFileBrowser.grid(row = 2, column = 0, padx = 2, pady = 2, sticky = EW)
        self.xstFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = 'namd_converter_mmtk_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
                
        return None        
        
    # Command hooks
    def validate(self):
                    
        self.pdbFile = self.pdbFileBrowser.getValue()
        # A PDB input file must have been set.
        if not self.pdbFile:
            raise Error('Please enter a PDB input file.')
        
        # A DCD input file must have been set.
        self.dcdFile = self.dcdFileBrowser.getValue()
        if not self.dcdFile:
            raise Error('Please enter a DCD input file.')

        self.xstFile = self.xstFileBrowser.getValue()
                    
        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'

        return True
        
    def apply(self):
        """
        This method starts the conversion.
        """

        converter = NAMDConverter(self.pdbFile, self.dcdFile, self.xstFile, self.outputFile)
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
    
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.dcdFileBrowser.getValue():
            self.outputFileBrowser.setValue(os.path.splitext(self.dcdFileBrowser.getValue())[0] + '.nc')
                
class PDBConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a PDB file to a MMTK NetCDF trajectory.
    
    @note: the conversion requires only the PDB file.
    """

    def __init__(self, parent):
        """
        The constructor.
        """

        ConverterDialog.__init__(self, parent, 'PDB --> MMTK NetCDF converter')
        
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
        
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the PDB input file.
        self.pdbFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "PDB input file",\
                                               tagName = 'pdb_converter_pdb_input_file',\
                                               contents = '',\
                                               filetypes = [('PDB', '*.pdb'), ('All files', '*')],\
                                               save = False)
        self.pdbFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.pdbFileBrowser.grid_columnconfigure(0, weight = 1)
        self.pdbFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        self.make_periodic_cb = ComboCheckbutton(settingsFrame,\
                                                 frameLabel = "Force periodicity",\
                                                 tagName = 'pdb_converter_make_periodic',\
                                                 onvalue = "yes",\
                                                 offvalue = "no")
        self.make_periodic_cb.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.make_periodic_cb.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = 'pdb_converter_mmtk_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 2, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
                
        return None        
        
    # Command hooks
    def validate(self):
                    
        self.pdbFile = self.pdbFileBrowser.getValue()
        # A PDB input file must have been set.
        if not self.pdbFile:
            raise Error('Please enter a PDB input file.')
        
        if self.make_periodic_cb.getValue() == "yes":
            self.make_periodic = True
        else:
            self.make_periodic = False
                            
        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'

        return True
        
    def apply(self):
        """
        This method starts the conversion.
        """

        converter = PDBConverter(self.pdbFile, self.outputFile, self.make_periodic)
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()

    
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.pdbFileBrowser.getValue():
            self.outputFileBrowser.setValue(os.path.splitext(self.pdbFileBrowser.getValue())[0] + '.nc')

class VASPBackConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a MMTK NetCDF trajectory to a VASP trajectory.
    
    @note: the conversion requires the MMTK NetCDF file.
    """

    def __init__(self, parent):
        """
        The constructor.
        """

        ConverterDialog.__init__(self, parent, 'MMTK NetCDF --> VASP converter')
        
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
        
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the MMTK input file.
        self.netcdfFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF input file",\
                                                  tagName = 'vasp_back_converter_mmtk_netcdf_input_file',\
                                                  contents = '',\
                                                  save = False)
        self.netcdfFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.netcdfFileBrowser.grid_columnconfigure(0, weight = 1)
        self.netcdfFileBrowser.entry.variable.trace('w', self.suggestOutputFilename)

        # The combo widget to browse the VASP CONTCAR file.
        self.contcarFileBrowser = ComboFileBrowser(settingsFrame,\
                                                   frameLabel = "CONTCAR output file",\
                                                   tagName = 'vasp_back_converter_contcar_output_file',\
                                                   contents = '',
                                                   save = True)
        self.contcarFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.contcarFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the VASP XDATCAR file.
        self.xdatcarFileBrowser = ComboFileBrowser(settingsFrame,\
                                                   frameLabel = "XDATCAR output file",\
                                                   tagName = 'vasp_back_converter_xdatcar_output_file',\
                                                   contents = '',
                                                   save = True)
        self.xdatcarFileBrowser.grid(row = 2, column = 0, padx = 2, pady = 2, sticky = EW)
        self.xdatcarFileBrowser.grid_columnconfigure(0, weight = 1)
        
        return None

    # Command hooks
    def validate(self):
        
        # The MMTK input file must have been set.
        self.netcdfFile = self.netcdfFileBrowser.getValue()
        if not self.netcdfFile:
            raise Error('Please enter a MMTK NetCDF input file.')
                    
        # A CONTCAR output file must have been set.
        self.contcarFile = self.contcarFileBrowser.getValue()
        if not self.contcarFile:
            raise Error('Please enter a CONTCAR or POSCAR output file.')
        
        # A XDATCAR output file must have been set.
        self.xdatcarFile = self.xdatcarFileBrowser.getValue()
        if not self.xdatcarFile:
            raise Error('Please enter a XDATCAR output file.')
                                
        return True
    
    def apply(self):        

        converter = VASPBackConverter(self.netcdfFile, self.contcarFile, self.xdatcarFile) 
        LogMessage('info', 'Conversion successful', ['gui'])
        self.cancel()
        
    def suggestOutputFilename(self, *dummy):
        """
        This method will propose a name for the output file based on the name of the loaded file.
        """

        if self.netcdfFileBrowser.getValue():
            name = self.netcdfFileBrowser.getValue()
            dirName = os.path.dirname(name)
            baseName = os.path.splitext(os.path.basename(name))[0]                        
            self.contcarFileBrowser.setValue(os.path.join(dirName, 'CONTCAR_'+baseName))
            self.xdatcarFileBrowser.setValue(os.path.join(dirName, 'XDATCAR_'+baseName))
            
class VASPConverterDialog(ConverterDialog):
    """Sets up a dialog used for the conversion from a VASP trajectory to a MMTK NetCDF trajectory.
    
    @note: the conversion requires the VASP CONTCAR or POSCAR and XDATCAR files.
    """

    def __init__(self, parent):
        """
        The constructor.
        """
        
        ConverterDialog.__init__(self, parent, 'VASP--> MMTK NetCDF converter')
        
    def body(self, master):
        """
        Create dialog body. Return widget that should have initial focus.
        """
                
        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the VASP XDATCAR file.
        self.xdatcarFileBrowser = ComboFileBrowser(settingsFrame,\
                                                   frameLabel = "XDATCAR input file",\
                                                   tagName = 'vasp_converter_xdatcar_input_file',\
                                                   contents = '',\
                                                   filetypes = [('XDATCAR files', 'XDATCAR*'),('All files', '*')],\
                                                   save = False)
        self.xdatcarFileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.xdatcarFileBrowser.grid_columnconfigure(0, weight = 1)
        
        # The combo widget to browse the VASP CONTCAR file.
        self.contcarFileBrowser = ComboFileBrowser(settingsFrame,\
                                                   frameLabel = "POSCAR/CONTCAR input file (required for version 4)",\
                                                   tagName = 'vasp_converter_poscar_contcar_input_file',\
                                                   contents = '',\
                                                   filetypes = [('POSCAR files', 'POSCAR*'), ('CONTCAR files', 'CONTCAR*'),('All files', '*')],\
                                                   save = False)
        self.contcarFileBrowser.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = EW)
        self.contcarFileBrowser.grid_columnconfigure(0, weight = 1)

        # The combo widget to the atomic contents of the system..
        self.atomSymbolsEntry = ComboStringEntry(settingsFrame,\
                                                  frameLabel = "Atom symbols (required for version 4)",\
                                                  tagName = 'vasp_converter_atom_symbols',\
                                                  contents = '')
        self.atomSymbolsEntry.grid(row = 2, column = 0, padx = 2, pady = 2, sticky = EW)
        self.atomSymbolsEntry.grid_columnconfigure(0, weight = 1)
            
        # The combo widget to select the time step in ps between two consecutive frames..
        self.timeStepEntry = ComboFloatEntry(settingsFrame,\
                                             frameLabel = "Time step (in ps) (required for version 5)",\
                                             tagName = 'vasp_converter_time_step',\
                                             contents = '0.001')
        self.timeStepEntry.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        self.timeStepEntry.grid_columnconfigure(0, weight = 1)

        # The combo widget to browse the MMTK output file.
        self.outputFileBrowser = ComboFileBrowser(settingsFrame,\
                                                  frameLabel = "MMTK NetCDF output file",\
                                                  tagName = 'vasp_converter_mmtk_netcdf_output_file',\
                                                  contents = '',\
                                                  save = True)
        self.outputFileBrowser.grid(row = 4, column = 0, padx = 2, pady = 2, sticky = EW)
        self.outputFileBrowser.grid_columnconfigure(0, weight = 1)
        
        return None

    # Command hooks
    def validate(self):
                
        # A XDATCAR input file must have been set.
        self.xdatcarFile = self.xdatcarFileBrowser.getValue()
        if not self.xdatcarFile:
            raise Error('Please enter a XDATCAR input file.')

        # A CONTCAR input file must have been set.
        self.contcarFile = self.contcarFileBrowser.getValue()
        
        try:
            self.atomSymbols = re.split('[;,\s]+', self.atomSymbolsEntry.getValue())
        except:
            self.atomSymbols = None
                                
        self.timeStep = self.timeStepEntry.getValue()

        # The MMTK output file must have been set.
        self.outputFile = self.outputFileBrowser.getValue()
        if not self.outputFile:
            raise Error('Please enter a MMTK NetCDF output file.')
        
        if self.outputFile[-3:] != '.nc':
            self.outputFile += '.nc'
                        
        return True
    
    def apply(self):        

        try:
            converter = VASPConverter(xdatcarFile = self.xdatcarFile, \
                                      contcarFile = self.contcarFile, \
                                      atomSymbols = self.atomSymbols, \
                                      timeStep = self.timeStep, \
                                      outputFile = self.outputFile)
        except:
            raise Error("An error occured during the conversion.")
            
        else:                        
            LogMessage('info', 'Conversion successful', ['gui'])
        
        self.cancel()
                    