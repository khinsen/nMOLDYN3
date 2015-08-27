"""This modules implements I{File-->Frame snapshot} dialog.

Classes:
    * PDBSnapshotGeneratorDialog: creates I{File-->Frame snapshot} dialog used to 
      extract a PDB file from one or several trajectory frame(s).
"""

# The python distribution modules
import os
import re
import sys
from tempfile import mkdtemp

# The Tcl/Tk modules
from tkFileDialog import askopenfilename
from Tkinter import *

# The ScientificPython modules
from Scientific.IO.FortranFormat import FortranFormat, FortranLine
from Scientific import N as Num

# The MMTK distribution modules
from MMTK.PDB import PDBOutputFile
from MMTK import Units
from MMTK.Trajectory import Trajectory

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import load_trajectory_file
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import *

atom_format = FortranFormat('A6,I5,1X,A4,A1,A4,A1,I4,A1,3X,3F8.3,2F6.2,6X,A4,2A2')

class PDBSnapshotGeneratorDialog(PortableToplevel):
    """Sets up a dialog used to export one or several trajectory frames into a PDB file. 
    
    @note: if a trajectory has been previously loaded in nMOLDYN this will be the one proposed for extraction 
        by default. Otherwise, the user can still choose a trajectory to visualize from the dialog.
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

        self.trajectory = load_trajectory_file(GVAR['current_traj'])
            
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
        self.fileBrowser = ComboFileBrowser(settingsFrame,\
                                            frameLabel = "MMTK trajectory file",\
                                            tagName = 'frame_snapshot_mmtk_trajectory_file',\
                                            contents = '',\
                                            save = False,\
                                            command = self.openTrajectory)
        self.fileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.fileBrowser.grid_columnconfigure(0, weight = 1)
        self.fileBrowser.entry.bind('<Return>', self.openTrajectory)
        
        # The combo widget for the first step setting.
        self.selectedStepEntry = ComboStringEntry(settingsFrame,\
                                                  frameLabel = 'Selected frames',\
                                                  tagName = 'frame_snapshot_selected_frames',\
                                                  contents = '1')
        self.selectedStepEntry.grid(row = 1, column = 0, sticky = EW, padx = 2, pady = 2)
        self.selectedStepEntry.grid_columnconfigure(0, weight = 1)

        self.outputStyle = ComboRadiobutton(settingsFrame,\
                                            frameLabel = "PDB Output Style",\
                                            tagName = "pdb_output_style",\
                                            contents = ["MMTK", "nMOLDYN"],\
                                            layout = (1,2))
        self.outputStyle.grid(row = 2, column = 0, sticky = EW, padx = 2, pady = 2)
        
        # The PDB output file browser.
        self.pdbFileBrowser = ComboFileBrowser(settingsFrame,\
                                               frameLabel = "PDB output file",\
                                               tagName = 'frame_snapshot_pdb_output_file',\
                                               contents = '',\
                                               save = True)
        self.pdbFileBrowser.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        self.pdbFileBrowser.grid_columnconfigure(0, weight = 1)
        
        # If a trajectory has been loaded in nMOLDYN for simulation purposes (e.g. MSD ...) then it is proposed by default for
        # visualization.
        if self.trajectory is not None:
            # The control variables are updated with the informations coming from the loaded trajectory.
            self.fileBrowser.setValue(self.trajectory.filename)
            self.selectedStepEntry.setValue('1')            
            self.pdbFileBrowser.setValue(os.path.splitext(self.trajectory.filename)[0] + '_Frame.pdb')
        else:
            self.pdbFileBrowser.setValue('Frame.pdb')
            
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

        self.update_idletasks()

        self.apply()
        
    def cancel(self, event=None):

        # Put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    # Command hooks
    def validate(self):
                    
        if not self.trajectory:
            raise Error('No MMTK trajectory file loaded for extraction.')
            
        if self.trajectory.filename != self.fileBrowser.getValue():
            raise Error('Mismatch between the loaded trajectory and the displayed MMTK trajectory file name.\
You should validate that file by pressing Return on its corresponding entry.')
            
        # A frame selection must have been set.
        if not self.selectedStepEntry.getValue():
            self.selectedFrames = '1'
            LogMessage('warning', 'No frame selected. Frame 1 will be extracted by default.', ['gui'])
        else:
            self.selectedFrames = self.selectedStepEntry.getValue().strip()

        try:
            self.selectedFrames = eval(self.selectedFrames)
            # The entry can be an integer.
            if isinstance(self.selectedFrames, int):
                self.selectedFrames = [self.selectedFrames]
                
            # The entry can be a list or a tuple.
            elif isinstance(self.selectedFrames, (list,tuple)):
                pass

        except:
            try:
                temp = []
                self.selectedFrames = [temp.extend(range(int(v[0]),int(v[1])+1,int(v[2]))) for v in re.findall('(\d+):(\d+):(\d+)',self.selectedFrames)]
                self.selectedFrames = temp
            except:            
                raise Error('Wrong format for frame selection.')
                    
        # Check that every selected step is valid.
        for s in self.selectedFrames:
            if not isinstance(s, int):
                raise Error('Wrong entry for frame selection.')

            if (s <= 0) or (s > len(self.trajectory)):
                raise Error('Some selected frame number are not within [1,%d].' % len(self.trajectory))

        # A PDB output file name must have been set.
        self.pdbFile = self.pdbFileBrowser.getValue()
        if not self.pdbFile:
            raise Error('Please enter a PDB output file.')
                                
        return True
        
    def apply(self):

        try:
            # pid = the pid number for the job
            pid = os.getpid()
            suffix = '.'.join(['',str(pid),'pdb','nmoldyn'])

            # This just creates a unique temporary directory.
            dirName = mkdtemp(suffix)

            outputFile = open(self.pdbFile,'w')
            
            for f in self.selectedFrames:
                
                pdbFrame = os.path.join(dirName,''.join(['Frame',str(f),'.pdb']))
                
                if self.outputStyle.getValue().lower() == "mmtk":
                    # This MMTK function write a PDB file given a trajectory and a selected frame index. 
                    PDBOutputFile(pdbFrame).write(self.trajectory.universe, self.trajectory.configuration[f-1])
                    
                else:
                                
                    pdbFile = open(pdbFrame, 'w')
                
                    resNum = 0
                    atNum = 0
                
                    for obj in self.trajectory.universe.objectList():
                    
                        resNum += 1
                
                        data = {}
                
                        for at in obj.atomList():
                        
                            atNum += 1
                    
                            try: occ = at.occupancy
                            except AttributeError: occ = 0.0
                    
                            try: temp = at.temperature_factor
                            except AttributeError: temp = 0.0

                            type = 'ATOM'
                                                        
                            data['name'] = at.name.upper()
                            p = at.position(self.trajectory.configuration[f-1])/Units.Ang
                            data['serial_number'] = atNum % 100000
                            data['residue_number'] = resNum
                            data['residue_name'] = at.topLevelChemicalObject().name
                            data['alternate'] = ""
                            data['occupancy'] = occ
                            data['temperature_factor'] = temp
                            data['element'] = at.type.symbol
                    
                            line = [type]
                     
                            line = line + [data.get('serial_number', 1),
                                           data.get('name'),
                                           data.get('alternate', ''),
                                           data.get('residue_name', '').rjust(3),
                                           data.get('chain_id', ''),
                                           data.get('residue_number', 1),
                                           data.get('insertion_code', ''),
                                           p[0], p[1], p[2],
                                           data.get('occupancy', 0.),
                                           data.get('temperature_factor', 0.),
                                           data.get('segment_id', ''),
                                           data.get('element', '').rjust(2),
                                           data.get('charge', '')]
                                        
                            pdbFile.write(str(FortranLine(line, atom_format)) + '\n')
                                        
                    pdbFile.close()
                    
                pdbFile = open(pdbFrame,'r')
                data = pdbFile.read()
                pdbFile.close()
                os.unlink(pdbFrame)
                outputFile.write('REMARK Frame %d\n' % (f,))
                outputFile.write(data)

            outputFile.close()

            LogMessage('info', 'Frame extraction successful', ['gui'])

        except:
            raise Error('Error when extracting PDB frame.')
                
    def openTrajectory(self, event = None):
        """
        Ths method is called when the user clicks on the 'Browse' button of the trajectory visualization dialog.
        It opens a file browser. After the file selection some of the dialog widgets are updated with the informations
        coming from the loaded trajectory.
        """

        # Case where the user enters a file name directly in the entry widget without using the browser.
        if event is not None:
            if event.widget == self.fileBrowser.entry:
                filename = self.fileBrowser.getValue()
            else:
                return
            
        else:
            # The name of the NetCDF file to load.
            filename = askopenfilename(parent = self,\
                                       filetypes = [('NetCDF file','*.nc')],\
                                       initialdir = PREFERENCES['trajfile_path'])

        # The file must exist.
        if filename:
            try:
                # The trajectory is loaded.                
                self.trajectory = Trajectory(None, filename, 'r')
            except IOError:
                raise Error('Can not read the trajectory.')
            else:
                # The control variables are updated with the informations about the loaded trajectory.
                self.fileBrowser.setValue(filename)
                self.selectedStepEntry.setValue('1')

        return 'break'
