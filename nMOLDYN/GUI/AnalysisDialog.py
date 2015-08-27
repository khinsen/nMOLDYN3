"""This modules implements I{Analysis-->selected analysis} dialog.

Classes:
    * AnalysisDialog: creates I{Analysis-->selected analysis} dialog used to setup and/orrun 
      an analysis.
"""

# The python distribution modules
import os
import sys

# The Tcl/Tk modules
from Tkinter import *
from tkMessageBox import askyesno

from Scientific import N as Num

# The MMTK modules.
from MMTK.Trajectory import Trajectory

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Analysis.Templates import *
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.IO import load_trajectory_file
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.AtomOrderDialog import AtomOrderDialog
from nMOLDYN.GUI.QVectorsDialog import QVectorsDialog
from nMOLDYN.GUI.PyroServerDialog import PyroServerDialog
from nMOLDYN.GUI.SelectionDialog import SelectionDialog
from nMOLDYN.GUI.Widgets import *

from nMOLDYN.__pkginfo__ import __version__ as NMOLDYN_VERSION
    
class AnalysisDialog(PortableToplevel):
    """Builds the dialog for nMOLDYN analysis.
    """

    def __init__(self, parent, analysis, chemicalObjectInfo):
        """The constructor.

        @param parent: the parent widget.

        @param analysis: the analysis to setup.
        @type analysis: a class object of one of the nMOLDYN.Analysis subclasses.

        @param chemicalObjectInfo:.
        @type chemicalObjectInfo: 
        """
                
        PortableToplevel.__init__(self, parent)
        self.transient(parent)

        self.parent = parent

        exec file(analysis) in None, self.__dict__
        
        self.trajectory = GVAR['current_traj']                

        self.chemicalObjectInfo = chemicalObjectInfo

        self.title(self.db_longname)

        statusFrame = LabelFrame(self, text = 'Progress (in %)', relief = GROOVE, bd = 2)
        self.statusBar = StatusBar(statusFrame)

        body = Frame(self)
        self.initial_focus = self.body(body)
        body.grid(row = 0, column = 0, sticky = EW)  

        self.buttonbox()        

        statusFrame.grid(row = 2, padx = 3, pady = 3, sticky = EW)
        statusFrame.grid_columnconfigure(0, weight = 1)

        self.statusBar.grid(row = 0, column = 0)

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.resizable(width = NO, height = NO)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50, parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    def body(self, master):
        """Create dialog body. Return widget that should have initial focus.
        """
        
        # The frame that will contain the widgets whose the analysis dialog is built up.
        widgetsFrame = LabelFrame(master, text = 'Setup', bd = 2, relief = GROOVE)
        widgetsFrame.grid(row = 0, column = 0, padx = 3, pady = 3, sticky = EW)
        widgetsFrame.grid_columnconfigure(0, weight = 1)

        self.widgets = {}

        for widget in self.db_parameternames:

            widget = widget.lower()
                        
            if widget == 'angminmax':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "HB angle min:max values (in deg)",\
                                                        tagName = 'angle_min_max',\
                                                        contents = "150.0:180.0")

            elif widget == 'armodelorder':
                self.widgets[widget] = ComboIntegerEntry(widgetsFrame,\
                                                         frameLabel = "Model order", \
                                                         tagName = 'model_order',\
                                                         contents = 50)

            elif widget == 'atomorder':
                
                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "Atom order",\
                                                   tagName = 'atom_selection',\
                                                   contents = "Select",\
                                                   withEntry = '')
                

                
                t = Trajectory(None, self.trajectory, 'r')

                self.widgets[widget].command = lambda s = self, n = sorted(set([at.name for at in t.universe.atomList()])) : AtomOrderDialog(s,n)
                
                t.close()

            elif widget == 'truebreakstep':
                self.widgets[widget] = ComboIntegerEntry(widgetsFrame,\
                                                         frameLabel = "True break step", \
                                                         tagName = 'true_break_step',\
                                                         contents = 1)

            elif widget == 'comselection':
                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "Center of mass selection",\
                                                   tagName = 'center_of_mass_selection',\
                                                   contents = "Select",\
                                                   withEntry = 'all')

                self.widgets[widget].command = lambda s = self, t = 'subset', c = self.chemicalObjectInfo, w = self.widgets[widget]: SelectionDialog(s,t,c,w) 

            elif widget == 'deuteration':
                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "Deuteration selection",\
                                                   tagName = 'deuteration_selection',\
                                                   contents = "Select",\
                                                   withEntry = 'no')

                self.widgets[widget].command = lambda s = self, t = 'deuteration', c = self.chemicalObjectInfo, w = self.widgets[widget] : SelectionDialog(s,t,c,w)

            elif widget == 'differentiation':
                
                t = Trajectory(None, GVAR['current_traj'], 'r')
                
                if 'velocities' in t.variables():
                    diffOrder = (0,1,2,3,4,5)
                else:
                    diffOrder = (1,2,3,4,5)
                    
                t.close()
                    
                self.widgets[widget] = ComboSpinbox(widgetsFrame,\
                                                    frameLabel = "Differentiation order",\
                                                    tagName = 'differentiation_order',\
                                                    contents = diffOrder)

            elif widget == 'direction':
                self.widgets[widget] = ComboRadiobutton(widgetsFrame,\
                                                        frameLabel = "Direction",\
                                                        tagName = 'direction',\
                                                        contents = ["X", "Y", "Z"],\
                                                        layout = (1,3))
                
            elif widget == 'disminmax':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "HB distance min:max values",\
                                                        tagName = 'dis_min_max',\
                                                        contents = "0.0:0.35")

            elif widget == 'distanceunits':
                self.widgets[widget] = ComboRadiobutton(widgetsFrame,\
                                                        frameLabel = "Distance units",\
                                                        tagName = 'distanceunits',\
                                                        contents = ["nm", "ang", "fm"],\
                                                        layout = (1,3))
                
            elif widget == 'filter':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Pass-Band filter (in tHz)",\
                                                        tagName = 'pass_band_filer',\
                                                        contents = "0.0:1000.0")

            elif widget == 'frequencyunits':
                self.widgets[widget] = ComboRadiobutton(widgetsFrame,\
                                                        frameLabel = "Frequency units",\
                                                        tagName = 'frequencyunits',\
                                                        contents = ["THz", "rad s^-1", "cm^-1", "meV", "ueV"],\
                                                        layout = (1,5))
                
            elif widget == 'group':
                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "Group selection",\
                                                   tagName = 'group_selection',\
                                                   contents = "Select",\
                                                   withEntry = 'all')

                self.widgets[widget].command = lambda s = self, t = 'group', c = self.chemicalObjectInfo, w = self.widgets[widget] : SelectionDialog(s,t,c,w)

            elif widget == 'normalize':
                self.widgets[widget] = ComboCheckbutton(widgetsFrame,\
                                                        frameLabel = "Normalize",\
                                                        tagName = 'normalize',\
                                                        onvalue = "yes",\
                                                        offvalue = "no")

            elif widget == 'output':
                # The basename of the input trajectory name.
                baseName = self.db_shortname.upper() + '_' + os.path.basename(self.trajectory)

                try:
                    outputFilename = os.path.join(PREFERENCES['outputfile_path'],baseName)

                except TypeError:
                    outputFilename = baseName

                self.widgets[widget] = ComboFileBrowser(widgetsFrame,\
                                                        frameLabel = "Output file",\
                                                        tagName= "output_file",\
                                                        contents = outputFilename,\
                                                        save = True,\
                                                        filetypes = [("NetCDF file", ".nc"),])

            elif widget == 'phivalues':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Phi values (in deg)",\
                                                        tagName = 'phi_values',\
                                                        contents = "-180.0:180.0:10.0")

            elif widget == 'projection':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Project displacement on",\
                                                        tagName = 'project_displacement_on',\
                                                        contents = "no")

            elif widget == 'pyroserver':                

                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "Pyro server",\
                                                   tagName = 'pyro_server',\
                                                   contents = "Select",\
                                                   withEntry = 'monoprocessor')

                self.widgets[widget].command = lambda s = self, w = self.widgets[widget] : PyroServerDialog(s,w)

            elif widget == 'qshellvalues':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Q values (in nm-1)",\
                                                        tagName = 'qshell_values',\
                                                        contents = "0.0:100.0:1.0")

            elif widget == 'qshellwidth':
                self.widgets[widget] = ComboFloatEntry(widgetsFrame,\
                                                       frameLabel = "Q shell width (in nm-1)",\
                                                       tagName = 'qshell_width',\
                                                       contents = 1.0)

            elif widget == 'qunits':
                self.widgets[widget] = ComboRadiobutton(widgetsFrame,\
                                                        frameLabel = "Q units",\
                                                        tagName = 'qunits',\
                                                        contents = ["nm^-1", "ang^-1"],\
                                                        layout = (1,2))

            elif widget == 'qvectors':
                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "Q vectors",\
                                                   tagName = 'qvectors',\
                                                   contents = "Setup",\
                                                   withEntry = {'qgeometry': 'spatial', \
                                                                'qshellwidth': 1.0, \
                                                                'qshellvalues': '0.0:10.0:1.0', \
                                                                'qvectorspershell': 50,\
                                                                'hkls' : None})

                self.widgets[widget].command = lambda s = self, w = self.widgets[widget]: QVectorsDialog(s,w) 
                
            elif widget == 'qvectorsdirection':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Q vectors direction",\
                                                        tagName = 'qvectors_direction',\
                                                        contents = "no")

            elif widget == 'qvectorsgenerator':
                self.widgets[widget] = ComboRadiobutton(widgetsFrame,\
                                                        frameLabel = "Q vectors generator",\
                                                        tagName = 'qvectors_generator',\
                                                        contents = ["3D isotropic", "2D isotropic", "anisotropic"],\
                                                        layout = (1,3))

            elif widget == 'qvectorspershell':
                self.widgets[widget] = ComboIntegerEntry(widgetsFrame,\
                                                         frameLabel = "Q vectors per shell",\
                                                         tagName= 'qvectors_per_shell',\
                                                         contents = 50)

            elif widget == 'referencedirection':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Reference direction",\
                                                        tagName = 'reference_direction',\
                                                        contents = "0.0,0.0,1.0")

            elif widget == 'referenceframe':
                self.widgets[widget] = ComboIntegerEntry(widgetsFrame,\
                                                         frameLabel = "Reference frame",\
                                                         tagName= 'reference_frame',\
                                                         contents = 1)

            elif widget == 'removetranslation':
                self.widgets[widget] = ComboCheckbutton(widgetsFrame,\
                                                        frameLabel = "Remove translation",\
                                                        tagName = 'remove_translation',\
                                                        onvalue = "yes",\
                                                        offvalue = "no")

            elif widget == 'resolution':
                self.widgets[widget] = ComboFloatEntry(widgetsFrame,\
                                                       frameLabel = "Resolution (FWHM meV)",\
                                                       tagName = 'resolution',\
                                                       contents = 10.0)

                self.updateResolution()                               

                self.widgets[widget].entry.variable.trace_vdelete('w', self.widgets[widget].entry.cbname)                

                self.widgets[widget].entry.variable.trace_variable('w', self.updateResolution)

                self.widgets[widget].entry.variable.trace_variable('w', self.widgets[widget].entry._callback)

            elif widget == 'rvalues':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Distances (in nm)",\
                                                        tagName = 'distances',\
                                                        contents = "0.0:1.0:0.1")

            elif widget == 'stepwiserbt':
                self.widgets[widget] = ComboCheckbutton(widgetsFrame,\
                                                        frameLabel = "Stepwise RBT",\
                                                        tagName = 'stepwise_rbt',\
                                                        onvalue = "yes",\
                                                        offvalue = "no")

            elif widget == 'storerbtdetails':
                self.widgets[widget] = ComboCheckbutton(widgetsFrame,\
                                                        frameLabel = "Store RBT details",\
                                                        tagName = 'store_rbt_details',\
                                                        onvalue = "yes",\
                                                        offvalue = "no")

            elif widget in ['subset','subset1','subset2']:
                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "%s selection" % widget.capitalize(),\
                                                   tagName = 'subset_selection',\
                                                   contents = "Select",\
                                                   withEntry = 'all')

                self.widgets[widget].command = lambda s = self, t = 'subset', c = self.chemicalObjectInfo, w = self.widgets[widget]: SelectionDialog(s,t,c,w) 

            elif widget == 'target':
                self.widgets[widget] = ComboButton(widgetsFrame,\
                                                   frameLabel = "Target selection",\
                                                   tagName = 'target_selection',\
                                                   contents = "Select",\
                                                   withEntry = 'all')

                if self.db_shortname == "GMFT":
                    self.widgets[widget].command = lambda s = self, t = 'subset', c = self.chemicalObjectInfo, w = self.widgets[widget]: SelectionDialog(s,t,c,w) 
                else:
                    self.widgets[widget].command = lambda s = self, t = 'group', c = self.chemicalObjectInfo, w = self.widgets[widget] : SelectionDialog(s,t,c,w)
                  
            elif widget == 'temperature':
                self.widgets[widget] = ComboFloatEntry(widgetsFrame,\
                                                       frameLabel = "Temperature (in K)",\
                                                       tagName = 'temperature',\
                                                       contents = 1.0)
                
            elif widget == 'thetavalues':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Theta values (in deg)",\
                                                        tagName = 'theta_values',\
                                                        contents = "0.0:180.0:10.0")

            elif widget == 'thickness':
                self.widgets[widget] = ComboFloatEntry(widgetsFrame,\
                                                       frameLabel = "Thickness (in nm)",\
                                                       tagName = 'thickness',\
                                                       contents = 0.05)
                
            elif widget == 'timeinfo':
                
                t = load_trajectory_file(self.trajectory)                                
                timeInfo = '%d:%d:%d' % (1, len(t), 1)
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Frame selection",\
                                                        tagName = 'frame_selection',\
                                                        contents = timeInfo)

                self.updateTimeInfo(t)

                self.widgets[widget].entry.variable.trace_variable('w', lambda: self.updateTimeInfo(t))

            elif widget == 'timeunits':
                self.widgets[widget] = ComboRadiobutton(widgetsFrame,\
                                                        frameLabel = "Time units",\
                                                        tagName = 'timeunits',\
                                                        contents = ["ps", "ns", "fs"],\
                                                        layout = (1,3))
                
            elif widget == 'trajectory':
                self.widgets[widget] = ComboLabel(widgetsFrame,\
                                                  frameLabel = "Trajectory file",\
                                                  tagName = 'tajectory_file',\
                                                  contents = self.trajectory)

            elif widget == 'weights':
                self.widgets[widget] = ComboRadiobutton(widgetsFrame,\
                                                        frameLabel = "Weights",\
                                                        tagName = 'weights',\
                                                        contents = ["equal", "mass", "coherent", "incoherent", "atomicNumber"],\
                                                        layout = (2,3))

            elif widget == 'wignerindexes':
                self.widgets[widget] = ComboStringEntry(widgetsFrame,\
                                                        frameLabel = "Wigner indexes",\
                                                        tagName = 'wigner_indexes',\
                                                        contents = "0,0,0")

            # The tagName for the combo widget is set.
            setattr(self.widgets[widget],'tagName','%s_%s' % (self.db_shortname.lower(),self.widgets[widget].tagName.lower()))

            # And displayed into the analysis dialog.
            self.widgets[widget].grid(column = 0, sticky = EW, padx = 2, pady = 2)
            self.widgets[widget].grid_columnconfigure(0, weight = 1)

        if hasattr(self, 'db_modulator'):
            for k, v in self.db_modulator.items():
                k = k.lower()
                if k == 'weights':
                    self.widgets[k].setValue(v)

                elif k == 'differentiation':
                    self.widgets[k].spinbox.config(values = v)

        return None

    def buttonbox(self):
        """Add standard button box.
        """

        # The frame that contains the 'Cancel' and 'OK' buttons.
        box = LabelFrame(self, text = 'Actions', bd = 2, relief = GROOVE)
        box.grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3, ipady = 5)
        box.grid_columnconfigure(0, weight = 1)

        # The 'Cancel' button.
        w = Button(box, text = 'Cancel', width = 8, command = self.cancel)
        w.grid(row = 0, column = 0, padx = 3, pady = 3, sticky = E)

        # The 'Save' button.
        w = Button(box, text = 'Save', width = 8, command = lambda : self.ok(runMode = 'save'))
        w.grid(row = 0, column = 1, padx = 3, pady = 3)

        # The 'Run' button.
        w = Button(box, text = 'Run', width = 8, command = lambda : self.ok(runMode = 'run'), default = ACTIVE)
        w.grid(row = 0, column = 2, padx = 3, pady = 3)

        # The 'Save and Run' button.
        w = Button(box, text = 'Save/Run', width = 8, command = lambda : self.ok(runMode = 'saveandrun'))
        w.grid(row = 0, column = 3, padx = 3, pady = 3)

        if self.db_estimable:
            # The 'Estimate' button.
            w = Button(box, text = 'Estimate', width = 8, command = lambda : self.ok(runMode = 'estimate'))
            w.grid(row = 0, column = 4, padx = 3, pady = 3)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

    # Standard button semantics.
    def ok(self, event = None, runMode = 'run'):

        if not self.validate(runMode):
            self.initial_focus.focus_set()
            return

        self.update_idletasks()

        self.apply(runMode)

    def cancel(self, event=None):

        # Put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    # Command hooks
    def validate(self, runMode):

        self.parameters = {'version' : NMOLDYN_VERSION}

        # Loop over the self.widgets dictionnary.
        for widgetName, widget in self.widgets.items():
        
            try:
                
                if widgetName == "qvectors":                                    
                    self.parameters[widgetName] = eval(widget.getValue(), {"__builtins__" : None})
                else:
                    self.parameters[widgetName] = widget.getValue()

            except:
                LogMessage('warning','Bad input for %s. Please try again.' % widgetName,['gui'])
                return False

        return True

    def apply(self, runMode):

        try:                            
            # The widgets are disabled during the analysis.
            for w in self.winfo_children():
                self.widgetsState(w, DISABLED)

            if runMode == 'estimate':
                self.estimateAnalysis()

            elif runMode == 'save':
                self.saveAnalysis()

            elif runMode == 'run':
                self.runAnalysis()

            elif runMode == 'saveandrun':
                self.saveAndRunAnalysis()

            else:
                raise

        except:
            # Do not print anything in case of an error because, a message should have been displayed deeper.
            pass

        for w in self.winfo_children():
            self.widgetsState(w, NORMAL)

    def estimateAnalysis(self):
        """Estimates the time taken by the analysis directly from the GUI.
        """

        self.parameters['estimate'] = 'yes'
        self.actualAnalysis = eval('%s_serial(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)

        t = self.actualAnalysis.runAnalysis()

        # The information dialog with the predicted analysis time.
        LogMessage('info','\nThe analysis should take:\n %(days)s days %(hours)s hours %(minutes)s minutes %(seconds)s seconds.\n' % t,['gui'])

    def saveAnalysis(self):
        """Saves a python script of the analysis that can be run independantly of the GUI.
        """
        
        # Dialog box to choose a name for the python run script.
        pythonScript = asksaveasfilename(filetypes = [("nMOLDYN autostart file", "*.py"), ("nMOLDYN input file", "*.nmi")],\
                                         title = "Browse python file",\
                                         initialdir = PREFERENCES['outputfile_path'])

        if not pythonScript:
            return
                
        self.parameters['estimate'] = 'no'
        
        if not self.parameters.has_key('pyroserver'):
            self.actualAnalysis = eval('%s_serial(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)
            return

        if self.parameters['pyroserver'].lower() == 'monoprocessor':
            self.actualAnalysis = eval('%s_serial(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)

        elif self.parameters['pyroserver'][:14].lower() == 'multiprocessor':
            self.actualAnalysis = eval('%s_parallel(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)

        elif self.parameters['pyroserver'][:7].lower() == 'cluster':
            self.actualAnalysis = eval('%s_parallel(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)
            
        else:
            raise Error('The pyro server was not set properly')        

        self.actualAnalysis.saveAnalysis(pythonScript)

    def runAnalysis(self):
        """Runs the analysis directly from the GUI.
        """

        self.parameters['estimate'] = 'no'

        if not self.parameters.has_key('pyroserver'):
            self.actualAnalysis = eval('%s_serial(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)

        else:
            if self.parameters['pyroserver'].lower() == 'monoprocessor':
                self.actualAnalysis = eval('%s_serial(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)

            elif self.parameters['pyroserver'][:14].lower() == 'multiprocessor':
                self.actualAnalysis = eval('%s_parallel(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)

            elif self.parameters['pyroserver'][:7].lower() == 'cluster':
                self.actualAnalysis = eval('%s_parallel(parameters = self.parameters, statusBar = self.statusBar)' % self.db_internalname)
            else:
                raise Error('The pyro server was not set properly')

        t = self.actualAnalysis.runAnalysis()
        self.parent.info.insert(contents = '\n\n\nPERFORMED ANALYSIS:\n\n')
        self.parent.info.insert(contents = self.actualAnalysis.information)

        message = 'Analysis run successfully in \n %(days)s days %(hours)s hours %(minutes)s minutes %(seconds)s seconds.\n\n' % t
        if self.actualAnalysis.toPlot is None:
            # The information dialog with the predicted analysis time.
            LogMessage('info',message,['gui'])
        else:
            # Try to import the matplotlib module.
            try:
                from nMOLDYN.GUI.PlotNetCDFVariableDialog import PlotNetCDFVariableDialog
            except:
                return    
            else:
                message += 'Do you want to plot the results right now ?'
                q = askyesno(title = 'Question', message = message)
                if q:
                    PlotNetCDFVariableDialog(self, **self.actualAnalysis.toPlot)

    def saveAndRunAnalysis(self):
        """
        This method is called when the user presses the 'Save and Run' button of an analysis dialog. 
        It saves a python script of the analysis and run the analysis directly from the GUI.
        """

        self.saveAnalysis()
        self.runAnalysis()

    def widgetsState(self, widget, state = DISABLED):
        """
        This method is launched when an analysis is run from the GUI. Il will disable all the widgets of the analysis dialog.
        """

        if widget.winfo_children():
            for w in widget.winfo_children():
                self.widgetsState(w, state)
        else:
            try:
                widget.config(state = state)

            except:
                return

    def updateTimeInfo(self, trajectory):
        """
        """

        timeInfo = self.widgets['timeinfo'].entry.getValue()

        updatedTimeInfo = 'Frame selection'

        # The value must be three semi colon separated integer.
        try:
            f, l, s = [int(v) for v in timeInfo.split(':')]

            # framemin  must be > 0.
            if f <= 0:
                raise

            # framemax must be >= framemin
            if l < f:
                raise

        # Otherwise, throws an error.
        except:
            updatedTimeInfo += ' (invalid frame selection)'

        else:


            f -= 1

            t = trajectory.time[f:l:s]

            if len(t) <= 1:
                dt = 1.0
            else:
                dt = t[1] - t[0]

            updatedFirst, updatedLast, updatedStep = [round(v,3) for v in [t[0],t[-1],dt]]

            updatedTimeInfo += ' (%s ps to %s ps step %s ps)' % (updatedFirst, updatedLast, updatedStep)

        self.widgets['timeinfo'].frameLabel.set(updatedTimeInfo)

    def updateResolution(self, *dummy):
        """
        """

        updatedResolution = "Resolution (FWHM in meV)"

        fwhmE = self.widgets['resolution'].getValue()

        try:

            fwhmT = 2.354820045**2/(2.0*Num.pi*1.5192669*fwhmE)

        except:
            updatedResolution += " <=> NaN"

        else:
            updatedResolution += " <==> %s ps" % round(fwhmT,3)

        self.widgets['resolution'].frameLabel.set(updatedResolution)
