"""This modules implements I{View-->Plot} dialog.

Classes:
    * SettingsDialog: sets up the settings dialog.
    * ASCIIToNetCDFConversionDialog: creates I{View-->Plot} dialog used to plot NetCDF variables.
"""

# The python distribution modules
import os
import sys

# The matplotlib modules
import pylab
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.colors import ColorConverter, LogNorm, Normalize, rgb2hex
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor

# The Tcl/Tk modules
import _tkinter
from tkFileDialog import askopenfilename, asksaveasfilename
from Tkinter import *
import tkColorChooser

# The ScientificPython modules
from Scientific.IO.NetCDF import _NetCDFFile, NetCDFFile
from Scientific import N as Num

# The MMTK modules.
from MMTK.Trajectory import Trajectory

# The nMOLDYN modules
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.Widgets import *

# Some modules variables
# The tuple of the possible interpolations.
interpolations = ['bessel', 'bilinear', 'bicubic', 'blackman', 'catrom', 'hamming', 'hermite', 'hanning',\
                  'gaussian', 'kaiser', 'lanczos', 'mitchell', 'nearest', 'spline16', 'spline36', 'quadric', 'sinc']

# The tuple of the possible color maps.
colorMaps = ['autumn', 'bone', 'cool', 'copper', 'flag', 'gray', 'hot', 'hsv', 'jet', 'pink', 'prism',\
             'spring', 'summer', 'winter', 'spectral']

# The tuple of the possible line styles.
lineStyles = ['-', '--', '-.', ':', 'None']

# The tuple of the possible marker styles.
markerStyles = ['+', '.', '<', '>', 'o', 'p', 's', 'v', 'x', '|', 'None']

axisScales = ['linear', 'log']

class SettingsDialog(PortableToplevel):
    """Sets up a dialog tp perform some settings on the plots.
    """

    def __init__(self, parent):
        """The constructor.

        @param parent: the parent widget.
        """

        PortableToplevel.__init__(self, parent)        
        self.transient(parent)

        self.parent = parent       

        self.widgets = {}
        self.settings = {}

        self.selectedPlot = None

        body = Frame(self)
        self.initial_focus = self.body(body)
        body.grid(row = 0, column = 0, sticky = EW)

        self.buttonbox()        

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.resizable(width = NO, height = NO)

        self.title = 'Settings'        

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
        settingsFrame.grid_rowconfigure(0, weight = 1)

        self.gFrame = LabelFrame(settingsFrame, text = 'Global', bd = 2, relief = GROOVE)
        self.gFrame.grid(row = 0, column = 0, sticky = "NEW", padx = 3, pady = 3)
        self.gFrame.grid_columnconfigure(0, weight = 1)
        self.gFrame.grid_rowconfigure(0, weight = 1)

        # All the attributes concerning the global and plot settings are set.
        self.initGlobalSettings()

        # The widget to change plot title.
        self.widgets['plotTitle'] = ComboStringEntry(self.gFrame,frameLabel = 'Plot title', contents = self.settings['plotTitle'])
        self.widgets['plotTitle'].grid(row = 1, column = 0, sticky = "NEW")
        self.widgets['plotTitle'].grid_columnconfigure(0, weight = 1)        
        self.widgets['plotTitle'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'plotTitle'))

        # The widget to change X-axis label.
        self.widgets['xLabel'] = ComboStringEntry(self.gFrame,frameLabel = 'X label', contents = self.settings['xLabel'])
        self.widgets['xLabel'].grid(row = 2, column = 0, sticky = "NEW")
        self.widgets['xLabel'].grid_columnconfigure(0, weight = 1)        
        self.widgets['xLabel'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'xLabel'))

        # The widget to change Y-axis label.
        self.widgets['yLabel'] = ComboStringEntry(self.gFrame,frameLabel = 'Y label', contents = self.settings['yLabel'])
        self.widgets['yLabel'].grid(row = 3, column = 0, sticky = "NEW")
        self.widgets['yLabel'].grid_columnconfigure(0, weight = 1)        
        self.widgets['yLabel'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'yLabel'))

        # The widget to change X-axis ticks contents.
        self.widgets['xTicks'] = ComboStringEntry(self.gFrame,frameLabel = 'X ticks', contents = self.settings['xTicks'])
        self.widgets['xTicks'].grid(row = 4, column = 0, sticky = "NEW")
        self.widgets['xTicks'].grid_columnconfigure(0, weight = 1)        
        self.widgets['xTicks'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'xTicks'))

        # The widget to change X-axis ticks size.
        self.widgets['xTicksSize'] = ComboIntegerEntry(self.gFrame, frameLabel = 'X ticks size', contents = self.settings['xTicksSize'])
        self.widgets['xTicksSize'].grid(row = 5, column = 0, sticky = "NEW")
        self.widgets['xTicksSize'].grid_columnconfigure(0, weight = 1)        
        self.widgets['xTicksSize'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'xTicksSize'))

        # The widget to change X-axis ticks angle.
        self.widgets['xTicksAngle'] = ComboIntegerEntry(self.gFrame, frameLabel = 'X ticks angle', contents = self.settings['xTicksAngle'])
        self.widgets['xTicksAngle'].grid(row = 6, column = 0, sticky = "NEW")
        self.widgets['xTicksAngle'].grid_columnconfigure(0, weight = 1)        
        self.widgets['xTicksAngle'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'xTicksAngle'))

        # The widget to change Y-axis ticks contents.
        self.widgets['yTicks'] = ComboStringEntry(self.gFrame,frameLabel = 'Y ticks', contents = self.settings['yTicks'])
        self.widgets['yTicks'].grid(row = 7, column = 0, sticky = "NEW")
        self.widgets['yTicks'].grid_columnconfigure(0, weight = 1)        
        self.widgets['yTicks'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'yTicks'))

        # The widget to change Y-axis ticks size.
        self.widgets['yTicksSize'] = ComboIntegerEntry(self.gFrame, frameLabel = 'Y ticks size', contents = self.settings['yTicksSize'])
        self.widgets['yTicksSize'].grid(row = 8, column = 0, sticky = "NEW")
        self.widgets['yTicksSize'].grid_columnconfigure(0, weight = 1)        
        self.widgets['yTicksSize'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'yTicksSize'))

        # The widget to change Y-axis ticks angle.
        self.widgets['yTicksAngle'] = ComboIntegerEntry(self.gFrame, frameLabel = 'Y ticks angle', contents = self.settings['yTicksAngle'])
        self.widgets['yTicksAngle'].grid(row = 9, column = 0, sticky = "NEW")
        self.widgets['yTicksAngle'].grid_columnconfigure(0, weight = 1)        
        self.widgets['yTicksAngle'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'yTicksAngle'))

        # The widget to change X-axis range.
        self.widgets['xRange'] = ComboStringEntry(self.gFrame,frameLabel = 'X range', contents = self.settings['xRange'])
        self.widgets['xRange'].grid(row = 10, column = 0, sticky = "NEW")
        self.widgets['xRange'].grid_columnconfigure(0, weight = 1)        
        self.widgets['xRange'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'xRange'))

        # The widget to change Y-axis range.
        self.widgets['yRange'] = ComboStringEntry(self.gFrame,frameLabel = 'Y range', contents = self.settings['yRange'])
        self.widgets['yRange'].grid(row = 11, column = 0, sticky = "NEW")
        self.widgets['yRange'].grid_columnconfigure(0, weight = 1)        
        self.widgets['yRange'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'yRange'))

        if len(self.parent.figure.gca().get_images()) == 0:
            v = axisScales.index(self.settings['xScale'])
            # The widget to change X-axis scale.
            self.widgets['xScale'] = ComboRadiobutton(self.gFrame,\
                                                      frameLabel = 'X scale',\
                                                      contents = axisScales,\
                                                      default = v,\
                                                      layout = (1,2))
            self.widgets['xScale'].grid(row = 12, column = 0, sticky = "NEW")
            self.widgets['xScale'].grid_columnconfigure(0, weight = 1)
            for r in self.widgets['xScale'].radio:
                r.bind('<ButtonRelease-1>', lambda event : self.changeSettings(event, 'xScale'))

            v = axisScales.index(self.settings['yScale'])
            # The widget to change Y-axis scale.
            self.widgets['yScale'] = ComboRadiobutton(self.gFrame,\
                                                      frameLabel = 'Y scale',\
                                                      contents = axisScales,\
                                                      default = v,\
                                                      layout = (1,2))
            self.widgets['yScale'].grid(row = 13, column = 0, sticky = "NEW")
            self.widgets['yScale'].grid_columnconfigure(0, weight = 1)        
            for r in self.widgets['yScale'].radio:
                r.bind('<ButtonRelease-1>', lambda event : self.changeSettings(event, 'yScale'))
                
        else:
            v = axisScales.index(self.settings['zScale'])
            # The widget to change Z-axis scale.
            self.widgets['zScale'] = ComboRadiobutton(self.gFrame,\
                                                      frameLabel = 'Z scale',\
                                                      contents = axisScales,\
                                                      default = v,\
                                                      layout = (1,2))
            self.widgets['zScale'].grid(row = 12, column = 0, sticky = "NEW")
            self.widgets['zScale'].grid_columnconfigure(0, weight = 1)        
            for r in self.widgets['zScale'].radio:
                r.bind('<ButtonRelease-1>', lambda event : self.changeSettings(event, 'zScale'))

        # The widget to change grid line width.
        self.widgets['gridWidth'] = ComboFloatEntry(self.gFrame, frameLabel = 'Grid width', contents = self.settings['gridWidth'])
        self.widgets['gridWidth'].grid(row = 15, column = 0, sticky = "NEW")
        self.widgets['gridWidth'].grid_columnconfigure(0, weight = 1)
        self.widgets['gridWidth'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'gridWidth'))

        # The widget to change grid style.
        self.widgets['gridStyle'] = ComboSpinbox(self.gFrame, frameLabel = 'Grid style', default = self.settings['gridStyle'], contents = tuple(lineStyles))
        self.widgets['gridStyle'].grid(row = 16, column = 0, sticky = "NEW")
        self.widgets['gridStyle'].grid_columnconfigure(0, weight = 1)
        self.widgets['gridStyle'].spinbox.config({'width' : 10, 'command' : lambda : self.changeSettings(None, 'gridStyle')})

        # The button to change grid color.
        self.widgets['gridColor'] = ComboButton(self.gFrame,\
                                                frameLabel = 'Grid color')
        self.widgets['gridColor'].grid(row = 17, column = 0, sticky = "NEW")
        self.widgets['gridColor'].grid_columnconfigure(0, weight = 1)
        self.widgets['gridColor'].button.config({'bg' : self.settings['gridColor'],\
                                                 'activebackground' : self.settings['gridColor'],\
                                                 'width' : 1,\
                                                 'command' : lambda : self.changeSettings(None, 'gridColor')})
        self.widgets['gridColor'].results = self.settings['gridColor']

        self.pFrame = LabelFrame(settingsFrame, text = 'Plot', bd = 2, relief = GROOVE)
        self.pFrame.grid(row = 0, column = 1, sticky = "NEW", padx = 3, pady = 3)
        self.pFrame.grid_columnconfigure(1, weight = 1)
        self.pFrame.grid_rowconfigure(0, weight = 1)

        # The listbox that will contain the plot list.
        self.plots = ComboListbox(self.pFrame, frameLabel = 'Available plots', contents = [])
        self.plots.grid(row = 0, column = 0, sticky = "NEW")
        self.plots.grid_columnconfigure(0, weight = 1)
        self.plots.grid_rowconfigure(0, weight = 1)
        self.plots.lb.config({'exportselection' : 0, 'width' : 28, 'height' : 4, 'selectmode' : SINGLE})
        self.plots.lb.bind('<ButtonRelease-1>', self.addPlotSettingsWidgets)

        if self.parent.figure.gca().get_images():
            self.plots.lb.insert(END, self.parent.figure.gca().get_images()[0].get_label())

        else:                
            for f in self.parent.figure.gca().get_lines():
                self.plots.lb.insert(END, f.get_label())

        self.plots.lb.selection_set(0)

        self.addPlotSettingsWidgets()

        return None

    def buttonbox(self):
        """
        Add standard button box.
        """

        # The frame that contains the 'Cancel' and 'OK' buttons.
        box = LabelFrame(self, text = 'Actions', bd = 2, relief = GROOVE)
        box.grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3)
        box.grid_columnconfigure(0, weight = 1)

        w = Button(box, text = "Cancel", width = 10, command = self.cancel)
        w.grid(row = 0, column = 0, sticky = E)

        w = Button(box, text = "OK", width = 10, command = self.ok, default = ACTIVE)
        w.grid(row = 0, column = 4, sticky = E)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

    # Standard button semantics.
    def ok(self, event = None):

        if not self.validate():
            self.initial_focus.focus_set()
            return

#        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):

        # Put focus back to the parent window
        self.parent.grab_set()
        self.parent.focus_set()
        self.destroy()

    # Command hooks
    def validate(self):

        try:
            self.storeSettings()

        except:
            LogMessage('warning','Bad input. Please try again.',['gui'])
            return False

        return True

    def apply(self):
        try:
            for k, v in self.settings.items():
                if k == 'plotTitle':
                    self.parent.figure.gca().set_title(v)
                elif k == 'xRange':
                    xMin, xMax = [float(vv) for vv in v.split('-')]
                    self.parent.figure.gca().set_xlim((xMin, xMax))
                elif k == 'yRange':
                    yMin, yMax = [float(vv) for vv in v.split('-')]
                    self.parent.figure.gca().set_ylim((yMin,yMax))
                elif k == 'xLabel':
                    self.parent.figure.gca().set_xlabel(v)
                elif k == 'yLabel':
                    self.parent.figure.gca().set_ylabel(v)
                elif k == 'xScale':
                    self.parent.figure.gca().set_xscale(v)
                elif k == 'yScale':
                    self.parent.figure.gca().set_yscale(v)
                elif k == 'zScale':
                    if v == 'linear':
                        self.parent.figure.gca().get_images()[self.selectedPlot].set_norm(Normalize())
                    else:
                        self.parent.figure.gca().get_images()[self.selectedPlot].set_norm(LogNorm())
                elif k == 'xTicks':
                    ticks, ticksLabels = [eval(t) for t in v.split(';')]
                    ticks = [float(t) for t in ticks]
                    ticksLabels = [str(t) for t in ticksLabels]

                    self.parent.figure.gca().set_xticks(ticks)
                    self.parent.figure.gca().set_xticklabels(ticksLabels)
                elif k == 'xTicksSize':
                    [l.set_size(int(v)) for l in self.parent.figure.gca().get_xticklabels()]
                elif k == 'xTicksAngle':
                    [l.set_rotation(int(v)) for l in self.parent.figure.gca().get_xticklabels()]
                elif k == 'yTicks':
                    ticks, ticksLabels = [eval(t) for t in v.split(';')]
                    ticks = [float(t) for t in ticks]
                    ticksLabels = [str(t) for t in ticksLabels]

                    self.parent.figure.gca().set_yticks(ticks)
                    self.parent.figure.gca().set_yticklabels(ticksLabels)
                elif k == 'yTicksSize':
                    [l.set_size(int(v)) for l in self.parent.figure.gca().get_yticklabels()]
                elif k == 'yTicksAngle':
                    [l.set_rotation(int(v)) for l in self.parent.figure.gca().get_yticklabels()]
                elif k == 'gridWidth':
                    self.parent.figure.gca().grid(True, linewidth = v)
                elif k == 'gridStyle':
                    self.parent.figure.gca().grid(True, linestyle = v)
                elif k == 'gridColor':
                    self.parent.figure.gca().grid(True, color = v)
                elif isinstance(k, int):
                    for kk, vv in v.items():
                        if kk == 'lineLabel':
                            self.parent.figure.gca().get_lines()[k].set_label(vv)
                            self.plots.lb.delete(k)
                            self.plots.lb.insert(k, vv)
                        elif kk == 'lineWidth':
                            self.parent.figure.gca().get_lines()[k].set_linewidth(vv)
                        elif kk == 'lineStyle':
                            self.parent.figure.gca().get_lines()[k].set_linestyle(vv)
                        elif kk == 'lineColor':
                            self.parent.figure.gca().get_lines()[k].set_color(vv)
                        elif kk == 'markerSize':
                            self.parent.figure.gca().get_lines()[k].set_markersize(vv)
                        elif kk == 'markerStyle':
                            self.parent.figure.gca().get_lines()[k].set_marker(vv)
                        elif kk == 'markerColor':
                            self.parent.figure.gca().get_lines()[k].set_markeredgecolor(vv)
                        elif kk == 'colorMap':
                            colorMap = eval('pylab.cm.' + vv)
                            self.parent.figure.gca().get_images()[k].set_cmap(colorMap)
                        elif kk == 'interpolation':
                            self.parent.figure.gca().get_images()[k].set_interpolation(vv)
                        elif kk == 'alpha':
                            self.parent.figure.gca().get_images()[k].set_alpha(vv)

        except:
            raise Error('Error when updating plots.')
                
        else:
            self.parent.figure.gca().legend(loc = 'best')
            self.parent.canvas.show()

    def initGlobalSettings(self):

        self.settings['plotTitle'] = self.parent.figure.gca().title.get_text()
        self.settings['xRange'] = '-'.join([str(v) for v in self.parent.figure.gca().get_xlim()])
        self.settings['yRange'] = '-'.join([str(v) for v in self.parent.figure.gca().get_ylim()])
        self.settings['xLabel'] = self.parent.figure.gca().xaxis.get_label().get_text()
        self.settings['yLabel'] = self.parent.figure.gca().yaxis.get_label().get_text()

        # The selected plot is a 2D plot.
        if len(self.parent.figure.gca().get_images()) == 0:
            self.settings['xScale'] = self.parent.figure.gca().get_xscale()
            self.settings['yScale'] = self.parent.figure.gca().get_yscale()
        else:
            self.settings['xScale'] = 'linear'
            self.settings['yScale'] = 'linear'
            self.settings['zScale'] = 'linear'

        self.settings['xTicks'] = str(list(self.parent.figure.gca().get_xticks()))
        self.settings['xTicks'] += ';' + str([v.get_text() for v in self.parent.figure.gca().get_xticklabels()])
        self.settings['xTicksSize'] = int(self.parent.figure.gca().get_xticklabels()[0].get_size())
        self.settings['xTicksAngle'] = int(self.parent.figure.gca().get_xticklabels()[0].get_rotation())
        self.settings['yTicks'] = str(list(self.parent.figure.gca().get_yticks()))
        self.settings['yTicks'] += ';' + str([v.get_text() for v in self.parent.figure.gca().get_yticklabels()])
        self.settings['yTicksSize'] = int(self.parent.figure.gca().get_yticklabels()[0].get_size())
        self.settings['yTicksAngle'] = int(self.parent.figure.gca().get_yticklabels()[0].get_rotation())
        self.settings['gridWidth'] = float(self.parent.figure.gca().get_xgridlines()[0].get_linewidth())
        self.settings['gridStyle'] = lineStyles.index(self.parent.figure.gca().get_xgridlines()[0].get_linestyle())
        self.settings['gridColor'] = rgb2hex(ColorConverter().to_rgb(self.parent.figure.gca().get_xgridlines()[0].get_color()))

    def initPlotSettings(self):

        if self.selectedPlot is None:
            return

        else:
            self.settings[self.selectedPlot] = {}

        # The selected plot is a 2D plot.
        if len(self.parent.figure.gca().get_images()) == 0:
            # The instance of the plot to set up.
            plot = self.parent.figure.get_axes()[0].get_lines()[self.selectedPlot]

            self.settings[self.selectedPlot]['lineLabel'] = plot.get_label()
            self.settings[self.selectedPlot]['lineWidth'] = int(float(plot.get_linewidth()))

            self.settings[self.selectedPlot]['lineStyle'] = plot.get_linestyle()
            self.settings[self.selectedPlot]['lineColor'] = rgb2hex(ColorConverter().to_rgb(plot.get_color()))

            self.settings[self.selectedPlot]['markerStyle'] = plot.get_marker()
            self.settings[self.selectedPlot]['markerSize'] = int(plot.get_markersize())
            self.settings[self.selectedPlot]['markerColor'] = rgb2hex(ColorConverter().to_rgb(plot.get_markeredgecolor()))

        else:
            # The instance of the plot to set up.
            plot = self.parent.figure.gca().get_images()[self.selectedPlot]

            self.settings[self.selectedPlot]['colorMap'] = plot.cmap.name
            self.settings[self.selectedPlot]['interpolation'] = plot._interpolation
            self.settings[self.selectedPlot]['alpha'] = plot.get_alpha()

    def storeSettings(self):

        self.settings['plotTitle'] = self.widgets['plotTitle'].getValue()
        self.settings['xRange'] = self.widgets['xRange'].getValue()
        self.settings['yRange'] = self.widgets['yRange'].getValue()
        self.settings['xLabel'] = self.widgets['xLabel'].getValue()
        self.settings['yLabel'] = self.widgets['yLabel'].getValue()

        self.settings['xTicks'] = self.widgets['xTicks'].getValue()
        self.settings['xTicksSize'] = self.widgets['xTicksSize'].getValue()
        self.settings['xTicksAngle'] = self.widgets['xTicksAngle'].getValue()
        self.settings['yTicks'] = self.widgets['yTicks'].getValue()
        self.settings['yTicksSize'] = self.widgets['yTicksSize'].getValue()
        self.settings['yTicksAngle'] = self.widgets['yTicksAngle'].getValue()

        self.settings['gridWidth'] = self.widgets['gridWidth'].getValue()
        self.settings['gridStyle'] = self.widgets['gridStyle'].getValue()
        self.settings['gridColor'] = self.widgets['gridColor'].getValue()

        # The selected plot is a 2D plot.
        if len(self.parent.figure.gca().get_images()) == 0:
            self.settings['xScale'] = self.widgets['xScale'].getValue()
            self.settings['yScale'] = self.widgets['yScale'].getValue()

            self.settings[self.selectedPlot]['lineLabel'] = self.widgets['lineLabel'].getValue()
            self.settings[self.selectedPlot]['lineWidth'] = self.widgets['lineWidth'].getValue()
            self.settings[self.selectedPlot]['lineStyle'] = self.widgets['lineStyle'].getValue()
            self.settings[self.selectedPlot]['lineColor'] = self.widgets['lineColor'].getValue()

            self.settings[self.selectedPlot]['markerStyle'] = self.widgets['markerStyle'].getValue()
            self.settings[self.selectedPlot]['markerSize'] = self.widgets['markerSize'].getValue()
            self.settings[self.selectedPlot]['markerColor'] = self.widgets['markerColor'].getValue()

        else:
            self.settings['xScale'] = 'linear'
            self.settings['yScale'] = 'linear'
            self.settings['zScale'] = self.widgets['zScale'].getValue()

            self.settings[self.selectedPlot]['colorMap'] = self.widgets['colorMap'].getValue()
            self.settings[self.selectedPlot]['interpolation'] = self.widgets['interpolation'].getValue()
            self.settings[self.selectedPlot]['alpha'] = self.widgets['alpha'].getValue()

    def selectColor(self, widget):

        defaultColor = widget.button.cget('bg')

        junk, selectedColor = tkColorChooser.askcolor(color = defaultColor, master = self)
        if isinstance(selectedColor, _tkinter.Tcl_Obj):
            selectedColor = selectedColor.string

        elif isinstance(selectedColor, str):
            pass

        else:
            selectedColor = defaultColor

        widget.button.configure(bg = selectedColor, activebackground = selectedColor)
        widget.results = selectedColor

        return selectedColor

    def changeSettings(self, event, widgetName):
        """
        Argument:
            - event: either a Tkinter event, either a Tkinter control variable value that has been traced for changes.
        """

        if not self.widgets.has_key(widgetName):
            return

        w = self.widgets[widgetName]
        if widgetName == 'plotTitle':
            self.parent.figure.gca().set_title(w.getValue())
        elif widgetName == 'xRange':
            xMin, xMax = [float(v) for v in w.getValue().split('-')]
            self.parent.figure.gca().set_xlim((xMin,xMax))
        elif widgetName == 'yRange':
            yMin, yMax = [float(v) for v in w.getValue().split('-')]
            self.parent.figure.gca().set_ylim((yMin,yMax))
        elif widgetName == 'xLabel':
            self.parent.figure.gca().set_xlabel(w.getValue())
        elif widgetName == 'yLabel':
            self.parent.figure.gca().set_ylabel(w.getValue())
        elif widgetName == 'xScale':
            self.parent.figure.gca().set_xscale(w.getValue())
        elif widgetName == 'yScale':
            self.parent.figure.gca().set_yscale(w.getValue())
        elif widgetName == 'zScale':
            if w.getValue() == 'linear':
                self.parent.figure.gca().get_images()[self.selectedPlot].set_norm(Normalize())
            else:
                self.parent.figure.gca().get_images()[self.selectedPlot].set_norm(LogNorm())
        elif widgetName == 'gridWidth':
            self.parent.figure.gca().grid(True, linewidth = w.getValue())
        elif widgetName == 'gridStyle':
            self.parent.figure.gca().grid(True, linestyle = w.getValue())
        elif widgetName == 'gridColor':
            color = self.selectColor(w)
            self.parent.figure.gca().grid(True, color = color)
        elif widgetName == 'xTicks':
            ticks, ticksLabels = [eval(t) for t in w.getValue().split(';')]
            ticks = [float(t) for t in ticks]
            ticksLabels = [str(t) for t in ticksLabels]
            self.parent.figure.gca().set_xticks(ticks)
            self.parent.figure.gca().set_xticklabels(ticksLabels)
        elif widgetName == 'xTicksSize':
            [l.set_size(int(w.getValue())) for l in self.parent.figure.gca().get_xticklabels()]
        elif widgetName == 'xTicksAngle':
            [l.set_rotation(int(w.getValue())) for l in self.parent.figure.gca().get_xticklabels()]
        elif widgetName == 'yTicks':
            ticks, ticksLabels = [eval(t) for t in w.getValue().split(';')]
            ticks = [float(t) for t in ticks]
            ticksLabels = [str(t) for t in ticksLabels]
            self.parent.figure.gca().set_yticks(ticks)
            self.parent.figure.gca().set_yticklabels(ticksLabels)
        elif widgetName == 'yTicksSize':
            [l.set_size(int(w.getValue())) for l in self.parent.figure.gca().get_yticklabels()]
        elif widgetName == 'yTicksAngle':
            [l.set_rotation(int(w.getValue())) for l in self.parent.figure.gca().get_yticklabels()]
        elif widgetName == 'lineLabel':
            self.parent.figure.gca().get_lines()[self.selectedPlot].set_label(w.getValue())
            self.plots.lb.delete(self.selectedPlot)
            self.plots.lb.insert(self.selectedPlot, w.getValue())
        elif widgetName == 'lineWidth':
            self.parent.figure.gca().get_lines()[self.selectedPlot].set_linewidth(w.getValue())
        elif widgetName == 'lineStyle':            
            self.parent.figure.gca().get_lines()[self.selectedPlot].set_linestyle(w.getValue())
        elif widgetName == 'lineColor':
            color = self.selectColor(w)
            self.parent.figure.gca().get_lines()[self.selectedPlot].set_color(color)
        elif widgetName == 'markerSize':
            self.parent.figure.gca().get_lines()[self.selectedPlot].set_markersize(w.getValue())
        elif widgetName == 'markerStyle':
            self.parent.figure.gca().get_lines()[self.selectedPlot].set_marker(w.getValue())
        elif widgetName == 'markerColor':
            color = self.selectColor(w)
            self.parent.figure.gca().get_lines()[self.selectedPlot].set_markeredgecolor(color)
        elif widgetName == 'alpha':
            self.parent.figure.gca().get_images()[self.selectedPlot].set_alpha(w.getValue())
        elif widgetName == 'colorMap':
            colorMap = eval('pylab.cm.' + w.getValue())
            self.parent.figure.gca().get_images()[self.selectedPlot].set_cmap(colorMap)
        elif widgetName == 'interpolation':        
            self.parent.figure.gca().get_images()[self.selectedPlot].set_interpolation(w.getValue())

        self.parent.figure.gca().legend(loc = 'best')
        self.parent.canvas.show()

    def removePreviousPlotSettings(self):
        """
        This method removes the previous plot settings widgets.
        """

        # Loop over all the plot settings widgets.
        for widget in self.pFrame.winfo_children():
            # Keep the plot settings list box
            if widget == self.plots:
                continue

            # Destroys all the other ones.
            else:
                widget.destroy()

    def addPlotSettingsWidgets(self, event = None):

        if event is not None:
            if not event.widget.curselection():
                return

        self.selectedPlot = int(self.plots.lb.curselection()[0])

        self.removePreviousPlotSettings()

        self.initPlotSettings()

        # The selected plot is a 2D plot.
        if len(self.parent.figure.gca().get_images()) == 0:

            self.widgets['lineLabel'] = ComboStringEntry(self.pFrame, frameLabel = 'Plot label', contents = self.settings[self.selectedPlot]['lineLabel'])
            self.widgets['lineLabel'].grid(row = 1, column = 0, sticky = "NEW")
            self.widgets['lineLabel'].grid_columnconfigure(0, weight = 1)        
            self.widgets['lineLabel'].grid_rowconfigure(1, weight = 1)
            self.widgets['lineLabel'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'lineLabel'))

            self.widgets['lineWidth'] = ComboIntegerEntry(self.pFrame, frameLabel = 'Line width', contents = self.settings[self.selectedPlot]['lineWidth'])
            self.widgets['lineWidth'].grid(row = 2, column = 0, sticky = "NEW")
            self.widgets['lineWidth'].grid_columnconfigure(0, weight = 1)
            self.widgets['lineWidth'].grid_rowconfigure(2, weight = 1)
            self.widgets['lineWidth'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'lineWidth'))

            v = lineStyles.index(self.settings[self.selectedPlot]['lineStyle'])
            self.widgets['lineStyle'] = ComboSpinbox(self.pFrame,\
                                                     frameLabel = 'Line style',\
                                                     default = v,\
                                                     contents = tuple(lineStyles))
            self.widgets['lineStyle'].grid(row = 3, column = 0, sticky = "NEW")
            self.widgets['lineStyle'].grid_columnconfigure(0, weight = 1)
            self.widgets['lineStyle'].grid_rowconfigure(3, weight = 1)
            self.widgets['lineStyle'].spinbox.config({'width' : 10, 'command' : lambda : self.changeSettings(None, 'lineStyle')})

            self.widgets['lineColor'] = ComboButton(self.pFrame, frameLabel = 'Line color')
            self.widgets['lineColor'].grid(row = 4, column = 0, sticky = "NEW")
            self.widgets['lineColor'].grid_columnconfigure(0, weight = 1)
            self.widgets['lineColor'].grid_rowconfigure(4, weight = 1)
            self.widgets['lineColor'].button.config({'bg' : self.settings[self.selectedPlot]['lineColor'],\
                                                     'activebackground' : self.settings[self.selectedPlot]['lineColor'],\
                                                     'width' : 1,\
                                                     'command' : lambda : self.changeSettings(None, 'lineColor')})
            self.widgets['lineColor'].results = self.settings[self.selectedPlot]['lineColor']

            v = markerStyles.index(self.settings[self.selectedPlot]['markerStyle'])
            self.widgets['markerStyle'] = ComboSpinbox(self.pFrame, frameLabel = 'Marker style',\
                                                       default = v,\
                                                       contents = tuple(markerStyles))
            self.widgets['markerStyle'].grid(row = 5, column = 0, sticky = "NEW")
            self.widgets['markerStyle'].grid_columnconfigure(0, weight = 1)
            self.widgets['markerStyle'].grid_rowconfigure(5, weight = 1)
            self.widgets['markerStyle'].spinbox.config({'width' : 10, 'command' : lambda : self.changeSettings(None, 'markerStyle')})

            self.widgets['markerSize'] = ComboIntegerEntry(self.pFrame, frameLabel = 'Marker size', contents = self.settings[self.selectedPlot]['markerSize'])
            self.widgets['markerSize'].grid(row = 6, column = 0, sticky = "NEW")
            self.widgets['markerSize'].grid_columnconfigure(0, weight = 1)
            self.widgets['markerSize'].grid_rowconfigure(6, weight = 1)
            self.widgets['markerSize'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'markerSize'))

            self.widgets['markerColor'] = ComboButton(self.pFrame, frameLabel = 'Marker color')
            self.widgets['markerColor'].grid(row = 7, column = 0, sticky = "NEW")
            self.widgets['markerColor'].grid_columnconfigure(0, weight = 1)
            self.widgets['markerColor'].grid_rowconfigure(7, weight = 1)
            self.widgets['markerColor'].button.config({'bg' : self.settings[self.selectedPlot]['markerColor'],\
                                                       'activebackground' : self.settings[self.selectedPlot]['markerColor'],\
                                                       'width' : 1,\
                                                       'command' : lambda event : self.changeSettings(event, 'markerColor')})
            self.widgets['markerColor'].results = self.settings[self.selectedPlot]['markerColor']

        # The selected plot is a 3D plot.
        else:

            v = colorMaps.index(self.settings[self.selectedPlot]['colorMap'])
            self.widgets['colorMap'] = ComboSpinbox(self.pFrame,\
                                                    frameLabel = 'Color map',\
                                                    default = v,\
                                                    contents = tuple(colorMaps))
            self.widgets['colorMap'].grid(row = 1, column = 0, sticky = "NEW")
            self.widgets['colorMap'].grid_columnconfigure(0, weight = 1)
            self.widgets['colorMap'].grid_rowconfigure(5, weight = 1)
            self.widgets['colorMap'].spinbox.config({'width' : 10,\
                                                     'command' : lambda : self.changeSettings(None, 'colorMap')})

            v = interpolations.index(self.settings[self.selectedPlot]['interpolation'])
            self.widgets['interpolation'] = ComboSpinbox(self.pFrame,\
                                                         frameLabel = 'Interpolation',\
                                                         default = v,\
                                                         contents = tuple(interpolations))
            self.widgets['interpolation'].grid(row = 2, column = 0, sticky = "NEW")
            self.widgets['interpolation'].grid_columnconfigure(0, weight = 1)
            self.widgets['interpolation'].grid_rowconfigure(5, weight = 1)
            self.widgets['interpolation'].spinbox.config({'width' : 10,\
                                                          'command' : lambda : self.changeSettings(None, 'interpolation')})

            self.widgets['alpha'] = ComboFloatEntry(self.pFrame, frameLabel = 'Alpha', contents = self.settings[self.selectedPlot]['alpha'])
            self.widgets['alpha'].grid(row = 3, column = 0, sticky = "NEW")
            self.widgets['alpha'].grid_columnconfigure(0, weight = 1)
            self.widgets['alpha'].grid_rowconfigure(0, weight = 1)
            self.widgets['alpha'].entry.bind('<FocusOut>', lambda event : self.changeSettings(event, 'alpha'))

class PlotNetCDFVariableDialog(PortableToplevel):
    """Sets up a dialog used to plot variables present in a NetCDF file.
    """

    def __init__(self, parent, title = None, netcdf = None, xVar = None, yVar = None, zVar = None):
        """The constructor.

        @param parent: the parent widget.

        @param title: a string specifying the title of the dialog.
        @type title: string

        @param netcdf: the name of a NetCDF file to plot (string) or an opened NetCDF trajectory file.
        @type netcdf: a string or a Scientific.IO.NetCDF._NetCDFFile object

        @param xVar: the NetCDF variable name of the X variable to plot.
        @type xVar: string

        @param yVar: the NetCDF variable name of the Y variable to plot.
        @type yVar:

        @param zVar: the NetCDF variable name of the Z variable to plot.
        @type zVar:
        """

        PortableToplevel.__init__(self, parent)

        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent       

        try:
            if isinstance(netcdf, _NetCDFFile):
                self.data = netcdf

            elif isinstance(netcdf, Trajectory):
                self.data = NetCDFFile(netcdf.filename, 'r')

            elif isinstance(netcdf, str):
                self.data = NetCDFFile(netcdf, 'r')

            else:
                raise

        except:
            self.data = None

        self.xVar = xVar
        self.yVar = yVar
        self.zVar = zVar

        body = Frame(self)
        self.initial_focus = self.body(body)
        body.grid(row = 0, column = 0, sticky = EW)

        self.buttonbox()        

        # If a NetCDF file has been loaded in nMOLDYN then it is proposed by default for plotting.
        if self.data is not None:
            # The control variables are updated with the informations coming from the loaded trajectory.
            self.fileBrowser.setValue(str(self.data))

            # The listbox storing the numerical 1D and 2D NetCDF data is filled up.
            self.displayVariables()

            if (self.xVar is not None) & (self.yVar is not None):
                if self.zVar is not None:
                    self.plotXYZ()
                else:
                    self.plotXY()

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", lambda : self.cancel(self))

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
                                            frameLabel = "NetCDF input file",\
                                            contents = '',\
                                            save = False,\
                                            command = self.openNetCDF)
        self.fileBrowser.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = EW)
        self.fileBrowser.grid_columnconfigure(0, weight = 1)
        self.fileBrowser.entry.bind('<Return>', self.openNetCDF)

        # The frame that will host the different listboxes of the dialog.
        listboxesFrame =  LabelFrame(settingsFrame, text = 'Variables', bd = 2, relief = GROOVE)
        listboxesFrame.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = W)

        # The combo listbox that will contain the X variables.
        self.xVarLb = ComboListbox(listboxesFrame, 'X', [])
        self.xVarLb.lb.config({'exportselection' : 0, 'width' : 22, 'height' : 4, 'selectmode' : SINGLE})
        self.xVarLb.config({'relief' : FLAT})
        self.xVarLb.grid(row = 0, column = 0)
        self.xVarLb.lb.bind('<ButtonRelease-1>', self.selectXVariable)

        # The combo listbox that will contain the Y variables.
        self.yVarLb = ComboListbox(listboxesFrame, 'Y', [])
        self.yVarLb.config({'relief' : FLAT})
        self.yVarLb.grid(row = 0, column = 1)
        self.yVarLb.lb.config({'exportselection' : 0, 'width' : 22, 'height' : 4, 'selectmode' : SINGLE})

        # The combo listbox that will contain the Z variables.
        self.zVarLb = ComboListbox(listboxesFrame, 'Z', [])
        self.zVarLb.config({'relief' : FLAT})
        self.zVarLb.grid(row = 0, column = 2)
        self.zVarLb.lb.config({'exportselection' : 0, 'width' : 22, 'height' : 4, 'selectmode' : SINGLE})

        # The frame that will host the pylab figure displaying the plot and the global settings widgets.
        plotFrame = LabelFrame(master, text = 'Plot', bd = 2, relief = GROOVE)
        plotFrame.grid(row = 3, column = 0, padx = 2, pady = 2, sticky = EW)
        plotFrame.grid_columnconfigure(0, weight = 1)

        # The subframe that will host the figure.
        figureFrame = Frame(plotFrame, relief = FLAT)
        figureFrame.grid(row = 0, column = 0, sticky = 'EW')
        figureFrame.grid_columnconfigure(0, weight = 1)

        # The figure.
        self.figure = Figure(dpi = 72)
        self.canvas = FigureCanvasTkAgg(self.figure, master = figureFrame)

        self.canvas.get_tk_widget().pack()
        self.canvas._tkcanvas.pack()
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, figureFrame)
        self.toolbar.update()

        return None

    def buttonbox(self):
        """
        Add standard button box.
        """

        # The frame that contains the action buttons.
        box = LabelFrame(self, text = 'Actions', bd = 2, relief = GROOVE)
        box.grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3)
        box.grid_columnconfigure(0, weight = 1)

        w = Button(box, text = "Cancel", width = 10, command = lambda : self.cancel(self))
        w.grid(row = 0, column = 0, sticky = E)

        w = Button(box, text = "Settings", width = 10, command = self.openSettingsDialog)
        w.grid(row = 0, column = 1, sticky = E)

        w = Button(box, text = "Reset", width = 10, command = self.resetPlots)
        w.grid(row = 0, column = 3, sticky = E)

        w = Button(box, text = "Export plot", width = 10, command = self.exportPlotDialog)
        w.grid(row = 0, column = 4, sticky = E)

        w = Button(box, text = "OK", width = 10, command = self.ok, default = ACTIVE)
        w.grid(row = 0, column = 5, sticky = E)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", lambda event : self.cancel(self, event))

    # Standard button semantics.
    def ok(self, event = None):
        
        if not self.validate():
            self.initial_focus.focus_set()
            return

        self.update_idletasks()

        self.apply()

    def cancel(self, dialog, event = None):
        
        dialog.parent.grab_set()
        dialog.parent.focus_set()
        dialog.destroy()

    # Command hooks
    def validate(self):

        if not self.xVarLb.lb.curselection():
            raise Error('Please select a X variable.')
        self.xVar = self.xVarLb.lb.get(int(self.xVarLb.lb.curselection()[0]))
        XVarName = self.data.variables[self.xVar].dimensions[0]
        dimX = self.data.dimensions[XVarName]

        if not self.yVarLb.lb.curselection():
            raise Error('Please select a Y variable.')
        self.yVar = self.yVarLb.lb.get(int(self.yVarLb.lb.curselection()[0]))
        YVarName = self.data.variables[self.yVar].dimensions[0]
        dimY = self.data.dimensions[YVarName]

        if not self.zVarLb.lb.curselection():
            self.zVar = None

            # The variable must have the same dimension than the selected Y variable to be retained.
            if dimX != dimY:
                raise Error('X (%d) and Y (%d) do not have the same dimension.' % (dimX,dimY))

        else:
            self.zVar = self.zVarLb.lb.get(int(self.zVarLb.lb.curselection()[0]))
            ZVarName = self.data.variables[self.zVar].dimensions
            dimZ1 = self.data.dimensions[ZVarName[0]]
            dimZ2 = self.data.dimensions[ZVarName[1]]
            if (dimX != dimZ1) or (dimY != dimZ2):
                raise Error('X (%d) and Y (%d) dimensions do not match Z (%d,%d) dimensions.' % (dimX,dimY,dimZ1,dimZ2))

        return True

    def apply(self):
        
        try:
            if self.zVar is None:
                self.plotXY()

            else:
                self.plotXYZ()
        except:
            raise Error('Error when plotting datas.')

    # The class specific methods.
    def openSettingsDialog(self):
        """This method will open the dialog to set up the global settings.
        """

        if not self.figure.gca().get_lines():
            LogMessage('warning', 'No plot defined.', ['gui'])
            return

        g = SettingsDialog(self)

    def resetPlots(self):
        """This method will clear up all the displayed plots.
        """

        # The figure is cleared up.
        self.figure.clear()

        if hasattr(self.figure,'xValues'):
            delattr(self.figure, 'xValues')

        if hasattr(self.figure,'yValues'):
            delattr(self.figure, 'yValues')

        if hasattr(self.figure,'zValues'):
            delattr(self.figure, 'zValues')

        # The canvas is updated.
        self.canvas.show()

    def exportPlotDialog(self):
        """This method pops up a dialog from which the plotted datas can be exported to an ASCII file.
        """

        if not self.figure.gca().get_lines():
            LogMessage('warning', 'No plot defined.', ['gui'])
            return

        self.exportDialog = PortableToplevel(self)
        
        self.exportDialog.parent = self

        self.exportDialog.transient(self)

        self.exportDialog.title('Export plot to ASCII')

        self.exportDialog.grab_set()

        # The frame that contains the action buttons.
        settingsFrame = LabelFrame(self.exportDialog, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)

        # Windows not resizable
        self.exportDialog.resizable(width = NO, height = NO)

        self.exportDialog.geometry("+%d+%d" % (self.winfo_rootx()+50, self.winfo_rooty()+50))

        if self.figure.gca().get_images():
            self.xIndexMinEntry = ComboIntegerEntry(settingsFrame, frameLabel = 'X Index min', contents = 1)
            self.xIndexMinEntry.grid(row = 0, column = 0, sticky = "EW")
            self.xIndexMinEntry.grid_columnconfigure(0, weight = 1)        

            self.xIndexMaxEntry = ComboIntegerEntry(settingsFrame, frameLabel = 'X Index max',\
                                                    contents = self.figure.gca().get_images()[0].get_array().shape[1])
            self.xIndexMaxEntry.grid(row = 1, column = 0, sticky = "EW")
            self.xIndexMaxEntry.grid_columnconfigure(0, weight = 1)        
            
            self.yIndexMinEntry = ComboIntegerEntry(settingsFrame, frameLabel = 'Y Index min', contents = 1)
            self.yIndexMinEntry.grid(row = 2, column = 0, sticky = "EW")
            self.yIndexMinEntry.grid_columnconfigure(0, weight = 1)        

            self.yIndexMaxEntry = ComboIntegerEntry(settingsFrame, frameLabel = 'Y Index max',\
                                                    contents = self.figure.gca().get_images()[0].get_array().shape[0])
            self.yIndexMaxEntry.grid(row = 3, column = 0, sticky = "EW")
            self.yIndexMaxEntry.grid_columnconfigure(0, weight = 1)        

        else:
            self.indexMinEntry = ComboIntegerEntry(settingsFrame, frameLabel = 'Index min', contents = 1)
            self.indexMinEntry.grid(row = 0, column = 0, sticky = "EW")
            self.indexMinEntry.grid_columnconfigure(0, weight = 1)        

            self.indexMaxEntry = ComboIntegerEntry(settingsFrame, frameLabel = 'Index max',\
                                                   contents = len(self.figure.gca().get_lines()[0].get_xdata()))
            self.indexMaxEntry.grid(row = 1, column = 0, sticky = "EW")
            self.indexMaxEntry.grid_columnconfigure(0, weight = 1)        

        # The filebrowser for the output ASCII file.
        self.exportedFileBrowser = ComboFileBrowser(settingsFrame,\
                                                    frameLabel = "ASCII output file",\
                                                    contents = '',\
                                                    save = True)
        self.exportedFileBrowser.grid(row = 4, column = 0, padx = 2, pady = 2, sticky = EW)
        self.exportedFileBrowser.grid_columnconfigure(0, weight = 1)
        
        # The frame that contains the action buttons.
        box = LabelFrame(self.exportDialog, text = 'Actions', bd = 2, relief = GROOVE)
        box.grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3)
        box.grid_columnconfigure(0, weight = 1)

        w = Button(box, text = "Cancel", width = 10, command = lambda : self.cancel(self.exportDialog))
        w.grid(row = 0, column = 0, sticky = E)

        w = Button(box, text = "OK", width = 10, command = self.exportPlot, default = ACTIVE)
        w.grid(row = 0, column = 1, sticky = E)

        self.exportDialog.protocol("WM_DELETE_WINDOW", lambda : self.cancel(self.exportDialog))
        self.exportDialog.bind("<Return>", self.exportPlot)
        self.exportDialog.bind("<Escape>", lambda event : self.cancel(self.exportDialog, event))
        
    def exportPlot(self, event = None):
        """This method exports plotted datas to an ASCII file.
        """

        filename = self.exportedFileBrowser.getValue()
        if not filename:
            LogMessage('warning', 'No output file defined.', ['gui'])
            return

        if filename[-4:] != '.dat':
            filename += '.dat'

        try:

            fFmt = '%20.8e '
            sFmt = '%20s '
            f = open(filename, 'w')

            if self.figure.gca().get_images():
                xIndexMin = self.xIndexMinEntry.getValue() - 1
                xIndexMax = self.xIndexMaxEntry.getValue()
                
                if xIndexMin < 0:
                    xIndexMin = 0
                if xIndexMax > self.figure.gca().get_images()[0].get_array().shape[1]:
                    xIndexMax = self.figure.gca().get_images()[0].get_array().shape[1]

                yIndexMin = self.yIndexMinEntry.getValue() - 1
                yIndexMax = self.yIndexMaxEntry.getValue()
                if yIndexMin < 0:
                    yIndexMin = 0
                yIndexMax = self.yIndexMaxEntry.getValue()
                if yIndexMax > self.figure.gca().get_images()[0].get_array().shape[0]:
                    yIndexMax = self.figure.gca().get_images()[0].get_array().shape[0]

                f.write('#datas = %s\n' % self.figure.gca().get_images()[0].get_label())
                f.write('#1st line = %s\n' % self.figure.gca().get_images()[0].get_axes().get_ylabel())
                f.write('#1st column = %s\n' % self.figure.gca().get_images()[0].get_axes().get_xlabel())
                                
                data = self.figure.gca().get_images()[0].get_array()

                # The first line is written with the Y values.
                f.write(sFmt % 'NaN')
                for v in range(xIndexMin,xIndexMax):
                    f.write(fFmt % self.figure.xValues[v])
                f.write('\n')                    

                # The 2D datas are written.
                for iy in range(yIndexMin,yIndexMax):
                    f.write(fFmt % self.figure.yValues[iy])
                    for ix in range(xIndexMin,xIndexMax):
                        v = self.figure.zValues[iy,ix]
                        try:
                            f.write(fFmt % v)
                        except:
                            f.write(sFmt % v)

                    f.write('\n')

            else:
                indexMin = self.indexMinEntry.getValue() - 1
                indexMax = self.indexMaxEntry.getValue()
                if indexMin < 0:
                    indexMin = 0
                if indexMax > len(self.figure.gca().get_lines()[0].get_xdata()):
                    indexMax = len(self.figure.gca().get_lines()[0].get_xdata())

                data = []
                comp = 1
                for plot in self.figure.gca().get_lines():
                    f.write('columns-%d = %s\n' % (comp,plot.get_axes().get_xlabel()))
                    f.write('columns-%d = %s\n' % (comp+1,plot.get_label()))
                    data.append(list(plot.get_xdata())[indexMin:indexMax])
                    data.append(list(plot.get_ydata())[indexMin:indexMax])
                    comp += 2
                data = [[r[col] for r in data] for col in range(len(data[0]))]

                f.write('\n')

                for d in data:
                    for dd in d:
                        try:
                            f.write(fFmt % dd)
                        except:
                            f.write(sFmt % dd)
                    f.write('\n')
            f.close()            
        except:
            raise Error('Error when exporting plotted data.')
        
        self.cancel(self.exportDialog)
            
    def selectXVariable(self, event):

        self.yVarLb.lb.selection_clear(0, END)
        self.zVarLb.lb.selection_clear(0, END)        

    def displayVariables(self):
        """
        This method display the numeric variables found in the NetCDF file into their appropriate listbox.
        """

        # The X, Y and Z variables listboxes are cleared up first.
        self.xVarLb.lb.delete(0, END)
        self.yVarLb.lb.delete(0, END)
        self.zVarLb.lb.delete(0, END)

        # Loop over the variables NetCDF dictionnary items.
        for varName in sorted(self.data.variables.keys()):
            variable = self.data.variables[varName]

            if len(variable.shape) > 2:
                continue

            # If the NetCDF variable is not numeric then loop to the next one.
            if variable.typecode() not in [Num.Float, Num.Float0, Num.Int, Num.Int32, Num.Character]:
                continue

            # If the variable is 1D then it is appended to the X and Y listboxes.
            if len(variable.shape) == 1:
                self.xVarLb.lb.insert(END, varName)
                self.yVarLb.lb.insert(END, varName)

            # If the variable is 2D then it is appended to the Z listbox.
            elif len(variable.shape) == 2:
                if variable.typecode() == Num.Character:
                    self.xVarLb.lb.insert(END, varName)
                    self.yVarLb.lb.insert(END, varName)

                else:
                    self.zVarLb.lb.insert(END, varName)                        

    def openNetCDF(self, event = None):
        """
        This method opens a NetCDF file and updates the dialog with the data read from that file.
        Arguments:
            -event: Tkinter event.
        """

        if event is not None:
            filename = self.fileBrowser.getValue()

        else:
            # The name of the NetCDF file to load.
            filename = askopenfilename(parent = self,\
                                       filetypes = [('NetCDF file','*.nc')],\
                                       initialdir = PREFERENCES['outputfile_path'])            

        # The file must exist.
        if filename:
            try:
                self.data = NetCDFFile(filename, 'r')

                self.fileBrowser.setValue(filename)

                # The X, Y and Z variables listboxes are filled up.
                self.displayVariables()

            except:
                LogMessage('warning','The file %s could not be read.' % filename,['gui'])
                
        return 'break'

    def plotXY(self):
        """
        This method display a 2D plot.
        """
        
        if (self.xVar is None) or (self.yVar is None):
            return

        if self.yVar is None:
            return

        if self.figure.gca().get_images():
            self.resetPlots()

        if self.data.variables[self.xVar].typecode() == Num.Character:
            temp = self.data.variables[self.xVar].getValue()
            xValues = Num.arange(temp.shape[0])
        else:
            xValues = self.data.variables[self.xVar].getValue()

        if self.data.variables[self.yVar].typecode() == Num.Character:
            temp = self.data.variables[self.yVar].getValue()
            yValues = Num.arange(temp.shape[0])
        else:
            yValues = self.data.variables[self.yVar].getValue()

        self.figure.gca().plot(xValues, yValues, label = self.yVar + ' vs ' + self.xVar)

        self.figure.gca().grid(True)

        if self.data.variables[self.xVar].typecode() == Num.Character:
            xTicks = xValues[:]
            xTicksLabels = ['%s' % temp[v,:].tostring().rstrip('\x00') for v in range(len(temp))]
            self.figure.gca().set_xticks(xTicks)
            self.figure.gca().set_xticklabels(xTicksLabels)

        if self.data.variables[self.yVar].typecode() == Num.Character:
            yTicks = yValues[:]
            yTicksLabels = ['%s' % temp[v,:].tostring().rstrip('\x00') for v in range(len(temp))]
            self.figure.gca().set_yticks(yTicks)
            self.figure.gca().set_yticklabels(yTicksLabels)

        # Set the title.
        try:
            self.figure.get_axes()[0].set_title(self.data.title)
        except AttributeError:
            self.figure.get_axes()[0].set_title('')

        # Sets the X-axis label.
        try:
            self.figure.get_axes()[0].set_xlabel('%s (%s)' % (self.xVar, self.data.variables[self.xVar].units))
        except AttributeError:
            self.figure.get_axes()[0].set_xlabel(self.xVar)

        # Sets the Y-axis label.
        try:
            self.figure.get_axes()[0].set_ylabel('%s (%s)' % (self.yVar, self.data.variables[self.yVar].units))
        except AttributeError:
            self.figure.get_axes()[0].set_ylabel(self.yVar)

        # Set the legend.
        self.figure.get_axes()[0].legend(loc = 'best')

        self.canvas.show()

    def displayPlotSlices(self, event):
        """
        This call back plot the orthogonal slices defined by the moving cursor of a 3D plot.
        """

        # The event triggering the call back must be the right mouse button.
        if event.button != 3:
            return

        # If the event does not occur into the plot then do nothing.
        if event.inaxes is None:
            return

        # The plot that will contain the horizontal slice.
        horizontalSlice = pylab.subplot(211)
        if self.data.variables[self.xVar].typecode() == Num.Character:
            temp = self.data.variables[self.xVar].getValue()
            xValue = Num.arange(len(temp))
            xTicks = ['%s' % temp[v,:].tostring().rstrip('\x00') for v in range(len(temp))]
            horizontalSlice.set_xticks(range(len(xTicks)))
            horizontalSlice.set_xticklabels(xTicks)
        else:
            xValue = self.data.variables[self.xVar].getValue()

        # The horizontal slice plot.    
        horizontalSlice.plot(xValue, self.data.variables[self.zVar][:,int(event.ydata)])

        # Sets the X-axis label.
        if hasattr(self.data.variables[self.xVar], 'units'):
            horizontalSlice.set_xlabel('%s (%s)' % (self.xVar, self.data.variables[self.xVar].units))
        else:
            horizontalSlice.set_xlabel('%s (Unknown unit)' % self.xVar)

        # Sets the Y-axis label.
        if self.data.variables[self.yVar].typecode() == Num.Character:
            value = self.data.variables[self.yVar].getValue()[int(event.ydata),:].tostring()
            horizontalSlice.set_ylabel('%s %s = %s' % (self.zVar,self.yVar, value), fontsize = 7)
        else:
            value = self.data.variables[self.yVar].getValue()[int(event.ydata)]
            horizontalSlice.set_ylabel('%s %s = %.3f' % (self.zVar,self.yVar, value), fontsize = 7)

        # The plot that will contain the vertical slice.
        verticalSlice = pylab.subplot(212)

        if self.data.variables[self.yVar].typecode() == Num.Character:
            temp = self.data.variables[self.yVar].getValue()
            xValue = Num.arange(len(temp))
            xTicks = ['%s' % temp[v,:].tostring().rstrip('\x00') for v in range(len(temp))]
            verticalSlice.set_xticks(range(len(xTicks)))
            verticalSlice.set_xticklabels(xTicks)
        else:
            xValue = self.data.variables[self.yVar].getValue()

        # The vertical slice plot.    
        verticalSlice.plot(xValue, self.data.variables[self.zVar][int(event.xdata),:])
        # The x label of the horizontal slice plot.
        verticalSlice.set_xlabel(self.yVar)
        if self.data.variables[self.xVar].typecode() == Num.Character:
            value = self.data.variables[self.xVar].getValue()[int(event.xdata),:].tostring()
            verticalSlice.set_ylabel('%s %s = %s' % (self.zVar,self.xVar, value), fontsize = 7)
        else:
            value = self.data.variables[self.xVar].getValue()[int(event.xdata)]
            verticalSlice.set_ylabel('%s %s = %.3f' % (self.zVar,self.xVar, value), fontsize = 7)

        pylab.show()

    def plotXYZ(self):
        """
        This method display a 2D plot.
        """

        if (self.xVar is None) or (self.yVar is None) or (self.zVar is None):
            return

        self.resetPlots()

        xValues = self.data.variables[self.xVar].getValue()
        yValues = self.data.variables[self.yVar].getValue()
        zValues = pylab.transpose(self.data.variables[self.zVar].getValue())

        self.figure.gca().imshow(zValues, aspect = 'auto', label = self.zVar + ' vs ' + self.xVar + ' and ' + self.yVar)
        self.figure.xValues = xValues
        self.figure.yValues = yValues
        self.figure.zValues = zValues

        # The image colorbar.
        self.figure.colorbar(self.figure.gca().get_images()[0])

        # Creation of a cursor to plot orthogonal slices of the 3D plot.
        cursor = Cursor(self.figure.gca(), useblit = True, color ='red', linewidth = 1)
        cursor.ax.set_xlim((0, zValues.shape[1]-1))
        cursor.ax.set_ylim((0, zValues.shape[0]-1))        

        if self.data.variables[self.xVar].typecode() == Num.Character:
            self.figure.gca().set_xticks(range(len(xValues)))
            xTickLabels = ['%s' % xValues[v,:].tostring().rstrip('\x00') for v in range(len(xValues))]
            self.figure.gca().set_xticklabels(xTickLabels)
        else:
            xStep = len(xValues)/10
            if xStep < 1:
                xStep = 1
            self.figure.gca().set_xticks(range(0, len(xValues), xStep))
            self.figure.gca().set_xticklabels(['%.3f' % xValues[v] for v in range(0, len(xValues), xStep)])

        if self.data.variables[self.yVar].typecode() == Num.Character:
            self.figure.gca().set_yticks(range(len(yValues)))
            yTickLabels = ['%s' % yValues[v,:].tostring().rstrip('\x00') for v in range(len(yValues))]
            self.figure.gca().set_yticklabels(yTickLabels)
        else:
            yStep = len(yValues)/10
            if yStep < 1:
                yStep = 1
            self.figure.gca().set_yticks(range(0, len(yValues), yStep))
            self.figure.gca().set_yticklabels(['%.3f' % yValues[v] for v in range(0, len(yValues), yStep)])

        # Set the title.
        try:
            self.figure.gca().set_title(self.data.title)

        except AttributeError:
            self.figure.gca().set_title('')

        # Sets the X-axis label.
        try:
            self.figure.gca().set_xlabel('%s (%s)' % (self.xVar, self.data.variables[self.xVar].units))

        except AttributeError:
            self.figure.gca().set_xlabel(self.xVar)

        # Sets the Y-axis label.
        try:
            self.figure.gca().set_ylabel('%s (%s)' % (self.yVar, self.data.variables[self.yVar].units))
        except AttributeError:
            self.figure.gca().set_ylabel(self.yVar)

        # Bind the plot to the callback |self.displayPlotSlices|.
        self.canvas.mpl_connect('button_press_event', self.displayPlotSlices)

        self.canvas.show()

