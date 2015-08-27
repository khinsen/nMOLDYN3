import Tkinter, Dialog, FileDialog
import os, string
from Scientific.TkWidgets import FloatEntry, IntEntry, StatusBar, ModalDialog
from Scientific.TkWidgets.TkPlotCanvas import PlotCanvas,PolyLine,PlotGraphics
from Scientific.TkWidgets.TkVisualizationCanvas import VisualizationCanvas, \
     VisualizationGraphics, PolyLine3D
from Scientific.DictWithDefault import DictWithDefault
from Scientific.IO import ArrayIO
import Numeric 
#from Tkinter import *
GUIError = 'GUIError'


class FilenameEntry(Tkinter.Frame):

    def __init__(self, master, browse_pattern = '*'):
	self.pattern = browse_pattern
	Tkinter.Frame.__init__(self, master)
	#Tkinter.Label(self, text=text).pack(side=Tkinter.LEFT)
	self.filename = Tkinter.StringVar()
	Tkinter.Button(self, text="Browse...",
		       command=self.browse).pack(side=Tkinter.RIGHT)
	entry = Tkinter.Entry(self, textvariable=self.filename)
	entry.pack(side=Tkinter.RIGHT, expand=1, fill=Tkinter.X)
	entry.icursor("end")

    def browse(self):
	file = FileDialog.LoadFileDialog(self).go(pattern=self.pattern) # [0]
	if file: self.filename.set(file)

    def get(self):
	filename =  self.filename.get()
	return filename


class FloatEntry(Tkinter.Frame):

    def __init__(self, master, text, init = None, lower=None, upper=None,
                 name = None):
	self.text = text
	self.lower = lower
	self.upper = upper
        if name is None:
            name = text
        self.name = name
        if len(text) > 0:
            label_text = string.upper(text[0])+text[1:]+':'
            xx=10
        else:
            label_text = ''
            xx=1
	Tkinter.Frame.__init__(self, master)
	Tkinter.Label(self, text=label_text,height=1,width=xx).pack(side=Tkinter.LEFT)
	self.value = Tkinter.DoubleVar()
	if init is not None:
	    self.value.set(init)
	self.entry = Tkinter.Entry(self, textvariable=self.value)
	self.entry.pack(side=Tkinter.RIGHT, anchor=Tkinter.E,
		   expand=1, fill=Tkinter.X)
	self.entry.icursor("end")

    def bind(self, sequence=None, func=None, add=None):
        self.entry.bind(sequence, func, add)

    def set(self, value):
	return self.value.set(value)

    def get(self):
	try:
	    value = self.value.get()
	except Tkinter.TclError:
	    Dialog.Dialog(self, title='Illegal value',
			  text='The value of "' + self.name +
			       '" must be a number.',
			  bitmap='warning', default=0,
			  strings = ('Cancel',))
	    raise ValueError
	range_check = 0
	if self.lower is not None and value < self.lower:
	    range_check = -1
	if self.upper is not None and value > self.upper:
	    range_check = 1
	if range_check != 0:
	    text = 'The value of "' + self.name + '" must not be '
	    if range_check < 0:
		text = text + 'smaller than ' + `self.lower` + '.'
	    else:
		text = text + 'larger than ' + `self.upper` + '.'
	    Dialog.Dialog(self, title='Value out of range', text=text,
			  bitmap='warning', default=0,
			  strings = ('Cancel',))
	    raise ValueError
	return value


class IntEntry(FloatEntry):

    def get(self):
        value = FloatEntry.get(self)
        ivalue = int(value)
        if ivalue != value:
	    Dialog.Dialog(self, title='Illegal value',
			  text='The value of "' + self.name +
			       '" must be an integer.',
			  bitmap='warning', default=0,
			  strings = ('Cancel',))
            raise ValueError
        return ivalue

class ButtonBar(Tkinter.Frame):

    def __init__(self, master, left_button_list, right_button_list):
	Tkinter.Frame.__init__(self, master, bd=2, relief=Tkinter.SUNKEN)
	for button, action in left_button_list:
	    Tkinter.Button(self, text=button,
			   command=action).pack(side=Tkinter.LEFT)
	for button, action in right_button_list:
	    Tkinter.Button(self, text=button,
			   command=action).pack(side=Tkinter.RIGHT)


class StatusBar(Tkinter.Frame):

    def __init__(self, master):
	Tkinter.Frame.__init__(self, master, bd=2, relief=Tkinter.RAISED)
	self.text = Tkinter.Label(self, text='')
	self.text.pack(side=Tkinter.LEFT, expand=Tkinter.YES)

    def set(self, text):
	self.text.configure(text = text)
	self.text.update_idletasks()
	self.master.configure(cursor='watch')
	self.update()
	self.update_idletasks()

    def clear(self):
	self.text.configure(text = '')
	self.text.update_idletasks()
	self.master.configure(cursor='top_left_arrow')
	self.update_idletasks()


class Tkwindow(Tkinter.Toplevel):

    def __init__(self, master):
	if master is None:
	    master = _root
	else:
	    master.addDependent(self)
	Tkinter.Toplevel.__init__(self, master)
	self.protocol("WM_DELETE_WINDOW", self.close)
	self._dependents = []

    def open(self):
	global _mainloop_started
	_open_windows.append(self)
	if not _mainloop_started:
	    _root.withdraw()
	    _root.mainloop()
	    _mainloop_started = 1

    def close(self):
	for window in self._dependents:
	    window.close()
	self.destroy()
	if self.master is not _root:
	    self.master.removeDependent(self)
	_open_windows.remove(self)
	if not _open_windows:
	    _root.destroy()

    def addDependent(self, window):
	self._dependents.append(window)

    def removeDependent(self, window):
	self._dependents.remove(window)

#
# Plotting via external program
#
plot_program = None
try:
    plot_program = os.environ['PLOTPROGRAM']
except KeyError: pass

def externalPlot(data):
    if plot_program is not None:
        filename = tempfile.mktemp()
        ArrayIO.writeArray(data, filename)
        if os.fork() == 0:
            pipe = os.popen(plot_program + ' ' + filename + \
                            ' 1> /dev/null 2>&1', 'w')
            pipe.close()
            os.unlink(filename)
            os._exit(0)

#
# Sort variables into categories for display
#
_category = {
    'energy': 'energy',
    'temperature': 'thermodynamic',
    'pressure': 'thermodynamic',
    'volume': 'thermodynamic',
    'configuration': 'configuration',
    'coordinate': 'auxiliary',
    'masses': 'masses',
    'time': 'time',
    }

def categorizeVariables(variable_names):
    categories = DictWithDefault([])
    for name in variable_names:
        words = string.split(name, '_')
        try:
            c = _category[words[-1]]
        except KeyError:
            c = ''
        if c:
            categories[c].append(name)
    for list in categories.values():
        list.sort()
    if categories.has_key('energy'):
        list = categories['energy']
        for variable in ['kinetic_energy', 'potential_energy']:
            if variable in list:
                list.remove(variable)
                list = [variable] + list
        categories['energy'] = list
    return categories

#
# Determine 'nice' plot range for a given value interval
#
def plotRange(lower, upper):
    range = upper-lower
    if range == 0.:
        return lower-0.5, upper+0.5
    log = Numeric.log10(range)
    power = Numeric.floor(log)
    fraction = log-power
    if fraction <= 0.05:
        power = power-1
    grid = 10.**power
    lower = lower - lower % grid
    mod = upper % grid
    if mod != 0:
        upper = upper - mod + grid
    return lower, upper

#
# Simple text display window
#
class TextViewer(Tkinter.Toplevel):

    def __init__(self, master, text, title=''):
        Toplevel.__init__(self, master)
	self.title(title)
	Label(self, text=text,
		      justify=LEFT).pack(side=TOP)
	button_bar = Frame(self)
	button_bar.pack(side=BOTTOM)
	Button(button_bar, text='Close',
	       command=self.destroy).pack(side=RIGHT)

#
# Multiple plot window
#
class PlotWindow(Tkinter.Toplevel):

    def __init__(self, master, title):
        Tkinter.Toplevel.__init__(self, master)
	self.title(title)
        self.plotlist = []
        self.protocol("WM_DELETE_WINDOW", self.close)

    def close(self):
        for plot in self.plotlist:
            self.master._unregisterPlot(plot)
        self.destroy()

    def registerPlot(self, plot):
        self.plotlist.append(plot)
        self.master._registerPlot(plot)
    
    def setSelection(self, plot):
        if self.master.selection is not None:
            i1, i2 = self.master.selection
            plot.select((self.master.time[i1], self.master.time[i2]))

#
# Standard plot window for independent quantities
#
class PlotViewer(PlotWindow):

    def __init__(self, master, title, time, trajectory, names):
        PlotWindow.__init__(self, master, title)
        self.time = time
	self.step = max(1, len(time)/300)
	self.setup(time, trajectory, names)

    def setup(self, time, trajectory, names):
        time = time[::self.step]
        for n in names:
            data = Numeric.zeros((len(time), 2), Numeric.Float)
            data[:, 0] = time
            data[:, 1] = getattr(trajectory, n)[::self.step]
	    self.plotBox(n, data)

    def plotBox(self, name, data, data_range=None):
	box = Tkinter.Frame(self, border=2, relief=Tkinter.SUNKEN)
	box.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=Tkinter.YES)
	frame = Tkinter.Frame(box, background='grey')
	frame.pack(side=Tkinter.TOP, fill=Tkinter.X, expand=Tkinter.NO)
	Tkinter.Label(frame, text=string.capitalize(string.join(
	                  string.split(name , '_'), ' ')),
	      background='grey').pack(side=Tkinter.LEFT)
        if data_range is None:
            min = Numeric.minimum.reduce(data[:,1])
            max = Numeric.maximum.reduce(data[:,1])
            min, max = plotRange(min, max)
        else:
            min, max = data_range
        plot_objects = []
        plot_data = data
        time = plot_data[:,0]
        jumps = Numeric.repeat(Numeric.arange(len(time)-1),
                               Numeric.less(time[1:], time[:-1]))+1
        for i in self.master.restarts:
            plot_objects.append(PolyLine([(self.time[i], min),
                                          (self.time[i], max)],
                                         color='black',
                                         stipple='gray25'))
        plot_objects.insert(0, PolyLine(plot_data, color = 'red'))
        plot = PlotCanvas(box, 400, 100, zoom=1,
                          select=self.master._selectRange)
        plot.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=Tkinter.YES)
        plot.draw(PlotGraphics(plot_objects),
                  'automatic', (min, max))
        plot.bind('<Double-Button-1>', lambda event, d=plot_data:
                                       externalPlot(d))
        self.registerPlot(plot)
        self.setSelection(plot)

#
# Special plot window for energy terms: uniform scale and
# display of sum.
#
class EnergyViewer(PlotViewer):

    def setup(self, time, trajectory, names):
        time = time[::self.step]
	sum = 0.
        range = 0.
        for n in names:
            data = getattr(trajectory, n)[::self.step]
	    sum = sum + data
            upper = Numeric.maximum.reduce(data)
            lower = Numeric.minimum.reduce(data)
            lower, upper = plotRange(lower, upper)
            range = max(range, upper-lower)
        upper = Numeric.maximum.reduce(sum)
        lower = Numeric.minimum.reduce(sum)
        lower, upper = plotRange(lower, upper)
        range = 0.5*max(range, upper-lower)
        for n in names:
            data = getattr(trajectory, n)[::self.step]
	    data = Numeric.transpose(Numeric.array([time, data]))
            mean = Numeric.add.reduce(data[:, 1])/len(data)
	    self.plotBox(n, data, (mean-range, mean+range))
        if len(names) > 1:
            data = Numeric.transpose(Numeric.array([time, sum]))
            mean = Numeric.add.reduce(data[:, 1])/len(data)
            self.plotBox('Sum', data, (mean-range, mean+range))
