"""This module implements all classes used for the generation of combo widgets. 

A combo widget can be defined as a label frame encapsulating several other Tkinter widgets.

Classes:
    * ValidatingEntry   : sets up a Tkinter Entry widget that checks its contents whenever it is changed.
    * IntegerEntry      : sets up a Tkinter Entry widget that checks that the input is an integer.
    * FloatEntry        : sets up a Tkinter Entry widget that checks that the input is a float.
    * StringEntry       : sets up a Tkinter Entry widget that checks that the input is a string matching a given pattern.
    * ComboLabelEntry   : sets up a combo widget made of a Tkinter Label widget embedded in a Tkinter LabelFrame widget.
    * ComboIntegerEntry : sets up a combo widget made of a Tkinter Label widget and IntegerEntry combo widget embedded 
      in a Tkinter LabelFrame widget.
    * ComboFloatEntry   : sets up a combo widget made of a Tkinter Label widget and FloatEntry embedded in a Tkinter 
      LabelFrame widget.
    * ComboStringEntry  : sets up a combo widget made of a Tkinter Label widget and StringEntry embedded in a Tkinter 
      LabelFrame widget.
    * ComboRadiobutton  : sets up a combo widget made of a Tkinter Label widget and a set of Tkinter RadioButton widgets 
      embedded in a Tkinter LabelFrame widget.
    * ComboCheckbutton  : sets up a combo widget made of a Tkinter Label widget and a Tkinter CheckButton widget embedded 
      in a Tkinter LabelFrame widget.
    * ComboOptionMenu   : sets up a combo widget made of a Tkinter Label widget and a Tkinter OptionMenu widget embedded 
      in a Tkinter LabelFrame widget.
    * ComboSpinbox      : sets up a combo widget made of a Tkinter Label widget and a Tkinter Spinbox widget embedded in 
      a Tkinter LabelFrame widget.
    * ComboScale        : sets up a combo widget made of a Tkinter Label widget and a Tkinter Scale widget embedded in 
      a Tkinter LabelFrame widget.
    * ComboButton       : sets up a combo widget made of a Tkinter Button widget and optionnally a Tkinter Entry widget embedded 
      in a Tkinter LabelFrame widget.
    * ComboFileBrowser  : sets up a combo widget made of a Tkinter Label widget, a StringEntry widget and a Tkinter Button 
      widget for file browsing embedded in a Tkinter LabelFrame widget.
    * ComboText         : sets up a combo widget made of a Tkinter Text and a vertical Tkinter Scrollbar widget widget embedded in a 
      Tkinter LabelFrame widget.
    * ComboListbox      : sets up a combo widget made of a Tkinter Listbox widget and a vertical Tkinter Scrollbar widget embedded 
      in a Tkinter LabelFrame widget.
    * StatusBar         : sets up a combo widget made of a Tkinter Scale widget embedded in a Tkinter Frame 
      widget.


The following classes were adapted from the implementation of Frederic Lundh:
    - ValidatingEntry
    - IntegerEntry
    - FloatEntry
"""

# The python distribution modules
import platform
import os
import re
import subprocess
import sys

# The Tcl/Tk modules
from Tkinter import *
from tkFileDialog import askopenfilename, asksaveasfilename
import tkMessageBox

# The nMOLDYN modules
from nMOLDYN.GlobalVariables import GVAR
from nMOLDYN.Preferences import PREFERENCES
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.GUI.HTMLReader.tkhtml import tkHTMLViewer
from nMOLDYN.GUI.Tags import tags

# The platform name.
PLATFORM = platform.system().upper()

# Package path

class PortableToplevel(Toplevel):

    if PLATFORM != 'WINDOWS':
        def transient(self, parent):
            return

class ValidatingEntry(Entry):
    """Base class for |IntegerEntry|, |FloatEntry| and |StringEntry| widgets
    """

    def __init__(self, master, contents = '', **kwargs):
        """The constructor.
        
        @param master: the parent widget.
        
        @param contents: the contents of the Tkinter Entry widget.
        @type contents: string
        """
        
        Entry.__init__(self, master, **kwargs)
        self.variable = StringVar()
        self.cbname = self.variable.trace_variable('w', self._callback)
        self.variable.set(contents)
        self._contents = contents

        self.config(textvariable = self.variable, bg = 'white')

        self.bind('<Return>', lambda event, v = self.variable.get() : self.checkValue(event, v))
        self.bind('<Tab>', lambda event, v = self.variable.get() : self.checkValue(event, v))

    def _callback(self, *dummy):
        contents = self.variable.get()
        newContents = self.checkValue(contents = contents)
        if newContents is not None:
            self._contents = newContents
        self.variable.set(self._contents)

    def checkValue(self, event = None, contents = None):
        pass

    def setValue(self, value):
        """Sets the value of the control variable linked to the Tkinter Entry widget.
        """

        self.variable.set(value)

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter Entry widget.
        """
        return self.variable.get()

class IntegerEntry(ValidatingEntry):
    """Sets up a Tkinter Entry widget that accepts only integer values.
    
    @note: subclass of nMOLDYN.GUI.Widgets.ValidatingEntry class.
    """

    def checkValue(self, event = None, contents = None):
        """Check that the contents of the Tkinter Entry widget is actually an integer.

        @param event: the keyboard event that triggers the checking of the contents of the Tkinter Entry widget.
        @type event: a Tkinter.Event object
        
        @param contents: the contents of the Tkinter Entry widget to check.
        @type contents: integer
        """

        if contents == '':
            return ''

        try:
            return int(contents)
        except:
            return None

    def setValue(self, value = None):
        """Sets the value of the control variable linked to the Tkinter Entry widget.
        """

        v = self.checkValue(contents = value)
        if v is not None:
            self.variable.set(v)

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter Entry widget if it is an integer otherwise
        throws an error.
        """
        try:
            return int(self.variable.get())
        
        except ValueError:
            raise Error('%s is not an integer' % self.variable.get())

class FloatEntry(ValidatingEntry):
    """Sets up a Tkinter Entry widget that accepts only float values.
    
    @note: subclass of nMOLDYN.GUI.Widgets.ValidatingEntry class.
    """

    def checkValue(self, event = None, contents = None):
        """Check that the contents of the Tkinter Entry widget is actually a float.
        
        @param event: the keyboard event that triggers the checking of the contents of the Tkinter Entry widget.
        @type event: a Tkinter.Event object
        
        @param contents: the contents of the Tkinter Entry widget to check.
        @type contents: float
        """

        if contents == '':
            return ''

        try:
            return float(contents)
        except:
            return None

    def setValue(self, value):
        """Sets the value of the control variable linked to the Tkinter Entry widget.
        """

        v = self.checkValue(contents = value)
        if v is not None:
            self.variable.set(v)

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter Entry widget if it is a float otherwise
        throws an error.
        """
        
        try:
            return float(self.variable.get())
        
        except ValueError:
            raise Error('%s is not a float' % self.variable.get())

class StringEntry(ValidatingEntry):
    """Sets up a Tkinter Entry widget that accepts only string values.
    
    @note: subclass of nMOLDYN.GUI.Widgets.ValidatingEntry class.
    """
    
    def __init__(self, master, contents = '', pattern = None, **kwargs):
        """The constructor.
        
        @param master: the parent widget.
        
        @param contents: the contents of the Tkinter Entry widget.
        @type contents: string
        
        @param pattern: if not None, a string specifying a pattern the contents of the Tkinter Entry widget has to
            fulfill.
        @type pattern:string
        """
        
        self.pattern = pattern
        ValidatingEntry.__init__(self, master, contents, **kwargs) 

    def checkValue(self, event = None, contents = None):
        """Check that the contents of the Tkinter Entry widget is actually a string that possibly match the 
        |self.pattern| pattern it is not None..

        @param event: the keyboard event that triggers the checking of the contents of the Tkinter Entry widget.
        @type event: a Tkinter.Event object
        
        @param contents: the contents of the Tkinter Entry widget to check.
        @type contents: string
        """
        
        if self.pattern is None:
            return contents

        if self.variable.get():
            if not re.findall(self.pattern, self.variable.get()):
                LogMessage('warning','%s does not math the pattern %s' % (self.variable.get(), self.pattern),['gui'])
                return None
            
            else:
                return contents

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter Entry widget if it is a string otherwise
        throws an error.
        """
        try:
            if self.variable.get():
                return str(self.variable.get())
            else:
                return None
        except:
            raise Error('%s is not a string' % self.variable.get())

class ComboFrame(LabelFrame):
    """Sets up a Tkinter LabelFrame widget.
    
    Some bindings are set by default:
        - Tkinter '<Motion>' Event will change the mouse cursor to 'question_arrow' symbol.
        - Tkinter '<Leave>' Event will change the mouse cursor to 'arrow' symbol when leaving the widget.
        - Tkinter '<Button-1>' Event will give some information about the widget when clicking on it.
    """

    def __init__(self, master, frameLabel = '', tagName = ''):
        """The constructor.
        """
        LabelFrame.__init__(self, master, relief = GROOVE, bd = 2)
        
        self.tagName = tagName
        
        if frameLabel:                        

            self.frameLabel = StringVar()
            self.frameLabel.set(frameLabel)
            
            # The label of the LabelFrame is a Label widget.
            self.label = Label(self, textvariable = self.frameLabel)
            
            self.config(labelwidget = self.label)
            
            # Some bindings for the Label widget of the LabelFrame.
            self.label.bind('<Motion>', lambda event : changeMouseCursor(event,'question_arrow'))
            self.label.bind('<Leave>', lambda event : changeMouseCursor(event,'arrow'))
            self.label.bind('<Button-1>', informationAboutWidget)

class ComboLabel(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget embedded in a Tkinter LabelFrame widget.
    """

    def __init__(self, master, frameLabel = '', tagName = '', contents = ''):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: the contents of the Tkinter Label widget.
        @type contents: string
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        self.labelVar = StringVar()
        self.labelVar.set(str(contents))

        # The label is set up.
        self.label = Label(self, textvariable = self.labelVar, bd = 2, relief = GROOVE, anchor = W, width = len(self.labelVar.get()))
        self.label.grid(row = 0, column = 0, sticky = E, padx = 3)

    def getValue(self):
        """Returns the value of the Tkinter Label widget.
        """
        return self.labelVar.get()

    def setValue(self, value):
        """Returns the value of the Tkinter Label widget.
        """
        return self.labelVar.set(value)

class ComboIntegerEntry(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a IntegerEntry widget embedded in a Tkinter LabelFrame 
    widget.
    """

    def __init__(self, master, frameLabel = '', tagName = '', contents = ''):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: the contents of the IntegerEntry widget.
        @type contents: integer
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        self.entry = IntegerEntry(self, contents, width = 10)
        self.entry.grid(row = 0, column = 0, sticky = E, padx = 3)

    def getValue(self):
        """Returns the value of the IntegerEntry widget.
        """
        return self.entry.getValue()

    def setValue(self, value):
        """Sets the value of the IntegerEntry widget to |value| (integer).
        
        @param value: the contents to insert in the IntegerEntry widget.
        @type value: integer
        """
        self.entry.setValue(value)

class ComboFloatEntry(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a FloatEntry widget embedded in a Tkinter LabelFrame 
    widget.
    """
    def __init__(self, master, frameLabel = '', tagName = '', contents = ''):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: the contents of the FloatEntry widget.
        @type contents: float
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        self.entry = FloatEntry(self, contents, width = 10)
        self.entry.grid(row = 0, column = 0, sticky = E, padx = 3)

    def getValue(self):
        """Returns the value of the FloatEntry widget.
        """
        return self.entry.getValue()

    def setValue(self, value):
        """Sets the value of the FloatEntry widget to |value| (float).
        
        @param value: the contents to insert in the FloatEntry widget.
        @type value: float
        """
        self.entry.setValue(value)

class ComboStringEntry(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a StringEntry widget embedded in a Tkinter LabelFrame 
    widget.
    """
    def __init__(self, master, frameLabel = '', tagName = '', contents = '', pattern = None):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: the contents of the StringEntry widget.
        @type contents: string

        @param pattern: if not None, a string specifying a pattern the contents of the Tkinter Entry widget has 
            to fulfill.
        @type pattern:string
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        self.entry = StringEntry(self, contents, pattern, width = 30)
        self.entry.grid(row = 0, column = 0, sticky = E, padx = 3)

    def getValue(self):
        """Returns the value of the StringEntry widget.
        """
        return self.entry.getValue()

    def setValue(self, value):
        """Sets the value of the StringEntry widget to |value| (string).
        
        @param value: the contents to insert in the StringEntry widget.
        @type value: string
        """
        self.entry.setValue(value)

class ComboRadiobutton(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a set of Tkinter RadioButton widgets 
    embedded in a Tkinter LabelFrame widget.
    """    

    def __init__(self, master, frameLabel = '', tagName = '', contents = [], default = 0, layout = None):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: a list (string) specifying the names of each Tkinter Radiobutton widget.
        @type contents: list

        @param default: an integer specifying which Tkinter Radiobutton widget will be switched on by default.
        @type default: integer

        @param layout: a tuple (integer) of the form (nrows,ncolumns) specifying the way the set of Tkinter Radiobutton widget will be displayed.
        @type layout: tuple
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        # At least one radiobutton must be provided.    
        if len(contents) == 0:
            raise Error('No radiobutton given to set up.')

        self.contents = contents
        self.radioVar = StringVar()
        self.radioVar.set(self.contents[default])

        radioFrame = Frame(self)
        radioFrame.grid(row = 0, column = 0, sticky = E, padx = 3)

        self.radio = []
        # If no |layout| argument is given, then build the radiobutton on a single row.
        if layout is None:
            row = 0
            for v in contents:
                self.radio.append(Radiobutton(radioFrame,\
                                              text = v,\
                                              variable = self.radioVar,\
                                              value = v,\
                                              width = 14,\
                                              justify = LEFT,\
                                              anchor = W))
                self.radio[-1].radio.grid(row = row, sticky = W)
                row += 1
        else:
            nrows, ncols = layout
            for r in range(nrows):
                for c in range(ncols):                    
                    ind = r*ncols + c
                    if ind < len(contents):                        
                        self.radio.append(Radiobutton(radioFrame,\
                                                      text = contents[r*ncols + c],\
                                                      variable = self.radioVar,\
                                                      value = contents[r*ncols + c],\
                                                      width = 14,\
                                                      justify = LEFT,\
                                                      anchor = W))
                        self.radio[-1].grid(row = r, column = c, sticky = W)

    def getValue(self):
        """Returns the value of the control variable linked to the set of Tkinter Radiobutton widgets.
        """
        return str(self.radioVar.get())

    def setValue(self, value):
        """Sets the value of the control variable linked to the set of Tkinter Radiobutton widgets.
        
        @param value: the value the control variable linked to the set of Tkinter Radiobutton widgets will take. 
            Must be an element of |self.contents|
        @type value: string
        """
        if value in self.contents:
            self.radioVar.set(value)
        else:
            self.radioVar.set(self.contents[0])

class ComboCheckbutton(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a set of Tkinter Checkbutton widget 
    embedded in a Tkinter LabelFrame widget.
    """    

    def __init__(self, master, frameLabel = '', tagName = '', contents = '', onvalue = True, offvalue = False, default = False):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: the name of the Tkinter Checkbutton widget.
        @type contents: string

        @param onvalue: the value of the Tkinter Checkbutton widget when switched on.
        @type onvalue: boolean, integer or string

        @param offvalue: the value of the Tkinter Checkbutton widget when switched off.
        @type offvalue: boolean, integer or string

        @param default: a boolean specifying whether the Tkinter Checkbutton widget is switched on (True) or switched off (False)
            by default.
        @type default: boolean
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName, )

        checkFrame = Frame(self)
        checkFrame.grid(row = 0, column = 0, sticky = E, padx = 3)

        self.offvalue = offvalue            
        self.onvalue = onvalue

        if isinstance(self.onvalue, str):
            self.checkVar = StringVar()

        elif isinstance(self.onvalue, bool):
            self.checkVar = BooleanVar()

        elif isinstance(self.onvalue, int):
            self.checkVar = IntVar()

        else:
            raise Error('%s can not be used to define the state of a checkbutton.' % onvalue)

        if onvalue.__class__ != offvalue.__class__:
            raise Error('%s and %s must be of the same class to setup properly the checkbutton.' % (offvalue, onvalue))

        if default:
            self.checkVar.set(self.onvalue)

        else:
            self.checkVar.set(self.offvalue)

        self.check = Checkbutton(checkFrame,\
                                 text = contents,\
                                 textvariable = self.checkVar,\
                                 variable = self.checkVar,\
                                 onvalue = onvalue,\
                                 offvalue = offvalue,\
                                 width = 8,\
                                 justify = LEFT)
        self.check.grid(sticky = E)

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter Checkbutton widget.
        """
        return self.checkVar.get()

    def setValue(self, value):
        """Sets the value of the control variable linked to the Tkinter Checkbutton widget.
        
        @param value: the value the control variable linked to the Tkinter Checkbutton widget will take. 
            Must be an element of |self.contents|
        @type value: string
        """
        try:
            self.checkVar.set(value)

        except:
            self.checkVar.set(self.offvalue)

class ComboOptionMenu(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a Tkinter OptionMenu widget embedded 
    in a Tkinter LabelFrame widget.
    """    

    def __init__(self, master, frameLabel = '', tagName = '', contents = [], default = 0):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: a list (string) specifying the names of for the entries each Tkinter OptionMenu widget.
        @type contents: list

        @param default: an integer specifying which entry of the Tkinter optionMenu widget will be displayed by default.
        @type default: integer
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        # At least one option must be provided.
        if len(contents) == 0:
            raise Error('No option given to set up this combowidget.')

        # The frame that will host the option menu widget is set up.
        optionFrame = Frame(self)
        optionFrame.grid(row = 0, column = 1, sticky = E, padx = 3)

        self.contents = contents
        self.optionVar = StringVar()
        self.optionVar.set(self.contents[default])

        self.option = OptionMenu(optionFrame, self.optionVar, *contents)
        self.option.grid(sticky = E)

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter OptionMenu widget.
        """
        return str(self.optionVar.get())

    def setValue(self, value):
        """Sets the value of the control variable linked to the Tkinter OptionMenu widget.
        
        @param value: the value the control variable linked to the Tkinter OptionMenu widget will take. 
            Must be an element of |self.contents|
        @type value: string
        """
        if value in self.contents:
            self.optionVar.set(value)
        else:
            self.optionVar.set(self.contents[0])

class ComboSpinbox(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a Tkinter Spinbox widget embedded 
    in a Tkinter LabelFrame widget.
    """    

    def __init__(self, master, frameLabel = '', tagName = '', contents = [], default = 0):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: a tuple (string) specifying the names of for the entries to insert in the Tkinter Spinbox widget.
        @type contents: tuple

        @param default: an integer specifying which entry of the Tkinter Spinbox widget will be displayed by default.
        @type default: integer
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        # At least one option must be provided.
        if len(contents) == 0:
            raise Error('No option given to set up this combowidget.')

        self.type = type(contents[0])
        for v in contents:
            if type(v) != self.type:
                self.type = str

        # The frame that will host the option menu widget is set up.
        spinboxFrame = Frame(self)
        spinboxFrame.grid(row = 0, column = 1, sticky = E, padx = 3)

        self.contents = contents
        self.default = default

        self.spinbox = Spinbox(spinboxFrame, values = self.contents, wrap = True, width = 4, bg = 'white')
        self.spinbox.grid(sticky = E)

        self.setValue(self.contents[self.default])

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter Spinbox widget.
        """
        return self.type(self.spinbox.get())

    def setValue(self, value):
        """Sets the value of the control variable linked to the Tkinter Spinbox widget.
        
        @param value: the value the control variable linked to the Tkinter Spinbox widget will take. 
            Must be an element of |self.contents|
        @type value: string
        """

        if value in self.contents:
            v = self.type(self.spinbox.get())
            while v != value:
                self.spinbox.invoke('buttonup')
                v = self.type(self.spinbox.get())

class ComboScale(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget and a Tkinter Scale widget embedded 
    in a Tkinter LabelFrame widget.
    """    

    def __init__(self, master, frameLabel = '', tagName = '', contents = []):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: a 3 elements list (integer) specifying respectively the first, the last and the resolution values
            for the Tkinter Scale widget.
        @type contents: list
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        # At least one option must be provided.
        if len(contents) == 0:
            raise Error('No values given to set up this combowidget.')

        self.from_, self.to, self.resolution = contents

        self.scale = Scale(self, from_ = self.from_, to = self.to, resolution = self.resolution, length = 200, orient = HORIZONTAL)
        self.scale.grid(row = 0, column = 0, sticky = E, padx = 3)

    def getValue(self):
        """Returns the value of the control variable linked to the Tkinter Scale widget.
        """
        return self.scale.get()

    def setValue(self, value):
        """Sets the value of the control variable linked to the Tkinter Scale widget.
        
        @param value: the value the Tkinter Scale widget will take. 
        @type value: integer
        """
        self.scale.set(value)

class ComboButton(ComboFrame):
    """Sets up a combo widget made of a Tkinter Button widget and optionnally a Tkinter Entry widget embedded 
    in a Tkinter LabelFrame widget.
    """    

    def __init__(self, master, frameLabel = '', tagName = '', contents = '', withEntry = None, command = None):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: string specifying the label of the Tkinter Button widget.
        @type contents: string

        @param withEntry: if not None, a Tkinter Entry widget will be placed beside the Tkinter Button
            with |withEntry| as contents.
        @type withEntry: boolean
        
        @param command: the function to call when clicking on the Tkinter Button widget.
        @type command: function
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        self.withEntry = withEntry

        self.results = None
        self.command = command

        if self.withEntry is not None: 
            # The entry widget is set up.
            self.entry = StringEntry(self, self.withEntry, width = 60)
            self.entry.grid(row = 0, column = 0, sticky = E)

        self.button = Button(self, text = contents, command = self.pushButton, width = 8)
        self.button.grid(row = 0, column = 1, sticky = E, padx = 3)

    def pushButton(self):
        """Executes the command linked to the Tkinter Button widget.
        """

        if self.command is not None:
            self.results = self.command()
            if self.withEntry is not None:
                self.entry.setValue(str(self.results.getValue()))

    def getValue(self):
        """Returns the output value of the |self.command| function.
        """

        if self.withEntry is not None:
            val = self.entry.getValue()
        else:
            # Case where |self.results| is instance of a ComboWidget having the getValue method.
            try:
                val = self.results.getValue()
            except:
                val = self.results

        return val

    def setValue(self):
        pass

class ComboFileBrowser(ComboFrame):
    """Sets up a combo widget made of a Tkinter Label widget, a StringEntry widget and a Tkinter Button 
    widget for file browsing embedded in a Tkinter LabelFrame widget.
    """    

    def __init__(self, master, frameLabel = '', tagName = '', contents = '', save = False, command = None, filetypes = []):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: string specifying the contents of the Tkinter Entry widget.
        @type contents: string

        @param save: boolean specifying whether the file browser is for saving (True) or for loading (False).
        @type save: boolean
        
        @param command: the function to call when browsing the file.
        @type command: function
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        self.save = save
        self.filetypes = filetypes
        self.command = command

        # The entry widget is set up.
        self.entry = StringEntry(self, contents, width = 60)
        self.entry.grid(row = 0, column = 0, sticky = E)

        self.button = Button(self, bitmap = "@"+os.path.join(GVAR['nmoldyn_path'], 'Images', 'openfile.xbm'), command = self.browse)
        self.button.grid(row = 0, column = 1, sticky = E)

    def browse(self):
        """Executes the command linked to the Tkinter Button widget.
        """

        # No command has been given when the instance of the ComboFileBrowser was created so just launch
        # a classical TkFileDialog browser and sets the entry of the Combo FileBrowser to the filename.
        if self.command is None:
            if self.save:
                filename = asksaveasfilename(parent = self,\
                                             title = 'Save file as',\
                                             filetypes = self.filetypes,\
                                             initialdir = PREFERENCES['outputfile_path'])
            else:
                filename = askopenfilename(parent = self,\
                                           title = 'Open file',\
                                           filetypes = self.filetypes,\
                                           initialdir = PREFERENCES['trajfile_path'])
            self.entry.setValue(filename)
        else:
            self.command()

    def getValue(self):
        """Returns the value of the StringEntry widget.
        """
        return self.entry.getValue()

    def setValue(self, value):
        """Sets the value of the StringEntry widget to |value| (string).
        
        @param value: the contents to insert in the StringEntry widget.
        @type value: string
        """
        self.entry.setValue(value)

class ComboText(ComboFrame):
    """Sets up a combo widget made of a Tkinter Text widget and a vertical Tkinter Scrollbar widget embedded in a Tkinter LabelFrame 
    widget.
    """
    
    def __init__(self, master = None, frameLabel = '', tagName = '', contents = None):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: string specifying the contents to insert in the Tkinter Text widget.
        @type contents: string
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        # The text widget and its associated scrollbar.
        sb = Scrollbar(self, orient = VERTICAL)
        sb.grid(row = 0, column = 1, sticky = NSEW)

        self.text = Text(self,\
                         state = DISABLED,\
                         yscrollcommand = sb.set,\
                         bg = 'white')
        self.text.grid(row = 0, column = 0, sticky = "NSEW")
        sb['command'] = self.text.yview

        self.insert(contents = contents)

    def cleanup(self):
        """Deletes the contents of the Tkinter Text widget.
        """
        
        self.text.config(state = NORMAL)
        self.text.delete('0.0', END)
        self.text.config(state = DISABLED)

    def insert(self, cursor = END, contents = None, tag = None):
        """Inserts |contents| text in the Tkinter Text widget.
        
        @param contents: a string specifying the text to insert in the Tkinter Text widget.
        @type contents: string
        
        @param tag: if not None, a string specifying the tag name to associate with the inserted text.
        @type tag: string
        """

        # If the contents is None, then do nothing.
        if contents is None:
            return

        # Activate the text widget for writing.
        self.text.config(state = NORMAL)

        # Insert the contents with a tag eventually.
        self.text.insert(cursor, contents, tag)

        # Desactivates the text widget for safety reason.
        self.text.config(state = DISABLED)

    def tag(self, pattern, tag):
        """Searches and tags all the occurences of a  given pattern.

        @param pattern: the pattern to search and tag.
        @type pattern: string

        @param tag: the tag to associate with |pattern|.
        @type tag: string
        """

        self.patternLength = IntVar()
        start = '1.0'
        # Loop over all the text widget contents to find all the coocurences of |pattern|.
        # Each time an occurence is found, it is tagged with the given tag |tag|.
        while 1:
            match1 = self.text.search(pattern, index = start, stopindex = END, count = self.patternLength)
            # The break condition of the infinite loop. No pattern has been found anymore.
            if not match1:
                break
            match2 = '%s+%sc' % (match1, self.patternLength.get())
            if not self.text.tag_names(match1):
                self.text.tag_add(tag, match1, match2)
            start = match2

class ComboListbox(ComboFrame):
    """Sets up a combo widget made of a Tkinter Listbox widget and a vertical Tkinter Scrollbar widget embedded in a Tkinter LabelFrame 
    widget.
    """

    def __init__(self, master = None, frameLabel = '', tagName = '', contents = []):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        
        @param frameLabel: the label for the Tkinter LabelFrame widget.
        @type frameLabel: string or Tkinter.StringVar object.
        
        @param tagLabel: the tag used for to document the widget. If set to '', the widget will not be documented.
        @type tagLabel: string.

        @param contents: a list (string) specifying the items to insert in the Tkinter Listbox widget.
        @type contents: list
        """

        # The comboframe is set up.
        ComboFrame.__init__(self, master, frameLabel, tagName)

        # The listbox that will display and its associated scrollbar.
        sb = Scrollbar(self, orient = VERTICAL)
        sb.grid(row = 0, column = 1, sticky = "NSEW")

        self.lb = Listbox(self, yscrollcommand = sb.set, bg = 'white')
        self.lb.grid(row = 0, column = 0, sticky = "NSEW")

        sb['command'] = self.lb.yview
        self.lb.bind('<ButtonRelease-1>', self.onListboxEntrySelection)
        self.lb.bind('<Motion>', self.focus)

        self.insert(contents)

    def focus(self, event):
        """Sets the focus to the Tkinter Listbox widget.
        
        @param event: the event triggering the focus.
        @type event: a Tkinter.Event object
        """
        self.lb.focus_set()

    def insert(self, contents = []):
        """Inserts a list of items (string) in the Tkinter Listbox widget.
        
        @param contents: a list (string) specifying the items to insert in the Tkinter Listbox widget.
        @type contents: list
        """

        # Some checkings.
        if not isinstance(contents, (list, tuple, dict)):
            raise Error('The object of type %s cannot be inserted into a listbox.' % contents.__class__.__name__)

        self.cleanup()

        # Loop over list to insert in the listbox.
        for c in contents:
            self.lb.insert(END, c)

    def cleanup(self):
        """Deletes all the items of the Tkinter Listbox widget.
        """

        self.lb.delete(0, END)

    def onListboxEntrySelection(self, event):
        """Updates the selection of the Tkinter Listbox widget when |event| Tkinter Event occurs.
        
        @param event: the event triggering the callback.
        @type event: a Tkinter.Event object
        """
        index = self.lb.nearest(event.y)
        if index == -1:
            return

        if not self.lb.selection_includes(index):
            self.lb.selection_clear(index)
        else:
            self.lb.selection_set(index)

class StatusBar(Frame):
    """Sets up a combo widget made of a Tkinter Scale widget embedded in a Tkinter Frame widget.
    """

    def __init__(self, master = None):
        """The constructor.
        
        @param master: the parent widget of the combo widget.
        """

        Frame.__init__(self, master, bd = 2, relief = RAISED)

        self.status = Scale(self,\
                            orient = HORIZONTAL,\
                            from_ = 0,\
                            to = 100,\
                            state = DISABLED,\
                            troughcolor = 'red')

        self.status.grid(row = 0, column = 0, padx = 3, pady = 3)

    def setValue(self, text):
        """Sets the text that will be displayed within the status bar.
        
        @param text: the value to display in the Status Bar.
        @type text: string or float
        """

        self.status.config(state = NORMAL)
        self.status.set(text)
        self.status.config(state = DISABLED)
        self.status.update()

def changeMouseCursor(event, cursor):
    """Changes the mouse cursor aspect.
    
    @param event: the event triggering the change of mouse cursor. 
    @type event: a Tkinter.Event object
    
    @param cursor: the Tkinter name of the mouse cursor.
    @type cursor: string
    """

    event.widget.config(cursor = cursor)
    
def informationAboutWidget(event):
    """Checks the label of the Tkinter LabelFrame of a combo widget and displays the corresponding informations
    documented in nMOLDYN users guide.

    @param event: the event triggering the request for information about the combo widget. 
    @type event: a Tkinter.Event object
    """

    tagName = event.widget.nametowidget(event.widget.winfo_parent()).tagName
        
    if tags.has_key(tagName):
        
        docStyle = PREFERENCES['documentation_style'] 
        
        if docStyle == 'html':
            url = os.path.join(GVAR['nmoldyn_path'], 'Doc', 'UsersGuide', 'HTML', tags[tagName]['html'])
            
            top = PortableToplevel()
            top.title('Information')

            top.wm_minsize(400, 200)

            viewer = tkHTMLViewer(top)
            viewer.display(url)

            top.wait_visibility()
            top.grab_set()
            top.focus_set()            

            top.wait_window(top)

        elif docStyle == 'pdf':
            # The name of the users guide pdf file.
            pdfFile = os.path.join(GVAR['nmoldyn_path'], 'Doc', 'UsersGuide', 'PDF', 'nMOLDYN_ug.pdf')

            if PLATFORM == 'WINDOWS':
                try:
                    subprocess.call([PREFERENCES['acroread_path'], '/a nameddest=%s' % tags[tagName]['pdf'],  pdfFile])
                except:
                    os.startfile(pdfFile)

            elif PLATFORM == 'DARWIN':
                try:
                    os.system('osascript -e \'tell application "Adobe Reader"\nopen "%s" options "nameddest=%s"\nend tell\'' % (pdfFile, tags[tagName]['pdf']))
                except:
                    subprocess.call(['open', pdfFile])

            else:
                try:
                    subprocess.call([PREFERENCES['acroread_path'], '/a nameddest=%s' % tags[tagName]['pdf'], pdfFile])
                except:
                    LogMessage('warning','Unable to read the file %s\nwith the acrobat reader located in %s.' % (pdfFile, PREFERENCES['acroread_path']),['gui'])
                    pass
                else:
                    return

    else:
        tkMessageBox.showinfo("nMOLDYN help", "Not documented.")

