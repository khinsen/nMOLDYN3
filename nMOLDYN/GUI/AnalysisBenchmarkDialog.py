"""This modules implements I{Help-->nMOLDYN benchmark} dialog.

Classes:
    * AnalysisBenchmarkDialog: creates I{Help-->nMOLDYN benchmark} dialog used to run the stability tests 
      between nMOLDYN v2.2.5 and the current version of nMOLDYN.
"""

# The python distribution modules
import copy
import os
import re

# The Tcl/Tk modules
from Tkinter import *
from tkMessageBox import showinfo

# The nMOLDYN modules
from nMOLDYN.__pkginfo__ import __version__ as NMOLDYN_VERSION
from nMOLDYN.Core.Error import Error
from nMOLDYN.GUI.Widgets import *
from nMOLDYN.Tests.StabilityTests import *

class AnalysisBenchmarkDialog(PortableToplevel):
    """Builds the dialog for nMOLDYN analysis benchmarks.
    
    These benchmarks will guarantee the stability of the current version with a reference one that is the 
    version 2.2.5 the last official release of nMOLDYN.
    
    @note: these benchmarks should be performed every time the code source of nMOLDYN has been touched for 
        any reason.
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

        # List of the test sets available.
        self.availableTestSets = sorted(availableTests + ['ALL'])

        # This is the tag counter to allow the user to click on each failing test id to get more informations.
        self.tagCounter = 1

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
        """Creates dialog body. Returns widget that should have initial focus.
        """

        settingsFrame = LabelFrame(master, text = 'Settings', bd = 2, relief = GROOVE)
        settingsFrame.grid(row = 0, column = 0, sticky = EW, padx = 3, pady = 3)
        settingsFrame.grid_columnconfigure(0, weight = 1)
        settingsFrame.grid_rowconfigure(0, weight = 1)

        # The combo listbox that will contain the list of the available tests.
        self.testListLb = ComboListbox(settingsFrame,'Available tests',[])
        self.testListLb.grid(row = 0, column = 0, padx = 2, pady = 2, sticky = 'NS')
        self.testListLb.grid_rowconfigure(0, weight = 1)
        self.testListLb.lb.config({'selectmode' : MULTIPLE})
        self.testListLb.lb.bind('<Button-3>',self.infoAboutTest)

        # The available test sets are inserted in the listboxes.
        for testSet in self.availableTestSets:
            self.testListLb.lb.insert(END, testSet)

        self.testResultsText = ComboText(settingsFrame,\
                                         frameLabel = 'Results for the selected benchmark(s)')
        self.testResultsText.grid(row = 0, column = 1, padx = 2, pady = 2, sticky = 'NS')
        self.testResultsText.config({'relief' : SUNKEN, 'bd' : 2})
        self.testResultsText.grid_rowconfigure(0, weight = 1)
        self.testResultsText.text.config({'wrap' : WORD})

        self.testResultsText.text.tag_config('12', foreground = 'black', font = ('Courier','12'))
        self.testResultsText.text.tag_config('g12', foreground = 'green', font = ('Courier','12'))
        self.testResultsText.text.tag_config('12B', foreground = 'black', font = ('Courier','12','bold'))

        return None
        
    def buttonbox(self):
        """Adds standard button box."""

        # The frame that contains the 'Cancel' and 'OK' buttons.
        box = LabelFrame(self, text = 'Actions', bd = 2, relief = GROOVE)
        box.grid(row = 1, column = 0, sticky = EW, padx = 3, pady = 3)
        box.grid_columnconfigure(0, weight = 1)

        w = Button(box, text = "Cancel", width=10, command = self.cancel)
        w.grid(row = 0, column = 0, sticky = E)
        w = Button(box, text = "Clear", width=10, command = self.clearResults)
        w.grid(row = 0, column = 1, sticky = E)
        w = Button(box, text = "OK", width=10, command = self.ok, default=ACTIVE)
        w.grid(row = 0, column = 2, sticky = E)
        
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

    # Standard button semantics.
    def ok(self, event = None):
        if not self.validate():
            self.initial_focus.focus_set()
            return

        self.update_idletasks()

        self.apply()
        
    def cancel(self, event = None):

        # Put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    # Command hooks
    def validate(self):
        
        # This is the list of the test sets selected by the |checkSimulationsDialog| dialog.
        self.selectedTestSets = [self.testListLb.lb.get(int(testSet)) for testSet in self.testListLb.lb.curselection()]

        if not self.selectedTestSets:
            raise Error('No Test selected.')
                                            
        return True

    def apply(self):

        # The tests results is switched to NORMAL state. Some informations can be written in.
        self.testResultsText.text.config(state = NORMAL)
                            
        if 'ALL' in self.selectedTestSets:
            self.selectedTestSets = availableTests

        comp = 1
        # Loop over each selected set of test (e.g. MSD, AT ...).
        for test in self.selectedTestSets:
                                
            t = eval(test+'Tests()')
            t.run()
                                                                            
            self.testResultsText.insert(END, t.summary)
                
            self.testResultsText.tag('#', 'g12')
                
            [self.testResultsText.tag(p, '12B') for p in re.findall('Summary.*\n', t.summary)]
                                        
            [self.testResultsText.tag(p, '12B') for p in re.findall('Number of.*\n', t.summary)]

            [self.testResultsText.tag(p, 'success%s' % comp) for p in re.findall('Successful tests:.*\n', t.summary)]
                    
            [self.testResultsText.tag(p, 'failure%s' % comp) for p in re.findall('Failing tests:.*\n', t.summary)]

            self.testResultsText.text.tag_config('success%s' % comp, foreground = 'blue', font = ('Courier','10','bold'))

            self.testResultsText.text.tag_config('failure%s' % comp, foreground = 'red', font = ('Courier','10','bold'))
                    
            [self.testResultsText.tag(p, '12B') for p in re.findall('Percentage.*\n', t.summary)]
                    
            [self.testResultsText.tag(p, '12B') for p in re.findall('Time to run.*\n', t.summary)]                    
                                    
            # This binding displays the contents of the successful tests when the user click on it with the mouse left button.
            self.testResultsText.text.tag_bind('success%s' % comp, '<Button-1>', lambda event, p = t.successes : self.displayTest(event, p))

            # This binding changes the mouse icon to a question when the user passes over a successful test.
            self.testResultsText.text.tag_bind('success%s' % comp, '<Motion>', lambda event : changeMouseCursor(event,'question_arrow'), '+')
                
            # This binding returns the mouse icon to its default value when the user passes out a successful test.
            self.testResultsText.text.tag_bind('success%s' % comp, '<Leave>', lambda event : changeMouseCursor(event,'arrow'), '+')

            # Some tag bindings for the failing tests.
            # This binding displays the contents of the failing tests when the user click on it with the mouse left button.
            self.testResultsText.text.tag_bind('failure%s' % comp, '<Button-1>', lambda event, p = t.failures : self.displayTest(event, p))
                
            # This binding changes the mouse icon to a question when the user passes over a failing test.
            self.testResultsText.text.tag_bind('failure%s' % comp,'<Motion>', lambda event : changeMouseCursor(event,'question_arrow'), '+')
                
            # This binding returns the mouse icon to its default value when the user passes out a failing test.
            self.testResultsText.text.tag_bind('failure%s' % comp,'<Leave>', lambda event : changeMouseCursor(event,'arrow'), '+')
                
            comp += 1
                           
        # The tests results is switched to DISABLED state. Nothing can be written in.
        self.testResultsText.text.config(state = DISABLED)
            
    def infoAboutTest(self, event):
        """This callback displays basic information about a selected test before running it by clicking on 
        the right mouse button.
        """

        selectedInd = self.testListLb.lb.nearest(event.y)
        selectedTest = self.testListLb.lb.get(selectedInd)
        t = eval(selectedTest+'Tests()')
        showinfo('Information',t.info)       
        del(t)
        
    def displayTest(self, event, parameters):
        """Displays in a text widget the nMOLDYN v2.2.5 input file of a selected failing test.
        
        Arguments:
            - event: the mouse left button click Tkinter event that will trigger the opening of the window.
            - testSet: the test set to which belong the failing test.
        """
                
        # The text index of the event.
        first = event.widget.index('@%d,%d' % (event.x, event.y))
        c = self.testResultsText.text.get(first)
                
        try:
            i = int(c)
            
            while c not in  [' ','\n']:
                first = first + '-1c'
                c = self.testResultsText.text.get(first)
                if c == '\n':
                    break

            first = first + '+1c'
            c = self.testResultsText.text.get(first)

            last = first
            while c not in  [' ','\n']:
                last = last + '+1c'
                c = self.testResultsText.text.get(last)
                if c == '\n':
                    break

            selectedTest =  int(self.testResultsText.text.get(first, last).strip())
            
        except:
            return        
        
        # Creation of a toplevel window that will display the reference and current
        # nMOLDYN input files of the failing test.
        self.detailsDialog = PortableToplevel(self)
        self.detailsDialog.transient(self.parent)
        
        # The text combo widget that will display the contents of the reference version of the input file.
        referenceContents = ComboText(self.detailsDialog,\
                                      frameLabel = 'nMOLDYN reference version parameters')
        referenceContents.grid(row = 0, column = 0)
        
        referenceContents.config({'width' : 100,\
                                  'height' : 20,\
                                  'relief' : SUNKEN,\
                                  'bd' : 2})
        referenceContents.text.config({'wrap' : WORD,\
                                       'font' : ('Courier', 10, 'normal'),\
                                       'width' : 80,\
                                       'height' : 10,\
                                       'bg' : 'salmon'})
        
        for k, v in parameters[selectedTest]['REF'].items():
            referenceContents.insert(END, '%s = %s\n' % (repr(k), repr(v)))
        referenceContents.text.config(state = DISABLED)

        # The text combo widget that will display the contents of the current version of the input file.
        currentContents = ComboText(self.detailsDialog,\
                                    frameLabel = 'nMOLDYN version %s parameters' % nMOLDYN_VERSION)
        currentContents.grid(row = 1, column = 0)
        
        currentContents.config({'width' : 100,\
                                'height' : 20,\
                                'relief' : SUNKEN,\
                                'bd' : 2})
        currentContents.text.config({'wrap' : WORD,\
                                     'width' : 80,\
                                     'height' : 10,\
                                     'bg' : '#ffff80',\
                                     'font' : ('Arial', 10, 'normal')})
        
        for k, v in parameters[selectedTest]['NEW'].items():
            currentContents.insert(END, '%s = %s\n' % (repr(k), repr(v)))
        currentContents.text.config(state = DISABLED)

        self.detailsDialog.grab_set()

        self.detailsDialog.bind("<Escape>", self.cancel1)
        self.detailsDialog.protocol("WM_DELETE_WINDOW", self.cancel1)

        self.detailsDialog.resizable(width = NO, height = NO)

        self.detailsDialog.geometry("+%d+%d" % (self.winfo_rootx()+50, self.winfo_rooty()+50))

        self.detailsDialog.title('Details about the test %s.' % selectedTest)

        self.detailsDialog.focus_set()

        self.detailsDialog.wait_window(self.detailsDialog)

    def cancel1(self, event = None):
        
        # Put focus back to the parent window
        self.grab_set()

        # Put focus back to the parent window
        self.focus_set()

        self.detailsDialog.destroy()

    def clearResults(self):
        """
        This methods clears up all the results about the evaluated test sets. 
        """

        self.testResultsText.cleanup()
