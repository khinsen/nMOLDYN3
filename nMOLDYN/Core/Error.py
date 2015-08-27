# The python distribution modules
import os
import sys
import traceback

# The nMOLDYN modules
from nMOLDYN.Core.Logger import LogMessage

class Error(Exception):
    """A subclass of Exception. Used to handles nMOLDYN internal error messages.
    """

    def __init__(self, extra = None):
        """
        The constructor.
                
        @param message: the error message.
        @type message: string        
        """
                
        # Store the information about the most recent exception caught by an except if any.
        self.type, self.value, self.traceback = sys.exc_info()
        
        # Store the extra message provided by the usert.
        self.extra = extra

        # The default handlers on which the error will be displayed.
        self.handlers = ['file', 'console']
                                                                                                
        msg = "\n"*2 + "!"*90 + "\n"
        msg += "\nnMOLDYN ERROR\n\n"
                
        if self.type is not None:
            msg += "\n\t".join(traceback.format_exception(self.type, self.value, self.traceback))
            
        if self.extra is not None:
            msg += "\nAdditional message:\n\t%s\n" % self.extra
            
        msg += '\nPlease look in the log file for more details.'
        msg += "\n\n" + "!"*90

        # Send the error message into the selected handlers.
        LogMessage('error', msg, self.handlers)
                                