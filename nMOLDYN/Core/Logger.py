"""This module implements the classes and procedures used to handle the nMOLDYN logging system.
Three kinds of loggers are availble in nMOLDYN: a console logger, a file logger and a gui logger.

Classes:
    * LogToGUI     : sets up the GUI logger.
    * LogToFile    : sets up the file logger.
    * LogToConsole : sets up the console logger.
    
Procedures:
    * LogMessage   : displays a logging message to the specified logger(s) with a specific logging level.
"""

# The python distribution modules
import logging
import os
import platform
from time import strftime

if platform.system().upper() == "WINDOWS":
    LOGFILE_PATH = os.environ["USERPROFILE"]
else:
    LOGFILE_PATH = os.environ["HOME"]

LOGGER_NAME = "nMOLDYN_%s.log" % strftime('%a%d%b%Y_%Hh%Mm%Ss')

# The nMOLDYN file logger is set up.
FILE_LOGGER = logging.getLogger('LOGFILE')
# With a minimum logging level.
FILE_LOGGER.setLevel(logging.DEBUG)
        
# The nMOLDYN console logger is set up.
CONSOLE_LOGGER = logging.getLogger('CONSOLE')
# With a minimum logging level. The message with a 'DEBUG' level will not be displayed on the console.
CONSOLE_LOGGER.setLevel(logging.INFO)
       
# The nMOLDYN GUI logger is set up.
GUI_LOGGER = logging.getLogger('GUI')
# With a minimum logging level. The message with a 'DEBUG' level will not be displayed on the console.
GUI_LOGGER.setLevel(logging.INFO)

def LogCreate(loggers):
    
    loggers = [l.lower() for l in loggers]
            
    for l in loggers:
        if l == "console":
            CONSOLE_LOGGER.addHandler(LogToConsole())
            CONSOLE_LOGGER.info(LOGGER_NAME)
            
        elif l == "file":
            logFileName = os.path.join(LOGFILE_PATH, LOGGER_NAME)
            FILE_LOGGER.addHandler(LogToFile(logFileName))
            
        elif l == "gui":
            GUI_LOGGER.addHandler(LogToGUI())

class LogToGUI(logging.Handler):
    """Sets up a GUI handler for the nMOLDYN logger.
    
    Emits the logging messages to a Tk dialog.

    @note: inherits from logging.Handler class that sets a generic handler.
    """
    
    def __init__(self):
        """The constructor.
        
        Sets the logger.
        """
        
        # Sets up a generic handler.
        logging.Handler.__init__(self)
        
        # Configure the log message for the GUI handler.
        self.setFormatter(logging.Formatter('nMOLDYN %(levelname)s %(filename)s %(lineno)s \n%(message)s'))

    def emit(self, record):
        """Emits the logging message in a tkMessageBox.
        
        @param record: the logging message.
        @type record: instance of LogRecord class.
        
        @note: the message box called will depend on the logging level:
            - tkMessageBox.showerror for 'ERROR' and 'CRITICAL' logging levels.
            - tkMessageBox.showwarning for 'WARNING' logging level.
            - tkMessageBox.showinfo for other logging levels.

        @note: this is the overloaded emit method of Handler objects.
        """

        if not GUI_LOGGER.handlers:
            return

        # Try to import the compiled tkinter module.
        try:
            import _tkinter
        except ImportError:    
            return
        else:
            import tkMessageBox
            del _tkinter        
                
        # Case of 'ERROR' and 'CRITICAL' logging levels.
        if record.levelname in ['ERROR','CRITICAL']:
            # An error message box is opened and display the message.
            tkMessageBox.showerror(record.levelname.capitalize(), self.format(record))

        # Case of 'WARNING' logging level.
        elif record.levelname == 'WARNING':
            # An warning message box is opened and display the message.
            tkMessageBox.showwarning(record.levelname.capitalize(), self.format(record))

        # Case of 'INFO' logging level.
        elif record.levelname in ['INFO','DEBUG']:
            # An warning message box is opened and display the message.
            tkMessageBox.showinfo(record.levelname.capitalize(), self.format(record))
                        
class LogToFile(logging.FileHandler):
    """Sets up a file logger.

    Emits the logging messages to the nMOLDYN log file.
    
    @note: inherits from logging.FileHandler class that sends logging output to a file.
    """
    
    def __init__(self, loggerName):
        """The constructor. Sets the logger.
        
        @param loggerName: the name of the file where all the logging messages will be emitted.
        @type loggerName: string        
        """

        # The file handler is set up using |filename| file.
        logging.FileHandler.__init__(self, loggerName)

    def emit(self, record):
        """Emits the logging message in a file.
        
        @param record: the logging message.
        @type record: instance of LogRecord class.
        
        @note: this is the overloaded emit method of Handler objects.
        """
        
        # The arrow used to decorate the messages. Its size depends on the logging level.
        arrow = '-'*(12 - len(record.levelname)) + '>'
        
        self.setFormatter(logging.Formatter('nMOLDYN %(levelname)s ' + arrow + ' %(message)s'))            
                    
        # The formatted record is written to the log file.
        self.stream.write(self.format(record)+'\n')
        
        # And flushed.
        self.stream.flush()

    def close(self):
        """Closes the file logger."""
        self.stream.close()

class LogToConsole(logging.StreamHandler):
    """Sets up a console logger.

    Emits the logging messages to the console.

    @note: inherits from logging.StreamHandler class that sends logging output to a file-like objects.
    """

    def __init__(self):
        """The constructor.
        
        Sets the logger.
        """
        
        # The stream handler is set up.
        logging.StreamHandler.__init__(self)
                
    def emit(self, record):
        """Emits the logging message to the console.
        
        @param record: the logging message.
        @type record: instance of LogRecord class.

        @note: this is the overloaded emit method of Handler objects.
        """

        # The arrow used to decorate the messages. Its size depends on the logging level.
        arrow = '-'*(12 - len(record.levelname)) + '>'

        self.setFormatter(logging.Formatter('nMOLDYN %(levelname)s ' + arrow + ' %(message)s'))

        # The formatted record is written to the console.
        self.stream.write(self.format(record)+'\n')
        
        # And flushed.
        self.stream.flush()
                
def LogMessage(level = 'debug', message = '', logger = None):
    """Displays a logging message.
    
    The logging message can be sent on one or simultaneoulsy several loggers being a console logger, a file logger or a gui logger.    
    
    @param level: the logging level of the message. This will change slightly the way the logging message will be 
        displayed.
    @type level: string being one of 'debug', 'info', 'warning', 'error' or 'critical'
    
    @param message: the message.
    @type message: string
    
    @param logger: a list containing 'gui' and/or 'file' and/or 'console' specifying on which logger(s) the logging
        message should emitted.
    @type logger: list
    """
    
    if logger is None:
        logger = ['gui', 'file', 'console']
            
    # Check whether the logging level is supported by logging module.
    try:
        loggingLevel = getattr(logging,level.upper())
        
    # If not, set it arbitrary to 'INFO' level.
    except AttributeError:
        loggingLevel = logging.INFO
    
    # A local copy of of |logger| argument.
    logger = [s.lower() for s in logger]
    
    for l in logger:
        try:
            if not NMOLDYN_LOGGERS[l].handlers:
                continue
            NMOLDYN_LOGGERS[l].log(loggingLevel, message)
        except:
            continue    
    
NMOLDYN_LOGGERS = {'console' : CONSOLE_LOGGER, 'file' : FILE_LOGGER, 'gui' : GUI_LOGGER}
        