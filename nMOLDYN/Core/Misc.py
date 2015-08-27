"""This modules implements the functions and procedures that can not be classified anywhere else in 
the library.

There should not be too much stuff here in order to not mess up the code.

Functions:
    * convertTime   : converts a duration in seconds in days, hours, minutes and seconds.
    * parseInterval : parse an interval of the form "xmin:xmax:dx".
"""

from numpy import linspace

from nMOLDYN.Core.Error import Error

def convertTime(time):
    """Converts a duration in seconds in days, hours, minutes and seconds.
        
    @param time: the time in seconds to convert.
    @type time: integer.
        
    @return: a dictionnary of the form {'days' : d, 'hours' : h, 'minutes' : m, 'seconds' : s}
        where d, h, m and s are integers resulting from the time conversion respectively in days, hours, minutes
        and seconds.
    @rtype: dict           
    """

    # If no time is provided, the results of the conversion is None.
    if time is None:
        converted = None

    # Otherwise do the conversion.
    else:
        converted = {'days' : 0, 'hours' : 0, 'minutes' : 0, 'seconds' : 0}

        # These are the formulas to convert a time in seconds to a time in days, hours, minutes and seconds.
        converted['days'], time = divmod(time, 86400)
        converted['hours'], time = divmod(time, 3600)
        converted['minutes'], converted['seconds'] = divmod(time, 60)

    return converted   

def parseInterval(interval):
    """This function convert a string of the form xmin:xmax:dx to an array
    
    @param interval: a string of the form 'mini:maxi:delta' to convert.
    @type interval: string
    
    @return: a list of the form [xmin, xmin + dx, xmin + 2*dx ...]
    @rtype: list
    """

    if not isinstance(interval, str):
        raise Error('Error when parsing %s parameter: wrong format for the string.' % pName)
    
    # The string must be three semi colon-separated floats. 
    try:
        xMin, xMax, dX = [float(v) for v in interval.split(':')]
                            
    except:
        raise Error('Error when parsing % string: must be of the form "xmin:xmax:dx".' % interval)
                    
    else:
        
        # The number of bins corresponding to the given interval.
        nBins = int((xMax - xMin)/dX)
        
        # The actual value of xMax is recomputed.
        xMax = xMin + nBins * dX
        
        values = linspace(xMin, xMax, nBins+1, True)
        
        values.sort()
                
        return values
    
    