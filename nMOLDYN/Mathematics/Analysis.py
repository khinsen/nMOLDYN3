"""This modules contains functions related to general mathematics (analysis, geometry, algebra ...).

Functions:
    * differentiate       : performs a numerical differentiation of 1D Numeric array
    * correlation         : performs the numerical correlation between two 1D Numeric arrays
    * convolution         : performs the numerical convolution between two 1D Numeric arrays
    * FFT                 : performs the FFT of a 1D Numeric array
    * invFFT              : performs the inverse FFT of a 1D Numeric array
    * gaussianWindow      : performs a Gaussian smoothing of 1D Numeric array.
    * factorial           : computes factorial (n) where n is an integer.
"""

# The ScientificPython modules
from Scientific import N
from Scientific.FFT import fft, inverse_fft

# The nMOLDYN modules
from nMOLDYN.Core.Error import Error

###############################################################################
# Analysis
###############################################################################

"""a2 = array used to perform order 2 numerical differentiation scheme."""
a2 = N.array([[   -3.,   4.,  -1.],
              [   -1.,   0.,   1.],
              [    1.,  -4.,   3.]], typecode = N.Float)

"""a3 = array used to perform order 3 numerical differentiation scheme."""
a3 = N.array([[  -11.,  18.,  -9.,   2.],
              [   -2.,  -3.,   6.,  -1.],
              [    1.,  -6.,   3.,   2.],
              [   -2.,   9., -18.,  11.]], typecode = N.Float)

"""a4 = array used to perform order 4 numerical differentiation scheme."""
a4 = N.array([[  -50.,  96., -72.,  32.,  -6.],
              [   -6., -20.,  36., -12.,   2.],
              [    2., -16.,   0.,  16.,  -2.],
              [   -2.,  12., -36.,  20.,   6.],
              [    6., -32.,  72., -96.,  50.]], typecode = N.Float)

"""a5 = N.array used to perform order 5 numerical differentiation scheme."""
a5 = N.array([[  -274.,  600., -600.,  400., -150.,  24.],
              [   -24., -130.,  240., -120.,   40.,  -6.],
              [     6.,  -60.,  -40.,  120.,  -30.,   4.],
              [    -4.,   30., -120.,   40.,   60.,  -6.],
              [     6.,  -40.,  120., -240.,  130.,  24.],
              [   -24.,  150., -400.,  600., -600., 274.]], typecode = N.Float)

def differentiate(inputSeries, order, dx):
    """Returns the numerical derivative of a signal.
    
    @param inputSeries: the signal to differentiate.
    @index inputSeries: NumPy array.
    
    @param order: the numerical differentiation order.
    @index order: integer in [1,5]

    @param dx: the differentiation step. Assumed to be constant over all the spectrum.
    @index dx: float
    
    @return: the differentiated signal.
    @rtype: NumPy array

    @see: M. Abramowitz, I.A. Stegun; 'Handbook of mathematical functions', Dover, New-York, 1972 p.914.
    """

    if len(inputSeries) == 0:
        raise Error('The inputSeries is empty.')

    if dx <= 0.0:
        raise Error('The time step for derivation must strictly positive.')

    # outputSeries is the output resulting from the differentiation
    outputSeries = N.zeros(len(inputSeries), typecode = N.Float)

    if order == 1:
        # Classical order 1
        # The formula used here is the double-sided differentiation where
        # x'(t) = (x(t+1) - x(t-1))/2dt
        # this formula triggers two special indexes cases 0 and len(x) - 1

        # An intermediate factor representing the denominator of the derivation formula
        fact = 1.0/(2.0*dx)

        # array of (x(ti+1) - x(ti)) with i ranging from 0 to len(x) - 2
        gj = inputSeries[1:] - inputSeries[:-1]

        # Special case for the first and last elements
        outputSeries[0] = N.add.reduce(a2[0,:]*inputSeries[:3])*fact
        outputSeries[-1] = N.add.reduce(a2[2,:]*inputSeries[-3:])*fact

        # gj[1:] = array of (x(ti+1) - x(ti)) with i ranging from 1 to len(x) - 2
        # gj[:-1] = array of (x(ti+1) - x(ti)) with i ranging from 0 to len(x) - 3
        # This is the double-sided differentiation formula
        outputSeries[1:-1] = (gj[1:]+gj[:-1])*fact

    elif order == 2:
        # Case of the order 2

        # An intermediate factor representing the denominator of the derivation formula
        fact = 1.0/(2.0*dx)

        # Special case for the first and last elements
        outputSeries[0]  = N.add.reduce(a2[0,:]*inputSeries[:3])*fact
        outputSeries[-1] = N.add.reduce(a2[2,:]*inputSeries[-3:])*fact

        # General case
        gj      = N.zeros((len(inputSeries)-2,3), typecode = N.Float)
        gj[:,0] = a2[1,0]*inputSeries[:-2]
        gj[:,1] = a2[1,1]*inputSeries[1:-1]
        gj[:,2] = a2[1,2]*inputSeries[2:]
        outputSeries[1:-1] = N.add.reduce(gj,-1)*fact

    elif order == 3:
        # Case of the order 3

        # An intermediate factor representing the denominator of the derivation formula
        fact = 1./(6.*dx)

        # Special case for the first and last elements
        outputSeries[0]  = N.add.reduce(a3[0,:]*inputSeries[:4])*fact
        outputSeries[1]  = N.add.reduce(a3[1,:]*inputSeries[:4])*fact
        outputSeries[-1] = N.add.reduce(a3[3,:]*inputSeries[-4:])*fact

        # General case
        gj      = N.zeros((len(inputSeries)-3,4), typecode = N.Float)
        gj[:,0] = a3[2,0]*inputSeries[:-3]
        gj[:,1] = a3[2,1]*inputSeries[1:-2]
        gj[:,2] = a3[2,2]*inputSeries[2:-1]
        gj[:,3] = a3[2,3]*inputSeries[3:]
        outputSeries[2:-1] = N.add.reduce(gj,-1)*fact

    elif order == 4:
        # Case of the order 4

        # An intermediate factor representing the denominator of the derivation formula
        fact = 1./(24.*dx)

        # Special case for the first and last elements
        outputSeries[0]  = N.add.reduce(a4[0,:]*inputSeries[:5])*fact
        outputSeries[1]  = N.add.reduce(a4[1,:]*inputSeries[:5])*fact
        outputSeries[-2] = N.add.reduce(a4[3,:]*inputSeries[-5:])*fact
        outputSeries[-1] = N.add.reduce(a4[4,:]*inputSeries[-5:])*fact

        # General case
        gj      = N.zeros((len(inputSeries)-4,5), typecode = N.Float)
        gj[:,0] = a4[2,0]*inputSeries[:-4]
        gj[:,1] = a4[2,1]*inputSeries[1:-3]
        gj[:,2] = a4[2,2]*inputSeries[2:-2]
        gj[:,3] = a4[2,3]*inputSeries[3:-1]
        gj[:,4] = a4[2,4]*inputSeries[4:]
        outputSeries[2:-2] = N.add.reduce(gj,-1)*fact

    elif order == 5:
        # Case of the order 5

        # An intermediate factor representing the denominator of the derivation formula
        fact = 1./(120.*dx)

        # Special case for the first and last elements
        outputSeries[0]  = N.add.reduce(a5[0,:]*inputSeries[:6])*fact
        outputSeries[1]  = N.add.reduce(a5[1,:]*inputSeries[:6])*fact
        outputSeries[2]  = N.add.reduce(a5[2,:]*inputSeries[:6])*fact
        outputSeries[-2] = N.add.reduce(a5[4,:]*inputSeries[-6:])*fact
        outputSeries[-1] = N.add.reduce(a5[5,:]*inputSeries[-6:])*fact

        # General case
        gj      = N.zeros((len(inputSeries)-5,6), typecode = N.Float)
        gj[:,0] = a5[3,0]*inputSeries[:-5]
        gj[:,1] = a5[3,1]*inputSeries[1:-4]
        gj[:,2] = a5[3,2]*inputSeries[2:-3]
        gj[:,3] = a5[3,3]*inputSeries[3:-2]
        gj[:,4] = a5[3,4]*inputSeries[4:-1]
        gj[:,5] = a5[3,5]*inputSeries[5:]
        outputSeries[3:-2] = N.add.reduce(gj,-1)*fact

    else:
        raise Error('Unknown differentiation order.')

    return outputSeries

# Just for testing purpose.
def myautocorrelation(inputSeries):
    
    nTimes = len(inputSeries)
    corr = N.zeros((nTimes,), typecode = N.Float)
    for i in range(nTimes):
        s = 0.0
        for j in range(nTimes-i):
            try:
                s += sum(inputSeries[j]*inputSeries[j+i])
            except:
                s += inputSeries[j]*inputSeries[j+i]
        corr[i] = s/float(nTimes-i)
                    
    return corr            

def correlation(inputSeries1, inputSeries2 = None):
    """Returns the numerical correlation between two signals.

    @param inputSeries1: the first signal.
    @type inputSeries1: NumPy array   
    
    @param inputSeries2: if not None, the second signal otherwise the correlation will be an autocorrelation.
    @type inputSeries2: NumPy array or None
    
    @return: the result of the numerical correlation.
    @rtype: NumPy array

    @note: if |inputSeries1| is a multidimensional array the correlation calculation is performed on
    the first dimension.

    @note: The correlation is computed using the FCA algorithm.
    """

    # The signal must not be empty.
    if len(inputSeries1) <= 0:
        raise Error('One or both time series are empty.')

    # The length of inputSeries1 is stored in inputSeries1Length
    inputSeries1Length = len(inputSeries1)

    # extendedLength = 2*len(inputSeries1)
    extendedLength = 2*inputSeries1Length

    # The FCA algorithm:

    # 1) computation of the FFT of inputSeries1 zero-padded until extendedLength
    # The computation is done along the 0-axis
    FFTSeries1 = fft(inputSeries1,extendedLength,0)        

    if inputSeries2 is None:
            # Autocorrelation case
        FFTSeries2 = FFTSeries1
    else:
        # 2) computation of the FFT of inputSeries2 zero-padded until extendedLength
        # The computation is  done along the 0-axis
        FFTSeries2 = fft(inputSeries2,extendedLength,0)

    # 3) Product between FFT(inputSeries1)* and FFT(inputSeries2)
    FFTSeries1 = N.conjugate(FFTSeries1)*FFTSeries2

    # 4) inverse FFT of the product
    # The computation is done along the 0-axis
    FFTSeries1 = inverse_fft(FFTSeries1,len(FFTSeries1),0)

    # This refers to (1/(N-m))*Sab in the published algorithm.
    # This is the correlation function defined for positive indexes only.
    if len(FFTSeries1.shape) == 1:
        corr = FFTSeries1.real[:inputSeries1Length] / (inputSeries1Length-N.arange(inputSeries1Length))
    else:
        corr = N.add.reduce(FFTSeries1.real[:inputSeries1Length],1) / (inputSeries1Length-N.arange(inputSeries1Length))

    return corr

def convolution(inputSeries1, inputSeries2):
    """Returns the numerical convolution between two signals.

    @param inputSeries1: the first signal   
    @type inputSeries1: NumPy array   
    
    @param inputSeries2: the second signal.
    @type inputSeries2: NumPy array
    
    @return: the result of the convolution.
    @rtype: NumPy array

    @note: if |inputSeries1| is a multidimensional array the convolution calculation is performed on
    the first dimension.

    @note: the convolution is computed using the convolve function of NumPy package.
    """
    
    # The signals must not be empty.
    if (len(inputSeries1) == 0) | (len(inputSeries2) == 0):
        raise Error('One or both inputSeries are empty.')

    # The convolution between inputSeries1 and inputSeries2.
    # A 1D Numeric array of dimension len(inputSeries1)*len(inputSeries2)-1
    return N.convolve(inputSeries1,inputSeries2)

def FFT(inputSeries):
    """Returns the FFT of a signal.

    @param inputSeries: the signal.
    @type inputSeries: NumPy array   

    @return: the FFT transformed signal.
    @rtype: NumPy array

    @note: the FFT is computed using the fft function of Scientific.FFT package.
    """

    # The signal must not be empty.
    if len(inputSeries) == 0:
        raise Error('The inputSeries is empty.')

    # The FFT of inputSeries
    return fft(inputSeries,len(inputSeries),0)

def invFFT(inputSeries):
    """Returns the inverse FFT of a signal.

    @param inputSeries: the signal.
    @type inputSeries: NumPy array   

    @return: the inverse FFT transformed signal.
    @rtype: NumPy array

    @note: the inverse FFT is computed using the inverse_fft function of Scientific.FFT package.
    """

    # The signal must not be empty.
    if len(inputSeries) == 0:
        raise Error('The inputSeries is empty.')

    # The inverse FFT of inputSeries
    return inverse_fft(inputSeries,len(inputSeries),0)

def gaussian(x, sigma = 1.0, mu = 0.0, normalize = False):
    """Returns a gaussian.
    """
                    
    gauss_fun = N.exp(-0.5*( (x - mu) / sigma)**2)
        
    if normalize:
    
        gauss_fun *= 1.0/(sigma*N.sqrt(2.0*N.pi))
        
    return gauss_fun

def gaussianWindow(inputSeries, gaussian):
    """Returns a gaussian smoothed signal.
    
    @param inputSeries: the signal to smooth.
    @type inputSeries: NumPy array   

    @param alpha: a float specifying the width of the smoothing gaussian.
    @type alpha: float

    @return: the smoothed signal.
    @rtype: NumPy array
    """

    # outputSeries is an array of length 2*len(series)-1 directionized to zero 
    # outputSeries is the smoothed version of inputSeries obtained by applying
    # a gaussian kernel to intputSeries
    outputSeries = N.zeros((2*len(inputSeries) - 2,), typecode = N.Float)
            
    # exp(...) is the gaussian window used to smooth the spectrum series 
    res = inputSeries*gaussian

    # The second half of outputSeries is filled using periodic conditions
    outputSeries[:len(inputSeries)] = res
    outputSeries[len(inputSeries):] = res[-2:0:-1]

    return outputSeries

def factorial(n):
    """Returns n!
    
    @param n: n.
    @type: integer
    
    @return: n!.
    @rtype: integer
    """
    
    # 0! = 1! = 1
    if n < 2:
        return 1
    else:
        # multiply is a ufunc of Numeric module
        return N.multiply.reduce(N.arange(2,n+1, typecode = N.Float))

def pgcd(n):
    """Computes the pgcd for a set of integers.
    
    @param n: n.
    @type: list of integers
    
    @return: pgcd([i1,i2,i3...]).
    @rtype: integer
    """   
    
    def _pgcd(a,b):
        while b: a, b = b, a%b
        return a
    
    p = _pgcd(n[0],n[1])
    
    for x in n[2:]:
        p = _pgcd(p, x)
        
    return p
    