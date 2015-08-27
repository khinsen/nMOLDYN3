# Autoregressive Model in arbitrary precision
#
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# last revision: 2002-12-3
#

from Scientific import N
from Scientific.Functions.Interpolation import InterpolatingFunction
from Scientific.Functions.Rational import RationalFunction

from gmpy import mpf

class MPAutoRegressiveModel:

    """Auto-regressive model for stochastic process

    This is an arbitrary-precision implementation that permits the
    calculation of memory functions at high orders. A standard precision
    model is used as input. The GMP library is used for high-precision
    computations.

    Constructor: MPAutoRegressiveModel(|model|, |precision|=100)

    Arguments:

    |model| -- a standard precision model

    |precision| -- the number of bits in the floating-point representation
    """

    def __init__(self, model, precision=0):
        self.precision = precision
        self.order = model.order
        self.delta_t = model.delta_t
        if model.coeff.typecode() == N.Complex64:
            self.coeff = [N.Complex(mpf(x.real, precision),
                                  mpf(x.imag, precision)) for x in model.coeff]
        else:
            self.coeff = [mpf(x, precision) for x in model.coeff]
        self.sigsq = mpf(model.sigsq, precision)
        self.sigma = mpf(model.sigma, precision)
        self.variance = mpf(model.variance, precision)
        self._poles = None

    def poles(self):
        if self._poles is None:
            c = map(lambda x: -x, self.coeff)
            c.append(mpf(1., self.precision))
            self._poles = Polynomial(c, self.precision).zeros()
        return self._poles

    def memoryFunctionZ(self):
        poles = self.poles()
        cpoles = [p.conjugate() for p in poles]
        coeff0 = conjugate(self.coeff[0])
        beta = self.order*[0]
        sum = 0.
        for i in range(self.order):
            pole = poles[i]
            prod = N.Complex(self.variance, 0.)
            for j in range(i):
                prod *= (pole-poles[j])
            for j in range(i+1, self.order):
                prod *= (pole-poles[j])
            for j in range(self.order):
                prod *= (pole-N.Complex(1., 0.)/cpoles[j])
            beta[i] = -((pole**(self.order-1))*self.sigsq/coeff0) / prod
            sum += beta[i]
        for i in range(self.order):
            beta[i] /= sum
        sum = 0.
        for i in range(self.order):
            sum += RationalFunction([beta[i]], [-poles[i], 1.])
        mz = (1./sum+Polynomial([1., -1.]))/self.delta_t**2
        if not isComplex(self.coeff[0]):
            mz.numerator.coeff = [c.real for c in mz.numerator.coeff]
            mz.denominator.coeff = [c.real for c in mz.denominator.coeff]
        return mz

    def memoryFunction(self, nsteps):
        mz = self.memoryFunctionZ()
        mem = mz.divide(nsteps-1)[0].coeff[:]
        mem.reverse()
        if len(mem) == nsteps+1:
            mem = mem[1:]
        if isComplex(self.coeff[0]):
            mem = N.array([complex(m.real, m.imag) for m in mem])
        else:
            mem = N.array([float(m) for m in mem])
        mem[0] = 2.*realPart(mem[0])
        time = self.delta_t*N.arange(nsteps)
        return InterpolatingFunction((time,), mem)


class Complex:

    def __init__(self, real, imag=0.):
        self.real = real
        self.imag = imag

    is_complex = 1

    def __repr__(self):
        return "Complex(%s, %s)" % (repr(self.real), repr(self.imag))

    def __str__(self):
        return "Complex(%s, %s)" % (str(self.real), str(self.imag))

    def conjugate(self):
        return Complex(self.real, -self.imag)

    def __int__(self):
        return int(self.real)

    def __coerce__(self, other):
        if hasattr(other, 'is_complex'):
            return self, other
        else:
            return self, Complex(other, 0.)

    def __add__(self, other):
        return Complex(self.real+other.real, self.imag+other.imag)
    __radd__ = __add__

    def __sub__(self, other):
        return Complex(self.real-other.real, self.imag-other.imag)
    def __rsub__(self, other):
        return Complex(other.real-self.real, other.imag-self.imag)

    def __neg__(self):
        return Complex(-self.real, -self.imag)

    def __mul__(self, other):
        return Complex(self.real*other.real-self.imag*other.imag,
                       self.real*other.imag+self.imag*other.real)
    __rmul__ = __mul__

    def __div__(self, other):
        oc = other.conjugate()
        prod = self*oc
        den = other*oc
        return Complex(prod.real/den.real, prod.imag/den.real)
    def __rdiv__(other, self):
        oc = other.conjugate()
        prod = self*oc
        den = other*oc
        return Complex(prod.real/den.real, prod.imag/den.real)

    def __pow__(self, n):
        p = Complex(1., 0.)
        for i in range(n):
            p = p*self
        return p

    def __abs__(self):
        return (self.real**2+self.imag**2).sqrt()

    def __cmp__(self, other):
        return cmp(self.real, other.real)**2 + cmp(self.imag, other.imag)**2

    def sqrt(self):
        if self.real == 0 and self.imag == 0:
            return self
        s = (0.5*(abs(self.real) + (self.real**2+self.imag**2).sqrt())).sqrt()
        d = 0.5*self.imag/s
        if self.real > 0:
            return Complex(s, d)
        elif self.imag >= 0:
            return Complex(d, s)
        else:
            return Complex(-d, -s)


class Polynomial:

    def __init__(self, coefficients, precision=0):
	self.coeff = coefficients
	self.dim = 1
        self.precision = precision

    is_polynomial = 1

    def __call__(self, x):
        p = 1.
        sum = 0.
        for c in self.coeff:
            sum += c*p
            p *= x
	return sum

    def __repr__(self):
        return "Polynomial(%s)" % repr(self.coeff)

    def __coerce__(self, other):
        if hasattr(other, 'is_polynomial'):
            return (self, other)
        elif hasattr(other, 'is_rational_function'):
            return None
        else:
            return (self, Polynomial([other]))

    def __add__(self, other):
        coeff = max(len(self.coeff), len(other.coeff))*[0]
        coeff[:len(self.coeff)] = self.coeff
        for i in range(len(other.coeff)):
            coeff[i] += other.coeff[i]
        return Polynomial(coeff)

    def __mul__(self, other):
        coeff = (len(self.coeff)+len(other.coeff)-1)*[0]
        for i in range(len(self.coeff)):
            for j in range(len(other.coeff)):
                coeff[i+j] += self.coeff[i]*other.coeff[j]
        return Polynomial(coeff)

    def __div__(self, other):
        if self.dim != 1 or other.dim != 1:
            raise ValueError, "not implemented"
        return RationalFunction(self, other)

    def __rdiv__(self, other):
        return RationalFunction(other, self)

    def zeros(self):
        order = len(self.coeff)-1
        coeff = self.coeff[:]
        roots = []
        for i in range(order):
            r = self._laguerre(coeff, Complex(mpf(0, self.precision),
                                              mpf(0, self.precision)))
            roots.append(r)
            if i == order-1: break
            rem = coeff[-1]
            for j in range(len(coeff)-2, -1, -1):
                temp = coeff[j]
                coeff[j] = rem
                rem = temp + rem*r
            coeff = coeff[:-1]
##          for i in range(order):
##              roots[i] = self._laguerre(self.coeff, roots[i])
        return roots

    def _laguerre(self, coeff, x, max_iter = 1000):
        order = len(coeff)-1
        for i in range(max_iter):
            p = coeff[-1]
            p1 = 0
            p2 = 0
            err = abs(p)
            absx = abs(x)
            for n in range(order-1, -1, -1):
                p2 = x*p2 + p1
                p1 = x*p1 + p
                p  = x*p  + coeff[n]
                err = absx*err + abs(p)
            if abs(p) <= mpf('1.e-30', self.precision)*err:
                return x
            g = p1/p
            h = g*g-p2/p
            sq = ((order-1)*(order*h-g*g)).sqrt()
            gp = g + sq
            gm = g - sq
            if abs(gp) < abs(gm):
                dx = order/gm
            else:
                dx = order/gp
            xnew = x - dx
            if xnew == x:
                return x
            if i % 10 == 9:
                x = x - (i/10)*dx
            else:
                x = xnew
        raise ValueError, "too many iterations"


class RationalFunction:

    def __init__(self, numerator, denominator=[1.]):
        if hasattr(numerator, 'is_polynomial'):
            self.numerator = numerator
        else:
            self.numerator = Polynomial(numerator)
        if hasattr(denominator, 'is_polynomial'):
            self.denominator = denominator
        else:
            self.denominator = Polynomial(denominator)
        self._normalize()

    is_rational_function = 1

    def __call__(self, value):
        return self.numerator(value)/self.denominator(value)

    def __repr__(self):
        return "RationalFunction(%s,%s)" % (repr(self.numerator.coeff),
                                            repr(self.denominator.coeff))

    def _normalize(self):
        self.numerator = self._truncate(self.numerator)
        self.denominator = self._truncate(self.denominator)
        n = 0
        while 1:
            if self.numerator.coeff[n] != 0:
                break
            if self.denominator.coeff[n] != 0:
                break
            n = n + 1
            if n == len(self.numerator.coeff):
                break
        if n > 0:
            self.numerator = Polynomial(self.numerator.coeff[n:])
            self.denominator = Polynomial(self.denominator.coeff[n:])
        factor = self.denominator.coeff[-1]
        if factor != 1.:
            for i in range(len(self.numerator.coeff)):
                self.numerator.coeff[i] /= factor
            for i in range(len(self.denominator.coeff)):
                self.denominator.coeff[i] /= factor

    def _truncate(self, poly):
        if poly.coeff[-1] != 0.:
            return poly
        coeff = poly.coeff
        while len(coeff) > 1 and coeff[-1] == 0.:
            coeff = coeff[:-1]
        return Polynomial(coeff)

    def __coerce__(self, other):
        if hasattr(other, 'is_rational_function'):
            return (self, other)
        elif hasattr(other, 'is_polynomial'):
            return (self, RationalFunction(other, [1.]))
        else:
            return (self, RationalFunction([other], [1.]))

    def __mul__(self, other):
        return RationalFunction(self.numerator*other.numerator,
                                self.denominator*other.denominator)
    __rmul__ = __mul__

    def __div__(self, other):
        return RationalFunction(self.numerator*other.denominator,
                                self.denominator*other.numerator)

    def __rdiv__(self, other):
        return RationalFunction(other.numerator*self.denominator,
                                other.denominator*self.numerator)

    def __add__(self, other):
        return RationalFunction(self.numerator*other.denominator+
                                self.denominator*other.numerator,
                                self.denominator*other.denominator)
    __radd__ = __add__

    def __sub__(self, other):
        return RationalFunction(self.numerator*other.denominator-
                                self.denominator*other.numerator,
                                self.denominator*other.denominator)

    def __rsub__(self, other):
        return RationalFunction(other.numerator*self.denominator-
                                other.denominator*self.numerator,
                                self.denominator*other.denominator)

    def divide(self, shift=0):
        num = self.numerator.coeff[:]
        if shift > 0:
            num = shift*[0.] + num
        den = self.denominator.coeff
        den_order = len(den)
        coeff = []
        while len(num) >= den_order:
            q = num[-1]/den[-1]
            coeff.append(q)
            for i in range(den_order):
                num[len(num)-den_order+i] -= q*den[i]
            num = num[:-1]
        if not coeff:
            coeff = [0]
        coeff.reverse()
        if not num:
            num = [0]
        return Polynomial(coeff), RationalFunction(num, den)

    def zeros(self):
        "Returns an array containing the zeros."
        return self.numerator.zeros()

    def poles(self):
        "Returns an array containing the poles."
        return self.denominator.zeros()


# Check if data is complex
def isComplex(x):
    try:
        x.imag
        return 1
    except (AttributeError, ValueError):
        return 0

# Return real part
def realPart(x):
    try:
        return x.real
    except (AttributeError, ValueError):
        return x

# Return conjugate
def conjugate(x):
    try:
        return x.conjugate()
    except AttributeError:
        return x
