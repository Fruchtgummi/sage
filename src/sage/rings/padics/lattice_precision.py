import _weakref as weakref
from sage.misc.misc import walltime

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity

from sage.rings.padics.precision_error import PrecisionError


# Global variables

DEFAULT_THRESOLD_DELETION = 50
STARTING_ADDITIONAL_PREC = 5



# Class pRational
#################

class pRational:
    """
    This class implements rational numbers as approximations
    of `p`-adic numbers.

    Only for internal use.
    """
    def __init__(self, p, x, exponent=0, valuation=None):
        """
        Construct the element ``x * p^exponent``

        INPUT:

        - ``p`` -- a prime number

        - ``x`` -- a rational number

        - ``exponent`` -- a relative integer (default: 0)

        - ``valuation`` -- an integer or None (default: ``None``),
          the ``p``-adic valuation of this element

        If not ``None``, this method trusts the given value to the
        attribute ``valuation``.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: pRational(2, 5)
            5
            sage: pRational(2, 5/3, 2)
            2^2 * 5/3
        """
        self.p = p
        if x in ZZ:
            self.x = ZZ(x)
        else:
            self.x = x
        self.exponent = exponent
        self._valuation = valuation

    def __repr__(self):
        """
        Return a string representation of this element

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: pRational(2, 5, 2)  # indirect doctest
            2^2 * 5
        """
        if self.exponent == 0:
            return str(self.x)
        else:
            return "%s^%s * %s" % (self.p, self.exponent, self.x)

    def reduce(self, prec):
        """
        Return this element reduced modulo ``p^prec``

        INPUT:

        - ``prec`` -- a relative integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 1234567); x
            1234567
            sage: x.reduce(12)
            1671

            sage: x = pRational(2, 1234/567); x
            1234/567
            sage: x.reduce(12)
            190
        """
        if prec is Infinity:
            return self
        x = self.x
        exp = self.exponent
        if x.parent() is ZZ:
            if prec > exp:
                x = x % (self.p ** (prec-exp))
            else:
                x = 0
        elif x.parent() is QQ:
            num = x.numerator()
            denom = x.denominator()
            valdenom = denom.valuation(self.p)
            denom /= self.p ** valdenom
            exp -= valdenom
            modulo = self.p ** (prec - exp)
            # probably we should use Newton iteration instead 
            # (but it is actually slower for now - Python implementation)
            _, inv, _ = denom.xgcd(modulo)
            x = (num*inv) % modulo
        if self.x == 0:
            val = Infinity
        else:
            val = self._valuation
        return self.__class__(self.p, x, exp, valuation=val)

    def reduce_relative(self, prec):
        """
        Return this element reduced modulo ``p^n`` where ``n = prec + val(x)``

        INPUT:

        - ``prec`` -- a nonnegative integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 1234567); x
            1234567
            sage: x.reduce_relative(12)
            1671

            sage: x = pRational(2, 1234/567); x
            1234/567
            sage: x.reduce_relative(12)
            190
        """
        v = self.valuation()
        if v is Infinity:
            return self
        return self.reduce(prec+v)

    def normalize(self):
        """
        Normalize this element, i.e. write it as ``p^v * u`` where 
        ``u`` is coprime to ``p``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.normalize(); x
            2^13 * 1929
        """
        if self.x == 0:
            self.exponent = 0
        else:
            val = self.valuation()
            exp = self.exponent
            self.x /= self.p ** (val-exp)
            if self.x in ZZ:
                self.x = ZZ(self.x)
            self.exponent = val

    def valuation(self):
        """
        Return the ``p``-adic valuation of this element

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.valuation()
            13
        """
        if self._valuation is None:
            valx = self.x.valuation(self.p)
            self._valuation = self.exponent + valx
        return self._valuation

    def is_p_power(self):
        """
        Return true if this element is a power of ``p``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 1024, 2); x
            2^2 * 1024
            sage: x.is_p_power()
            True

            sage: y = pRational(2, 123456, 7); y
            2^7 * 123456
            sage: y.is_p_power()
            False
        """
        self.normalize()
        return self.x == 1

    def is_zero(self):
        """
        Return true if this element vanishes

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.is_zero()
            False

            sage: (x-x).is_zero()
            True
        """
        return self.x == 0

    def __add__(self, other):
        """
        Return the sum of ``self`` and ``other``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x + y
            2^7 * 28635808
        """
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if self._valuation is None or other._valuation is None:
            val = None
        elif self._valuation < other._valuation:
            val = self._valuation
        elif self._valuation > other._valuation:
            val = other._valuation
        else:
            val = None
        if sexp < oexp:
            return self.__class__(p, self.x + other.x * p**(oexp-sexp), sexp, valuation=val)
        else:
            return self.__class__(p, self.x * p**(sexp-oexp) + other.x, oexp, valuation=val)

    def __sub__(self, other):
        """
        Return the subtraction of ``self`` by ``other``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x - y
            2^7 * -28388896
        """
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if self._valuation is None or other._valuation is None:
            val = None
        elif self._valuation < other._valuation:
            val = self._valuation
        elif self._valuation > other._valuation:
            val = other._valuation
        else:
            val = None
        if sexp < oexp:
            return self.__class__(p, self.x - other.x * p**(oexp-sexp), sexp, valuation=val)
        else:
            return self.__class__(p, self.x * p**(sexp-oexp) - other.x, oexp, valuation=val)

    def __neg__(self):
        """
        Return the opposite of this element

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: -x
            2^7 * -123456
        """
        return self.__class__(self.p, -self.x, self.exponent, valuation=self._valuation)

    def __mul__(self, other):
        """
        Return the product of ``self`` and ``other``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x * y
            2^19 * 110000654016
        """
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = self._valuation + other._valuation
        return self.__class__(self.p, self.x * other.x, self.exponent + other.exponent, valuation=val)

    def __div__(self, other):
        """
        Return the quotient of ``self`` by ``other``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011
            sage: x / y
            2^-5 * 123456/891011
        """
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = self._valuation - other._valuation
        return self.__class__(self.p, self.x / other.x, self.exponent - other.exponent, valuation=val)

    def __lshift__(self, n):
        """
        Return the product of this element by ``p^n``

        INPUT:

        - ``n`` -- a relative integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x << 10
            2^17 * 123456
        """
        if self._valuation is None:
            val = None
        else:
            val = self._valuation + n
        return self.__class__(self.p, self.x, self.exponent + n, valuation=val)

    def __rshift__(self, n):
        """
        Return the quotient of this element by ``p^n``

        INPUT:

        - ``n`` -- a relative integer

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x >> 10
            2^-3 * 123456
        """
        return self << (-n)

    def unit_part(self):
        """
        Return the unit part of this element, that is the part ``u``
        in the writing ``u * p^v`` with ``u`` coprime to ``p``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.unit_part()
            1929
        """
        if self.is_zero():
            raise ValueError("the unit part of zero is not defined")
        p = self.p
        val = self.valuation()
        x = self.x / (p ** (val-self.exponent))
        return self.__class__(p, x, 0, valuation=0)

    def xgcd(self,other):
        """
        Return the gcd of ``self`` and ``other`` together with two
        element ``u`` and ``v`` such that ``u*self + v*other = gcd``

        The ``gcd`` is normalized so that it is a power of ``p``

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: y = pRational(2, 891011, 12); y
            2^12 * 891011

            sage: d, u, v = x.xgcd(y)
            sage: d
            2^7 * 32
            sage: d.normalize(); d
            2^12 * 1

            sage: u*x + v*y
            2^7 * 32
        """
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if sexp < oexp:
            a = ZZ(self.x)
            b = ZZ(other.x * (p ** (oexp-sexp)))
            exp = sexp
        else:
            a = ZZ(self.x * (p ** (sexp-oexp)))
            b = ZZ(other.x)
            exp = oexp
        d, u, v = a.xgcd(b)
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = min(self._valuation, other._valuation)
        d = self.__class__(p, d, exp, valuation=val)
        u = self.__class__(p, u)
        v = self.__class__(p, v)
        return d, u, v

    def value(self):
        """
        Return this element as a rational number

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456, 7); x
            2^7 * 123456
            sage: x.value()
            15802368
        """
        return (self.p ** self.exponent) * self.x

    def list(self, prec):
        """
        Return the list of the digits of this element (written in radix 
        ``p``) up to position ``prec``.

        The first zeros are omitted.

        TESTS::

            sage: from sage.rings.padics.lattice_precision import pRational
            sage: x = pRational(2, 123456); x
            123456
            sage: x.list(5)
            []
            sage: x.list(20)
            [1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]

            sage: y = pRational(2, 123/456); y
            41/152
            sage: y.list(10)
            [1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1]
        """
        if self.x not in ZZ:
            self = self.reduce(prec)
        val = self.valuation()
        p = self.p
        x = ZZ(self.x * p**(self.exponent - val))
        l = [ ]; i = val
        while i < prec:
            x, digit = x.quo_rem(p)
            l.append(digit)
            i += 1
        return l


# Helper function
#################

def list_of_padics(elements):
    """
    Convert a list of p-adic composed elements (as polynomials, matrices)
    to a list of weak refererences of their p-adic coefficients.

    This is an helper function for the methods :meth:`precision_lattice`

    TESTS::

        sage: from sage.rings.padics.lattice_precision import list_of_padics
        sage: R = ZpLC(2)
        sage: M = random_matrix(R,2,2)
        sage: list_of_padics(M)
        [<weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
         <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
         <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
         <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>]
    """
    from sage.rings.padics.padic_lattice_element import pAdicLatticeElement
    if isinstance(elements, pAdicLatticeElement):
        return [ weakref.ref(elements) ]
    if not isinstance(elements, list):
        elements = list(elements)
    ans = [ ]
    for x in elements:
        ans += list_of_padics(x)
    return ans

def format_history(tme, status, timings):
    """
    Return a formated output for the history.

    This is an helper function for the methods :meth:`history`.

    TESTS::

        sage: from sage.rings.padics.lattice_precision import format_history
        sage: format_history(1.23456789, ['o','o','o','o','o','o','~','o','o'], true)
        '1.234568s  oooooo~oo'
        sage: format_history(1.23456789, ['o','o','o','o','o','o','~','o','o'], false)
        'oooooo~oo'

        sage: format_history(12.3456789, ['o','o','o','o','o','o','~','o','o'], true)
        '  >= 10s   oooooo~oo'
        sage: format_history(10^(-10), ['o','o','o','o','o','o','~','o','o'], true)
        '   ---     oooooo~oo'
        sage: format_history(-1, ['o','o','o','o','o','o','~','o','o'], true)
        ' Timings   oooooo~oo'
    """
    status = ''.join(status)
    if timings:
        if tme < 0:
            s = " Timings "
        elif tme < 0.000001:
            s = "   ---   "
        elif tme >= 10:
            s = "  >= 10s "
        else:
            s = "%.6fs" % tme
        return s + "  " + status
    else:
        return status


class DifferentialPrecisionGeneric(UniqueRepresentation, SageObject):
    """
    A generic class for precision objects obtained by automatic
    differentiation
    """
    def __init__(self, p, type, label):
        self._p = p
        self._label = label
        self._type = type
        self._elements = [ ]   # Probably better to use a double chained list
        self._matrix = { }
        self._marked_for_deletion = [ ]
        self._approx_zero = pRational(p, ZZ(0))
        self._thresold_deletion = DEFAULT_THRESOLD_DELETION
        # History
        self._history_init = None
        self._history = None

    def label(self):
        """
        Return the label of the parent to which this precision lattice
        corresponds.

        EXAMPLE::

            sage: R = ZpLC(2, label="mylabel")
            sage: R.precision().label()
            'mylabel'
        """
        return self._label

    def _repr_(self):
        """
        Return a string representation of this precision lattice

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: R.precision()
            Precision lattice on ... objects

        If a label has been specified, it is included in the representation

            sage: R = ZpLC(2, label="mylabel")
            sage: R.precision()
            Precision lattice on 0 object (label: mylabel)
        """
        n = len(self._elements)
        if self._label is None:
            if n > 1:
                return "Precision %s on %s objects" % (self._type, len(self._elements))
            else:
                return "Precision %s on %s object" % (self._type, len(self._elements))
        else:
            if n > 1:
                return "Precision %s on %s objects (label: %s)" % (self._type, len(self._elements), self._label)
            else:
                return "Precision %s on %s object (label: %s)" % (self._type, len(self._elements), self._label)

    def thresold_deletion(self, thresold=None):
        """
        Return or set the thresold for column deletion

        When a variable dies, the ambiant space in which the precision
        module lives can be reduced (by projection onto the hyperplane
        defined by the dead variable).
        However this reduction has a cost because it leads to re-echonize
        a part of the matrix that encodes the precision. The size of this
        part is roughly the distance between the last column and the one
        corresponding to the dead variable.

        The thresold deletion is the maximal distance until which the
        above reduction is performed. After the thresold, the column of
        the dead variable is kept in this matrix as if the variable were
        not destroyed.

        INPUT:

        - ``thresold`` -- a non negative integer, ``Infinity`` or ``None`` 
          (default: ``None``): if ``None``, return the current thresold;
          otherwise set the thresold to the given value

        NOTE::

        Setting the thresold to ``0`` disables the dimension reduction.

        Setting the thresold to ``Infinity``

        EXAMPLES::

            sage: R = ZpLC(2, label='thresold_deletion')
            sage: prec = R.precision()
            sage: prec.thresold_deletion()
            50

            sage: prec.thresold_deletion(20)
            sage: prec.thresold_deletion()
            20

            sage: prec.thresold_deletion(-2)
            Traceback (most recent call last):
            ...
            ValueError: The thresold must be a nonnegative integer or Infinity

        TESTS::

            sage: L = [ R.random_element() for _ in range(100) ]
            sage: prec.ambiant_dimension()
            100

            sage: del L[:10]
            sage: s = sum(L)
            sage: prec.ambiant_dimension()
            102

            sage: del L[-10:]
            sage: s = sum(L)
            sage: prec.ambiant_dimension()
            93
        """
        if thresold is None:
            return self._thresold_deletion
        if thresold is Infinity or (thresold in ZZ and thresold >= 0):
            self._thresold_deletion = thresold
        else:
            raise ValueError("The thresold must be a nonnegative integer or Infinity")

    def prime(self):
        """
        Return the underlying prime number attached to this precision lattice.

        EXAMPLE::

            sage: R = ZpLC(2, label="mylabel")
            sage: R.precision().prime()
            2
        """
        return self._p

    def _index(self, ref):
        """
        Return the index of the element whose reference is ``ref``
        """
        return self._elements.index(ref)

    def ambiant_dimension(self):
        """
        Return the dimension of the vector space in which the precision
        module/lattice lives

        EXAMPLES:

            sage: R = ZpLC(2, label='ambient_dim')
            sage: prec = R.precision()

            sage: x,y = R(1,10), R(1,5)
            sage: prec.ambiant_dimension()
            2
            sage: prec.dimension()
            2

            sage: u = x + y
            sage: prec.ambiant_dimension()
            3
            sage: prec.dimension()
            3

        In the case of ``ZpLC`` (lattice-cap precision), it is always
        equal to the dimension of the lattice.

        In the case of ``ZpLF`` (lattice-float precision), the precision
        object is not necessarily a lattice and then may have smaller
        dimension::

            sage: R = ZpLF(2, label='ambient_dim')
            sage: prec = R.precision()

            sage: x,y = R(1,10), R(1,5)
            sage: prec.ambiant_dimension()
            2
            sage: prec.dimension()
            2

            sage: u = x + y
            sage: prec.ambiant_dimension()
            3
            sage: prec.dimension()
            2
        """
        return len(self._matrix)

    def dimension(self):
        """
        Return the dimension of this precision module

        This method must be implemented in subclasses
        """
        raise NotImplementedError("implement this function is subclasses")

    def new_element(self, *args):
        """
        Insert a new element in this precision module

        This method must be implemented in subclasses
        """
        raise NotImplementedError("implement this function is subclasses")

    def mark_for_deletion(self, ref):
        """
        Mark for deletion an element of this precision module.

        This method must be implemented in subclasses
        """
        raise NotImplementedError("implement this function is subclasses")

    def del_elements(self, thresold=None):
        """
        Delete all elements marked for deletion up to the given thresold

        This method must be implemented in subclasses
        """
        raise NotImplementedError("implement this function in subclasses")

    def precision_absolute(self, x):
        """
        Return the absolute precision of the given element

        This method must be implemented in subclasses
        """
        raise NotImplementedError("implement this function in subclasses")

    def precision_lattice(self, elements=None):
        """
        Return a lattice modeling the precision on the given set of elements
        or, if not given, on the whole set of elements tracked by the precision
        module

        This method must be implemented in subclasses
        """
        raise NotImplementedError("implement this function in subclasses")

    def number_of_diffused_digits(self, elements=None):
        """
        Return the number of diffused digits of precision within a 
        subset of elements.

        NOTE:

        A diffused digit of precision is a known digit which is not
        located on a single variable but only atppears on a suitable
        linear combinaison of variables.

        The number of diffused digits of precision quantifies the
        quality of the approximation of the lattice precision by a
        jagged precision (that is a precision which is split over 
        all variables).

        We refer to [] for a detail exposition of the notion of
        diffused digits.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: prec = R.precision()
            sage: x = R(1,10); y = R(1,5)
            sage: u = x + y
            sage: v = x - y

            sage: prec.number_of_diffused_digits([x,y])
            0
            sage: prec.number_of_diffused_digits([u,v])
            6

        The elements `u` and `v` are known at absolute precision `O(2^5)`.
        However, the sum `u + v = 2x` is known at precision `O(2^11)`, that
        is with `6` more digits.
        That is where the `6` diffused digits of precision comes from.

        Here is another example with matrices::

            sage: M = matrix(R, 2, 2, [R(3,5),R(7,5),R(1,5),R(11,1)])
            sage: N = M^10

        The next syntax provides as easy way to select an interesting
        subset of variables (the selected subset consists of the four
        entries of the matrix ``N``)::

            sage: prec.number_of_diffused_digits(N)
            17
        """
        M = self.precision_lattice(elements)
        n = M.nrows()
        if M.ncols() > n:
            return Infinity
        p = self._p
        diffused = 0
        for j in range(n):
            val = minval = M[j,j].valuation(p)
            for i in range(j):
                v = M[i,j].valuation(p)
                if v < minval: minval = v
            diffused += val - minval
        return diffused

    def number_of_tracked_elements(self, dead=True):
        """
        Return the number of tracked elements through this precision
        lattice

        INPUT:

        - ``dead`` -- a boolean (default: ``True``); whether dead
          elements for which the corresponding column is still not
          erased should be counted or not

        EXAMPLES::

            sage: R = ZpLC(2, label='count')
            sage: prec = R.precision()
            sage: x = R(1,10); y = R(1,5)
            sage: prec.number_of_tracked_elements()
            2

            sage: u = x + y
            sage: v = x - y
            sage: prec.number_of_tracked_elements()
            4

            sage: del x; del y
            sage: prec.number_of_tracked_elements()
            4
            sage: prec.number_of_tracked_elements(dead=False)
            2

            sage: prec.del_elements()
            sage: prec.number_of_tracked_elements()
            2
        """
        if dead:
            return len(self._elements)
        else:
            count = 0
            for x_ref in self._elements:
                if x_ref() is not None: count += 1
            return count

    def tracked_elements(self, values=True, dead=True):
        """
        Return the list of tracked elements

        INPUT:

        - ``values`` -- a boolean (default: ``True``); if false,
          the method returns a list of weak references on tracked
          elements instead

        - ``dead`` -- a boolean (default: ``True``); whether dead
          elements for which the corresponding column is still not
          erased should be listed or not

        EXAMPLES::

            sage: R = ZpLC(2, label='tracked')
            sage: prec = R.precision()
            sage: x = R(1,10); y = R(1,5)
            sage: prec.tracked_elements()
            [1 + O(2^10), 1 + O(2^5)]
            sage: prec.tracked_elements(values=False)
            [<weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; dead>]
            sage: prec.tracked_elements(values=False, dead=False)
            [<weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>]

            sage: u = x + y
            sage: v = x - y
            sage: prec.tracked_elements()
            [1 + O(2^10), 1 + O(2^5), 2 + O(2^5), O(2^5)]
            sage: prec.tracked_elements(values=False)
            [<weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; dead>]

            sage: del x; del y
            sage: prec.tracked_elements()
            [None, None, 2 + O(2^5), O(2^5), None]
            sage: prec.tracked_elements(values=False)
            [<weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; to 'pAdicLatticeCapElement' at 0x...>,
             <weakref at 0x...; dead>]
        """
        if values:
            if dead:
                return [ ref() for ref in self._elements ]
            else:
                return [ ref() for ref in self._elements if ref() is not None ]
        else:
            if dead:
                return list(self._elements)
            else:
                return [ ref for ref in self._elements if ref() is not None ]


    # History

    def history_enable(self):
        """
        Enable history.

        We refer to the documentation of the method :meth:`history` for 
        a complete documentation (including examples) about history.

        TESTS::

            sage: R = ZpLC(2, label='history_en')
            sage: prec = R.precision()

            sage: print(prec.history())  # history is disabled by default
            Traceback (most recent call last):
            ...
            ValueError: History is not tracked

            sage: prec.history_enable()
            sage: print(prec.history())
             Timings
               ---     

        .. SEEALSO::

            :meth:`history`, :meth:`history_disable`, :meth:`history_clear`
        """
        if self._history is None:
            self._history_init = ( len(self._elements), list(self._marked_for_deletion) )
            self._history = [ ]

    def history_disable(self):
        """
        Disable history.

        We refer to the documentation of the method :meth:`history` for 
        a complete documentation (including examples) about history.

        TESTS::

            sage: R = ZpLC(2, label='history_dis')
            sage: prec = R.precision()

            sage: print(prec.history())  # history is disabled by default
            Traceback (most recent call last):
            ...
            ValueError: History is not tracked

            sage: prec.history_enable()
            sage: print(prec.history())
             Timings
               ---     

            sage: prec.history_disable()
            sage: print(prec.history())
            Traceback (most recent call last):
            ...
            ValueError: History is not tracked

        .. SEEALSO::

            :meth:`history`, :meth:`history_enable`, :meth:`history_clear`
        """
        self._history = self._history_init = None

    def history_clear(self):
        """
        Clear history

        We refer to the documentation of the method :meth:`history` for 
        a complete documentation (including examples) about history.

        TESTS::

            sage: R = ZpLC(2, label='history_clear')
            sage: prec = R.precision()
            sage: prec.history_enable()

            sage: x = R(1,10); y = R(1,5)
            sage: x,y = x+y, x-y
            sage: print(prec.history())
             Timings
               ...     oooo
               ...     ~~oo

        When we clear history, only the last line is kept::

            sage: prec.history_clear()
            sage: print(prec.history())
             Timings   ~~oo
               ---     ~~oo

            sage: del x;

            sage: print(prec.history())
             Timings   ~~oo
               ...     ~~~o

        .. SEEALSO::

            :meth:`history`, :meth:`history_enable`, :meth:`history_disable`
        """
        if self._history is None:
            raise ValueError("History is not tracked")
        self._history_init = ( len(self._elements), list(self._marked_for_deletion) )
        self._history = [ ]

    def history(self, compact=True, separate_reduce=False, timings=True, output_type='asciiart'):
        """
        Show history

        The history records creations and deletions of elements attached 
        to this precision lattice, together with many timings.

        INPUT:

        - ``compact`` -- a boolean (default: ``True``); if true, all 
          consecutive operations of the same type appear on a single row

        - ``separate_reduce`` -- a boolean (default: ``False``); specify
          whether partial/full Hermite reduction should be displayed
          separatedly

        - ``timings`` -- a boolean (default: ``True``); specify whether
          timings should be displayed

        - ``output_type`` -- only ``asciiart`` is implemented for now.

        IMPORTANT NOTE:

        History is disabled by default.
        It should then be enabled (through a call to the method :meth:`history_enable`)
        before use.

        EXAMPLES::

            sage: R = ZpLC(2, label='history_en')
            sage: prec = R.precision()

        We first enable history::

            sage: prec.history_enable()

        At the beginning, the history is of course empty::

            sage: print(prec.history())
             Timings     
               ---

        Now we start creating and deleting elements::

            sage: L = [ R.random_element() for _ in range(20) ]
            sage: for p in range(20):
            ....:    if is_prime(p): L[p] = None
            sage: prec.del_elements()

            sage: print(prec.history())
             Timings    
               ...   oooooooooooooooooooo
               ...   oo~~o~o~ooo~o~ooo~o~
               ...   oooooooooooo

        The legend is the following::
        - the symbol ``o`` represents a tracked element
        - the symbol ``~`` represents an element which is marked for deletion

        On the history, we see:
        - 1st line: twenty new elements were created
          (this corresponds to the affectation of the list ``L``)
        - 2nd line: elements at prime positions were marked for deletion
          (this corresponds to the ``for`` loop)
        - 3rd line: the above elements are indeed deleted
          (this corresponds to the call of the method :meth:`del_elements`

        Here are some variants::

            sage: print(prec.history(timings=False))
            oooooooooooooooooooo
            oo~~o~o~ooo~o~ooo~o~
            oooooooooooo

            sage: print(prec.history(separate_reduce=True))
             Timings
               ...   oooooooooooooooooooo
               ...   oo~~o~o~ooo~o~ooo~o~
               ...   oo~~o~o~ooo~ooooo
               ...   oo~~o~o~ooo~orrrr 
               ...   oo~~o~o~oooooooo
               ...   oo~~o~o~ooorrrrr
               ...   oo~~o~ooooooooo
               ...   oo~~o~orrrrrrrr
               ...   oo~~oooooooooo
               ...   oo~~orrrrrrrrr
               ...   oo~oooooooooo
               ...   oo~rrrrrrrrrr
               ...   oooooooooooo
               ...   oorrrrrrrrrr
               ---   oooooooooooo

        The symbol ``r`` represents a column of the precision matrix which is
        currently under partial Hermite reduction

        Timings for automatic reduction do not appear because they are included
        in the timings for deletion.

        The symbol ``R`` is used to symbolize a column which is under full 
        Hermite reduction. Note that full Hermite reduction are never performed 
        automatically but needs to be called by hand::

            sage: prec.reduce()
            sage: print(prec.history(separate_reduce=True))
             Timings
               ...
               ...   RRRRRRRRRRRR
               ---   oooooooooooo

        .. SEEALSO::

            :meth:`history_enable`, :meth:`history_disable`, :meth:`history_clear`
        """
        if self._history is None:
            raise ValueError("History is not tracked")
        total_time = 0
        if output_type == 'asciiart':
            # Legend:
            #  o : tracked element
            #  ~ : element marked for deletion
            #  r : partial reduction
            #  R : full Hermite reduction
            (n, mark) = self._history_init
            status = n*['o']
            for index in mark:
                status[index] = '~'
            hist = [ format_history(-1, status, timings) ]
            oldevent = ''; total_time = 0
            for (event, index, tme) in self._history:
                if event == 'partial reduce' or event == 'full reduce':
                    if separate_reduce:
                        if total_time > 0:
                            hist.append(format_history(total_time, status, timings))
                        if event == 'partial reduce': code = 'r'
                        else: code = 'R'
                        status_red = status[:index] + (len(status) - index) * [code]
                        hist.append(format_history(tme, status_red, timings))
                        total_time = 0
                        oldevent = ''
                    else:
                        total_time += tme
                    continue
                if not compact or event != oldevent:
                    if total_time > 0:
                        hist.append(format_history(total_time, status, timings))
                    total_time = 0
                    oldevent = event
                total_time += tme
                if event == 'add':
                    if index is None:
                        status.append('o')
                    else:
                        status = status[:index] + ['o'] + status[index:]
                elif event == 'mark':
                    status[index] = '~'
                elif event == 'del':
                    del status[index]
            if total_time > 0 or oldevent == '':
                hist.append(format_history(total_time, status, timings))
            return '\n'.join(hist)
        else:
            raise NotImplementedError

    def timings(self, action=None):
        """
        Return cumulated timings (grouped by actions) since the last 
        time history has been cleared.

        INPUT:

        - ``action`` -- ``None`` (the default), ``add``, ``mark``, ``del``,
          ``partial reduce`` or ``full reduce``; if not None, return the 
          cumulated timing corresponding to this action; otherwise, return
          a dictionary

        Here are the meanings of the keywords above:
        - ``add``: time spent in adding new colunmns to the precision matrix
          (corresponding to the creation of new elements)
        - ``mark``: time spent in marking elements for deletion
        - ``del``: time spent in deleting columns of the precision matrix
          and re-echelonizing the matrix
        - ``partial reduce``: time spent in partial Hermite reduction
        - ``full reduce``: time spent in full Hermite reduction

        EXAMPLES::

            sage: R = ZpLC(2, label='timings')
            sage: prec = R.precision()
            sage: prec.history_enable()
            sage: M = random_matrix(R,5,5)
            sage: N = M^10
            sage: prec.timings()    # somewhat random
            {'add': 1.0530245304107666,
             'del': 0.24358701705932617,
             'mark': 0.0013289451599121094,
             'partial reduce': 0.21604204177856445
             'full reduce': 0}

        TESTS::

            sage: prec.history_clear()
            sage: prec.timings()
            {'add': 0, 'del': 0, 'full reduce': 0, 'mark': 0, 'partial reduce': 0}
        """
        if self._history is None:
            raise ValueError("History is not tracked")
        tme_by_event = { 'add': 0, 'del': 0, 'mark': 0, 'partial reduce': 0, 'full reduce': 0 }
        for (event, _, tme) in self._history:
            tme_by_event[event] += tme
        if action is None:
            return tme_by_event
        if tme_by_event.has_key(action):
            return tme_by_event[action]
        else:
            raise ValueError("invalid event")


# class PrecisionLattice
########################

class PrecisionLattice(DifferentialPrecisionGeneric):
    """
    A class for handling precision lattices which are used to
    track precision in the ZpLC model.

    The precision lattice is stored as a triangular matrix whose
    rows are generators of the lattice.
    """
    def __init__(self, p, label):
        DifferentialPrecisionGeneric.__init__(self, p, 'lattice', label)
        self._absolute_precisions = { }
        self._capped = { }

    def _index(self, ref):
        """
        Return the index of the element whose reference is ``ref``
        """
        return len(self._matrix[ref]) - 1

    def dimension(self):
        """
        Return the dimension of this lattice

        EXAMPLES::

            sage: R = ZpLC(5, label='dimension')
            sage: prec = R.precision()
            sage: prec.dimension()
            0

            sage: x = R(1,10)
            sage: prec.dimension()
            1
        """
        return len(self._matrix)

    def reduce(self, index=0, partial=False):
        """
        Reduce the size of the entries above the diagonal of the precision matrix

        INPUT:

        - ``index`` -- an integer, the starting row for which the reduction
          is performed

        - ``partial`` -- a boolean (default: False) specifying whether a
          partial or a full Hermite reduction should be performed

        NOTE:

        The partial reduction has cost `O(m^2)` where `m` is the number of 
        rows that need to be reduced (that is the difference between the 
        total number of rows and ``index``).

        The full Hermite reduction has cost `O(m^3)`.

        NOTE:

        The software ensures that the precision lattice is always 
        partially reduced.
        Calling the function manually with the argument ``partial=True``
        should then just do nothing.

        TESTS::

            sage: R = ZpLC(2)
            sage: x = R.random_element()
            sage: del x
            sage: R.precision().del_elements()   # indirect doctest
        """
        n = len(self._elements)
        if index >= n-1:
            return
        if partial:
            # Partial reduction
            # Cost: O(m^2) with m = n-index
            tme = walltime()
            diffval = (n-index) * [0]
            for j in range(n-1, index, -1):
                col = self._matrix[self._elements[j]]
                prec = col[j].valuation() - diffval[j-index]
                for i in range(index,j):
                    col[i] = col[i].reduce(prec)
                    col[i].normalize()  # seems to be faster then
                    dval = col[i].valuation() - prec
                    if dval < diffval[i-index]:
                        diffval[i-index] = dval
            # We update history
            if self._history is not None:
                self._history.append(('partial reduce', index, walltime(tme)))
        else:
            # Full Hermite reduction
            # Cost: O(m^3) with m = n-index
            tme = walltime()
            for j in range(index+1, n):
                # In what follows, we assume that col[j] is a power of p
                col = self._matrix[self._elements[j]]
                valpivot = col[j].valuation()
                for i in range(index, j):
                    reduced = col[i].reduce(valpivot)
                    scalar = (col[i] - reduced) >> valpivot
                    if scalar.is_zero(): continue
                    col[i] = reduced
                    col[i].normalize()
                    for j2 in range(j+1, n):
                        col2 = self._matrix[self._elements[j2]]
                        col2[i] -= scalar*col2[i]
                        col2[i].normalize()
            # We update history
            if self._history is not None:
                self._history.append(('full reduce', index, walltime(tme)))

    def new_element(self, x, dx, bigoh, dx_mode='linear_combinaison', capped=False):
        """
        Update the lattice when a new element is created.

        This function is not meant to be called manually.
        It is automatically called by the parent when a new
        element is created.

        INPUT:

        - ``x`` -- the newly created element

        - ``dx`` -- a dictionary representing the differential of ``x``

        - ``dx_mode`` -- a string, either ``linear_combinaison`` (the default)
          or ``values``

        - ``capped`` -- a boolean, whether this element has been capped 
          according to the parent's cap

        If ``dx_mode`` is ``linear_combinaison``, the dictionary ``dx`` 
        encodes the expression of the differential of ``x``. 
        For example, if ``x`` was defined as ``x = y*z`` then:

        .. MATH::

            dx = y dz + z dy

        and the corresponding dictionary is ``{z: y, y: z}`` (except
        that the keys are not the elements themselves but weak references
        to them).

        If ``dx_mode`` is ``values``, the dictionary ``dx`` directly
        specifies the entries that have to stored in the precision lattice.
        This mode is only used for multiple conversion between different
        parents (see :meth:`multiple_conversion`).

        TESTS::

            sage: R = ZpLC(2)
            sage: x = R.random_element()
            sage: y = R.random_element()
            sage: z = x*y    # indirect doctest
        """
        # First we delete some elements marked for deletion
        if self._marked_for_deletion:
            self.del_elements(thresold=self._thresold_deletion)

        # Then we add the new element
        tme = walltime()
        p = self._p
        n = len(self._elements)
        x_ref = weakref.ref(x, self.mark_for_deletion)
        self._elements.append(x_ref)
        col = n * [self._approx_zero]
        if dx_mode == 'linear_combinaison':
            for elt,scalar in dx:
                ref = weakref.ref(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                c = self._matrix[ref]
                for i in range(len(c)):
                    col[i] += scalar * c[i]
        elif dx_mode == 'values':
            for elt,scalar in dx:
                ref = weakref.ref(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                i = self._index(ref)
                col[i] = scalar
        else:
            raise ValueError("dx_mode must be either 'linear_combinaison' or 'values'")
        for i in range(n):
            col[i] = col[i].reduce(bigoh)
        col.append(pRational(p, ZZ(1), bigoh))
        self._matrix[x_ref] = col
        self._capped[x_ref] = capped

        # We update history
        if self._history is not None:
            self._history.append(('add', None, walltime(tme)))

    def mark_for_deletion(self, ref):
        """
        Mark an element for deletion.

        This function is not meant to be called manually.
        It is automatically called by the garbage collection when 
        an element is collected.

        INPUT:

        - ``ref`` -- a weak reference to the destroyed element

        NOTE::

        This method does not update the precision lattice.
        The actual update is performed when the method :meth:`del_elements`
        is called. This is automatically done at the creation of a new
        element but can be done manually as well.

        EXAMPLES::

            sage: R = ZpLC(2, label='markdel')
            sage: prec = R.precision()
            sage: x = R(1,10)
            sage: prec
            Precision lattice on 1 object (label: markdel)
            sage: del x   # indirect doctest: x is here marked for deletion
            sage: prec
            Precision lattice on 1 object (label: markdel)
            sage: prec.del_elements()       # x is indeed deleted
            sage: prec
            Precision lattice on 0 object (label: markdel)
        """
        tme = walltime()
        try:
            index = self._index(ref)
        except (IndexError, KeyError):
            return
        self._marked_for_deletion.append(index)
        if self._history is not None:
            self._history.append(('mark', index, walltime(tme)))

    def del_elements(self, thresold=None):
        """
        Erase columns of the lattice precision matrix corresponding to
        elements which are marked for deletion and echelonize the matrix
        in order to keep it upper triangular.

        INPUT:

        - ``thresold`` -- an integer or ``None`` (default: ``None``):
          a column whose distance to the right at greater than the
          thresold is not erased

        EXAMPLES::

            sage: R = ZpLC(2, label='delelts')
            sage: prec = R.precision()

            sage: x = R(1,10)
            sage: prec
            Precision lattice on 1 object (label: delelts)
            sage: prec.precision_lattice()
            [1024]

            sage: del x
            sage: prec
            Precision lattice on 1 object (label: delelts)
            sage: prec.precision_lattice()
            [1024]

            sage: prec.del_elements()
            sage: prec
            Precision lattice on 0 object (label: delelts)
            sage: prec.precision_lattice()
            []
        """
        p = self._p
        n = len(self._elements)

        self._marked_for_deletion.sort(reverse=True)
        count = 0
        for index in self._marked_for_deletion:
            if thresold is not None and index < n - thresold: break
            n -= 1; count += 1

            tme = walltime()
            ref = self._elements[index]
            del self._elements[index]
            del self._matrix[ref]
            capped = self._capped[ref]
            del self._capped[ref]

            # Now, we echelonize
            for i in range(index,n):
                ref = self._elements[i]
                col = self._matrix[ref]
                if col[i].valuation() < col[i+1].valuation():
                    self._capped[ref], capped = capped, capped or self._capped[ref]
                else:
                    capped = capped or self._capped[ref]
                      
                d, u, v = col[i].xgcd(col[i+1])
                up, vp = col[i+1]/d, col[i]/d
                col[i] = d
                del col[i+1]
                for j in range(i+1,n):
                    col = self._matrix[self._elements[j]]
                    col[i], col[i+1] = u*col[i] + v*col[i+1], up*col[i] - vp*col[i+1]

            # We update history
            if self._history is not None:
                self._history.append(('del', index, walltime(tme)))

            # And we reduce a bit
            # (we do not perform a complete reduction because it is costly)
            self.reduce(index, partial=True)

        del self._marked_for_deletion[:count]

    def lift_to_precision(self, x, prec):
        """
        Lift the specified element to the specified precision

        INPUT:

        - ``x`` -- the element whose precision has to be lifted

        - ``prec`` -- the new precision

        NOTE:

        The new precision lattice is computed as the intersection
        of the current precision lattice with the subspace

        ..MATH::

            p^{prec} \Z_p dx \oplus \bigoplus_{y \neq x} \Q_p dy

        This function may change at the same time the precision of 
        other elements having the same parent.

        NOTE:

        This function is not meant to be called directly.
        You should prefer call the method :meth:`lift_to_precision`
        of ``x`` instead.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1,10); x
            1 + O(2^10)
            sage: y = R(1,5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)

            sage: prec = R.precision()
            sage: prec.lift_to_precision(z, 12)
            sage: z
            2 + O(2^12)
            sage: y
            1 + O(2^10)
        """
        ref = weakref.ref(x)
        col = self._matrix[ref]
        n = len(self._elements)

        rows_by_val = { }
        for i in range(len(col)):
            v = col[i].valuation()
            if v >= prec: continue
            if rows_by_val.has_key(v):
                rows_by_val[v].append(i)
            else:
                rows_by_val[v] = [i]
        vals = rows_by_val.keys()
        vals.sort()
        vals.append(prec)

        for t in range(len(vals)-1):
            v, w = vals[t], vals[t+1]
            rows = rows_by_val[v]
            piv = max(rows)
            alpha = col[piv].unit_part()
            for i in rows:
                if i == piv: continue
                # We clear the entry on the i-th line
                beta = col[i].unit_part()
                for j in range(piv,n):
                    col_cur = self._matrix[self._elements[j]]
                    col_cur[i] = alpha*col_cur[i] - beta*col_cur[piv]
            # We rescale the piv-th line
            for j in range(piv,n):
                col_cur = self._matrix[self._elements[j]]
                col_cur[piv] = col_cur[piv] << (w-v)
            # Now the entry on the piv-th line has valuation w
            # We update the dictionary accordingly
            if w < prec:
                rows_by_val[w].append(piv)

        # We clear the cached absolute precisions
        self._absolute_precisions = { }

    def _compute_precision_absolute(self, ref):
        """
        Compute the absolute precision of the given element and cache it

        For internal use.
        """
        col = self._matrix[ref]
        absprec = Infinity
        capped = False
        for i in range(len(col)):
            v = col[i].valuation()
            if v < absprec:
                absprec = v
                capped = self._capped[self._elements[i]]
            elif v == absprec:
                capped = capped and self._capped[self._elements[i]]
        self._absolute_precisions[ref] = [absprec, capped]

    def precision_absolute(self, x):
        """
        Return the absolute precision of the given element

        INPUT:

        - ``x`` -- the element whose absolute precision is requested

        NOTE:

        The absolute precision is obtained by projecting the precision
        lattice onto the line of coordinate ``dx``

        NOTE:

        This function is not meant to be called directly.
        You should prefer call the method :meth:`precision_absolute`
        of ``x`` instead.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1,10); x
            1 + O(2^10)
            sage: y = R(1,5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)
            sage: z.precision_absolute()  # indirect doctest
            5
        """
        ref = weakref.ref(x)
        if not self._absolute_precisions.has_key(ref):
            self._compute_precision_absolute(ref)
        return self._absolute_precisions[ref][0]

    def is_precision_capped(self, x):
        """
        Return whether the absolute precision on the given 
        results from a cap coming from the parent

        INPUT:

        - ``x`` -- the element

        This function is not meant to be called directly.
        You should prefer call the method :meth:`is_precision_capped`
        of ``x`` instead.

        EXAMPLES::

            sage: R = ZpLC(2)
            sage: x = R(1,10); x
            1 + O(2^10)
            sage: x.is_precision_capped()  # indirect doctest
            False

            sage: y = x-x; y
            O(2^40)
            sage: y.is_precision_capped()  # indirect doctest
            True
        """
        ref = weakref.ref(x)
        if not self._absolute_precisions.has_key(ref):
            self._compute_precision_absolute(ref)
        return self._absolute_precisions[ref][1]

    def precision_lattice(self, elements=None):
        """
        Return a matrix representing the precision lattice on a
        subset of elements.

        INPUT:

        - ``elements`` -- a list of elements or ``None`` (default: ``None``)

        - ``echelon`` -- a boolean (default: ``True``); specify whether
          the result should be in echelon form

        EXAMPLES::

            sage: R = ZpLC(2, label='preclattice')
            sage: prec = R.precision()
            sage: x = R(1,10); y = R(1,5)
            sage: u = x + y
            sage: v = x - y
            sage: prec.precision_lattice()
            [         1024             0          1024          1024]
            [            0            32            32 1099511627744]
            [            0             0       2097152             0]
            [            0             0             0 1099511627776]
            sage: prec.precision_lattice([u,v])
            [  32 2016]
            [   0 2048]

        Here is another example with matrices::

            sage: M = matrix(R, 2, 2, [R(3,5),R(7,5),R(1,5),R(11,1)])
            sage: N = M^10
            sage: prec.precision_lattice()
            23 x 23 dense matrix over Integer Ring (use the '.str()' method to see the entries)

        The next syntax provides as easy way to select an interesting
        subset of variables (the selected subset consists of the four
        entries of the matrix ``N``)::

            sage: prec.precision_lattice(N)
            [  2048    512  28160 230400]
            [     0   2048  14336 258048]
            [     0      0  65536  65536]
            [     0      0      0 262144]

        We can give a list of matrices as well::

            sage: prec.precision_lattice([M,N])
            [       32         0         0         0 226115584  96788480  52174848  82804736]
            [        0        32         0         0  52174848 121765888  11829248  28516352]
            [        0         0        32         0  96788480  42762240 121765888 199614464]
            [        0         0         0         2   5175296  12475904   1782272   4045824]
            [        0         0         0         0 268435456         0         0         0]
            [        0         0         0         0         0 268435456         0         0]
            [        0         0         0         0         0         0 268435456         0]
            [        0         0         0         0         0         0         0 268435456]
        """
        if elements is None:
            elements = self._elements
        else:
            elements = list_of_padics(elements)
        n = len(self._elements)
        rows = [ ]; val = 0
        for ref in elements:
            col = self._matrix[ref]
            row = [ x.value() for x in col ]
            valcol = min([ x.valuation() for x in col ])
            if valcol < val: val = valcol
            row += (n-len(row)) * [ZZ(0)]
            rows.append(row)
        from sage.matrix.constructor import matrix
        M = matrix(rows).transpose()
        if val < 0:
            M *= self._p ** (-val)
        M = M.change_ring(ZZ)
        M.echelonize()
        n = len(elements)
        M = M.submatrix(0,0,n,n)
        if val < 0:
            M *= self._p ** val
        return M


# class PrecisionModule
#######################

class PrecisionModule(DifferentialPrecisionGeneric):
    """
    A class for handling precision modules which are used to
    track precision in the ZpLF model.

    The precision module (which is not necessarily a lattice)
    is stored as a matrix whose rows are generators.
    """
    def __init__(self, p, label, prec):
        DifferentialPrecisionGeneric.__init__(self, p, 'module', label)
        self._absolute_precisions = { }
        self._zero_cap = prec
        self._internal_prec = prec + STARTING_ADDITIONAL_PREC
        self._count = 0
        self._thresold = 1

    def internal_prec(self):
        """
        Return the relative precision at which computations has handled
        internally

        It is slightly greater than the actual precision and increases
        a bit (at a logarithmic rate) when new elements are created 
        and/or computed

        EXAMPLES::

            sage: R = ZpLF(5, prec=20, label='internal_prec')
            sage: prec = R.precision()

            sage: prec.internal_prec()
            25

            sage: L = [ R.random_element() for _ in range(50) ]
            sage: prec.internal_prec()
            28
        """
        return self._internal_prec

    def dimension(self):
        """
        Return the dimension of this precision module

        EXAMPLES:

        In general, the dimension increases by 1 when a new
        element with a given precision is created::

            sage: R = ZpLF(2, label='dimension')
            sage: prec = R.precision()

            sage: prec.dimension()
            0
            sage: x = R.random_element(prec=10)
            sage: prec.dimension()
            1
            sage: y = R.random_element(prec=10)
            sage: prec.dimension()
            2

        However in general it does not increase while
        doing computations::

            sage: u = x + y
            sage: v = x^2 + 3*y + x*y + y^3
            sage: prec.dimension()
            2

        Of course, it may also decrease when a sufficient 
        number of variables are collected::

            sage: del x, y, u
            sage: prec.del_elements()
            sage: prec.dimension()
            1

            sage: del v
            sage: prec.del_elements()
            sage: prec.dimension()
            0
        """
        if len(self._elements) == 0:
            return 0
        return len(self._matrix[self._elements[-1]])

    def is_lattice(self):
        """
        Return ``True`` if this precision module is a lattice
        (i.e. has maximal dimension)

        EXAMPLES::

            sage: R = ZpLF(2, label='is_lattice')
            sage: prec = R.precision()

            sage: x = R(1,10)
            sage: y = R(1,5)
            sage: prec.is_lattice()
            True

            sage: u = x + y
            sage: prec.is_lattice()
            False

        .. SEEALSO::

            :meth:`dimension`
        """
        return self.dimension() == len(self._elements)

    def new_element(self, x, dx, bigoh, dx_mode='linear_combinaison'):
        """
        Update the lattice when a new element is created.

        This function is not meant to be called manually.
        It is automatically called by the parent when a new
        element is created.

        INPUT:

        - ``x`` -- the newly created element

        - ``dx`` -- a dictionary representing the differential of ``x``

        - ``dx_mode`` -- a string, either ``linear_combinaison`` (the default)
          or ``values``

        - ``capped`` -- a boolean, whether this element has been capped 
          according to the parent's cap

        If ``dx_mode`` is ``linear_combinaison``, the dictionary ``dx`` 
        encodes the expression of the differential of ``x``. 
        For example, if ``x`` was defined as ``x = y*z`` then:

        .. MATH::

            dx = y dz + z dy

        and the corresponding dictionary is ``{z: y, y: z}`` (except
        that the keys are not the elements themselves but weak references
        to them).

        If ``dx_mode`` is ``values``, the dictionary ``dx`` directly
        specifies the entries that have to stored in the precision module.
        This mode is only used for multiple conversion between different
        parents (see :meth:`multiple_conversion`).

        TESTS::

            sage: R = ZpLF(2)
            sage: x = R.random_element()
            sage: y = R.random_element()
            sage: z = x*y    # indirect doctest
        """
        # First we delete some elements marked for deletion
        if self._marked_for_deletion:
            self.del_elements(thresold=self._thresold_deletion)

        # We first increase the internal prec
        self._count += 1
        if self._count > self._thresold:
            self._internal_prec += 1
            self._thresold *= self._p

        tme = walltime()
        p = self._p
        n = self.dimension()
        x_ref = weakref.ref(x, self.mark_for_deletion)
        col = n * [self._approx_zero]
        if dx_mode == 'linear_combinaison':
            expected_vals = n * [ Infinity ]
            for elt,scalar in dx:
                ref = weakref.ref(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                c = self._matrix[ref]
                for i in range(len(c)):
                    summand = scalar * c[i]
                    expected_vals[i] = min(expected_vals[i], summand.valuation())
                    col[i] += summand
            for i in range(n):
                if col[i].valuation() >= expected_vals[i] + self._zero_cap:
                    col[i] = self._approx_zero
        elif dx_mode == 'values':
            for elt,scalar in dx:
                ref = weakref.ref(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                i = self._index(ref)
                col[i] = scalar
        else:
            raise ValueError("dx_mode must be either 'linear_combinaison' or 'values'")

        for i in range(n):
            col[i] = col[i].reduce_relative(self._internal_prec)
        if bigoh is not None:
            col.append(pRational(p, ZZ(1), bigoh))

        self._elements.append(x_ref)
        self._matrix[x_ref] = col

        # We update history
        if self._history is not None:
            self._history.append(('add', None, walltime(tme)))


    def mark_for_deletion(self, ref):
        """
        Mark an element for deletion.

        This function is not meant to be called manually.
        It is automatically called by the garbage collection when 
        an element is collected.

        INPUT:

        - ``ref`` -- a weak reference to the destroyed element

        NOTE::

        This method may do not update the precision module.
        The actual update is performed when the method :meth:`del_elements`
        is called. This is automatically done at the creation of a new
        element but can be done manually as well.

        EXAMPLES::

            sage: R = ZpLF(2, label='markdel')
            sage: prec = R.precision()
            sage: x = R(1,10)
            sage: prec
            Precision module on 1 object (label: markdel)
            sage: del x   # indirect doctest: x is here marked for deletion
            sage: prec
            Precision module on 1 object (label: markdel)
            sage: prec.del_elements()       # x is indeed deleted
            sage: prec
            Precision module on 0 object (label: markdel)
        """
        tme = walltime()
        try:
            index = self._index(ref)
        except (IndexError, KeyError):
            return
        if index == 0:
            length_before = 0
        else:
            length_before = len(self._matrix[self._elements[index-1]])
        length = len(self._matrix[ref])
        if length > length_before:
            self._marked_for_deletion.append(index)
            if self._history is not None:
                self._history.append(('mark', index, walltime(tme)))
        else:
            # if the column is not a pivot, we erase it without delay
            del self._elements[index]
            def decr(i):
                if i < index: return i
                else: return i-1
            self._marked_for_deletion = [ decr(i) for i in self._marked_for_deletion ]
            if self._history is not None:
                self._history.append(('del', index, walltime(tme)))


    def del_elements(self, thresold=None):
        """
        Erase columns of the lattice precision matrix corresponding to
        elements which are marked for deletion and reduce the matrix
        in order to keep it in echelon form.

        INPUT:

        - ``thresold`` -- an integer or ``None`` (default: ``None``):
          a column whose distance to the right at greater than the
          thresold is not erased

        EXAMPLES::

            sage: R = ZpLF(2, label='delelts')
            sage: prec = R.precision()

            sage: x = R(1,10)
            sage: prec
            Precision module on 1 object (label: delelts)
            sage: prec.precision_lattice()
            [1024]

            sage: del x
            sage: prec
            Precision module on 1 object (label: delelts)
            sage: prec.precision_lattice()
            [1024]

            sage: prec.del_elements()
            sage: prec
            Precision module on 0 object (label: delelts)
            sage: prec.precision_lattice()
            []
        """
        p = self._p
        n = len(self._elements)

        self._marked_for_deletion.sort(reverse=True)
        count = 0
        for index in self._marked_for_deletion:
            if thresold is not None and index < n - thresold: break
            n -= 1; count += 1

            tme = walltime()

            length = len(self._matrix[self._elements[index]])
            del self._matrix[self._elements[index]]
            del self._elements[index]
            start = index
            while start < n:
                i = start
                val = Infinity
                end = n
                while i < n:
                    col = self._matrix[self._elements[i]]
                    if len(col) > length: 
                        end = i
                        break
                    v = col[-1].valuation()
                    if v < val:
                        val = v
                        piv = i
                    i += 1
                if val < Infinity:
                    # another pivot has been found, we place it in front
                    self._elements[start], self._elements[piv] = self._elements[piv], self._elements[start]
                    break

                # No pivot was found. We re-echelonize
                for i in range(start, end):
                    del self._matrix[self._elements[i]][-1]
                if end == n: break
                # col is the column of index "end"
                # its size is (length + 1)
                d, u, v = col[length-1].xgcd(col[length])
                up, vp = col[length]/d, col[length-1]/d
                col[length-1] = d.reduce_relative(self._internal_prec)
                del col[length]
                start = end + 1
                for j in range(start, n):
                    col = self._matrix[self._elements[j]]
                    a1 = u*col[length-1]; a2 = v*col[length]; a = a1 + a2
                    b1 = up*col[length-1]; b2 = vp*col[length]; b = b1 + b2
                    if a.valuation() > min(a1.valuation(), a2.valuation()) + self._zero_cap:
                        col[length-1] = self._approx_zero
                    else:
                        col[length-1] = a.reduce_relative(self._internal_prec)
                    if b.valuation() > min(b1.valuation(), b2.valuation()) + self._zero_cap:
                        col[length] = self._approx_zero
                    else:
                        col[length] = b.reduce_relative(self._internal_prec)
                length += 1

            # We update history
            if self._history is not None:
                self._history.append(('del', index, walltime(tme)))

        del self._marked_for_deletion[:count]


    def _compute_precision_absolute(self, ref):
        """
        Compute the absolute precision of the given element and cache it

        For internal use.
        """
        col = self._matrix[ref]
        if len(col) == 0:
            self._absolute_precisions[ref] = Infinity
        else:
            self._absolute_precisions[ref] = min( [ c.valuation() for c in col ] )

    def precision_absolute(self, x):
        """
        Return the absolute precision of the given element

        INPUT:

        - ``x`` -- the element whose absolute precision is requested

        NOTE:

        The absolute precision is obtained by projecting the precision
        module onto the line of coordinate ``dx``

        NOTE:

        This function is not meant to be called directly.
        You should prefer call the method :meth:`precision_absolute`
        of ``x`` instead.

        EXAMPLES::

            sage: R = ZpLF(2)
            sage: prec = R.precision()

            sage: x = R(1,10); x
            1 + O(2^10)
            sage: y = R(1,5); y
            1 + O(2^5)
            sage: z = x + y; z
            2 + O(2^5)
            sage: z.precision_absolute()  # indirect doctest
            5

        In some cases, the absolute precision returned by this function
        may be infinite::

            sage: y = R(1)
            sage: prec.precision_absolute(y)
            +Infinity

        However calling the method :meth:`absolute_precision` of the
        element itself reintroduces a cap::

            sage: y.precision_absolute()
            20
        """
        ref = weakref.ref(x)
        if not self._absolute_precisions.has_key(ref):
            self._compute_precision_absolute(ref)
        return self._absolute_precisions[ref]

    def precision_lattice(self, elements=None):
        """
        Return a matrix representing the precision lattice on a
        subset of elements.

        INPUT:

        - ``elements`` -- a list of elements or ``None`` (default: ``None``)

        EXAMPLES::

            sage: R = ZpLF(2, label='preclattice')
            sage: prec = R.precision()
            sage: x = R(1,10); y = R(1,5)
            sage: prec.precision_lattice()
            [1024    0]
            [   0   32]

            sage: u = x + y
            sage: v = x - y
            sage: prec.precision_lattice([u,v])
            [  32 2016]
            [   0 2048]

        If the precision module does not project to a lattice,
        an error is raised

            sage: prec.precision_lattice([x,y,u,v])
            Traceback (most recent call last):
            ...
            PrecisionError: the differential is not surjective

        Here is another example with matrices::

            sage: M = matrix(R, 2, 2, [R(3,5),R(7,5),R(1,5),R(11,1)])
            sage: N = M^10

        The next syntax provides as easy way to select an interesting
        subset of variables (the selected subset consists of the four
        entries of the matrix ``N``)::

            sage: prec.precision_lattice(N)
            [  2048    512  28160 230400]
            [     0   2048  14336 258048]
            [     0      0  65536  65536]
            [     0      0      0 262144]
        """
        if elements is None:
            elements = self._elements
        else:
            elements = list_of_padics(elements)
        n = len(self._elements)
        rows = [ ]; val = 0
        for ref in elements:
            col = self._matrix[ref]
            row = [ x.value() for x in col ]
            valcol = min([ x.valuation() for x in col ])
            if valcol < val: val = valcol
            row += (n-len(row)) * [ZZ(0)]
            rows.append(row)
        from sage.matrix.constructor import matrix
        M = matrix(rows).transpose()
        if val < 0:
            M *= self._p ** (-val)
        M = M.change_ring(ZZ)
        M.echelonize()
        n = len(elements)
        if len(M.pivots()) < n:
            raise PrecisionError("the differential is not surjective")
        for i in range(n):
            v = M[i,i].valuation(self._p)
            M[i,i] = self._p ** v
        M.echelonize()
        M = M.submatrix(0,0,n,n)
        if val < 0:
            M *= self._p ** val
        return M
