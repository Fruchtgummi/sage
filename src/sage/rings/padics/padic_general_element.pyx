r"""
Elements of general extensions of p-adic rings, realized as polynomials over
the base ring.

AUTHORS:

- Julian Rueth (2013-01-08): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.padics.padic_ext_element cimport pAdicExtElement
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.integer import Integer
from sage.structure.element import is_Element
from sage.structure.element cimport RingElement, ModuleElement, CommutativeRingElement, Element
from sage.rings.padics.precision_error import PrecisionError

cdef class pAdicGeneralElement(pAdicExtElement):
    r"""
    Shared functionality of all elements of general extensions of p-adic rings.

    EXAMPLES::

        sage: K = ZpFM(3,10)
        sage: R.<x> = K[]
        sage: M.<x> = pAdicExtension(K, x - 1)
        sage: x
        1 + O(3^10)
        sage: M.zero()
        O(3^10)

    Elements can be initialized from lists of elements of the base field::

        sage: M([K.zero(),K.one(),K.one()])
        2 + O(3^10)

    Elements can be initialized from elements of the base field::

        sage: M(K.one())
        1 + O(3^10)
        sage: M(K(3))
        3 + O(3^10)
        sage: M(K(3^10))
        O(3^10)

    """
    def __init__(self, parent, x, absprec = None):
        """
        Initializes ``self`` from ``x``.

        INPUT:

            - ``parent`` -- the parent of ``self``

            - ``x`` -- an element from which ``self`` can be initialized

            - ``absprec`` -- integer or ``None`` (default: None), reduce the
              precision of ``self`` to that value if not ``None``

        AUTHORS:

        - Julian Rueth (2013-01-08): initial version

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: M(None)
            O(3^10)
            sage: M(None,3)
            O(3^10)

        """
        pAdicExtElement.__init__(self, parent)

        if x is None:
            self._element = parent.implementation_ring().zero()
        elif is_Element(x) and x.parent() is self.parent():
            self._element = x._element
        elif is_Element(x) and x.parent() is self.parent().implementation_ring():
            self._element = x
        elif is_Element(x) and x.parent() is self.parent().base_ring():
            self._element = parent._implementation_ring()[2](parent._polynomial_ring(x))
        elif is_Element(x) and x.parent() is self.parent().residue_field():
            self._element = parent.implementation_ring()(x)
        elif isinstance(x,list):
            if all([c in self.parent().base_ring() for c in x]):
                self._element = parent._implementation_ring()[2](parent._polynomial_ring(x))
            else:
                raise NotImplementedError("initialization from %s"%(x))
        else:
            raise NotImplementedError("initialization to %s from %s (%s)"%(self.parent(),x,x.parent()))

        assert self._element.parent() is parent.implementation_ring()

        if absprec is not None:
            self._element = self._element.add_bigoh(absprec)

    def precision_absolute(self):
        """
        The precision to which ``self`` is accurate.

        OUTPUT:

        An integer or infinity.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x.precision_absolute()
            10
            sage: M(0,3).precision_absolute()
            10
            sage: M.zero().precision_absolute()
            10

        """
        return self._element.precision_absolute()

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Compute the sum of ``self`` and ``right``.

        INPUT:

            - ``right`` -- an element of the same parent as ``self``

        OUTPUT:

        An element in the same parent as ``self``, the sum of ``self`` and
        ``right``.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x + x
            2 + O(3^10)

        """
        ret = self.parent()(None)
        ret._element = self._element + right._element
        return ret

    cpdef ModuleElement _neg_(self):
        """
        Compute ``-self``.

        OUTPUT:

        An element of the same parent as ``self``.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: -x
            2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10)

        """
        ret = self.parent()(None)
        ret._element = -self._element
        return ret

    def is_unit(self):
        """
        Returns whether this is a unit in its ring.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x.is_unit()
            True
            sage: (x-1).is_unit()
            False
            sage: M(None).is_unit()
            False

        """
        return self._element.is_unit()

    def __invert__(self):
        """
        Compute a multiplicative inverse of this element.

        OUTPUT:

        An element in the parent of ``self``. Raises an ``NotImplementedError``
        if ``self`` is not a unit.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: ~x
            1 + O(3^10)
            sage: ~(x-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot invert non-units

        """
        if not self.is_unit():
            raise NotImplementedError("cannot invert non-units")
        ret = self.parent()(None)
        ret._element = ~self._element
        return ret

    cpdef RingElement _div_(self, RingElement other):
        """
        Return the result of ``self`` divided by ``other``.

        INPUT:

            - ``other`` -- an element in the same ring as ``self``

        OUTPUT:

        An element in the same ring as ``self``. Raises a
        ``NotImplementedError`` if ``other`` is not invertible in the parent of
        ``other``.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x/x
            1 + O(3^10)
            sage: x/(x-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot invert non-units

        """
        return self * ~other

    def __nonzero__(self):
        """
        Returns whether this element is non-zero.

        OUTPUT:

        Returns False if this element is either an exact zero or an inexact
        zero.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x.is_zero() # indirect doctest
            False
            sage: (x-1).is_zero()
            True
            sage: M(0,3).is_zero()
            True

        """
        return not self._element.is_zero()

    cpdef RingElement _mul_(self, RingElement other):
        """
        Compute the multiplication of ``self`` by ``other``.

        INPUT:

            - ``other`` -- an element in the same ring as ``self``

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x*x
            1 + O(3^10)
            sage: x*3
            3 + O(3^10)

        """
        ret = self.parent()(None)
        ret._element = self._element*other._element
        return ret

    def __lshift__(self, shift):
        """
        Multiply this element by ``pi^shift`` with ``pi`` being a fixed
        uniformizer of the extension.

        INPUT:

            - ``shift`` -- an integer, if the element is not defined in a
              field, then ``-shift`` must not exceed its valuation.

        OUTPUT:

        An element in the same ring.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x << 3
            3^3 + O(3^10)

        """
        ret = self.parent()(None)
        ret._element = self._element << shift
        return ret

    def __rshift__(self, shift):
        """
        Multiply this element by ``pi^-shift`` with ``pi`` being a fixed
        uniformizer of the extension.

        INPUT:

            - ``shift`` -- an integer, if the element is not defined in a
              field, then ``shift`` must not exceed its valuation. TODO: not true see below

        OUTPUT:

        An element in the same ring.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x >> 3
            O(3^10)

        """
        return self<<-shift

    def is_zero(self, absprec=None):
        """
        Returns whether this element is zero.

        INPUT:

            - ``absprec`` -- an integer, ``None`` or infinity (default: None),
              if ``absprec`` is ``None``, then this method returns whether this
              element is zero to the precision it is defined to, otherwise
              returns whether the element zero if considered to precision
              ``absprec``.

        OUTPUT:

        Raises an error if the element is not defined to the degree given by
        ``absprec``.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x.is_zero()
            False
            sage: (x<<3).is_zero(2)
            True
            sage: x.is_zero(infinity)
            False

        """
        return self._element.is_zero(absprec)

    def unit_part(self):
        """
        Returns an element ``u`` such that the element is of the form
        ``u*pi^v`` with ``pi`` being a fixed uniformizer of the extension.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: (x<<3).unit_part()
            1 + O(3^10)

        """
        ret = self.parent()(None)
        ret._element = self._element.unit_part()
        return ret

    def lift_to_precision(self, absprec=None):
        """
        Returns a lift of this element which is defined at least to precision
        ``absprec``.

        INPUT:

            - ``absprec`` -- an integer or infinity

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: M(0,3).lift_to_precision(5)
            O(3^10)

        """
        ret = self.parent()(None)
        ret._element = self._element.lift_to_precision(absprec)
        return ret

    def is_equal_to(self, right, absprec=None):
        """
        Returns whether the element is equal to ``right``.

        INPUT:

            - ``right`` -- an element that coerces to the parent ring

            - ``absprec`` -- an integer, infinity, or ``None`` (default:
              ``None``), if ``None``, return whether the elements are equal
              when considered to the smallest of their precisions, otherwise,
              return whether the elements are equal when considered to
              precision ``absprec``.

        OUTPUT:

        Raises an error if ``absprec`` exceeds the precision of the elements.

        EXAMPLES::

            sage: K = ZpFM(3,10)
            sage: R.<x> = K[]
            sage: M.<x> = pAdicExtension(K, x - 1)
            sage: x.is_equal_to(-x)
            False
            sage: (3*x).is_equal_to(-3*x, 1)
            True

        """
        return self._element.is_equal_to(right._element, absprec)

    def _repr_(self):
        return repr(self._element)
        polynomial = self.parent()._from_implementation_ring(self)
        if self.parent().degree() == 1:
            assert polynomial.degree()<=0
            ret = repr(polynomial[0])
        else:
            ret = " + ".join(['('+repr(c)+')'+("*"+polynomial.variable_name()+'^'+str(i) if i else '') for i,c in enumerate(polynomial.coeffs())])

        return "%s [approximate]\n %s [exact]"%(ret,repr(self._element))

    cdef long valuation_c(self):
        return (<pAdicGenericElement>self._element).valuation_c()

    cpdef bint _is_inexact_zero(self) except -1:
        return (<pAdicGenericElement>self._element)._is_inexact_zero()

    cpdef bint _is_exact_zero(self) except -1:
        return (<pAdicGenericElement>self._element)._is_exact_zero()

    cdef int _cmp_units(left, pAdicGenericElement right) except -2:
        return (<pAdicGenericElement>left._element)._cmp_units(<pAdicGenericElement>right._element)

    def residue(self, absprec=1):
        return self._element.residue(absprec)

    def precision_relative(self):
        return self._element.precision_relative()
