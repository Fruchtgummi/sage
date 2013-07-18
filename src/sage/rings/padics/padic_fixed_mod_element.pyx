"""
`p`-Adic Fixed-Mod Elements

Elements of `p`-Adic Rings with Fixed Modulus

AUTHORS:

- David Roe
- Genya Zaytman: documentation
- David Harvey: doctests
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/libs/linkages/padics/mpz.pxi"
include "FM_template.pxi"

from sage.libs.pari.gen cimport PariInstance
cdef PariInstance P = sage.libs.pari.all.pari
from sage.rings.finite_rings.integer_mod import Mod

cdef class PowComputer_(PowComputer_base):
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field):
        _prec_type = 'fixed-mod'
        PowComputer_base.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field)

cdef class pAdicFixedModElement(FMElement):
    r"""
    INPUT:

    - ``parent`` -- a ``pAdicRingFixedMod`` object.

    - ``x`` -- input data to be converted into the parent.

    - ``absprec`` -- ignored; for compatibility with other `p`-adic rings

    - ``relprec`` -- ignored; for compatibility with other `p`-adic rings

    NOTES::

    The following types are currently supported for x:

    - Integers
    - Rationals -- denominator must be relatively prime to `p`
    - FixedMod `p`-adics
    - Elements of ``IntegerModRing(p^k)`` for ``k` less than or equal to the
      modulus

    The following types should be supported eventually:

    - Finite precision p-adics
    - Lazy p-adics
    - Elements of local extensions of THIS `p`-adic ring that actually lie in `\ZZ_p`

    EXAMPLES::

        sage: R = Zp(5, 20, 'fixed-mod', 'terse')

    Construct from integers::

        sage: R(3)
        3 + O(5^20)
        sage: R(75)
        75 + O(5^20)
        sage: R(0)
        0 + O(5^20)

        sage: R(-1)
        95367431640624 + O(5^20)
        sage: R(-5)
        95367431640620 + O(5^20)

    Construct from rationals::

        sage: R(1/2)
        47683715820313 + O(5^20)
        sage: R(-7875/874)
        9493096742250 + O(5^20)
        sage: R(15/425)
        Traceback (most recent call last):
        ...
        ValueError: p divides denominator

    Construct from IntegerMod::

        sage: R(Integers(125)(3))
        3 + O(5^20)
        sage: R(Integers(5)(3))
        3 + O(5^20)
        sage: R(Integers(5^30)(3))
        3 + O(5^20)
        sage: R(Integers(5^30)(1+5^23))
        1 + O(5^20)
        sage: R(Integers(49)(3))
        Traceback (most recent call last):
        ...
        TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

        sage: R(Integers(48)(3))
        Traceback (most recent call last):
        ...
        TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

    Some other conversions::

        sage: R(R(5))
        5 + O(5^20)

    # todo: doctests for converting from other types of p-adic rings
    """
    def lift(self):
        r"""
        Return an integer congruent to this element modulo the precision.

        .. WARNING::

            Since fixed modulus elements don't track their precision,
            the result may not be correct modulo
            `i^{\mbox{prec_cap}}` if the element was defined by
            constructions that lost precision.

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.lift()
            8
            sage: type(a.lift())
            <type 'sage.rings.integer.Integer'>
        """
        return self.lift_c()

    cdef lift_c(self):
        """
        Returns an integer congruent to this element modulo the precision.

        .. WARNING::

            Since fixed modulus elements don't track their precision,
            the result may not be correct modulo
            `i^{\mbox{prec_cap}}` if the element was defined by
            constructions that lost precision.

        EXAMPLES::

            sage: R = ZpFM(7,4); a = R(8); a.lift() # indirect doctest
            8
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def _pari_(self):
        """
        Conversion to PARI.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: pari(R(1777)) #indirect doctest
            2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Converts this element to an equivalent pari element.

        EXAMPLES::

            sage: R = ZpFM(5, 10); a = R(17); pari(a) #indirect doctest
            2 + 3*5 + O(5^10)
            sage: pari(R(0))
            O(5^10)
            sage: pari(R(0,5))
            O(5^10)
        """
        # holder is defined in the linkage file
        cdef long val = mpz_remove(holder.value, self.value, self.prime_pow.prime.value)
        return P.new_gen_from_padic(val, self.prime_pow.prec_cap - val,
                                    self.prime_pow.prime.value,
                                    self.prime_pow.pow_mpz_t_tmp(self.prime_pow.prec_cap - val)[0],
                                    holder.value)

    def _integer_(self, Z=None):
        """
        Returns an integer congruent to this element modulo the precision.

        .. WARNING::

            Since fixed modulus elements don't track their precision,
            the result may not be correct modulo
            `p^{\mbox{prec_cap}}` if the element was defined by
            constructions that lost precision.

        EXAMPLES::

            sage: R = ZpFM(5); R(-1)._integer_()
            95367431640624
        """
        return self.lift_c()

    def residue(self, absprec=1):
        r"""
        Reduce this element mod `p^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` - an integer (default: 1)

        OUTPUT:

        element of Z/(p^prec Z) -- self reduced mod p^prec

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.residue(1)
            1
        """
        cdef Integer selfvalue, modulus
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_sgn((<Integer>absprec).value) < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
        cdef long aprec = mpz_get_ui((<Integer>absprec).value)
        modulus = PY_NEW(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        selfvalue = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        return Mod(selfvalue, modulus)

    def multiplicative_order(self):
        r"""
        Returns the minimum possible multiplicative order of this element.

        OUTPUT:

        an integer -- the multiplicative order of this element.  This is the
        minimum multiplicative order of all elements of `\ZZ_p` lifting this
        element to infinite precision.

        EXAMPLES::

            sage: R = ZpFM(7, 6)
            sage: R(1/3)
            5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + 4*7^5 + O(7^6)
            sage: R(1/3).multiplicative_order()
            +Infinity
            sage: R(7).multiplicative_order()
            +Infinity
            sage: R(1).multiplicative_order()
            1
            sage: R(-1).multiplicative_order()
            2
            sage: R.teichmuller(3).multiplicative_order()
            6
        """
        cdef mpz_t tmp
        cdef Integer ans
        if mpz_divisible_p(self.value, self.prime_pow.prime.value):
            return infinity
        if mpz_cmp_ui(self.value, 1) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 1)
            return ans
        mpz_init(tmp)
        mpz_sub_ui(tmp, self.prime_pow.pow_mpz_t_top()[0], 1)
        if mpz_cmp(self.value, tmp) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(tmp, self.value, self.prime_pow.prime.value, self.prime_pow.pow_mpz_t_top()[0])
        if mpz_cmp(tmp, self.value) == 0:
            mpz_clear(tmp)
            return self.residue(1).multiplicative_order()
        else:
            mpz_clear(tmp)
            return infinity

    def padded_list(self, n, list_mode = 'simple'):
        """
        Returns a list of coefficients of p starting with $p^0$ up to
        $p^n$ exclusive (padded with zeros if needed)

        INPUT::

            - self -- a p-adic element
            - n -- an integer

        OUTPUT::

            - list -- the list of coefficients of self

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(2*7+7**2); a.padded_list(5)
            [0, 2, 1, 0, 0]

        NOTE::

            For elements with positive valuation, this function will
            return a list with leading 0s, unlike for field elements.

            The slice operators throw an error if asked for a slice
            above the precision, while this function works
        """
        if list_mode == 'simple' or list_mode == 'smallest':
            zero = Integer(0)
        else:
            zero = self.parent()(0)
        L = self.list()
        return L[:n] + [zero] * (n - len(L))

    def precision_absolute(self):
        """
        Returns the absolute precision of self.

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- the absolute precision of self

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_absolute()
            4
        """
        return self.parent().precision_cap()

    def precision_relative(self):
        r"""
        Returns the relative precision of self

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- the relative precision of self

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_relative()
            3
            sage: a = R(0); a.precision_relative()
            0
        """
        cdef unsigned long diff
        cdef Integer ans
        ans = PY_NEW(Integer)
        diff = self.prime_pow.prec_cap - self.valuation_c()
        mpz_set_si(ans.value, diff)
        return ans

    def residue(self, absprec=1):
        r"""
        Reduces ``self`` modulo `p^\mathrm{absprec}`.

        INPUT:

        - ``absprec`` - a non-negative integer (default: 1)

        OUTPUT:

        ``self`` reduced modulo `p^\mathrm{absprec}` as an element of
        `\mathbb{Z}/p^\mathrm{absprec}\mathbb{Z}`

        EXAMPLES::

            sage: R = Zp(7,4,'fixed-mod')
            sage: a = R(8)
            sage: a.residue(1)
            1
            sage: a.residue(2)
            8

        TESTS::

            sage: R = Zp(7,4,'fixed-mod')
            sage: a = R(8)
            sage: a.residue(0)
            0
            sage: a.residue(-1)
            Traceback (most recent call last):
            ...
            ValueError: Cannot reduce modulo a negative power of p.
            sage: a.residue(5)
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision known in order to compute residue.

        """
        cdef Integer selfvalue, modulus
        selfvalue = PY_NEW(Integer)
        modulus = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        cdef unsigned long aprec
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if absprec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif absprec < 0:
            raise ValueError, "Cannot reduce modulo a negative power of p."
        else:
            aprec = mpz_get_ui((<Integer>absprec).value)
        if aprec > self.prime_pow.prec_cap:
            sig_on()
            mpz_pow_ui(modulus.value, self.prime_pow.prime.value, aprec)
            sig_off()
        else:
            mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        return Mod(selfvalue, modulus)

    #def square_root(self):
    #    r"""
    #    Returns the square root of this p-adic number

    #    INPUT:
    #        self -- a p-adic element
    #    OUTPUT:
    #        p-adic element -- the square root of this p-adic number

    #        The square root chosen is the one whose reduction mod p is in
    #        the range [0, p/2).

    #        Note that because this is a fixed modulus ring, garbage digits
    #        may be introduced, if either
    #        (a) the valuation of the input is positive, or
    #        (b) p = 2.

    #        If no square root exists, a ValueError is raised.
    #        (This may be changed later to return an element of an extension
    #        field.)

    #    EXAMPLES:
    #        sage: R = Zp(3,20,'fixed-mod')
    #        sage: R(0).square_root()
    #            O(3^20)
    #        sage: R(1).square_root()
    #            1 + O(3^20)
    #        sage: R(2).square_root()
    #        Traceback (most recent call last):
    #        ...
    #        ValueError: element is not a square
    #        sage: R(4).square_root() == R(-2)
    #            True
    #        sage: R(9).square_root()
    #            3 + O(3^20)
    #        sage: R2 = Zp(2,20,'fixed-mod')
    #        sage: R2(0).square_root()
    #            O(2^20)
    #        sage: R2(1).square_root()
    #            1 + O(2^20)
    #        sage: R2(4).square_root()
    #            2 + O(2^20)
    #        sage: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
    #            True
    #        sage: R2(17).square_root()
    #            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^20)
    #        sage: R3 = Zp(5,20,'fixed-mod', 'terse')
    #        sage: R3(0).square_root()
    #            0 + O(5^20)
    #        sage: R3(1).square_root()
    #            1 + O(5^20)
    #        sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
    #            True
    #    """
    #    #todo: make more efficient
    #    try:
    #        # use pari
    #        return self.parent()(pari(self).sqrt())
    #    except PariError:
    #        # todo: should eventually change to return an element of
    #        # an extension field
    #        raise ValueError, "element is not a square"

    def valuation(self, prime=None):
        """
        Returns the valuation of self.

        If self is zero, the valuation returned is the precision of the ring.

        INPUT::

            - self -- a p-adic element

        OUTPUT::

            - integer -- the valuation of self.

        EXAMPLES::

            sage: R = Zp(17, 4,'fixed-mod')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'fixed-mod')
            sage: R(0).valuation()
            4
            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.valuation_c())
        return ans

    cdef long valuation_c(self):
        """
        Returns the valuation of self.

        EXAMPLES::

            sage: R = ZpFM(5, 5); R(0).valuation() #indirect doctest
            5
        """
        if mpz_sgn(self.value) == 0:
            return self.prime_pow.prec_cap
        cdef mpz_t tmp
        cdef long ans
        mpz_init(tmp)
        ans = mpz_remove(tmp, self.value, self.prime_pow.prime.value)
        mpz_clear(tmp)
        return ans

    cpdef val_unit(self):
        """
        Returns a 2-tuple, the first element set to the valuation of
        self, and the second to the unit part of self.

        If self == 0, then the unit part is O(p^self.parent().precision_cap()).

        EXAMPLES::

            sage: R = ZpFM(5,5)
            sage: a = R(75); b = a - a
            sage: a.val_unit()
            (2, 3 + O(5^5))
            sage: b.val_unit()
            (5, O(5^5))
        """
        cdef Integer val
        cdef pAdicFixedModElement unit
        if mpz_sgn(self.value) == 0:
            return (self.parent().precision_cap(), self)
        val = PY_NEW(Integer)
        unit = self._new_c()
        mpz_set_ui(val.value, mpz_remove(unit.value, self.value, self.prime_pow.prime.value))
        return (val, unit)

def make_pAdicFixedModElement(parent, value):
    """
    Unpickles a fixed modulus element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_fixed_mod_element import make_pAdicFixedModElement
        sage: R = ZpFM(5)
        sage: a = make_pAdicFixedModElement(R, 17*25); a
        2*5^2 + 3*5^3 + O(5^20)
    """
    return unpickle_fme_v2(pAdicFixedModElement, parent, value)
