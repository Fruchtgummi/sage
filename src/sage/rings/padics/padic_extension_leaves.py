"""
p-Adic Extension Leaves

The final classes for extensions of Zp and Qp (ie classes that are not
just designed to be inherited from).

AUTHORS:

- David Roe

- Julian Rueth (2012-10-24): added two step extensions
"""
#*****************************************************************************
#       Copyright (C) 2007-2009 David Roe <roed@math.harvard.edu>
#                     2012      Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from pow_computer_ext import PowComputer_ext_maker
from sage.libs.ntl.ntl_ZZ_pX import ntl_ZZ_pX

from unramified_extension_generic import UnramifiedExtensionGeneric
from eisenstein_extension_generic import EisensteinExtensionGeneric
from two_step_extension_generic import TwoStepExtensionGeneric
from general_extension_generic import GeneralExtensionGeneric

from generic_nodes import pAdicCappedRelativeRingGeneric, \
                          pAdicCappedRelativeFieldGeneric, \
                          pAdicCappedAbsoluteRingGeneric, \
                          pAdicFixedModRingGeneric

from padic_ZZ_pX_FM_element import pAdicZZpXFMElement
from padic_ZZ_pX_CR_element import pAdicZZpXCRElement
from padic_ZZ_pX_CA_element import pAdicZZpXCAElement

from padic_laurent_FM_element import pAdicLaurentFMElement
from padic_laurent_CR_element import pAdicLaurentCRElement
from padic_laurent_CA_element import pAdicLaurentCAElement

from padic_general_FM_element import pAdicGeneralFMElement
from padic_general_CR_element import pAdicGeneralCRElement

class UnramifiedExtensionRingCappedRelative(UnramifiedExtensionGeneric, pAdicCappedRelativeRingGeneric):
    """
    TESTS::

        sage: R.<a> = ZqCR(27,10000); R == loads(dumps(R))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        """
        A capped relative representation of Zq.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Zp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R.<a> = ZqCR(27,10000); R #indirect doctest
            Unramified Extension of 3-adic Ring with capped relative precision 10000 in a defined by (1 + O(3^10000))*x^3 + (O(3^10000))*x^2 + (2 + O(3^10000))*x + (1 + O(3^10000))

            sage: R.<a> = ZqCR(next_prime(10^30)^3, 3); R.prime()
            1000000000000000000000000000057
        """
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**prec)
        if prec <= 30:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), prec, prec, prec, False, ntl_poly, "small", "u")
        else:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), 30, prec, prec, False, ntl_poly, "big", "u")
        self._shift_seed = None
        self._pre_poly = prepoly
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXCRElement)

class UnramifiedExtensionFieldCappedRelative(UnramifiedExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    """
    TESTS::

        sage: R.<a> = QqCR(27,10000); R == loads(dumps(R))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        """
        A representation of Qq.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Qp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R.<a> = Qq(27,10000); R #indirect doctest
            Unramified Extension of 3-adic Field with capped relative precision 10000 in a defined by (1 + O(3^10000))*x^3 + (O(3^10000))*x^2 + (2 + O(3^10000))*x + (1 + O(3^10000))

            sage: R.<a> = Qq(next_prime(10^30)^3, 3); R.prime()
            1000000000000000000000000000057
        """
        # Currently doesn't support polynomials with non-integral coefficients
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**prec)
        if prec <= 30:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), prec, prec, prec, True, ntl_poly, "small", "u")
        else:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), 30, prec, prec, True, ntl_poly, "big", "u")
        self._shift_seed = None
        self._pre_poly = prepoly
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXCRElement)

class UnramifiedExtensionRingCappedAbsolute(UnramifiedExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    """
    TESTS::

        sage: R.<a> = ZqCA(27,10000); R == loads(dumps(R))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        """
        A capped absolute representation of Zq.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Zp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R.<a> = ZqCA(27,10000); R #indirect doctest
            Unramified Extension of 3-adic Ring with capped absolute precision 10000 in a defined by (1 + O(3^10000))*x^3 + (2 + O(3^10000))*x + (1 + O(3^10000))

            sage: R.<a> = ZqCA(next_prime(10^30)^3, 3); R.prime()
            1000000000000000000000000000057
        """
        # Currently doesn't support polynomials with non-integral coefficients
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**prec)
        if prec <= 30:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), prec, prec, prec, True, ntl_poly, "small", "u")
        else:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), 30, prec, prec, True, ntl_poly, "big", "u")
        self._shift_seed = None
        self._pre_poly = prepoly
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXCAElement)

class UnramifiedExtensionRingFixedMod(UnramifiedExtensionGeneric, pAdicFixedModRingGeneric):
    """
    TESTS::

        sage: R.<a> = ZqFM(27,10000); R == loads(dumps(R))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        """
        A fixed modulus representation of Zq.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Qp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R.<a> = ZqFM(27,10000); R #indirect doctest
            Unramified Extension of 3-adic Ring of fixed modulus 3^10000 in a defined by (1 + O(3^10000))*x^3 + (2 + O(3^10000))*x + (1 + O(3^10000))

            sage: R.<a> = ZqFM(next_prime(10^30)^3, 3); R.prime()
            1000000000000000000000000000057
        """
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**prec)
        self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), max(min(prec - 1, 30), 1), prec, prec, False, ntl_poly, "FM", "u")
        self._shift_seed = None
        self._pre_poly = prepoly
        UnramifiedExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXFMElement)

class EisensteinExtensionRingCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeRingGeneric):
    """
    TESTS::

        sage: R = Zp(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
        sage: W.<w> = R.ext(f); W == loads(dumps(W))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        """
        A capped relative representation of an eisenstein extension of Zp.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Zp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R = Zp(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
            sage: W.<w> = R.ext(f); W #indirect doctest
            Eisenstein Extension of 3-adic Ring with capped relative precision 10000 in w defined by (1 + O(3^10000))*x^3 + (O(3^10001))*x^2 + (3^2 + O(3^10001))*x + (-3 + O(3^10001))
            sage: W.precision_cap()
            30000

            sage: R.<p> = Zp(next_prime(10^30), 3, print_pos=False); S.<x> = ZZ[]; f = x^3 + p^2*x - p
            sage: W.<w> = R.ext(f); W.prime()
            1000000000000000000000000000057
            sage: W.precision_cap()
            9
        """
        unram_prec = (prec + poly.degree() - 1) // poly.degree()
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**unram_prec)
        shift_seed = self._compute_shift_seed(prepoly, poly.base_ring())
        shift_poly = ntl_ZZ_pX([a.lift() for a in shift_seed.list()], shift_seed.base_ring().prime()**unram_prec)
        if unram_prec <= 30:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), unram_prec, unram_prec, prec, False, ntl_poly, "small", "e", shift_poly)
        else:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), 30, unram_prec, prec, False, ntl_poly, "big", "e", shift_poly)
        self._shift_seed = shift_seed
        self._pre_poly = prepoly
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXCRElement)

class EisensteinExtensionFieldCappedRelative(EisensteinExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    """
    TESTS::

        sage: R = Qp(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
        sage: W.<w> = R.ext(f); W == loads(dumps(W))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        """
        A capped relative representation of an eisenstein extension of Qp.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Qp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R = Qp(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
            sage: W.<w> = R.ext(f); W #indirect doctest
            Eisenstein Extension of 3-adic Field with capped relative precision 10000 in w defined by (1 + O(3^10000))*x^3 + (O(3^10001))*x^2 + (3^2 + O(3^10001))*x + (-3 + O(3^10001))
            sage: W.precision_cap()
            30000

            sage: R.<p> = Qp(next_prime(10^30), 3, print_pos=False); S.<x> = ZZ[]; f = x^3 + p^2*x - p
            sage: W.<w> = R.ext(f); W.prime()
            1000000000000000000000000000057
            sage: W.precision_cap()
            9
        """
        # Currently doesn't support polynomials with non-integral coefficients
        unram_prec = (prec + poly.degree() - 1) // poly.degree()
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**unram_prec)
        shift_seed = self._compute_shift_seed(prepoly, poly.base_ring())
        shift_poly = ntl_ZZ_pX([a.lift() for a in shift_seed.list()], shift_seed.base_ring().prime()**unram_prec)
        if unram_prec <= 30:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), unram_prec, unram_prec, prec, True, ntl_poly, "small", "e", shift_poly)
        else:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), 30, unram_prec, prec, True, ntl_poly, "big", "e", shift_poly)
        self._shift_seed = shift_seed
        self._pre_poly = prepoly
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXCRElement)

class EisensteinExtensionRingCappedAbsolute(EisensteinExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    """
    TESTS::

        sage: R = ZpCA(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
        sage: W.<w> = R.ext(f); W == loads(dumps(W))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        """
        A capped absolute representation of an eisenstein extension of Zp.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Zp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R = ZpCA(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
            sage: W.<w> = R.ext(f); W #indirect doctest
            Eisenstein Extension of 3-adic Ring with capped absolute precision 10000 in w defined by (1 + O(3^10000))*x^3 + (3^2 + O(3^10000))*x + (-3 + 3^10000 + O(3^10000))
            sage: W.precision_cap()
            30000

            sage: R.<p> = ZpCA(next_prime(10^30), 3, print_pos=False); S.<x> = ZZ[]; f = x^3 + p^2*x - p
            sage: W.<w> = R.ext(f); W.prime()
            1000000000000000000000000000057
            sage: W.precision_cap()
            6
        """
        unram_prec = (prec + poly.degree() - 1) // poly.degree()
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**unram_prec)
        shift_seed = self._compute_shift_seed(prepoly, poly.base_ring())
        shift_poly = ntl_ZZ_pX([a.lift() for a in shift_seed.list()], shift_seed.base_ring().prime()**unram_prec)
        if unram_prec <= 30:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), unram_prec, unram_prec, prec, False, ntl_poly, "small", "e", shift_poly)
        else:
            self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), 30, unram_prec, prec, False, ntl_poly, "big", "e", shift_poly)
        self._shift_seed = shift_seed
        self._pre_poly = prepoly
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXCAElement)

class EisensteinExtensionRingFixedMod(EisensteinExtensionGeneric, pAdicFixedModRingGeneric):
    """
    TESTS::

        sage: R = ZpFM(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
        sage: W.<w> = R.ext(f); W == loads(dumps(W))
        True
    """
    def __init__(self, prepoly, poly, prec, halt, print_mode, names, shift_seed):
        """
        A fixed modulus representation of an eisenstein extension of Zp.

        INPUTS::

            - prepoly -- The original polynomial defining the
              extension.  This could be a polynomial with integer
              coefficients, for example, while poly has coefficients
              in Zp.

            - poly -- The polynomial with coefficients in
              self.base_ring() defining this extension.

            - prec -- The precision cap of this ring.

            - halt -- unused

            - print_mode -- A dictionary of print options.

            - names -- a 4-tuple, (variable_name, residue_name, unramified_subextension_variable_name, uniformizer_name)

        EXAMPLES::

            sage: R = ZpFM(3, 10000, print_pos=False); S.<x> = ZZ[]; f = x^3 + 9*x - 3
            sage: W.<w> = R.ext(f, check=False); W #indirect doctest
            Eisenstein Extension of 3-adic Ring of fixed modulus 3^10000 in w defined by (1 + O(3^10000))*x^3 + (3^2 + O(3^10000))*x + (-3 + 3^10000 + O(3^10000))
            sage: W.precision_cap()
            30000

            sage: R.<p> = ZpFM(next_prime(10^30), 3, print_pos=False); S.<x> = ZZ[]; f = x^3 + p^2*x - p
            sage: W.<w> = R.ext(f); W.prime()
            1000000000000000000000000000057
            sage: W.precision_cap()
            9
        """
        unram_prec = (prec + poly.degree() - 1) // poly.degree()
        ntl_poly = ntl_ZZ_pX([a.lift() for a in poly.list()], poly.base_ring().prime()**unram_prec)
        shift_seed = self._compute_shift_seed(prepoly, poly.base_ring())
        shift_poly = ntl_ZZ_pX([a.lift() for a in shift_seed.list()], shift_seed.base_ring().prime()**unram_prec)
        # deal with prec not a multiple of e better.
        self.prime_pow = PowComputer_ext_maker(poly.base_ring().prime(), max(min(unram_prec - 1, 30), 1), unram_prec, prec, False, ntl_poly, "FM", "e", shift_poly)
        self._shift_seed = shift_seed
        self._pre_poly = prepoly
        EisensteinExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicZZpXFMElement)

class TwoStepExtensionRingCappedRelative(TwoStepExtensionGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        TwoStepExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicLaurentCRElement)

class TwoStepExtensionFieldCappedRelative(TwoStepExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        TwoStepExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicLaurentCRElement)

class TwoStepExtensionRingCappedAbsolute(TwoStepExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        TwoStepExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicLaurentCAElement)

class TwoStepExtensionRingFixedMod(TwoStepExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        TwoStepExtensionGeneric.__init__(self, poly, prec, print_mode, names, pAdicLaurentFMElement)

class GeneralExtensionFieldCappedRelative(GeneralExtensionGeneric, pAdicCappedRelativeFieldGeneric):
    def __init__(self, prepoly, poly, prec, halt, print_mode, names):
        self.prime_pow = None # general extensions do not use the pow_computer yet
        GeneralExtensionGeneric.__init__(self, prepoly, poly, prec, print_mode, names, pAdicGeneralCRElement)

class GeneralExtensionRingCappedRelative(GeneralExtensionGeneric, pAdicCappedRelativeRingGeneric):
    def __init__(self, prepoly, poly, abs_ring, to_abs_ring, to_abs_ring_base, from_abs_ring, prec, halt, print_mode, names):
        self.prime_pow = None # general extensions do not use the pow_computer yet
        GeneralExtensionGeneric.__init__(self, poly, abs_ring, to_abs_ring, to_abs_ring_base, from_abs_ring, prec, print_mode, names, pAdicGeneralCRElement)

class GeneralExtensionRingCappedAbsolute(GeneralExtensionGeneric, pAdicCappedAbsoluteRingGeneric):
    def __init__(self, prepoly, poly, abs_ring, to_abs_ring, to_abs_ring_base, from_abs_ring, prec, halt, print_mode, names):
        self.prime_pow = None # general extensions do not use the pow_computer yet
        GeneralExtensionGeneric.__init__(self, poly, abs_ring, to_abs_ring, to_abs_ring_base, from_abs_ring, prec, print_mode, names, pAdicGeneralCRElement)

class GeneralExtensionRingFixedMod(GeneralExtensionGeneric, pAdicFixedModRingGeneric):
    def __init__(self, prepoly, poly, abs_ring, to_abs_ring, to_abs_ring_base, from_abs_ring, prec, halt, print_mode, names):
        self.prime_pow = None # general extensions do not use the pow_computer yet
        GeneralExtensionGeneric.__init__(self, poly, abs_ring, to_abs_ring, to_abs_ring_base, from_abs_ring, prec, print_mode, names, pAdicGeneralFMElement)
