"""
Two step extensions of `\mathbb{Q}_p` and `\mathbb{Z}_p`

This file implements the shared functionality for two step extensions, i.e., Eisenstein extensions of unramified extensions.

AUTHORS:

- Julian Rueth (2012-10-22): initial version

"""
#*****************************************************************************
#       Copyright (C) 2012 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from padic_extension_generic import pAdicExtensionGeneric
from sage.rings.finite_rings.constructor import GF

class TwoStepExtensionGeneric(pAdicExtensionGeneric):
    """
    An Eisenstein extension of an unramified extension of `\mathbb{Q}_p` or `\mathbb{Z}_p`.
    """
    def __init__(self, poly, upoly, epoly, prec, print_mode, names, element_class):
        """
        Initializes ``self``.

        INPUT:

            - ``poly`` -- a tuple of ``upoly`` and ``epoly`` used to define ``self``

            - ``upoly`` -- a polynomial over the `p`-adic base ring which
              defines the unramified step of the extension

            - ``epoly`` -- a polynomial over a polynomial ring in the variable
              of the unramified step of the extension which defines the
              Eisenstein step of the extension

            - ``prec`` -- a positive integer, the precision to use when
              representing elements in ``self``

            - ``print_mode`` -- a dictionary of print options

            - ``names`` -- a tuple ``((u,a),ubar,u,a)`` where ``u`` is the
              variable of the unramified step of the extension, ``a`` is the
              variable of the Eisenstein step of the extension, and ``u0`` is
              the variable of the residue field extension.

            - ``element_class`` -- the class which implements the elements of
              ``self``

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import TwoStepExtensionRingFixedMod
            sage: K = ZpFM(3,10)
            sage: Ru.<u> = K[]
            sage: upoly = u^2 + 3*u + 4
            sage: Ra.<a> = Ru[]
            sage: epoly = a^3 - 9*u*a^2 + 3*u
            sage: M = TwoStepExtensionRingFixedMod((upoly, epoly), upoly, epoly, 30, None, {}, (('u','a'),'u0','u','a')) # indirect doctest

        """
        # do we need this?
        #fixed_prec = min([min([c.precision_absolute() for c in upoly.coeffs() if not c.is_zero()]),min([min([c.precision_absolute() for c in d]) for d in epoly.coeffs() if not c.is_zero()])])
        #if fixed_prec != upoly.base_ring().precision_cap():
        #    print "WARNING: dropping precision of base ring from %s to %s"%(upoly.base_ring().precision_cap(),fixed_prec)
        #    functor,param = upoly.base_ring().construction()
        #    functor.prec = fixed_prec
        #    upoly = upoly.change_ring(functor(param))
        #    epoly = epoly.change_ring(upoly.parent())
        #    prec = fixed_prec*epoly.degree()
        #    poly = (upoly,epoly)

        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        self._upoly = upoly
        from sage.rings.padics.factory import is_unramified
        assert is_unramified(upoly),"%s is not unramified"%upoly
        self._epoly = epoly
        self._prime = upoly.base_ring().prime()
        self._inertia_subring = upoly.base_ring().extension(upoly,names=(names[2],),res_name=names[1])
        from sage.rings.padics.factory import is_eisenstein
        assert is_eisenstein(epoly.change_ring(self._inertia_subring)),"%s is not Eisenstein"%epoly.change_ring(self._inertia_subring)
        self._res_field = GF(self.prime()**upoly.degree(), name = names[1], modulus = upoly.change_ring(upoly.base_ring().residue_field()))
        self._populate_coercion_lists_(coerce_list=[self._inertia_subring],element_constructor=element_class)

    def hom(self, im_gens):
        if len(im_gens)!=2:
            raise ValueError
        from sage.categories.morphism import SetMorphism
        from sage.categories.rings import Rings
        from sage.categories.homset import Hom
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(self.base_ring(),names=('T',))
        S = PolynomialRing(im_gens[0].parent(),names=('S',))
        return SetMorphism(Hom(self,im_gens[0].parent(),Rings()), lambda x:S([R(c.vector(base=self.base_ring()))(im_gens[0]) for c in  x.vector(base=self.inertia_subring())])(im_gens[1]))

    def _repr_(self):
        """
        Returns a print representation of this extension.

        EXAMPLES::

            sage: K = ZpCR(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M
            Eisenstein extension of unramified extension of 3-adic Ring with capped relative precision 10 in ('u', 'a') defined by ((1 + O(3^10))*u^2 + (3 + O(3^11))*u + (1 + 3 + O(3^10)), ((1 + O(3^10)))*a^3 + ((2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + O(3^12))*u)*a^2 + (3 + O(3^11))*u)

        """
        return "Eisenstein extension of unramified extension of %s in %s defined by %s"%(self.ground_ring(), self.variable_name(), self.defining_polynomial())

    def ramification_index(self):
        """
        Return the absolute ramification index of ``self``, i.e., the valuation
        of ``p`` in ``self``.

        EXAMPLES::

            sage: K = ZpCR(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M.ramification_index()
            3

        """
        return self._epoly.degree()

    def degree(self):
        """
        Returns the degree of this ring over the padic base ring.

        OUTPUT:

        A positive integer, the degree of the Eisenstein step times the degree
        of the unramified step.

        EXAMPLES::

            sage: K = ZpCR(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M.degree()
            6

        """
        return self._epoly.degree() * self._upoly.degree()

    def gen(self, n=0):
        """
        Returns a generator for self as an extension of its ground ring.

        EXAMPLES::

            sage: K = ZpCR(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M.gen()
            u + O(a^30)
            sage: M.gen(1)
            a + O(a^31)

        """
        if n == 0:
            return self([self.residue_field().gen()] + [self.residue_field().zero()]*(self.precision_cap()-1))
        if n == 1:
            return self([self.residue_field().zero(),self.residue_field().one()] + [self.residue_field().zero()]*(self.precision_cap()-1))
        raise IndexError, "only two generators"

    def gens(self):
        """
        Returns a list of generators of self.

        OUTPUT:

        Returns a list consisting of the generator of the Eisenstein extension
        followed by the generator of the unramified extension.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M.gens()
            [u + O(a^30), a + O(a^31)]

        """
        return [self.gen(0),self.gen(1)]

    def inertia_subring(self):
        """
        Returns the inertia subring of ``self``.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M.inertia_subring()
            Unramified Extension of 3-adic Field with capped relative precision 10 in u defined by (1 + O(3^10))*u^2 + (3 + O(3^10))*u + (1 + 3 + O(3^10))

        """
        return self._inertia_subring

    def residue_class_field(self):
        """
        Returns the residue class field.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M.residue_class_field()
            Finite Field in u0 of size 3^2

        """
        return self._res_field

    def prime(self):
        """
        Returns the characteristic of the residue field.

        EXAMPLES::

            sage: K = Qp(3,10)
            sage: Ru.<u> = K[]
            sage: Ra.<a> = Ru[]
            sage: M.<u,a> = pAdicExtension(K, (u^2 + 3*u + 4,a^3 - 9*u*a^2 + 3*u))
            sage: M.prime()
            3

        """
        return self._prime

    def uniformizer(self):
        return self.gens()[1]

    def uniformizer_pow(self, n):
        return self.uniformizer()**n
