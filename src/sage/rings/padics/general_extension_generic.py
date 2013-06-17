"""
General extensions of `\mathbb{Q}_p` and `\mathbb{Z}_p`

This file implements the shared functionality for general extensions of p-adic rings.

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

from sage.misc.cachefunc import cached_method

from padic_extension_generic import pAdicExtensionGeneric

class GeneralExtensionGeneric(pAdicExtensionGeneric):
    """
    An extension of a p-adic ring, defined by an irreducible polynomial.
    """
    def __init__(self, prepoly, poly, prec, print_mode, names, element_class):
        """
        Initializes ``self``.

        INPUT:

            TODO

            - ``poly`` -- an irreducible polynomial over a p-adic ring

            - ``abs_ring`` -- a ring isomorphic to this ring which will be used
              to perform all calculations

            - ``to_abs_ring`` -- the image of the generator of this ring in ``abs_ring``

            - ``to_abs_ring_base`` -- a ring homomorphism from the base ring to
              ``abs_ring`` which, together with ``to_abs_ring`` defines an
              isomorphism

            - ``from_abs_ring`` -- a tuple consisting of the images of the
              generators of ``abs_ring`` in this ring (as a polynomial or
              rational function over the base ring)

            - ``prec`` -- a positive integer, the precision to use when
              representing elements this ring

            - ``print_mode`` -- a dictionary of print options

            - ``names`` -- a tuple ``(x,ubar)`` where ``x`` is a root of
              ``poly`` and ``ubar`` is the variable of the residue field.

            - ``element_class`` -- the class which implements the elements of
              this ring

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3')) # indirect doctest

        """
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        self._poly = poly
        self._polynomial_ring = poly.parent().change_var(names[0])
        self._prime = poly.base_ring().prime()

    @cached_method
    def _implementation_ring(self):
        """
        Return a ring isomorphic to this ring which implements the arithmetic
        in this ring.

        OUTPUT:

        A triple ``(R,from,to)`` where ``R`` is a ring, ``from`` is a morphism
        from ``R`` to this ring, and ``to`` is a morphism from this ring to
        ``R``.

        """
        return GeneralExtensionGeneric.__implementation_ring(self.base_ring(), self._polynomial_ring, self._poly)

    @cached_method
    def __implementation_ring(base, polynomial_ring, modulus)
        """
        Return a ring isomorphic to the extension of ``modulus`` over ``base``.

        INPUT:

        - ``base`` -- the base ring of ``modulus``

        - ``polynomial_ring`` -- the ring used to represent element in the
          extension

        - ``modulus`` -- a polynomial over ``base``

        OUTPUT:

        A triple as in :meth:`_implementation_ring`.

        """
        if any([c.valuation()<0 for c in modulus.list()]):
            raise NotImplementedError("minimal polynomial is not integral after normalization")

        from sage.rings.padics.factory import GenericExtensionFactory
        if GenericExtensionFactory.is_trivial(modulus):
            return GeneralExtensionGeneric.__implementation_ring_trivial(self.base_ring(), self._polynomial_ring, self._poly)
        elif not isinstance(base, TwoStepExtensionGeneric):
            # reduce to the case where base is an Eisenstein over unramified extension
            base2, base2_to_base, base_to_base2 = base.two_step_ring()
            modulus_over2 = modulus.map_coefficients(base_to_base2, base2)
            polynomial_ring2 = polynomial_ring.change_ring(base2)

            implementation_ring, from_implementation_ring, to_implementation_ring = GeneralExtensionGeneric.__implementation_ring(base2, polynomial_ring2, modulus_over2)

            from_implementation_ring = lambda f:from_implementation_ring(f).map_coefficients(base2_to_base, base)
            to_implementation_ring = lambda f: to_implementation_ring(f.map_coefficients(base_to_base2, base2))
            return implementation_ring, from_implementation_ring, to_implementation_ring
        else:
            if GenericExtensionFactory.is_eisenstein(modulus):
                return GeneralExtensionGeneric.__implementation_ring_eisenstein(base, polynomial_ring, modulus)
            elif GeneralExtensionGeneric.is_unramified(modulus):
                return GeneralExtensionGeneric.__implementation_ring_unramified(base, polynomial_ring, modulus)
            elif GeneralExtensionGeneric.is_totally_ramified(modulus):
                return GeneralExtensionGeneric.__implementation_ring_totally_ramified(base, polynomial_ring, modulus)
            else:
                return GeneralExtensionGeneric.__implementation_ring_general(base, polynomial_ring, modulus)

    @staticmethod
    def __implementation_ring_trivial(base, polynomial_ring, modulus):
        assert modulus.degree() == 1
        assert modulus.is_monic()
        implementation_ring = base
        from_implementation_ring = polynomial_ring
        to_implementation_ring = lambda f:f(-modulus[0])
        return base, from_implementation_ring, to_implementation_ring

    @staticmethod
    def __implementation_ring_totally_ramified(base, polynomial_ring, modulus):
        from padic_valuation import pAdicValuation
        v = pAdicValuation(base)
        is_totally_ramified, ramification_steps = v.is_totally_ramified(modulus, include_steps=True, assume_squarefree=True)
        assert is_totally_ramified
        assert any([v(v.phi()).denominator() == modulus.degree() for v in ramification_steps])
        slopes = [v(v.phi()) for v in ramification_steps]
        keys = [v.phi() for v in ramification_steps]
        numerators = [(slope.numerator()*modulus.degree()//slope.denominator(), key) for slope,key in zip(slopes,keys)]
        numerators.append((-modulus.degree(),base(base.prime())**-1))
        bfs = {0:[]}
        while not 1 in bfs:
            for onum in bfs.keys():
                for num, key in numerators:
                    nnum = onum+num
                    if nnum not in bfs or len(bfs[onum])+1 < len(bfs[nnum]):
                        bfs[nnum] = bfs[onum]+[key]

        uniformizer = modulus.parent().quo(modulus).one()
        for key in bfs[1]:
            uniformizer*=key

        epoly = uniformizer.minpoly()

        assert epoly.degree() == modulus.degree(), epoly
        assert epoly[0].valuation() == 1, "%s is not Eisenstein"%epoly

        implementation_ring = base.extension(epoly,names=epoly.parent().variable_name())

        from_implementation_ring = uniformizer.lift()(self.gen())
        from_implementation_ring = lambda f: f.polynomial()(from_implementation_ring)

        to_implementation_ring = modulus.change_ring(implementation_ring).any_root()
        to_implementation_ring = lambda f:f(to_implementation_ring)

        return implementation_ring, from_implementation_ring, to_implementation_ring

    @staticmethod
    def __implementation_ring_eisenstein(base, polynomial_ring, modulus):
        if not isinstance(base, TwoStepExtensionGeneric):
        else:
            # create the bivariate ring
            # compute the resultant
            raise NotImplementedError()

    def __implementation_ring_unramified(base, polynomial_ring, modulus):
        # extend the maximal unramified subextension of base and put its Eisenstein poly on top
        raise NotImplementedError()

    #def _mod_minpoly(primitive_element, total_degree, base):
    #    powers = [primitive_element**j for j in range(total_degree+1)]
    #    powers = [pow.lift() for pow in powers]
    #    powers = [pow.list() for pow in powers]
    #    assert all([c.is_zero() for c in pow[total_degree:]])
    #    powers = [pow[:total_degree] for pow in powers]
    #    powers = [pow + [powers[0][0].parent().zero()]*(total_degree-len(pow)) for pow in powers]
    #    from sage.misc.flatten import flatten
    #    powers = [flatten([c.vector(base=base) for c in pow]) for pow in powers]
    #    from sage.matrix.constructor import matrix
    #    A = matrix(powers[:-1])
    #    b = -A.matrix_space().row_space()(powers[-1])
    #    x = A.solve_left(b)
    #    assert x*A == b, x*A-b

    #    x = x.list()
    #    assert len(x) == total_degree, total_degree
    #    R = primitive_element.parent().modulus().parent().change_ring(base)
    #    minpoly = R(x+[1])
    #    assert sum([c*primitive_element**i for i,c in enumerate(minpoly.coeffs())]).is_zero()
    #    return minpoly

    def implementation_ring(self):
        """
        Return the ring implementing the arithmetic in this ring.

        """
        return self._implementation_ring()[0]

    def to_implementation_ring(self, element):
        return element._element

    def from_implementation_ring(self, element):
        return self(element)

#    def _to_implementation_ring(self, element):
#        """
#        Represent ``element`` in the implementing p-adic ring.
#
#        INPUT:
#
#            - ``element`` -- a polynomial with coefficients in the base ring
#
#        OUTPUT:
#
#            An element in the ``abs_ring`` implementing this extension.
#
#        EXAMPLES::
#
#            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
#            sage: R = ZpFM(3,10)
#            sage: Rx.<x> = R[]
#            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
#            sage: M._to_implementation_ring(x)
#            1 + O(3^10)
#
#        """
#        if not element.parent() is self._polynomial_ring:
#            raise TypeError("%s must be a polynomial in %s"%(element,self._polynomial_ring))
#        return element.map_coefficients(self._implementation_ring()[1][0], self._implementation_ring()[0])(self._implementation_ring()[1][1])
#
#    def _from_implementation_ring(self, element):
#        """
#        Take an element from the implementing p-adic ring and bring it back
#        into this ring.
#
#        INPUT:
#
#            - ``element`` -- an element in the implementing p-adic ring
#
#        OUTPUT:
#
#            A polynomial with coefficients in the base ring.
#
#        EXAMPLES::
#
#            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
#            sage: R = ZpFM(3,10)
#            sage: Rx.<x> = R[]
#            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
#            sage: M._from_implementation_ring(R.one())
#            1 + O(3^10)
#
#        """
#        from sage.rings.padics.two_step_extension_generic import TwoStepExtensionGeneric
#        if self._implementation_ring is self.base_ring():
#            return self([element])
#        elif isinstance(self._implementation_ring,TwoStepExtensionGeneric):
#            valpart = (self.__from_implementation_ring[1]**element.valuation()).change_variable_name(self.variable_name())
#            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
#            R=PolynomialRing(self._polynomial_ring,'S')
#            ret = R([self._polynomial_ring(list(c.matrix().row(0)))(self.__from_implementation_ring[0]) for c in element.list(self._implementation_ring.inertia_subring())])(self.__from_implementation_ring[1].change_variable_name(self.variable_name()))
#            ret=ret.parent().quo(self.defining_polynomial())(ret)
#            ret=ret.list()[:1+ret.lift().degree()]
#            ret=self._polynomial_ring(ret)
#            return self(ret)
#        else:
#            raise NotImplementedError("Not implemented from %s to %s"%(self._implementation_ring,self))

    def _repr_(self):
        """
        Returns a print representation of this extension.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            General extension of 3-adic Ring of fixed modulus 3^10 in x defined by (1 + O(3^10))*x + (2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10))

        """
        type = "General"
        if self.is_trivial():
            type = "Trivial"
        elif self.is_eisenstein():
            type = "Eisenstein"
        elif self.is_unramified():
            type = "Unramified"
        elif self.is_totally_ramified():
            type = "Totally ramified"
        return "%s extension of %s in %s defined by %s"%(type, self.ground_ring(), self.variable_name(), self.defining_polynomial())

    def is_trivial(self):
        return self._poly.degree() == 1

    def is_eisenstein(self):
        from sage.rings.padics.factory import GenericExtensionFactory
        return GenericExtensionFactory.is_eisenstein(self._poly)

    def is_unramified(self):
        from sage.rings.padics.factory import GenericExtensionFactory
        return GenericExtensionFactory.is_unramified(self._poly)

    def is_totally_ramified(self):
        F = self._poly.map_coefficients(lambda c:c.residue(), self.base_ring().residue_class_field())
        F = F.factor()
        assert len(F)==1
        return F[0][0].degree()==1

    def ramification_index(self):
        """
        Return the absolute ramification index of ``self``, i.e., the valuation
        of ``p`` in ``self``.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            sage: M.ramification_index()
            1

        """
        return self.implementation_ring().ramification_index()

    def degree(self):
        """
        Returns the degree of this ring over its base ring.

        OUTPUT:

        A positive integer, the degree of the defining polynomial.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            sage: M.degree()
            1

        """
        return self.defining_polynomial().degree()

    def gen(self, n=0):
        """
        Returns a generator for self as an extension of its ground ring.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            sage: M.gen()
            1 + O(3^10)

        """
        if n == 0:
            return self([self.base_ring().zero(), self.base_ring().one()])
        raise IndexError, "only one generator"

    def gens(self):
        """
        Returns a list of generators of self.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            sage: M.gens()
            [1 + O(3^10)]

        """
        return [self.gen(0)]

    def inertia_subring(self):
        """
        Returns the inertia subring of ``self``.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            sage: M.inertia_subring()
            3-adic Ring of fixed modulus 3^10

        """
        return self._implementation_ring.inertia_subring()

    def residue_class_field(self):
        """
        Returns the residue class field.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            sage: M.residue_class_field()
            Finite Field of size 3

        """
        return self._implementation_ring()[0].residue_class_field()

    def prime(self):
        """
        Returns the characteristic of the residue field.

        EXAMPLES::

            sage: from sage.rings.padics.padic_extension_leaves import GeneralExtensionRingFixedMod
            sage: R = ZpFM(3,10)
            sage: Rx.<x> = R[]
            sage: M = GeneralExtensionRingFixedMod(x - 1, R, R(1), Hom(R,R).identity(), (), 10, None, {}, ('x','u0','u','3'))
            sage: M.prime()
            3

        """
        return self._prime

    def uniformizer(self):
        ret = self(None)
        ret._element = self._implementation_ring.uniformizer()
        return ret

    def _uniformizer_print(self):
        return "not_used"
