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

        EXAMPLES::

            sage: K = Qp(3)
            sage: R.<x> = K[]
            sage: L.<a> = K.extension(x + 1)
            sage: L._implementation_ring()

        """
        return GeneralExtensionGeneric.__implementation_ring(self.base_ring(), self._polynomial_ring, self._poly)

    @staticmethod
    def __implementation_ring(base, polynomial_ring, modulus):
        """
        Return a ring isomorphic to the extension of ``modulus`` over ``base``.

        INPUT:

        - ``base`` -- the base ring of ``modulus``

        - ``polynomial_ring`` -- the ring used to represent element in the
          extension

        - ``modulus`` -- a polynomial over ``base``

        OUTPUT:

        A triple as in :meth:`_implementation_ring`.

        EXAMPLES::

            sage: K = Qp(3)
            sage: R.<x> = K[]
            sage: L.<x> = K.extension(x + 1)
            sage: L.implementation_ring() # indirect doctest
            3-adic Field with capped relative precision 20

        """
        GeneralExtensionGeneric.__implementation_ring_check_args(base, polynomial_ring, modulus)

        from sage.rings.padics.factory import GenericExtensionFactory
        if GenericExtensionFactory.is_trivial(modulus):
            return GeneralExtensionGeneric.__implementation_ring_trivial(base, polynomial_ring, modulus)
        elif not base.is_eisenstein_over_unramified():
            # reduce to the case where base is an Eisenstein over unramified extension
            base2, base2_to_base, base_to_base2 = base.make_eisenstein_over_unramified()
            modulus_over2 = modulus.map_coefficients(base_to_base2, base2)
            polynomial_ring2 = polynomial_ring.change_ring(base2)

            implementation_ring, from_implementation_ring, to_implementation_ring = GeneralExtensionGeneric.__implementation_ring(base2, polynomial_ring2, modulus_over2)

            from_implementation_ring = lambda f:from_implementation_ring(f).map_coefficients(base2_to_base, base)
            to_implementation_ring = lambda f: to_implementation_ring(f.map_coefficients(base_to_base2, base2))
            return implementation_ring, from_implementation_ring, to_implementation_ring
        else:
            if GenericExtensionFactory.is_eisenstein(modulus):
                return GeneralExtensionGeneric.__implementation_ring_eisenstein(base, polynomial_ring, modulus)
            elif GenericExtensionFactory.is_unramified(modulus):
                return GeneralExtensionGeneric.__implementation_ring_unramified(base, polynomial_ring, modulus)
            elif GenericExtensionFactory.is_totally_ramified(modulus):
                return GeneralExtensionGeneric.__implementation_ring_totally_ramified(base, polynomial_ring, modulus)
            else:
                return GeneralExtensionGeneric.__implementation_ring_general(base, polynomial_ring, modulus)

    @staticmethod
    def __implementation_ring_check_args(base, polynomial_ring, modulus):
        """
        Helper method for :meth:`__implementation_ring` to check the validity
        of the parameters.

        Input is the same as in :meth:`__implementation_ring`.

        EXAMPLES:

            sage: K = Qp(3)
            sage: R.<x> = K[]
            sage: L.<a> = K.extension(x + 1/3) # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: modulus must be integral

        """
        if any([c.valuation()<0 for c in modulus.list()]):
            raise NotImplementedError("modulus must be integral")
        if not modulus.is_monic():
            raise NotImplementedError("modulus must be monic")
        if modulus.base_ring() is not base:
            raise ValueError("modulus must be defined over base")
        if polynomial_ring.base_ring() is not base:
            raise ValueError("polynomial_ring must be over base")

    @staticmethod
    def __implementation_ring_trivial(base, polynomial_ring, modulus):
        """
        Helper method for :meth:`__implementation_ring`` if ``modulus`` defines
        a trivial extension.

        Input and output are the same as in :meth:`__implementation_ring`.

        EXAMPLES::

            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a + 1)
            sage: L.implementation_ring() # indirect doctest
            3-adic Field with capped relative precision 20

        """
        GeneralExtensionGeneric.__implementation_ring_check_args(base, polynomial_ring, modulus)

        if modulus.degree() != 1 or not modulus.is_monic():
            raise ValueError("modulus must be a monic polynomial of degree 1")

        implementation_ring = base
        from_implementation_ring = polynomial_ring
        to_implementation_ring = lambda f:f(-modulus[0])
        return base, from_implementation_ring, to_implementation_ring

    @staticmethod
    def __implementation_ring_totally_ramified(base, polynomial_ring, modulus):
        """
        Helper method for :meth:`__implementation_ring`` if ``modulus`` defines
        a totally ramified extension but ``modulus`` is not Eisenstein.

        Input and output are the same as in :meth:`__implementation_ring`.

        ALGORITHM:

        TODO: why is this guaranteed to work?
        TODO: is this in some sense optimal?

        Uses :meth:`pAdicValuation.is_totally_ramified` to compute MacLane
        approximants. A product of the key polynomials of these defines a
        uniformizer.

        Empirically, this product gives a better minimal polynomial if it uses
        the minimal number of factors possible.  Currently, these factors are
        determined using brute force; there is probably a better way.

        EXAMPLES::

            sage: K = Qp(3,5)
            sage: R.<x> = K[]
            sage: L.<a> = K.extension(x^2 + 3*x + 9) # indirect doctest
            sage: L.implementation_ring()
            Eisenstein Extension of 3-adic Field with capped relative precision 5 in x defined by (1 + O(3^3))*x^2 + (2*3 + 2*3^2 + 2*3^3 + O(3^4))*x + (3 + O(3^3))

        """
        GeneralExtensionGeneric.__implementation_ring_check_args(base, polynomial_ring, modulus)

        from factory import GenericExtensionFactory
        if GenericExtensionFactory.is_eisenstein(modulus):
            raise ValueError("modulus must not be Eisenstein")
        if not base.is_field():
            raise NotImplementedError("currently only implemented over fields")

        from padic_valuation import pAdicValuation
        v = pAdicValuation(base)
        is_totally_ramified, ramification_steps = v.is_totally_ramified(modulus, include_steps=True, assume_squarefree=True)
        assert is_totally_ramified, "MacLane approximants terminated at: %s"%ramification_steps
        #TODO: why is this true
        assert any([v(v.phi()).denominator() == modulus.degree() for v in ramification_steps])

        # find a product of keys which gives a uniformizer
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

        # a uniformizer and its minpoly
        uniformizer = modulus.parent().quo(modulus).one()
        for key in bfs[1]:
            uniformizer*=key
        epoly = uniformizer.minpoly()

        assert epoly.degree() == modulus.degree(), epoly
        assert epoly[0].valuation() == 1, "%s is not Eisenstein"%epoly

        # now we reduced to the case of an Eisenstein extension; solve that
        # case recursively and combine the pieces
        implementation_ring = base.extension(epoly,names=epoly.parent().variable_name())

        # to get from the Eisenstein extension back to this ring, we note that
        # the uniformizer from above is the uniformizer in the Eisenstein
        # extension (can we get more precision here?)
        image_of_uniformizer = uniformizer.lift()(polynomial_ring.gen())
        from_implementation_ring = lambda f: f.polynomial()(image_of_uniformizer)

        # to go the other way, we send the generator of this ring to any root
        # in the implementation ring (is there a better way?)
        image_of_gen = modulus.change_ring(implementation_ring).any_root()
        to_implementation_ring = lambda f:f(image_of_gen)

        return implementation_ring, from_implementation_ring, to_implementation_ring

    @staticmethod
    def __implementation_ring_eisenstein(base, polynomial_ring, modulus):
        """
        Helper method for :meth:`__implementation_ring`` if ``modulus`` is
        Eisenstein.

        Input and output are the same as in :meth:`__implementation_ring`.

        ALGORITHM:

        ``base`` is an Eisenstein extension of an unramified extension;
        computing the resultant of the Eisenstein polynomial and modulus gives
        a single Eisenstein polynomial which defines the total extension over
        the unramified extension.

        EXAMPLES::

            sage: K = Qp(3)
            sage: R.<x> = K[]
            sage: L.<x> = K.extension(x + 3)
            sage: R.<a> = L[]
            sage: M.<a> = L.extension(a^2 + x)
            sage: M.implementation_ring() # indirect doctest

        """
        GeneralExtensionGeneric.__implementation_ring_check_args(base, polynomial_ring, modulus)

        if not GenericExtensionFactory.is_eisenstein(modulus):
            raise ValueError("modulus must be Eisenstein")

        unramified_base = base.maximal_unramified_subextension()
        raise NotImplementedError("create bivariate ring")
        raise NotImplementedError("bring both Eisenstein polynomials into that ring")
        epoly = f.resultant(g, var)

        implementation_ring = unramified_base.extension(epoly, names=epoly.variable_name())

        # to get back to the Eisenstein over Eisenstein over unramified
        # extension, we ignore one of the Eisenstein extensions and simply not
        # reduce modulo the middle extension - this is probably not what we
        # want in the long run
        from_implementation_ring = lambda f: polynomial_ring(f.polynomial().list())

        # we find the image of the uniformizer of base in implementation_ring
        # by factoring its minpoly
        base_uniformizer = base.eisenstein_polynomial().change_ring(implementation_ring).any_root()
        to_implementation_ring = lambda f: f.map_coefficients(lambda c:c.polynomial(unramified_base.gen(), base_uniformizer), implementation_ring)(implementation_ring.uniformizer())

        return implementation_ring, from_implementation_ring, to_implementation_ring

    def __implementation_ring_unramified(base, polynomial_ring, modulus):
        """
        Helper method for :meth:`__implementation_ring`` if ``modulus`` is
        unramified.

        Input and output are the same as in :meth:`__implementation_ring`.

        ALGORITHM:

        ``base`` is an Eisenstein extension of an unramified extension; since
        unramified extensions are uniquely determined by their degree, we
        create a new unramified extension of the appropriate degree, and put
        the original Eisenstein extension on top of that extension.

        EXAMPLES::

            sage: K = Qp(2)
            sage: R.<u> = K[]
            sage: L.<u> = K.extension(u^2 + u + 1)
            sage: R.<a> = L[]
            sage: M.<a> = L.extension(a^2 + 2*u + 2)
            sage: R.<v> = M[]
            sage: N.<v> = M.extension(TODO unramified poly)
            sage: N.implementation_ring() # indirect doctest

        """
        GeneralExtensionGeneric.__implementation_ring_check_args(base, polynomial_ring, modulus)

        if not GenericExtensionFactory.is_unramified(modulus):
            raise ValueError("modulus must be unramified")

        unramified_degree = base.maximal_unramified_subextension().degree() * modulus.degree()
        from sage.rings.finite_rings.constructor import FiniteField as GF
        unramified_modulus = GF(base.prime()**unramified_degree).modulus().change_ring(ZZ).change_ring(base.ground_ring_of_tower())
        unramified_base = base.ground_ring_of_tower().extension(unramified_modulus, names=modulus.parent().variable_name())

        old_unramified_gen = base.maximal_unramified_subextension().modulus().change_ring(unramified_base).any_root()
        old_eisenstein = base.eisenstein_polynomial().map_coefficients(lambda c:c.polynomial()(old_unramified_gen), unramified_base)

        implementation_ring = unramified_base.extension(old_eisenstein, names=old_eisenstein.variable_name())

        base_to_implementation_ring = lambda f: f.bivariate_polynomial(old_unramified_gen, implementation_ring.uniformizer())
        gen_in_implementation_ring = modulus.map_coefficients(base_to_implementation_ring, implementation_ring).any_root()
        to_implementation_ring = lambda f: f.map_coefficients(base_to_implementation_ring, implementation_ring)(gen_in_implementation_ring)

        #from_implementation_ring = lambda f: f.bivariate_polynomial(???,base.uniformizer())

        return implementation_ring, from_implementation_ring, to_implementation_ring

    @staticmethod
    def __implementation_ring_general(base, polynomial_ring, modulus):
        """
        Helper method for :meth:`__implementation_ring`` if ``modulus`` has an
        unramified and a totally ramfieid part.

        Input and output are the same as in :meth:`__implementation_ring`.

        ALGORITHM:

        we extend ``base`` with the unramified part of ``modulus``; the factors
        of ``modulus`` then define totally ramified extensions, handled by
        :meth:`__implementation_ring_totally_ramified`.

        EXAMPLES::

            sage: K = Qp(2)
            sage: R.<u> = K[]
            sage: L.<u> = K.extension(u^2 + u + 1)
            sage: R.<a> = L[]
            sage: M.<a> = L.extension(a^2 + 2*u + 2)
            sage: R.<v> = M[]
            sage: N.<v> = M.extension(general poly)
            sage: N.implementation_ring() # indirect doctest

        """
        raise NotImplementedError

    def implementation_ring(self):
        """
        Return the ring implementing the arithmetic in this ring.

        """
        return self._implementation_ring()[0]

    def to_implementation_ring(self, element):
        """
        Return ``element`` as an element of :meth:`implementation_ring`.

        INPUT:

        - ``element`` -- an element of this ring

        """
        return element._element

    def from_implementation_ring(self, element):
        """
        Return ``element`` as an element of this ring.

        INPUT:

        - ``element`` -- an elemenf of :meth:`implementation_ring`

        """
        return self(element)

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

    def is_eisenstein_over_unramified(self):
        return False
