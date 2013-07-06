"""
Valuations on polynomial rings based on `\phi`-adic expansions

This file implements a base class for discrete valuations on polynomial rings,
defined by a `\phi`-adic expansion.

AUTHORS:

- Julian Rueth (15-04-2013): initial version

REFERENCES:

.. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

.. [ML1936'] MacLane, S. (1936). A construction for absolute values in
polynomial rings. Transactions of the American Mathematical Society, 40(3),
363-395.

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from discrete_valuation import DiscreteValuation
from sage.misc.abstract_method import abstract_method

def _lift_to_maximal_precision(c):
    """
    Lift ``c`` to maximal precision if the parent is not exact.

    EXAMPLES::

        sage: from sage.rings.padics.developing_valuation import _lift_to_maximal_precision
        sage: R = Zp(2,5)
        sage: x = R(1,2); x
        1 + O(2^2)
        sage: _lift_to_maximal_precision(x)
        1 + O(2^5)

        sage: x = 1
        sage: _lift_to_maximal_precision(x)
        1

    """
    return c if c.parent().is_exact() else c.lift_to_precision()

class DevelopingValuation(DiscreteValuation):
    """
    An abstract base class for a discrete valuation of polynomials defined over
    the polynomial ring ``domain`` by the `\phi`-adic development.

    INPUT:

    - ``domain`` -- a polynomial ring

    - ``phi`` -- a monic element of ``domain``

    EXAMPLES::

        sage: from sage.rings.padics.developing_valuation import DevelopingValuation
        sage: R = Zp(2,5)
        sage: S.<x> = R[]
        sage: DevelopingValuation(S, x)
        `(1 + O(2^5))*x`-adic valuation of Univariate Polynomial Ring in x over 2-adic Ring with capped relative precision 5

    """
    def __init__(self, domain, phi):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.developing_valuation import DevelopingValuation
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = DevelopingValuation(S, x)
            sage: type(v)
            <class 'sage.rings.padics.developing_valuation.DevelopingValuation'>

        """
        if phi.parent() is not domain:
            raise ValueError("phi must be in the domain of the valuation")
        if phi.is_constant():
            raise ValueError("phi must not be constant")
        if not phi.leading_coefficient().is_one():
            raise ValueError("phi must be monic")

        DiscreteValuation.__init__(self, domain)

        self._phi = phi

    def phi(self):
        """
        Return the polynomial `\phi`, the key polynomial of this valuation.

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.phi()
            (1 + O(2^5))*x

        """
        return self._phi

    def effective_degree(self, f):
        """
        Return the effective degree of ``f`` with respect to this valuation.

        The effective degree of `f` is the largest `i` such that the valuation
        of `f` and the valuation of `f_i\phi^i` in the development `f=\sum_j
        f_j\phi^j` coincide.

        INPUT:

        - ``f`` -- a non-zero polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.effective_degree(x)
            1
            sage: v.effective_degree(2*x + 1)
            0

        """
        # defined on p.497 of [ML1936']
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_zero():
            raise ValueError("the effective degree is only defined for non-zero polynomials")

        v = self(f)
        return [i for i,w in enumerate(self.valuations(f)) if w == v][-1]

    def is_equivalence_unit(self, f):
        """
        Return whether ``f`` is an equivalence unit, i.e., an element of
        :meth:`effective_degree` zero.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_unit(x)
            False
            sage: v.is_equivalence_unit(S.zero())
            False
            sage: v.is_equivalence_unit(2*x + 1)
            True

        """
        # defined on p.497 of [ML1936']
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")

        if f.is_zero():
            return False
        return self.effective_degree(f) == 0

    def equivalence_reciprocal(self, f):
        """
        Return an equivalence reciprocal of ``f``.

        An equivalence reciprocal of `f` is a polynomial `h` such that `f\cdot
        h` is equivalent to 1 modulo this valuation.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation which is an
          equivalence unit

        EXAMPLES::

            sage: R = Zp(3,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = 3*x + 2
            sage: h = v.equivalence_reciprocal(f); h
            (2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5))
            sage: v.is_equivalent(f*h, 1)
            True

        In an extended valuation over an extension field::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.extension(x^2 + x + u, 1)
            sage: f = 2*x + u
            sage: h = v.equivalence_reciprocal(f); h
            (u*2 + (u + 1)*2^2 + 2^3 + u*2^4 + 2^5 + O(2^6))*x + (u + 1) + (u + 1)*2 + 2^2 + u*2^3 + 2^4 + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        Extending the valuation once more::

            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: h = v.equivalence_reciprocal(f); h
            ((u + 1)*2^3 + 2^5 + 2^6 + (u + 1)*2^7 + O(2^8))*x^3 + (2^2 + u*2^4 + (u + 1)*2^5 + (u + 1)*2^6 + O(2^7))*x^2 + (u*2 + u*2^2 + u*2^3 + u*2^4 + 2^5 + O(2^6))*x + (u + 1) + (u + 1)*2 + (u + 1)*2^2 + (u + 1)*2^3 + u*2^4 + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        .. SEEALSO::

            :meth:`is_equivalence_unit`

        """
        # defined on p.497 of [ML1936']
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if not self.is_equivalence_unit(f):
            raise ValueError("f must be an equivalence unit")

        e0 = self.coefficients(f)[0]
        one,g,h = self.phi().xgcd(e0)
        assert one.is_one()

        assert self(f*h) == 0
        assert self(f*h - 1) > 0

        return h

    def extension(self, phi, mu, check=True):
        """
        Return the inductive valuation which extends this valuation by mapping
        ``phi`` to ``mu``.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation; this must be
          a key polynomial, see :meth:`is_key` for properties of key
          polynomials.

        - ``mu`` -- a rational number, the valuation of ``phi`` in the extended
          valuation

        - ``check`` -- whether or not to check the correctness of the
          parameters

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.extension(x^2 + x + u, 1)
            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: v
            [ Gauss valuation induced by 2-adic valuation, v((1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)) = 1, v((1 + O(2^5))*x^4 + (2^2 + O(2^6))*x^3 + (1 + (u + 1)*2 + O(2^5))*x^2 + ((u + 1)*2^2 + O(2^6))*x + (u + 1) + (u + 1)*2 + (u + 1)*2^2 + (u + 1)*2^3 + (u + 1)*2^4 + O(2^5)) = 3 ]
            sage: v.residue_field()
            Univariate Quotient Polynomial Ring in u2 over Univariate Quotient Polynomial Ring in u1 over Finite Field in u0 of size 2^2 with modulus u1^2 + u1 + u0 with modulus u2^2 + u1*u2 + u1

        .. SEEALSO::

            :class:`AugmentedValuation`

        """
        from augmented_valuation import AugmentedValuation
        return AugmentedValuation(self, phi, mu, check)

    def is_key(self, phi, explain=False):
        """
        Return whether ``phi`` is a key polynomial for this valuation.

        A key polynomial must satisfy the following conditions:

        - it must be monic
        - it must be equivalence-irreducible (see :meth:`is_equivalence_irreducible`)
        - it must be minimal (see :meth:`is_minimal`)

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation

        - ``explain`` -- a boolean (default: ``False``), if ``True``, return a
          string explaining why ``phi`` is not a key polynomial

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_key(x)
            True
            sage: v.is_key(2*x, explain = True)
            (False, 'phi must be monic')
            sage: v.is_key(x^2, explain = True)
            (False, 'phi must be equivalence irreducible')

            sage: w = v.extension(x, 1)
            sage: w.is_key(x + 1, explain = True)
            (False, 'phi must be minimal')

        """
        if phi.parent() is not self.domain():
            raise ValueError("phi must be in the domain of the valuation")

        reason = None

        if not phi.is_monic():
            reason = "phi must be monic"
        elif not self.is_minimal(phi):
            reason = "phi must be minimal"
        elif not self.is_equivalence_irreducible(phi):
            reason = "phi must be equivalence irreducible"

        if explain:
            return reason is None, reason
        else:
            return reason is None

    @abstract_method
    def is_commensurable_inductive(self):
        """
        Return whether this valuation is a commensurable inductive valuation
        over the discrete valuation of the base ring of the polynomial ring, as
        defined in section 4 of [ML1936].

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_commensurable_inductive()
            True
            sage: w = v.extension(x, 1)
            sage: w.is_commensurable_inductive()
            True

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

        """
        pass

    def is_minimal(self, f):
        """
        Return whether the polynomial ``f`` is minimal with respect to this
        valuation, as defined in definition 4.1 of [ML1936].

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_minimal(x + 1)
            True
            sage: w = v.extension(x, 1)
            sage: w.is_minimal(x + 1)
            False

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.
        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_constant():
            raise ValueError("f must not be constant")

        if self.is_commensurable_inductive():
            # use the characterization of theorem 9.4 in [ML1936]
            if not f.is_monic():
                raise NotImplementedError("is_minimal() only implemented for monic polynomials")
            return self.valuations(f)[-1] == self(f)

        raise NotImplementedError("is_minimal() only implemented for commensurable inductive values")

    def is_equivalence_irreducible(self, f):
        """
        Return whether the polynomial ``f`` is equivalence irreducible, i.e.,
        whether its :meth:`equivalence_decomposition` is irreducible.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_irreducible(x)
            True
            sage: v.is_equivalence_irreducible(x^2)
            False
            sage: v.is_equivalence_irreducible(x^2 + 2)
            False

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_constant():
            raise ValueError("f must not be constant")

        if self.is_commensurable_inductive():
            # use the characterization of Theorem 13.1 in [ML1936]
            if not f.is_monic():
                raise NotImplementedError("is_equivalence_irreducible() only implemented for monic polynomials")
            F = self.equivalence_decomposition(f, lift_to_keys=True) # OUCH: do something smarter here - there has to be a way to make this work with lift_to_keys=False
            return len(F)==1 and F[0][1]==1 and F.unit().degree()==0

        raise NotImplementedError("is_equivalence_irreducible() only implemented for inductive values")

    def equivalence_decomposition(self, f, lift_to_keys=True):
        """
        Return an equivalence decomposition of ``f``, i.e., a polynomial
        `g(x)=e(x)\prod_i \phi_i(x)` with `e(x)` an equivalence unit (see
        :meth:`is_equivalence_unit()`) and the `\phi_i` key polynomials (see
        :meth:`is_key`) such that ``f`` :meth:`is_equivalent` to `g`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        - ``lift_to_keys`` -- a boolean (default: ``True``), if ``False``, only
          return the factorization of `Rf` in the residue field (see
          description of the algorithm below), and do not lift it to an
          equivalence decomposition

        ALGORITHM:

        We use the algorithm described in Theorem 4.4 of [ML1936']. After
        removing all factors `\phi` from a polynomial `f`, there is an
        equivalence unit `R` such that `Rf` has valuation zero. Now `Rf` can be
        factored as `\prod_i \alpha_i` over the :meth:`residue_field`. Lifting
        all `\alpha_i` to key polynomials `\phi_i` gives `Rf=\prod_i R_i f_i`
        for suitable equivalence units `R_i` (see :meth:`lift_to_key`). Taking
        `R'` an :meth:`equivalence_reciprocal` of `R`, we have `f` equivalent
        to `(R'\prod_i R_i)\prod_i\phi_i`.

        EXAMPLES::

            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_decomposition(S.zero())
            Traceback (most recent call last):
            ...
            ValueError: equivalence decomposition of zero is not defined
            sage: v.equivalence_decomposition(S.one())
            1 + O(2^10)
            sage: v.equivalence_decomposition(x^2+2)
            ((1 + O(2^10))*x)^2
            sage: v.equivalence_decomposition(x^2+1)
            ((1 + O(2^10))*x + 1 + O(2^10))^2

        A polynomial that is an equivalence unit, is returned as the unit part
        of a :class:`sage.structure.factorization.Factorization`, leading to a unit
        non-minimal degree::

            sage: w = v.extension(x, 1)
            sage: F = w.equivalence_decomposition(x^2+1); F
            (1 + O(2^10))*x^2 + 1 + O(2^10)
            sage: F.unit()
            (1 + O(2^10))*x^2 + 1 + O(2^10)

        However, if the polynomial has a non-unit factor, then the unit might
        be replaced by a factor of lower degree::

            sage: f = x * (x^2 + 1)
            sage: F = w.equivalence_decomposition(f); F
            (1 + O(2^10))*x
            sage: F.unit()
            1 + O(2^10)

        Examples over an iterated unramified extension:

            sage: v = v.extension(x^2 + x + u, 1)
            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)

            sage: v.equivalence_decomposition(x)
            (1 + O(2^10))*x
            sage: F = v.equivalence_decomposition( v.phi() )
            sage: len(F)
            1
            sage: F = v.equivalence_decomposition( v.phi() * (x^4 + 4*x^3 + (7 + 2*u)*x^2 + (8 + 4*u)*x + 1023 + 3*u) )
            sage: len(F)
            2

        REFERENCES:

        .. [ML1936'] MacLane, S. (1936). A construction for absolute values in
        polynomial rings. Transactions of the American Mathematical Society, 40(3),
        363-395.

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_zero():
            raise ValueError("equivalence decomposition of zero is not defined")

        if not self.domain().base_ring().is_field():
            raise NotImplementedError("equivalence decomposition only implemented over fields")

        from sage.structure.factorization import Factorization
        if self.is_equivalence_unit(f):
            return Factorization([],unit=f)

        if not self.is_commensurable_inductive():
            raise NotImplementedError("only implemented for inductive valuations")

        f0 = f # used to check correctness of the output

        phi_divides = 0
        while self.valuations(f)[0] > self(f):
            f = f-self.coefficients(f)[0]
            assert self.phi().divides(f)
            f,_ = f.quo_rem(self.phi())
            phi_divides += 1

        R = self.equivalence_unit(-self(f))
        F = self.reduce(f*R)
        F = F.factor()
        unit = F.unit()
        F = list(F)

        # used to speed up is_equivalence_irreducible()
        if not lift_to_keys:
            for i,(g,e) in enumerate(F):
                if g == self.residue_ring().gen():
                    F[i] = (g, e+phi_divides)
                    break
            else:
                F.append((self.residue_ring().gen(), phi_divides))
            return Factorization(F)

        unit = self.lift( self.residue_ring()(unit) )

        from sage.misc.all import prod
        unit *= self.lift(self.residue_ring()(prod([ psi.leading_coefficient()**e for psi,e in F ])))
        F = [(self.lift_to_key(psi/psi.leading_coefficient()),e) for psi,e in F]

        unit *= self.equivalence_reciprocal(R) * prod([self.equivalence_unit(-self(g))**e for g,e in F])

        if phi_divides:
            for i,(g,e) in enumerate(F):
                if g == self.phi():
                    F[i] = (self.phi(),e+phi_divides)
                    break
            else:
                F.append((self.phi(),phi_divides))

        ret = Factorization(F, unit=unit)
        # assert self.is_equivalent(ret.prod(), f0) -- this might fail because of leading zeros
        assert self((ret.prod() - f0).map_coefficients(lambda c:_lift_to_maximal_precision(c)))
        assert self.is_equivalence_unit(ret.unit())
        return ret

    def minimal_representative(self, f):
        """
        Return a minimal representative for ``f``, i.e., a pair `e, a`
        such that ``f`` :meth:`is_equivalent`` to `e a`, `e` is
        an equivalence unit and `a` is minimal and monic.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        A factorization which has `e` as its unit and `a` as its unique factor.

        ALGORITHM:

        We use the algorithm described in the proof of Lemma 4.1 of [ML1936'].
        In the expansion `f=\sum_i f_i\phi^i` take `e=f_i` for the largest `i`
        with `f_i\phi^i` minimal (see :meth:`effective_degree`).
        Let `h` be the :meth:`equivalence_reciprocal` of `e` and take `a` given
        by the terms of minimal valuation in the expansion of `e f`.

        EXAMPLES::

            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x

            sage: v = v.extension(x, 1)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x + 2 + O(2^11)
            sage: f = x^3 + 6*x + 4
            sage: F = v.minimal_representative(f); F
            (2 + 2^2 + O(2^11)) * ((1 + O(2^10))*x + 2 + 2^2 + 2^4 + 2^6 + 2^8 + 2^10 + O(2^11))
            sage: v.is_minimal(F[0][0])
            True
            sage: v.is_equivalent(F[0][0], f)
            True

        REFERENCES:

        .. [ML1936'] MacLane, S. (1936). A construction for absolute values in
        polynomial rings. Transactions of the American Mathematical Society, 40(3),
        363-395.

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if f.is_zero():
            raise ValueError("the minimal representative of zero is not defined")

        if not self.is_commensurable_inductive():
            raise NotImplemented("only implemented for inductive valuations")

        f0 = f
        e = self.coefficients(f)[self.effective_degree(f)]
        f *= self.equivalence_reciprocal(e).map_coefficients(lambda c:_lift_to_maximal_precision(c))

        coeffs = [c if v == self(f) else c.parent().zero() for v,c in zip(self.valuations(f),self.coefficients(f))]
        coeffs[self.effective_degree(f0)] = self.domain().base_ring().one()
        ret = sum([c*self._phi**i for i,c in enumerate(coeffs)])
        assert self.effective_degree(ret) == self.effective_degree(f0)
        assert ret.is_monic(), coeffs
        assert self.is_minimal(ret)
        from sage.structure.factorization import Factorization
        ret = Factorization([(ret,1)],unit=e)
        # assert self.is_equivalent(ret.prod(), f0) -- this might fail because of leading zeros
        assert self((ret.prod() - f0).map_coefficients(lambda c:_lift_to_maximal_precision(c)))
        return ret

    def _normalize_leading_coefficients(self, f):
        """
        This method removes leading zero coefficients from ``f`` when
        appropriate.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        ``f`` with leading zero coefficients removed.

        .. NOTE::

            When ``f`` has leading zero coefficients, one could argue that we
            should never strip these but they often arise naturally, e.g., when
            when working with expressions as ``g-g`` or ``(g+c)-g``. We strip
            such coefficients if they are zero to sufficient precision. To be
            precise, if their precision exceeds the valuation of any other
            coefficient.
            It is not clear that this is the right way to do this.

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: g = x
            sage: v.coefficients(g-g) # indirect doctest
            []
            sage: (g-g).list()
            [0, O(2^5)]
            sage: f = x*R(0,1) + R(1,2); f
            (O(2))*x + (1 + O(2^2))
            sage: v.coefficients(f) # indirect doctest
            [(1 + O(2^2))]
            sage: f = x*R(0,1) + R(2,2); f
            (O(2))*x + (2 + O(2^2))
            sage: v.coefficients(f) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: f must not have leading zero coefficients

        """
        if len(f.list()) > f.degree()+1:
            # f has leading zero coefficients
            if f.is_zero():
                f= f.parent().zero()
            elif min([self.constant_valuation()(c) for c in f.list()[f.degree()+1:]]) > max([self.constant_valuation()(c) for c in f.list()[:f.degree()+1]]):
                f= self.domain()(f.list()[:f.degree()+1])
            else:
                raise ValueError("f must not have leading zero coefficients")

        return f

    def coefficients(self, f):
        """
        Return the `\phi`-adic expansion of ``f``.

        INPUT:

        - ``f`` -- a monic polynomial in the domain of this valuation

        OUTPUT:

        A list `[f_0,f_1,\dots]` of polynomials in the domain of this valuation
        such that `f=\sum_i f_i\phi^i`

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v.coefficients(f) # note that these constants are in the polynomial ring
            [(1 + 2 + O(2^5)), (2 + O(2^6)), (1 + O(2^5))]
            sage: v = v.extension( x^2 + x + 1, 1)
            sage: v.coefficients(f)
            [(1 + O(2^5))*x + (2 + O(2^5)), (1 + O(2^5))]

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        f = self._normalize_leading_coefficients(f)

        ret = []
        while f.degree() >= 0:
            f,r = f.quo_rem(self._phi)
            ret.append(r)
        return ret

    def newton_polygon(self, f):
        """
        Return the newton polygon the `\phi`-adic development of ``f``.

        INPUT::

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v.newton_polygon(f)
            Newton Polygon with vertices [(0, 0), (2, 0)]

            sage: v = v.extension( x^2 + x + 1, 1)
            sage: v.newton_polygon(f)
            Newton Polygon with vertices [(0, 0), (1, 1)]
            sage: v.newton_polygon( f * v.phi()^3 )
            Newton Polygon with vertices [(0, +Infinity), (3, 3), (4, 4)]

        .. SEEALSO::

            :class:`newton_polygon.NewtonPolygon`

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")

        from newton_polygon import NewtonPolygon
        return NewtonPolygon(self.valuations(f))

    def _call_(self, f):
        """
        Evaluate this valuation at ``f``.

        INPUT::

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3
            sage: v(f)
            0

            sage: v = v.extension( x^2 + x + 1, 1)
            sage: v(f)
            0
            sage: v(f * v.phi()^3 )
            3
            sage: v(S.zero())
            +Infinity

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation %s but is in %s"%(self.domain(),f.parent()))

        if f.is_zero():
            from sage.rings.all import infinity
            return infinity

        return min(self.valuations(f))

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: from sage.rings.padics.developing_valuation import DevelopingValuation
            sage: DevelopingValuation(S, x)
            `(1 + O(2^5))*x`-adic valuation of Univariate Polynomial Ring in x over 2-adic Field with capped relative precision 5

        """
        return "`%s`-adic valuation of %s"%(self._phi, self.domain())

    def residue_ring(self):
        """
        Return the residue ring of this valuation, i.e., a polynomial ring over
        the :meth:`residue_field`

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.residue_ring()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)

        """
        return self.domain().change_ring(self.residue_field())

    def mac_lane_step(self, G, assume_squarefree=False):
        assert not G.is_constant()
        R = G.parent()
        if R is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity

        if self(G) is infinity:
            raise ValueError("G must not have valuation infinity")

        if self.is_key(G):
            return [self.extension(G, infinity)]

        F = self.equivalence_decomposition(G)

        ret = []
        for phi,e in F:
            if G == phi: continue
            w = self.extension(phi, self(phi), check=False)
            NP = w.newton_polygon(G).principal_part()
            assert len(NP)
            for i in range(len(NP.slopes())):
                slope = NP.slopes()[i]
                side = NP.sides()[i]
                new_mu = self(phi) - slope
                base = self
                if phi.degree() == base.phi().degree():
                    assert new_mu > self(phi)
                    from gauss_valuation import GaussValuation
                    if not isinstance(base, GaussValuation):
                        base = base._base_valuation

                new_leaf = base.extension(phi, new_mu)
                assert slope is -infinity or 0 in new_leaf.newton_polygon(G).slopes()
                ret.append(new_leaf)

        assert ret
        return ret
