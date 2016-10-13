"""
The discrete valuation of a `p`-adic ring

This file makes it possible to use `p`-adics, integers, and rationals in the
general discrete valuation framework.

AUTHORS:

- Julian Rueth (2013-03-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the mac_lane/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.getcwd())
    sys.path.append(os.path.dirname(os.getcwd()))

from discrete_valuation import DiscreteValuation, DiscretePseudoValuation
from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method
from sage.misc.fast_methods import WithEqualityById

class PadicValuationFactory(UniqueFactory):
    """
    Create a ``prime``-adic valuation on ``R``.

    INPUT:

    - ``R`` -- a ring

    - ``prime`` -- a prime or ``None`` (default: ``None``), if ``None``, try
      to automatically determine the appropriate prime for ``R``

    EXAMPLES::

        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

    ``prime`` may be ``None`` for `p`-adic rings::

        sage: pAdicValuation(Qp(2))
        2-adic valuation
        sage: pAdicValuation(ZZ)
        Traceback (most recent call last):
        ...
        ValueError: prime must be specified for this ring

    TESTS:

    A case that took very long due to the hashing of number fields::

        sage: from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: Delta1= x^12 + 20*x^11 + 154*x^10 + 664*x^9 + 1873*x^8 + 3808*x^7 + 5980*x^6 + 7560*x^5 + 7799*x^4 + 6508*x^3 + 4290*x^2 + 2224*x + 887 
        sage: K.<theta>=NumberField(x^9-2)
        sage: vK=pAdicValuation(QQ,2).extension(K)
        sage: G=Delta1.change_ring(K)
        sage: v0=GaussValuation(G.parent(),vK)
        sage: V=v0.mac_lane_step(G)
        sage: while len(V) == 1: V=V[0].mac_lane_step(G)
        sage: v=V[0]
        sage: while v(v.phi()) < 50: v=v.mac_lane_step(G)[0]
        sage: F=v.phi()
        sage: L.<alpha>=K.extension(F)
        sage: vK.mac_lane_approximants(F)
        [[ Gauss valuation induced by theta-adic valuation, v(x + 1) = 9/4 ]]
        sage: vL = vK.extension(L)

    """
    def create_key_and_extra_args(self, R, prime=None, check=True):
        from sage.rings.all import ZZ, QQ
        from sage.rings.number_field.number_field import is_NumberField
        from sage.rings.padics.padic_generic import pAdicGeneric
        if R is ZZ or R is QQ:
            if prime is None:
                raise ValueError("prime must be specified for this ring")
            return (R, prime), {}
        elif isinstance(R, pAdicGeneric):
            if prime is None:
                prime = R.prime()
            if isinstance(prime, DiscretePseudoValuation):
                if prime(R.prime())>0:
                    prime = R.prime()
            if prime != R.prime():
                raise NotImplementedError("p-adic valuation not implemented for this ring")
            return (R, prime), {}
        elif is_NumberField(R):
            if prime is None:
                raise ValueError("prime must be specified for this ring")
            if not isinstance(prime, DiscretePseudoValuation):
                raise ValueError("prime must be a discrete (pseudo) valuation on a polynomial ring over the base of the number field")
            normalized_primes = prime.constant_valuation().mac_lane_approximants(R.relative_polynomial())
            normalized_prime = [v for v in normalized_primes if v(prime.phi()) == prime(prime.phi())]
            if len(normalized_prime) == 0:
                raise ValueError("prime does not define a valuation on this number field.")
            if len(normalized_prime) > 1:
                raise ValueError("prime does not define a unique valuation on this number field.")
            return (id(R), normalized_prime[0]), {"R":R}
        else:
            raise NotImplementedError("p-adic valuation not implemented for this ring")

    def create_object(self, version, key, **extra_args):
        from sage.rings.all import ZZ, QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        R, prime = key
        if R is ZZ or R is QQ:
            return pAdicValuation_int(R, prime)
        elif isinstance(R, pAdicGeneric):
            return pAdicValuation_padic(R)
        else:
            R = extra_args['R']
            return pAdicValuation_number_field(R, prime)

pAdicValuation = PadicValuationFactory("pAdicValuation")

class pAdicValuation_base(DiscreteValuation):
    """
    Common base class for `p`-adic valuations on integral domains.

    INPUT:

    - ``ring`` -- an integral domain whose elements provide a method
      ``valuation`` which takes ``prime`` as a parameter

    - ``prime`` -- a prime

    - ``check`` -- a boolean (default: ``True``) whether to perform checks on
      the parameters

    EXAMPLES::

        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: pAdicValuation(QQ, 5)
        5-adic valuation

     For `p`-adic rings, ``prime`` has to match the prime of the ring.

        sage: v = pAdicValuation(Zp(3), 2); v
        Traceback (most recent call last):
        ...
        NotImplementedError: p-adic valuation not implemented for this ring

    TESTS::

        sage: TestSuite(pAdicValuation(ZZ, 3)).run(skip="_test_category")
        sage: TestSuite(pAdicValuation(QQ, 5)).run(skip="_test_category")
        sage: TestSuite(pAdicValuation(Zp(5), 5)).run(skip="_test_category")

    """
    def __init__(self, ring, prime):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.padic_valuation import pAdicValuation_int # optional: integrated
            sage: isinstance(pAdicValuation(ZZ, 3), pAdicValuation_int)
            True

        """
        DiscreteValuation.__init__(self, ring)
        self._prime = prime

    def prime(self):
        """
        Return the prime `p` of this `p`-adic valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).prime()
            3

        """
        return self._prime

    def value_group(self):
        r"""
        Return the value group \ZZ of this valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).value_group()
            DiscreteValueGroup(1)

        """
        from discrete_value_group import DiscreteValueGroup
        return DiscreteValueGroup(1)

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        return "%s-adic valuation"%(self.prime())

    def _call_(self, x):
        """
        Evaluate this valuation at ``x``.

        INPUT::

        - ``x`` --  an element in the domain of this valuation

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3)(9)
            2

        """
        if x.parent() is not self.domain():
            raise ValueError("x must be in the domain of the valuation")

        return x.valuation(self.prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).reduce(4)
            1

        """
        if x.parent() is not self.domain():
            raise ValueError("x must be in the domain of the valuation")

        return self.residue_field()(x)

    def lift(self, x):
        """
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: xbar = v.reduce(4)
            sage: v.lift(xbar)
            1

        """
        if x.parent() is not self.residue_field():
            raise ValueError("x must be in the residue field of the valuation")

        return self.domain()(x)

    def uniformizer(self):
        """
        Return a uniformizer of this `p`-adic valuation, i.e., `p` as an
        element of the domain.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.uniformizer()
            3

        """
        return self.domain()(self._prime)

    def residue_field(self):
        """
        Return the residue field of the ring of integers of this valuation.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.residue_field()
            Finite Field of size 3

        """
        from sage.rings.all import GF
        return GF(self._prime)

    def shift(self, c, v):
        """
        Multiply ``c`` by a ``v``-th power of the uniformizer.

        INPUT:

        - ``c`` -- an element of the domain of this valuation

        - ``v`` -- an integer

        OUTPUT:

        If the resulting element has negative valation, then it will be brought
        to the fraction field, otherwise it is an element of the domain.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.shift(2,3)
            54
            sage: v.shift(2,-3)
            2/27

        """
        from sage.rings.all import ZZ
        if c.parent() is not self.domain():
            raise ValueError("c must be in the domain of the valuation")
        if not v in ZZ:
            raise ValueError("v must be an integer")
        v = ZZ(v)

        c = self.domain().fraction_field()(c)
        ret = c*self.uniformizer()**v

        if self(c) >= -v: ret = self.domain()(ret)
        return ret

    def is_unramified(self, G, include_steps=False, assume_squarefree=False):
        """
        Return whether ``G`` defines an unramified extension of the `p`-adic
        completion of the domain this valuation.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``); whether or not
          to include the approximate valuations that were used to determine
          the result in the return value.

        - ``assume_squarefree`` -- a boolean (default: ``False``); whether or
          not to assume that ``G`` is square-free over the `p`-adic completion
          of the domain of this valuation. Setting this to ``True`` can
          significantly improve the performance.

        EXAMPLES:

        We consider an extension as unramified if its ramification index is 1.
        Hence, a trivial extension is unramified::

            sage: pass


        NOTE: It is not enough to check whether ``G`` is irreducible in
        reduction for this.

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.E() > 1:
                ret = False
                break
            if v.F() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = True
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next)>1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def is_totally_ramified(self, G, include_steps=False, assume_squarefree=False):
        """
        Return whether ``G`` defines a totally ramified extension.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``), include the the
          valuations produced during the process of checking whether ``G`` is
          totally ramified in the return value

        - ``assume_squarefree`` -- a boolean (default: ``False``), whether to
          assume the input to be squarefree

        ALGORITHM:

        This is simplified version of :meth:`mac_lane_approximants`.

        EXAMPLES::

            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v.is_totally_ramified(G)
            False
            sage: G = x + 1
            sage: v.is_totally_ramified(G)
            True
            sage: G = x^2 + 2
            sage: v.is_totally_ramified(G)
            False
            sage: G = x^2 + 5
            sage: v.is_totally_ramified(G)
            True
            sage: v.is_totally_ramified(G, include_steps=True)
            (True, [Gauss valuation induced by 5-adic valuation, [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x) = 1/2 ]])

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.F() > 1:
                ret = False
                break
            if v.E() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = False
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next)>1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def montes_factorization(self, G, assume_squarefree=False):
        """
        Factor ``G`` by computing :meth:`mac_lane_approximants` to arbitrary
        precision.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``assume_squarefree`` -- a boolean (default: ``False``), whether to
          assume the input to be squarefree

        EXAMPLES::

            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v.montes_factorization(G) # long time
            ((1 + O(5^4))*x + 2 + 5 + 2*5^2 + 5^3 + O(5^4)) * ((1 + O(5^4))*x + 3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4))

        REFERENCES:

        .. [GMN2008] Jordi Guardia, Jesus Montes, Enric Nart (2008). Newton
        polygons of higher order in algebraic number theory. arXiv:0807.2620
        [math.NT]

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not G.is_monic():
            raise ValueError("G must be monic")
        if not all([self(c)>=0 for c in G.coefficients()]):
            raise ValueError("G must be integral")

        from sage.rings.all import infinity
        # W contains approximate factors of G
        W = self.mac_lane_approximants(G, precision_cap=infinity,assume_squarefree=assume_squarefree)
        ret = [w.phi() for w in W]

        # the polynomials in ret give "a" factorization; there is no guarantee
        # that every factorization is of this form... TODO
        from sage.structure.factorization import Factorization
        return Factorization([ (g,1) for g in ret ], simplify=False)

    def mac_lane_approximants(self, G, precision_cap=None, assume_squarefree=False):
        """
        Compute the extensions `w` of this valuation to a polynomial ring over
        the domain which have `w(G)=\infty`.

        INPUT:

        - ``G`` -- a monic polynomial over the domain of this valuation

        - ``precision_cap`` -- a rational, ``infinity`` or ``None`` (default:
          ``None``), stop computation as soon as no new branches can show up in
          the tree of valuations and the valuation of the key polynomials is at
          least ``precision_cap`` (if not ``None``).

        - ``assume_squarefree`` -- a boolean (default: ``False``), whether to
          assume the input to be squarefree

        EXAMPLES::

            sage: k=Qp(2,10)
            sage: v = pAdicValuation(k)

            sage: R.<x>=k[]
            sage: G = x
            sage: v.mac_lane_approximants(G)
            [Gauss valuation induced by 2-adic valuation]
            sage: v.mac_lane_approximants(G, precision_cap = infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = +Infinity ]]

            sage: G = x^2 + 1
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x + (1 + O(2^10))) = 1/2 ]]
            sage: v.mac_lane_approximants(G, precision_cap = infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x + (1 + O(2^10))) = 1/2, v((1 + O(2^10))*x^2 + (1 + O(2^10))) = +Infinity ]]

            sage: G = x^4 + 2*x^3 + 2*x^2 - 2*x + 2
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = 1/4 ]]
            sage: v.mac_lane_approximants(G,infinity)
            [[ Gauss valuation induced by 2-adic valuation, v((1 + O(2^10))*x) = 1/4, v((1 + O(2^10))*x^4 + (2 + O(2^11))*x^3 + (2 + O(2^11))*x^2 + (2 + 2^2 + 2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^8 + 2^9 + 2^10 + O(2^11))*x + (2 + O(2^11))) = +Infinity ]]

        The factorization of primes in the Gaussian integers can be read off
        the Mac Lane approximants::

            sage: v0 = pAdicValuation(QQ, 2)
            sage: R.<x> = QQ[]
            sage: G = x^2 + 1
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 1/2 ]]

            sage: v0 = pAdicValuation(QQ, 3)
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 3-adic valuation, v(x^2 + 1) = +Infinity ]]

            sage: v0 = pAdicValuation(QQ, 5)
            sage: v0.mac_lane_approximants(G)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 2) = 1 ], [ Gauss valuation induced by 5-adic valuation, v(x + 3) = 1 ]]
            sage: v0.mac_lane_approximants(G, precision_cap = 10) # long time
            [[ Gauss valuation induced by 5-adic valuation, v(x + 25670807) = 11 ], [ Gauss valuation induced by 5-adic valuation, v(x + 23157318) = 11 ]]

        The same example over the 5-adic numbers. In the quadratic extension
        `\QQ[x]/(x^2+1)`, 5 factors `-(x - 2)(x + 2)`, this behaviour can be
        read off the Mac Lane approximants::

            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v1,v2 = v.mac_lane_approximants(G); v1,v2
            ([ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + O(5^4))) = 1 ], [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + 3 + O(5^4)) = 1 ])
            sage: w1, w2 = v.mac_lane_approximants(G,precision_cap=2); w1,w2
            ([ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + 5 + 2*5^2 + O(5^4))) = 3 ], [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + 3 + 3*5 + 2*5^2 + O(5^4)) = 3 ])

        Note how the latter give a better approximation to the factors of `x^2 + 1`::

            sage: v1.phi() * v2.phi() - G
            (5 + O(5^4))*x + 5 + O(5^4)
            sage: w1.phi() * w2.phi() - G
            (5^3 + O(5^4))*x + 5^3 + O(5^4)

        In this example, the process stops with a factorization of `x^2 + 1`::

            sage: v.mac_lane_approximants(G, precision_cap=infinity)
            [[ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (2 + 5 + 2*5^2 + 5^3 + O(5^4))) = +Infinity ],
             [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x + (3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4))) = +Infinity ]]

        This obviously cannot happen over the rationals where we only get an
        approximate factorization::

            sage: v = pAdicValuation(QQ, 5)
            sage: R.<x>=QQ[]
            sage: G = x^2 + 1
            sage: v.mac_lane_approximants(G)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 2) = 1 ], [ Gauss valuation induced by 5-adic valuation, v(x + 3) = 1 ]]
            sage: v.mac_lane_approximants(G, precision_cap=5)
            [[ Gauss valuation induced by 5-adic valuation, v(x + 14557) = 6 ], [ Gauss valuation induced by 5-adic valuation, v(x + 32318) = 7 ]]

        TESTS:

        Initial versions ran into problems with the trivial residue field
        extensions in this case::

            sage: K = Qp(3,20)
            sage: R.<T> = K[]

            sage: alpha = T^3/4
            sage: G = 3^3*T^3*(alpha^4 - alpha)^2 - (4*alpha^3 - 1)^3
            sage: G = G/G.leading_coefficient()
            sage: pAdicValuation(K).mac_lane_approximants(G) # long time
            [[ Gauss valuation induced by 3-adic valuation, v((1 + O(3^20))*T + (2 + O(3^20))) = 1/9, v((1 + O(3^20))*T^9 + (2*3 + 2*3^2 + O(3^21))*T^8 + (3 + 3^5 + O(3^21))*T^7 + (2*3 + 2*3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^6 + O(3^21))*T^6 + (2*3 + 2*3^2 + 3^4 + 3^6 + 2*3^7 + O(3^21))*T^5 + (3 + 3^2 + 3^3 + 2*3^6 + 2*3^7 + 3^8 + O(3^21))*T^4 + (2*3 + 2*3^2 + 3^3 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + O(3^21))*T^3 + (2*3 + 2*3^2 + 3^3 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 2*3^8 + O(3^21))*T^2 + (3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^7 + 3^8 + O(3^21))*T + (2 + 2*3 + 2*3^2 + 2*3^4 + 2*3^5 + 3^7 + O(3^20))) = 55/27 ]]

        A similar example::

          sage: R.<x> = QQ[]
          sage: v = pAdicValuation(QQ, 3)
          sage: G = (x^3 + 3)^3 - 81
          sage: v.mac_lane_approximants(G)
          [[ Gauss valuation induced by 3-adic valuation, v(x) = 1/3, v(x^3 + 3*x + 3) = 13/9 ]]

        Another problematic case::

            sage: R.<x>=QQ[] 
            sage: Delta=x^12 + 20*x^11 + 154*x^10 + 664*x^9 + 1873*x^8 + 3808*x^7 + 5980*x^6 + 7560*x^5 + 7799*x^4 + 6508*x^3 + 4290*x^2 + 2224*x + 887 
            sage: K.<theta>=NumberField(x^6+108) 
            sage: K.is_galois()
            True
            sage: vK=pAdicValuation(QQ,2).extension(K)
            sage: vK(2) 
            3 
            sage: vK(theta) 
            1 
            sage: G=Delta.change_ring(K) 
            sage: V=vK.mac_lane_approximants(G) 

        OUTPUT:

        A list of valuations over the parent of ``G``.

        .. NOTE::

            For an irreducible ``G`` over `K`, these extensions correspond to
            the extensions of this valuation to the field `K[x]/(G(x))`, with
            `K` the number field corresponding to the domain of this valuation,
            see Theorem 2.1 in [ML1936]

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

        """
        R = G.parent()
        if R.base_ring() is not self.domain():
            raise ValueError("G must be defined over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")
        if not G.is_monic():
            raise ValueError("G must be monic")

        from sage.rings.all import infinity
        from gauss_valuation import GaussValuation

        leaves = [ GaussValuation(R, self)]
        while True:
            ef = [ v.E()*v.F() for v in leaves]
            if sum(ef) == G.degree():
                if precision_cap is None or all([v(v.phi())>precision_cap for v in leaves]):
                    return leaves

            expandables = []
            new_leaves = []
            for v in leaves:
                if v(G) is infinity:
                    new_leaves.append(v)
                else:
                    expandables.append(v)
            leaves = new_leaves

            if not expandables:
                return leaves

            for v in expandables:
                leaves.extend(v.mac_lane_step(G))

    def change_ring(self, base_ring):
        return pAdicValuation(base_ring)

    def extension(self, L, algorithm="mac_lane"):
        K = self.domain()
        if L is K:
            return self
        if L.base_ring() is not K:
            raise ValueError("L must be a simple finite extension of %s"%K)

        if algorithm == "ideal":
            I = L.ideal(self.prime()).factor()
            if len(I) > 1:
                raise ValueError("extension to %s is not unique"%L)
            I = I[0]
            return pAdicValuation(L, I, check=False)
        elif algorithm == "mac_lane":
            W = self.mac_lane_approximants(L.defining_polynomial())
            if len(W) > 1:
                raise ValueError("extension to %s is not unique"%L)
            return pAdicValuation(L, W[0], check=False)
        else: raise ValueError

class pAdicValuation_padic(pAdicValuation_base):
    """
    The `p`-adic valuation of a `p`-adic ring.

    INPUT:

    - ``R`` -- a `p`-adic ring

    EXAMPLES::

        sage: pAdicValuation(Qp(2)) #indirect doctest
        2-adic valuation

    """
    def __init__(self, R):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.padic_valuation import padicValuation_padic # optional: integrated
            sage: isinstance(pAdicValuation(Qp(2)), pAdicValuation_padic)
            True

        """
        pAdicValuation_base.__init__(self, R, R.prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element of the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: R = Zp(3)
            sage: pAdicValuation(Zp(3)).reduce(R(4))
            1

        """
        if x.parent() is not self.domain():
            raise ValueError("x must be in the domain of the valuation")

        return self.residue_field()(x.residue())

    def lift(self, x):
        """
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = pAdicValuation(R)
            sage: xbar = v.reduce(R(4))
            sage: v.lift(xbar)
            1 + O(3^20)

        """
        if x.parent() is not self.residue_field():
            raise ValueError("x must be in the residue field of the valuation")

        return self.domain()(x).lift_to_precision()

    def uniformizer(self):
        """
        Return a uniformizer of this valuation, i.e., `p` as an element of the
        domain.

        EXAMPLES::

            sage: v = pAdicValuation(Zp(3))
            sage: v.uniformizer()
            3 + O(3^21)

        """
        return self.domain().uniformizer()

    def residue_field(self):
        """
        Return the residue field of the ring of integers of this valuation.

        EXAMPLES::

            sage: v = pAdicValuation(Zp(3))
            sage: v.residue_field()
            Finite Field of size 3

        """
        return self.domain().residue_field()

    def shift(self, c, v):
        """
        Multiply ``c`` by a ``v``th power of the uniformizer.

        INPUT:

        - ``c`` -- an element of the domain of this valuation

        - ``v`` -- an integer

        OUTPUT:

        If the resulting element has negative valation, then it will be brought
        to the fraction field, otherwise it is an element of the domain.

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = pAdicValuation(Zp(3))
            sage: v.shift(R(2),3)
            2*3^3 + O(3^23)

        """
        from sage.rings.all import ZZ
        if c.parent() is not self.domain():
            raise ValueError("c must be in the domain of the valuation")
        if not v in ZZ:
            raise ValueError("v must be an integer")
        v = ZZ(v)

        return c<<v

class pAdicValuation_int(pAdicValuation_base):
    pass

class pAdicValuation_number_field(pAdicValuation_base):
    def __init__(self, R, valuation):
        p = valuation.residue_field().characteristic()
        assert(p>0)
        pAdicValuation_base.__init__(self, R, valuation.uniformizer()(R.gen()))
        self._valuation = valuation

    def _mac_lane_step(self):
        E,F = self._valuation.E(),self._valuation.F()
        self._valuation = self._valuation.mac_lane_step(self.domain().relative_polynomial())
        assert len(self._valuation)== 1
        self._valuation = self._valuation[0]
        assert E == self._valuation.E()
        assert F == self._valuation.F()

    @cached_method
    def residue_field(self):
        return self._valuation.residue_field()

    def _repr_(self):
        return "%s-adic valuation"%(self._prime)

    def _call_(self, x):
        if x.parent() is not self.domain():
            raise ValueError("x must be in the domain of the valuation")

        if x.is_zero():
            from sage.rings.all import infinity
            return infinity
        
        f = self._valuation.domain()(x.list())
        while True:
            vals = self._valuation.valuations(f)
            ret = min(vals)
            argmin = [1 for t in vals if t == ret]
            if len(argmin)>1:
                self._mac_lane_step()
                continue
            return ret*self._valuation.E()

    def reduce(self, x):
        if x.parent() is not self.domain():
            raise ValueError("x must be in the domain of the valuation")

        if self(x)>0:
            return self.residue_field().zero()
        f = self._valuation.domain()(x.list())
        while True:
            ret = self._valuation.reduce(f)
            if ret.degree():
                self._mac_lane_step()
                continue
            return ret[0]

    def lift(self, x):
        if x.parent() is not self.residue_field():
            raise ValueError("x must be in the residue field of the valuation")
        if x.is_one():
            return self.domain().one()
        if x.is_zero():
            return self.domain().zero()

        return self._valuation.lift(x)(self.domain().gen())

# Fix doctests so they work in standalone mode
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.path.dirname(os.getcwd()))