r"""
Finite Fields

Sage supports arithmetic in finite prime and extension fields.
Several implementation for prime fields are implemented natively in
Sage for several sizes of primes `p`. These implementations
are


-  ``sage.rings.finite_rings.integer_mod.IntegerMod_int``,

-  ``sage.rings.finite_rings.integer_mod.IntegerMod_int64``, and

-  ``sage.rings.finite_rings.integer_mod.IntegerMod_gmp``.


Small extension fields of cardinality `< 2^{16}` are
implemented using tables of Zech logs via the Givaro C++ library
(``sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro``).
While this representation is very fast it is limited to finite
fields of small cardinality. Larger finite extension fields of
order `q >= 2^{16}` are internally represented as
polynomials over smaller finite prime fields. If the
characteristic of such a field is 2 then NTL is used internally to
represent the field
(``sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e``).
In all other case the PARI C library is used
(``sage.rings.finite_rings.finite_field_ext_pari.FiniteField_ext_pari``).

However, this distinction is internal only and the user usually
does not have to worry about it because consistency across all
implementations is aimed for. In all extension field
implementations the user may either specify a minimal polynomial or
leave the choice to Sage.

For small finite fields the default choice are Conway polynomials.

The Conway polynomial `C_n` is the lexicographically first
monic irreducible, primitive polynomial of degree `n` over
`GF(p)` with the property that for a root `\alpha`
of `C_n` we have that
`\beta=
\alpha^{(p^n - 1)/(p^m - 1)}` is a root of
`C_m` for all `m` dividing `n`. Sage
contains a database of Conway polynomials which also can be queried
independently of finite field construction.

While Sage supports basic arithmetic in finite fields some more
advanced features for computing with finite fields are still not
implemented. For instance, Sage does not calculate embeddings of
finite fields yet.

EXAMPLES::

    sage: k = GF(5); type(k)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>

::

    sage: k = GF(5^2,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>

::

    sage: k = GF(2^16,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>

::

    sage: k = GF(3^16,'c'); type(k)
    <class 'sage.rings.finite_rings.finite_field_ext_pari.FiniteField_ext_pari_with_category'>

Finite Fields support iteration, starting with 0.

::

    sage: k = GF(9, 'a')
    sage: for i,x in enumerate(k):  print i,x
    0 0
    1 a
    2 a + 1
    3 2*a + 1
    4 2
    5 2*a
    6 2*a + 2
    7 a + 2
    8 1
    sage: for a in GF(5):
    ...    print a
    0
    1
    2
    3
    4

We output the base rings of several finite fields.

::

    sage: k = GF(3); type(k)
    <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

::

    sage: k = GF(9,'alpha'); type(k)
    <class 'sage.rings.finite_rings.finite_field_givaro.FiniteField_givaro_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

::

    sage: k = GF(3^40,'b'); type(k)
    <class 'sage.rings.finite_rings.finite_field_ext_pari.FiniteField_ext_pari_with_category'>
    sage: k.base_ring()
    Finite Field of size 3

Further examples::

    sage: GF(2).is_field()
    True
    sage: GF(next_prime(10^20)).is_field()
    True
    sage: GF(19^20,'a').is_field()
    True
    sage: GF(8,'a').is_field()
    True

AUTHORS:

- William Stein: initial version

- Robert Bradshaw: prime field implementation

- Martin Albrecht: Givaro and ntl.GF2E implementations
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import random

from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.structure.parent_gens import normalize_names

import sage.rings.arith as arith
import sage.rings.integer as integer

import sage.rings.polynomial.polynomial_element as polynomial_element
import sage.rings.polynomial.multi_polynomial_element as multi_polynomial_element

# We don't late import this because this means trouble with the Givaro library
# On a Macbook Pro OSX 10.5.8, this manifests as a Bus Error on exiting Sage.
# TODO: figure out why
from finite_field_givaro import FiniteField_givaro

import sage.interfaces.gap
import sage.databases.conway

from sage.structure.factory import UniqueFactory
from sage.structure.sage_object import SageObject
import weakref

class FiniteFieldFactory(UniqueFactory):
    """
    Return the globally unique finite field of given order with
    generator labeled by the given name and possibly with given
    modulus.

    INPUT:


    -  ``order`` - int

    -  ``name`` - string; must be specified unless
       - ``order`` is prime
       - ``modulus`` is 'conway', in which case a name is generated
          from the degree and the ``prefix`` string.

    -  ``modulus`` - (optional) either a defining polynomial for the
       field, i.e., generator of the field will be a root of this
       polynomial; or a string:

           - 'conway': force the use of a pseudo-Conway polynomial.
            If a conway polynomial is found in the database, that is used.
            Otherwise, a polynomial is used with the same algebraic properties
            but without the lexicographic constraint implying uniqueness.
           - 'random': use a random irreducible polynomial;
           - 'default': a Conway polynomial is used if in the database. Otherwise
            a sparse polynomial is used for binary fields and a
            random polynomial is used for other characteristics.
           - 'first_lexicographic': uses the first irreducible polynomial.

       Other options might be available depending on the
       implementation.

    -  ``elem_cache`` - cache all elements to avoid
       creation time (default: order < 500)

    -  ``check_irreducible`` - verify that the polynomial
       modulus is irreducible

    - ``prefix`` (default: 'z') - a string used to generate names
      for subfields and pushouts of this field (mainly for use with
      pseudo-Conway polynomials).

    - ``proof`` -- bool (default: True); if True use provable
      primality test; otherwise only use pseudoprimality test.

    -  ``args`` - additional parameters passed to finite
       field implementations

    -  ``kwds`` - additional keyword parameters passed to
       finite field implementations


    ALIAS: You can also use GF instead of FiniteField - they are
    identical.

    EXAMPLES::

        sage: k.<a> = FiniteField(9); k
        Finite Field in a of size 3^2
        sage: parent(a)
        Finite Field in a of size 3^2
        sage: charpoly(a, 'y')
        y^2 + 2*y + 2
        sage: L = GF(27); L
        Finite Field in z3 of size 3^3

    We illustrate the proof flag.  The following example would hang
    for a very long time if we didn't use proof=False.  (NOTE: Magma
    only supports proof=False for making finite fields, so falsely
    appears to be faster than Sage -- see Trac 10975.)::

        sage: k = FiniteField(10^1000 + 453, proof=False)
        sage: k = FiniteField((10^1000 + 453)^2, 'a', proof=False)      # long time -- about 5 seconds

    ::

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x +1 )
        sage: f = K.modulus(); f
        x^5 + 4*x + 1
        sage: type(f)
         <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>

    The modulus must be irreducible::

        sage: K.<a> = GF(5**5, name='a', modulus=x^5 - x )
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not.

    You can't accidentally fool the constructor into thinking the modulus
    is irreducible when it isn't mod p, since it actually tests
    irreducibility modulo p.  Also, the modulus has to be of the right degree.

    ::

        sage: F.<x> = QQ[]
        sage: factor(x^5 + 2)
        x^5 + 2
        sage: K.<a> = GF(5**5, name='a', modulus=x^5 + 2 )
        Traceback (most recent call last):
        ...
        ValueError: finite field modulus must be irreducible but it is not.
        sage: K.<a> = GF(5**5, name='a', modulus=x^3 + 3*x + 3)
        Traceback (most recent call last):
        ...
        ValueError: The degree of the modulus does not correspond to the
        cardinality of the field.

    If you wish to live dangerously, you can tell the constructor not
    to test irreducibility using check_irreducible=False, but this can
    easily lead to crashes and hangs - so do not do it unless you know
    that the modulus really is irreducible and has the correct degree!

    ::

        sage: F.<x> = GF(5)[]
        sage: K.<a> = GF(5**2, name='a', modulus=x^2 + 2, check_irreducible=False)

    ::

        sage: L = GF(3**2, name='a', modulus=QQ[x](x - 1), check_irreducible=False)
        sage: L.list()  # random
        [0, a, 1, 2, 1, 2, 1, 2, 1]

    The order of a finite field must be a prime power::

        sage: GF(1)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be > 1.
        sage: GF(100)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be a prime power.

    Finite fields with explicit random modulus are not cached::

        sage: k.<a> = GF(5**10, modulus='random')
        sage: n.<a> = GF(5**10, modulus='random')
        sage: n is k
        False
        sage: GF(5**10, 'a') is GF(5**10, 'a')
        True

    We check that various ways of creating the same finite field yield
    the same object, which is cached.

    ::

        sage: K = GF(7, 'a')
        sage: L = GF(7, 'b')
        sage: K is L
        True
        sage: K = GF(4,'a'); K.modulus()
        x^2 + x + 1
        sage: L = GF(4,'a', K.modulus())
        sage: K is L
        True
        sage: M = GF(4,'a', K.modulus().change_variable_name('y'))
        sage: K is M
        True

    You may print finite field elements as integers. This
    currently only works if the order of field is `<2^{16}`,
    though.

    ::

        sage: k.<a> = GF(2^8, repr='int')
        sage: a
        2

    """
    def create_key_and_extra_args(self, order, name=None, modulus=None, names=None, prefix='z', impl=None, proof=None, **kwds):
        """
        EXAMPLES::

            sage: GF.create_key_and_extra_args(9, 'a')
            ((9, ('a',), 'conwayz', None, '{}'), {})
            sage: GF.create_key_and_extra_args(9, 'a', foo='value')
            ((9, ('a',), 'conwayz', None, "{'foo': 'value'}"), {'foo': 'value'})
        """
        from sage.structure.proof.all import WithProof, arithmetic
        if proof is None: proof = arithmetic()
        with WithProof('arithmetic', proof):
            order = int(order)
            if order <= 1:
                raise ValueError("the order of a finite field must be > 1.")

            if arith.is_prime(order):
                name = None
                modulus = None
                p = integer.Integer(order)
                n = integer.Integer(1)
            elif arith.is_prime_power(order):
                if not names is None: name = names

                p, n = arith.factor(order)[0]

                if modulus is None or modulus == "default":
                    if exists_conway_polynomial(p, n) or name is None:
                        modulus = "conway"
                    else:
                        if p == 2:
                            modulus = "minimal_weight"
                        else:
                            modulus = "random"
                elif p != 2 and modulus == 'first_lexicographic':
                    R = FiniteField(p)['x']
                    x = R.gen()
                    for f in R.polynomials(of_degree=n):
                        if f.is_irreducible():
                            f = f.monic() # probably not necessary...
                            break
                    modulus = f
                    
                if modulus == "conway":
                    modulus += prefix
                    if name is None:
                        name = prefix + str(n)
                elif modulus == "random":
                    modulus += str(random.randint(0, 1<<128))
                name = normalize_names(1,name)

                if isinstance(modulus, (list, tuple)):
                    modulus = FiniteField(p)['x'](modulus)
                # some classes use 'random' as the modulus to
                # generate a random modulus, but we don't want
                # to cache it
                elif sage.rings.polynomial.polynomial_element.is_Polynomial(modulus):
                    modulus = modulus.change_variable_name('x')
                elif not isinstance(modulus, str):
                    raise ValueError("Modulus parameter not understood.")
            else:  # Neither a prime, nor a prime power
                raise ValueError("the order of a finite field must be a prime power.")

            return (order, name, modulus, impl, str(kwds), p, n, proof), kwds

    def create_object(self, version, key, check_irreducible=True, elem_cache=None,
                      names=None, **kwds):
        """
        EXAMPLES::

            sage: K = GF(19) # indirect doctest
            sage: TestSuite(K).run()
        """
        # IMPORTANT!  If you add a new class to the list of classes
        # that get cached by this factor object, then you *must* add
        # the following method to that class in order to fully support
        # pickling:
        #
        #     def __reduce__(self):   # and include good doctests, please!
        #         return self._factory_data[0].reduce_data(self)
        #
        # This is not in the base class for finite fields, since some finite
        # fields need not be created using this factory object, e.g., residue
        # class fields.

        if len(key) == 5:
            # for backward compatibility of pickles (see trac 10975).
            order, name, modulus, impl, _ = key
            p, n = arith.factor(order)[0]
            proof = True
        else:
            order, name, modulus, impl, _, p, n, proof = key

        if isinstance(modulus, str) and modulus.startswith("random"):
            modulus = "random"

        if elem_cache is None:
            elem_cache = order < 500

        if n == 1 and (impl is None or impl == 'modn'):
            from finite_field_prime_modn import FiniteField_prime_modn
            # Using a check option here is probably a worthwhile
            # compromise since this constructor is simple and used a
            # huge amount.
            K = FiniteField_prime_modn(order, check=False, **kwds)
        else:
            # We have to do this with block so that the finite field
            # constructors below will use the proof flag that was
            # passed in when checking for primality, factoring, etc.
            # Otherwise, we would have to complicate all of their
            # constructors with check options (like above).
            from sage.structure.proof.all import WithProof
            with WithProof('arithmetic', proof):
                if check_irreducible and polynomial_element.is_Polynomial(modulus):
                    if modulus.parent().base_ring().characteristic() == 0:
                        modulus = modulus.change_ring(FiniteField(p))
                    if not modulus.is_irreducible():
                        raise ValueError("finite field modulus must be irreducible but it is not.")
                    if modulus.degree() != n:
                        raise ValueError("The degree of the modulus does not correspond to the cardinality of the field.")
                if name is None:
                    raise TypeError, "you must specify the generator name or use Conway polynomials"
                if order < zech_log_bound and (impl is None or impl =='givaro'):  
                    # DO *NOT* use for prime subfield, since that would lead to
                    # a circular reference in the call to ParentWithGens in the
                    # __init__ method.
                    K = FiniteField_givaro(order, name, modulus, cache=elem_cache,**kwds)
                else:
                    if order % 2 == 0 and (impl is None or impl == 'ntl'):
                        from finite_field_ntl_gf2e import FiniteField_ntl_gf2e
                        K = FiniteField_ntl_gf2e(order, name, modulus, **kwds)
                    else:
                        from finite_field_ext_pari import FiniteField_ext_pari
                        K = FiniteField_ext_pari(order, name, modulus, **kwds)

        return K

    def other_keys(self, key, K):
        """
        EXAMPLES::

            sage: key, extra = GF.create_key_and_extra_args(9, 'a'); key
            (9, ('a',), 'conway', None, '{}', 3, 2, True)
            sage: K = GF.create_object(0, key); K
            Finite Field in a of size 3^2
            sage: GF.other_keys(key, K)
            [(9, ('a',), x^2 + 2*x + 2, None, '{}', 3, 2, True),
             (9, ('a',), x^2 + 2*x + 2, 'givaro', '{}', 3, 2, True)]
        """
        if len(key) == 5: # backward compat
            order, name, modulus, impl, _ = key
            p, n = arith.factor(order)[0]
            proof = True
        else:
            order, name, modulus, impl, _, p, n, proof = key

        from sage.structure.proof.all import WithProof
        with WithProof('arithmetic', proof):
            if K.degree() > 1:
                modulus = K.modulus().change_variable_name('x')
            new_keys = [(order, name, modulus, impl, _, p, n, proof)]
            from finite_field_prime_modn import FiniteField_prime_modn
            if isinstance(K, FiniteField_prime_modn):
                impl = 'modn'
            elif isinstance(K, FiniteField_givaro):
                impl = 'givaro'
            else:
                from finite_field_ntl_gf2e import FiniteField_ntl_gf2e
                from finite_field_ext_pari import FiniteField_ext_pari
                if isinstance(K, FiniteField_ntl_gf2e):
                    impl = 'ntl'
                elif isinstance(K, FiniteField_ext_pari):
                    impl = 'pari'
            new_keys.append( (order, name, modulus, impl, _, p, n, proof) )
            return new_keys


GF = FiniteField = FiniteFieldFactory("FiniteField")


def is_PrimeFiniteField(x):
    """
    Returns True if x is a prime finite field.

    EXAMPLES::

        sage: from sage.rings.finite_rings.constructor import is_PrimeFiniteField
        sage: is_PrimeFiniteField(QQ)
        False
        sage: is_PrimeFiniteField(GF(7))
        True
        sage: is_PrimeFiniteField(GF(7^2,'a'))
        False
        sage: is_PrimeFiniteField(GF(next_prime(10^90,proof=False)))
        True
    """
    from finite_field_prime_modn import FiniteField_prime_modn
    from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic

    return isinstance(x, FiniteField_prime_modn) or \
           (isinstance(x, FiniteField_generic) and x.degree() == 1)

##################################################################

def conway_polynomial(p, n):
    r"""
    Return the Conway polynomial of degree n over GF(p), which is
    loaded from a table.

    If the requested polynomial is not known, this function raises a
    RuntimeError exception.

    INPUT:


    -  ``p`` - int

    -  ``n`` - int


    OUTPUT:


    -  ``Polynomial`` - a polynomial over the prime finite
       field GF(p).


    .. note::

       The first time this function is called a table is read from
       disk, which takes a fraction of a second. Subsequent calls do
       not require reloading the table.

    See also the ``ConwayPolynomials()`` object, which is a
    table of Conway polynomials. For example, if
    ``c=ConwayPolynomials()``, then
    ``c.primes()`` is a list of all primes for which the
    polynomials are known, and for a given prime `p`,
    ``c.degree(p)`` is a list of all degrees for which the
    Conway polynomials are known.

    EXAMPLES::

        sage: conway_polynomial(2,5)
        x^5 + x^2 + 1
        sage: conway_polynomial(101,5)
        x^5 + 2*x + 99
        sage: conway_polynomial(97,101)
        Traceback (most recent call last):
        ...
        RuntimeError: requested conway polynomial not in database.
    """
    (p,n)=(int(p),int(n))
    R = FiniteField(p)['x']
    try:
        return R(sage.databases.conway.ConwayPolynomials()[p][n])
    except KeyError:
        raise RuntimeError("requested conway polynomial not in database.")

def exists_conway_polynomial(p, n):
    r"""
    Return True if the Conway polynomial over `F_p` of degree
    `n` is in the database and False otherwise.

    If the Conway polynomial is in the database, to obtain it use the
    command ``conway_polynomial(p,n)``.

    EXAMPLES::

        sage: exists_conway_polynomial(2,3)
        True
        sage: exists_conway_polynomial(2,-1)
        False
        sage: exists_conway_polynomial(97,200)
        False
        sage: exists_conway_polynomial(6,6)
        False
    """
    return sage.databases.conway.ConwayPolynomials().has_polynomial(p,n)

pseudo_conway_poly = {}

class PseudoConwayPolyTree(SageObject):
    """
    An object holding references to pseudo-Conway polynomials for divisors of the given degree, so that they aren't garbage collected.
    """
    def __init__(self, p, n, nodes_dict, f):
        """
        INPUT::

        - p -- The prime for which this defines an extension of GF(p).
        - n -- The degree of the extension.
        - nodes_dict -- (dict or bool).  A dictionary of PseudoConwayPolyTrees, indexed by prime divisors of n.  The entry for q corresponds to the pseudo-Conway extension of degree n//q.
                        OR a boolean.  If True, then this polynomial is actually in the Conway polynomials database, and no references to other PCPTs are stored.  If False, then n is prime and no references are stored (since there is no compatiblity condition).
        - f -- The polynomial defining this extension.

        EXAMPLES::

            sage: from sage.rings.finite_rings.constructor import find_pseudo_conway_polynomial_tree
            sage: PCPT = find_pseudo_conway_polynomial_tree(2, 6, False)
            sage: PCPT.get_pseudo_conway_poly(3)
            x^3 + x + 1
        """
        self.p = p
        self.n = n
        self.nodes_dict = nodes_dict
        self.f = f
        
    def get_pseudo_conway_poly(self, m):
        """
        Returns the pseudo-Conway polynomial associated to a divisor of the degree of this tree.

        EXAMPLES::
        
            sage: from sage.rings.finite_rings.constructor import find_pseudo_conway_polynomial_tree
            sage: PCPT = find_pseudo_conway_polynomial_tree(2, 6, False)
            sage: PCPT.get_pseudo_conway_poly(3)
            x^3 + x + 1
            sage: PCPT.get_pseudo_conway_poly(4)
            Traceback (most recent call last):
            ...
            RuntimeError: 4 does not divide 6
        """
        if m == self.n:
            return self.f
        if not m.divides(self.n):
            raise RuntimeError, "%s does not divide %s"%(m, self.n)
        try:
            return pseudo_conway_poly[self.p][m]().f
        except (KeyError, AttributeError):
            return conway_polynomial(self.p, m)
            
    def check_consistency(self):
        """
        Checks that the pseudo-Conway polynomials dividing the degree of this tree satisfy the required compatibility conditions.

        EXAMPLES::

            sage: from sage.rings.finite_rings.constructor import find_pseudo_conway_polynomial_tree
            sage: PCPT = find_pseudo_conway_polynomial_tree(2, 6, False)
            sage: PCPT.check_consistency()
            sage: PCPT = find_pseudo_conway_polynomial_tree(2, 60, False) # long
            sage: PCPT.check_consistency() # long
        """
        p, n = self.p, self.n
        K = GF(p**n, modulus = self.f, names='a')
        a = K.gen()
        for m in n.divisors():
            assert (a**((p**n-1)//(p**m-1))).minimal_polynomial() == self.get_pseudo_conway_poly(m)

def find_pseudo_conway_polynomial_tree(p, n, use_database=True):
    """
    Returns an object holding references to a set of consistent pseudo-Conway polynomials for degrees dividing n.

    A Conway polynomial `f \in \Bold{F}_p` of degree `n` satisfies the following three conditions:
    
    - `f` is irreducible.
    - In the quotient field `\Bold{F}_p[x]/(f)`, the indeterminant `x` is a multiplicative generator.
    - In this same quotient field, the minimal polynomial of `x^{\frac{p^n-1}{p^m-1}}` is the Conway polynomial of degree `m` for every divisor `m` of `n`.
    - `f` is lexicographically least among all such polynomials, under a certain ordering.

    The final condition is needed only in order to make the Conway polynomial unique.  We define a pseudo-Conway polynomial to be one satisfying the first three conditions.

    INPUT:

    - `p` -- a prime.
    - `n` -- an integer greater than 1.
    - `use_database` -- whether to use the Conway polynomials database if the Conway polynomial for `p, n` exists within it (versus computing a pseudo-Conway polynomial)

    EXAMPLES::

        sage: from sage.rings.finite_rings.constructor import find_pseudo_conway_polynomial_tree
        sage: PCPT = find_pseudo_conway_polynomial_tree(2, 12, False)
        sage: PCPT.f
        x^12 + x^10 + x^9 + x^8 + x^4 + x^2 + 1
    """
    if not pseudo_conway_poly.has_key(p):
        pseudo_conway_poly[p] = {}
    pdict = pseudo_conway_poly[p]
    n = sage.rings.integer.Integer(n)
    if pdict.has_key(n):
        pcp = pdict[n]()
        if pcp is not None:
            return pcp
    if use_database and exists_conway_polynomial(p, n):
        ans = PseudoConwayPolyTree(p, n, True, GF(p)['x'](conway_polynomial(p, n)))
        pdict[n] = weakref.ref(ans)
        return ans
    if n == 1 or n.is_prime():
        ans = PseudoConwayPolyTree(p, n, False, compute_pseudo_conway_polynomial(p, n, None))
    nodes = {}
    for q in n.prime_factors():
        nodes[q] = find_pseudo_conway_polynomial_tree(p, n//q, use_database)
    ans = PseudoConwayPolyTree(p, n, nodes, compute_pseudo_conway_polynomial(p, n, nodes))
    pdict[n] = weakref.ref(ans)
    return ans

def _find_pow_of_frobenius(p, x, y, field_degree):
    """
    Finds the power of Frobenius which yields `x` when applied to `y`.

    INPUT:

    - `p` -- a prime.
    - ``field_degree`` -- a positive integer
    - `x` -- an element of a field `K` of characeristic p so that the multiplicative order of `x` is ``p^field_degree``.
    - `y` -- an element of `K` with the same minimal polynomial.
    
    OUTPUT:
    
    An element `i` of ``Integers(field_degree)`` so that `x = y^{p^i}`

    EXAMPLES::

        sage: from sage.rings.finite_rings.constructor import _find_pow_of_frobenius
        sage: K.<a> = GF(3^14)
        sage: x = K.multiplicative_generator()
        sage: y = x^27
        sage: _find_pow_of_frobenius(3, x, y, 14)
        11
    """
    assert x.minimal_polynomial() == y.minimal_polynomial()
    assert x.multiplicative_order() == p**field_degree - 1
    assert y.multiplicative_order() == p**field_degree - 1
    from integer_mod import mod
    for i in range(field_degree):
        if x == y: break
        y = y**p
    else:
        raise RuntimeError, "No appropriate power of Frobenius found"
    return mod(i, field_degree)

def _crt_non_coprime(running, a):
    """
    Extends the ``crt`` method of IntegerMod to the case of non-relatively prime modulus.

    EXAMPLES::

        sage: from sage.rings.finite_rings.constructor import _crt_non_coprime
        sage: a = _crt_non_coprime(mod(14, 18), mod(20,30)); a
        50
        sage: a.modulus()
        90
        sage: _crt_non_coprime(mod(13, 18), mod(20,30))
        Traceback (most recent call last):
        ...
        AssertionError
    """
    g = running.modulus().gcd(a.modulus())
    if g == 1:
        return running.crt(a)
    else:
        assert running % g == a % g
        running_modulus = running.modulus()
        a_modulus = a.modulus()
        for qq in g.prime_divisors():
            a_val_unit = a_modulus.val_unit(qq)
            running_val_unit = running_modulus.val_unit(qq)
            if a_val_unit[0] > running_val_unit[0]:
                running_modulus = running_val_unit[1]
            else:
                a_modulus = a_val_unit[1]
        return (running % running_modulus).crt(a % a_modulus)

def _frobenius_shift(K, generators, check_only=False):
    """
    Given a field `K` of degree `n` over ``GF(p)`` and a dictionary holding, for each divisor `q` of `n`, an element with minimal polynomial a pseudo-Conway polynomial of degree `n/q`, modifies these generators into a compatibile system.

    Such a system of generators is said to be compatible if for each pair of prime divisors `q_1` and `q_2` and each common divisor `m` of `n/q_1` and `n/q_2`, the equality

    ``generators[q1]^((p^(n/q1)-1)/(p^m-1)) == generators[q2]^((p^(n/q2)-1)/(p^m-1))``

    holds.

    INPUT:

    - K -- a finite field of degree n    
    - generators -- a dictionary, indexed by prime divisors `q` of `n`, whose entries are elements of `K` satisfying the `n/q` pseudo-Conway polynomial.
    - check_only -- if True, just asserts that the generators given form a compatible system.

    EXAMPLES::

        sage: R.<x> = GF(2)[]
        sage: f30 = x^30 + x^28 + x^27 + x^25 + x^24 + x^20 + x^19 + x^18 + x^16 + x^15 + x^12 + x^10 + x^7 + x^2 + 1
        sage: f20 = x^20 + x^19 + x^15 + x^13 + x^12 + x^11 + x^9 + x^8 + x^7 + x^4 + x^2 + x + 1
        sage: f12 = x^12 + x^10 + x^9 + x^8 + x^4 + x^2 + 1
        sage: K.<a> = GF(2^60, modulus='first_lexicographic')
        sage: x30 = f30.any_root(K)
        sage: x20 = f20.any_root(K)
        sage: x12 = f12.any_root(K)
        sage: generators = {2: x30, 3: x20, 5: x12}
        sage: from sage.rings.finite_rings.constructor import _frobenius_shift, _find_pow_of_frobenius
        sage: _frobenius_shift(K, generators)
        sage: _find_pow_of_frobenius(2, x30, generators[2], 30)
        0
        sage: _find_pow_of_frobenius(2, x20, generators[3], 20)
        13
        sage: _find_pow_of_frobenius(2, x12, generators[5], 12)
        8
    """
    p = K.characteristic()
    n = K.degree()
    if len(generators) == 1:
        return generators
    compatible = {}
    from integer_mod import mod
    for m in n.divisors():
        compatible[m] = {}
    for q, x in generators.iteritems():
        for m in (n//q).divisors():
            compatible[m][q] = x**((p**(n//q)-1)//(p**m-1))
    if check_only:
        for m in n.divisors():
            try:
                q, x = compatible[m].popitem()
            except KeyError:
                break
            for qq, xx in compatible[m].iteritems():
                assert x == xx
        return
    crt = {}
    qlist = sorted(generators.keys())
    for j in range(1, len(qlist)):
        for i in range(j):
            crt[(i, j)] = []
    for m in n.divisors():
        mqlist = sorted(compatible[m].keys())
        for k in range(1,len(mqlist)):
            j = qlist.index(mqlist[k])
            i = qlist.index(mqlist[k-1])
            crt[(i,j)].append(_find_pow_of_frobenius(p, compatible[m][qlist[j]], compatible[m][qlist[i]], m))
    from integer_mod import mod
    pairs = crt.keys()
    for i, j in pairs:
        L = crt[(i,j)]
        #print qlist[i], qlist[j]
        #print " ".join(["(%s, %s)"%(mm, mm.modulus()) for mm in L])
        running = mod(0,1)
        for a in L:
            running = _crt_non_coprime(running, a)
        #print "(%s, %s)"%(running, running.modulus())
        crt[(i,j)] = [(mod(running, q**(running.modulus().valuation(q))), running.modulus().valuation(q)) for q in qlist]
        crt[(j,i)] = [(-a, level) for a, level in crt[(i,j)]]
    # Let x_j be the power of Frobenius we apply to generators[qlist[j]], for 0 < j < len(qlist)
    # We have some direct conditions on the x_j: x_j reduces to each entry in crt[(0,j)].
    # But we also have the equations x_j - x_i reduces to each entry in crt[(i,j)].
    # We solve for x_j one prime at a time.  For each prime, we have an equations of the form
    # x_j - x_i = c_ij.  The modulus of the currently known value of x_j, x_i and c_ij will all be powers 
    # (possibly 0, possibly different) of the same prime.

    # We can set x_0=0 everywhere, can get an initial setting of x_j from the c_0j.
    # We go through prime by prime.
    import bisect
    frob_powers=[mod(0,1) for q in qlist]
    def find_leveller(qindex, level, x, xleveled, searched, i):
        searched[i] = True
        crt_possibles = []
        for j in range(1,len(qlist)):
            if i==j: continue
            if crt[(i,j)][qindex][1] >= level:
                if xleveled[j]:
                    return [j]
                elif not searched.has_key(j):
                    crt_possibles.append(j)
        for j in crt_possibles:
            path = find_leveller(qindex, level, x, xleveled, searched, j)
            if path is not None:
                path.append(j)
                return path
        return None
    def propogate_levelling(qindex, level, x, xleveled, i):
        for j in range(1, len(qlist)):
            if i==j: continue
            if not xleveled[j] and crt[(i,j)][qindex][1] >= level:
                newxj = x[i][0] + crt[(i,j)][qindex][0]
                x[j] = (newxj, min(x[i][1], crt[(i,j)][qindex][1]))
                xleveled[j] = True
                propogate_levelling(qindex, level, x, xleveled, j)
                
    for qindex in range(len(qlist)):
        q = qlist[qindex]
        x = [0] + [crt[(0,j)][qindex] for j in range(1,len(qlist))] # we include the initial 0 to match up our indexing with crt.                
        # We first check that our equations are consistent and determine which powers of q occur as moduli.
        levels = []
        for j in range(2, len(qlist)):
            for i in range(j):
                # we need crt[(0,j)] = crt[(0,i)] + crt[(i,j)]
                if i != 0: assert x[j][0] == x[i][0] + crt[(i,j)][qindex][0]
                level = crt[(i,j)][qindex][1]
                if level > 0:
                    ins = bisect.bisect_left(levels,level)
                    if ins == len(levels):
                        levels.append(level)
                    elif levels[ins] != level:
                        levels.insert(ins, level)
        for level in levels:
            xleveled = [0] + [x[i][1] >= level for i in range(1,len(qlist))]
            while True:
                try:
                    i = xleveled.index(False, 1)
                    searched = {}
                    levelling_path = find_leveller(qindex, level, x, xleveled, searched, i)
                    if levelling_path is None:
                        # Any lift will work, since there are no constraints
                        x[i] = (mod(x[i][0].lift(), q**level), level)
                        xleveled[i] = True
                        propogate_levelling(qindex, level, x, xleveled, i)
                    else:
                        levelling_path.append(i)
                        for m in range(1,len(path)):
                            # this point on the path may have already been leveled in a previous propogation
                            if not xleveled[path[m]]:
                                newx = x[path[m-1]][0] + crt[(path[m-1],path[m])][qindex][0]
                                x[path[m]] = (newx, min(x[path[m-1]][1], crt[(path[m-1],path[m])][qindex][1]))
                                xleveled[path[m]] = True
                                propogate_levelling(qindex, level, x, xleveled, path[m])
                except ValueError:
                    break
        for j in range(1,len(qlist)):
            frob_powers[j] = frob_powers[j].crt(x[j][0])
    for j in range(1, len(qlist)):
        #print "frob_power[%s] = mod(%s, %s)"%(qlist[j], frob_powers[j], frob_powers[j].modulus())
        generators[qlist[j]] = generators[qlist[j]]**(p**(-frob_powers[j]).lift())
    _frobenius_shift(K, generators, check_only=True)

def compute_pseudo_conway_polynomial(p, n, nodes):
    """
    Computes a pseudo-Conway polynomial over the prime `p` of degree `n`.

    A Conway polynomial `f \in \Bold{F}_p` of degree `n` satisfies the following three conditions:
    
    - `f` is irreducible.
    - In the quotient field `\Bold{F}_p[x]/(f)`, the indeterminant `x` is a multiplicative generator.
    - In this same quotient field, the minimal polynomial of `x^{\frac{p^n-1}{p^m-1}}` is the Conway polynomial of degree `m` for every divisor `m` of `n`.
    - `f` is lexicographically least among all such polynomials, under a certain ordering.

    The final condition is needed only in order to make the Conway polynomial unique.  We define a pseudo-Conway polynomial to be one satisfying the first three conditions.

    INPUT:
    
    - `p` -- a prime
    - `n` -- a positive integer
    - ``nodes`` -- None (in the case that `n` is prime) or a dictionary of PseudoConwayPolyTree objects, indexed by prime divisors `q` of `n` with entries corresponding to pseudo-Conway polynomials of degree `n/q`.

    OUTPUT:

    A pseudo-Conway polynomial over the prime `p` of degree `n`.

    ALGORITHM:
    
    Uses an algorithm described in [HL99]_, modified to find pseudo-Conway polynomials rather than Conway polynomials.  The major modification was the addition of the function _frobenius_shift.

    REFERENCE:

    .. [HL99] Heath L. and Loehr N. (1999). 
       New algorithms for generating Conway polynomials over finite fields.
       Proceedings of the tenth annual ACM-SIAM symposium on Discrete algorithms, pp. 429-437.

    EXAMPLES::

        sage: from sage.rings.finite_rings.constructor import compute_pseudo_conway_polynomial, find_pseudo_conway_polynomial_tree
        sage: compute_pseudo_conway_polynomial(next_prime(10000),11,None) # long because the arithmetic in FiniteField_ext_pari fields is slow.
        x^11 + x + 7
        sage: PCPT30 = find_pseudo_conway_polynomial_tree(2, 30, False)
        sage: PCPT20 = find_pseudo_conway_polynomial_tree(2, 20, False)
        sage: PCPT12 = find_pseudo_conway_polynomial_tree(2, 12, False)
        sage: compute_pseudo_conway_polynomial(2, 60, {2:PCPT30,3:PCPT20,5:PCPT12})
        x^60 + x^59 + x^56 + x^53 + x^49 + x^45 + x^43 + x^42 + x^41 + x^40 + x^39 + x^37 + x^36 + x^34 + x^31 + x^30 + x^27 + x^26 + x^25 + x^23 + x^21 + x^19 + x^17 + x^15 + x^14 + x^11 + x^9 + x^8 + x^7 + x^6 + x^5 + x^2 + 1
    """
    k = GF(p)
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    R = PolynomialRing(k, names='x')
    if n == 1:
        return R([-k.multiplicative_generator(), 1])
    if nodes is None: 
        # n is prime.  We require only that the chosen generator (x) be a multiplicative generator for the units.
        if p in [2,3,5,7,11]: # here we have Cunningham data
            from sage.rings.factorint import factor_cunningham
            F = [(p**n-1)//a for a, e in factor_cunningham(n)]
        else:
            F = (p**n-1).prime_factors() # could be VERY expensive
        x = R.gen()
        from sage.rings.arith import power_mod
        for f in R.polynomials(of_degree=n):
            # By using this iterator our chosen polynomial will be fairly sparse.
            if f.is_irreducible():
                if not f.is_monic():
                    f = f.monic()
                for m in F:
                    if power_mod(x, m, f) == 1:
                        break
                else:
                    return f
    # Our strategy is as follows:
    # Work in an arbitrary field K of order p**n.
    # We're looking for a multiplicative generator of K whose 
    K = GF(p**n, modulus="first_lexicographic", names='a')
    r = p**n - 1
    xi = {}
    for q in nodes.iterkeys():
        #print nodes[q].f.degree(), n//q
        xi[q] = nodes[q].f.any_root(K, -n//q, True)
    _frobenius_shift(K, xi)
            
    q, x = xi.popitem()
    v = p**(n//q) - 1
    for q, xitem in xi.iteritems():
        fi = p**(n//q) - 1
        g, alpha, beta = v.xgcd(fi)
        x = x**beta * xitem**alpha
        v = v.lcm(fi)
    g = r // v
    #print g, n, g.divides(p**n-1)
    z = x.nth_root(g, cunningham=(p in [2,3,5,7,11]))
    groot1 = K.multiplicative_generator()**v
    while z.multiplicative_order() != r:
        z *= groot1
    #print z
    return z.minimal_polynomial()

zech_log_bound = 2**16
