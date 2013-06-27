r"""
Function Field Morphisms

AUTHORS:

- William Stein (2010): initial version

- Julian Rueth (2011-09-14): refactored class hierarchy

EXAMPLES::

    sage: K.<x> = FunctionField(QQ); R.<y> = K[]
    sage: K.hom(1/x)
    Morphism of function fields defined by x |--> 1/x
    sage: L.<y> = K.extension(y^2-x)
    sage: K.hom(y)
    Morphism of function fields defined by x |--> y
    sage: L.hom([y,x])
    Morphism of function fields defined by y |--> y,  x |--> x
    sage: L.hom([x,y])
    Traceback (most recent call last):
    ...
    ValueError: invalid morphism
"""
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Julian Rueth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism
from sage.categories.map import Map
from sage.rings.morphism import RingHomomorphism

class FunctionFieldDerivation(Map):
    r"""
    A base class for derivations of function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: d = K.derivation()
        sage: isinstance(d, sage.rings.function_field.maps.FunctionFieldDerivation)
        True

    """
    def __init__(self, K):
        r"""
        Initialize a derivation from ``K`` to ``K``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation() # indirect doctest

        """
        self.__field = K
        from sage.categories.homset import Hom
        Map.__init__(self, Hom(K,K))

    def _repr_type(self):
        r"""
        Return the type of this map (a derivation), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d._repr_type()
            'Derivation'

        """
        return "Derivation"

    def is_injective(self):
        r"""
        Return whether this derivation is injective.

        OUTPUT:

        Returns ``False`` since derivations are never injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d.is_injective()
            False

        """
        return False

    def domain(self):
        r"""
        Return the domain of this derivation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.derivation().domain()
            Rational function field in x over Rational Field

        """
        return self.__field

    def codomain(self):
        r"""
        Return the codomain of this derivation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: K.derivation().codomain()
            Rational function field in x over Rational Field

        """
        return self.__field

class FunctionFieldDerivation_rational(FunctionFieldDerivation):
    r"""
    A derivation of rational function fields over the constant base field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: d = K.derivation()
        sage: isinstance(d, sage.rings.function_field.maps.FunctionFieldDerivation_rational)
        True

    """
    def __init__(self, K, u):
        r"""
        Initialize a derivation of ``K`` which sends the generator of ``K`` to
        ``u``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation() # indirect doctest

        """
        FunctionFieldDerivation.__init__(self, K)
        self._u = u

    def _call_(self, x):
        r"""
        Compute the derivation of ``x``.

        INPUT:

        - ``x`` -- an element of the rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d(x) # indirect doctest
            1
            sage: d(x^3)
            3*x^2

        """
        f,g = x.numerator(),x.denominator()

        if not f.gcd(g).is_one():
            raise NotImplementedError("derivations only implemented for rational functions with coprime numerator and denominator.")

        numerator = f.derivative()*g - f*g.derivative()
        if numerator.is_zero():
            return self.codomain().zero()
        else:
            return self._u * self.codomain()( numerator / g**2 )

class FunctionFieldDerivation_polymod(FunctionFieldDerivation):
    r"""
    A derivation of function fields over the constant base field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: d = L.derivation()
        sage: isinstance(d, sage.rings.function_field.maps.FunctionFieldDerivation_polymod)
        True

    """
    def __init__(self, L, d):
        r"""
        Initialize a derivation of ``L`` which extends the derivation ``d`` on
        its base field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()

        """
        FunctionFieldDerivation.__init__(self, L)
        self._d = d

    def _call_(self, x):
        r"""
        Compute the derivative of ``x``.

        INPUT:

        - ``x`` -- an element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d(x) # indirect doctest
            1
            sage: d(y)
            (-1/2/-x)*y
            sage: d(y^2)
            1

        """
        if x.is_zero():
            return self.codomain().zero()

        f = x.minpoly()
        if not f.gcd(f.derivative()).is_one():
            raise NotImplementedError("derivations only for separable elements.")
        return - f.map_coefficients(lambda c:self._d(c))(x) / f.derivative()(x)

class FunctionFieldIsomorphism(Morphism):
    r"""
    A base class for isomorphisms between function fields and
    vector spaces.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space()
        sage: isinstance(f, sage.rings.function_field.maps.FunctionFieldIsomorphism)
        True
    """
    def _repr_type(self):
        """
        Return the type of this map (an isomorphism), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def is_injective(self):
        """
        Return True, since this isomorphism is injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return True, since this isomorphism is surjective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_surjective()
            True
        """
        return True

class MapVectorSpaceToFunctionField(FunctionFieldIsomorphism):
    r"""
    An isomorphism from a vector space to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); f
        Isomorphism morphism:
          From: Vector space of dimension 2 over Rational function field in x over Rational Field
          To:   Function field in y defined by y^2 - x*y + 4*x^3
    """
    def __init__(self, V, K):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(f)
            <class 'sage.rings.function_field.maps.MapVectorSpaceToFunctionField'>
        """
        self._V = V
        self._K = K
        self._R = K.polynomial_ring()
        from sage.categories.homset import Hom
        FunctionFieldIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f(x*V.0 + (1/x^3)*V.1)         # indirect doctest
            1/x^3*y + x
        """
        f = self._R(self._V(v).list())
        return self._K(f)

    def domain(self):
        """
        Return the vector space which is the domain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.domain()
            Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        return self._V

    def codomain(self):
        """
        Return the function field which is the codomain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.codomain()
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self._K

class MapFunctionFieldToVectorSpace(FunctionFieldIsomorphism):
    """
    An isomorphism from a function field to a vector space.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); t
        Isomorphism morphism:
          From: Function field in y defined by y^2 - x*y + 4*x^3
          To:   Vector space of dimension 2 over Rational function field in x over Rational Field
    """
    def __init__(self, K, V):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(t)
            <class 'sage.rings.function_field.maps.MapFunctionFieldToVectorSpace'>
        """
        self._V = V
        self._K = K
        self._zero = K.base_ring()(0)
        self._n = K.degree()
        from sage.categories.homset import Hom
        FunctionFieldIsomorphism.__init__(self, Hom(K, V))

    def domain(self):
        """
        Return the function field which is the domain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t.domain()
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self._K

    def codomain(self):
        """
        Return the vector space which is the domain of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t.codomain()
            Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        return self._V

    def _repr_type(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t(x + (1/x^3)*y)                       # indirect doctest
            (x, 1/x^3)
        """
        y = self._K(x)
        v = y.list()
        w = v + [self._zero]*(self._n - len(v))
        return self._V(w)

class FunctionFieldMorphism(RingHomomorphism):
    r"""
    Base class for morphisms between function fields.
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: f = K.hom(1/x); f
            Morphism of function fields defined by x |--> 1/x
            sage: isinstance(f, sage.rings.function_field.maps.FunctionFieldMorphism)
            True
        """
        RingHomomorphism.__init__(self, parent)

        self._im_gen = im_gen
        self._base_morphism = base_morphism

    def is_injective(self):
        """
        Returns True since homomorphisms of fields are injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: f = K.hom(1/x); f
            Morphism of function fields defined by x |--> 1/x
            sage: f.is_injective()
            True
        """
        return True

    def __repr__(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f.__repr__()
            'Morphism of function fields defined by y |--> 2*y'
        """
        return "Morphism of function fields defined by %s"%self._short_repr()

    def _short_repr(self):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f._short_repr()
            'y |--> 2*y'
        """
        a = '%s |--> %s'%(self.domain().gen(), self._im_gen)
        if self._base_morphism is not None:
            a += ',  ' + self._base_morphism._short_repr()
        return a

class FunctionFieldMorphism_polymod(FunctionFieldMorphism):
    """
    Morphism from a finite extension of a function field to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: f = L.hom(-y); f
        Morphism of function fields defined by y |--> -y
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2); f
            Morphism of function fields defined by y |--> 2*y
            sage: type(f)
            <class 'sage.rings.function_field.maps.FunctionFieldMorphism_polymod'>
            sage: factor(L.polynomial())
            y^3 + 6*x^3 + x
            sage: f(y).charpoly('y')
            y^3 + 6*x^3 + x
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, base_morphism)
        # Verify that the morphism is valid:
        R = self.codomain()['X']
        v = parent.domain().polynomial().list()
        if base_morphism is not None:
            v = [base_morphism(a) for a in v]
        f = R(v)
        if f(im_gen):
            raise ValueError, "invalid morphism"

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f(y/x + x^2/(x+1))            # indirect doctest
            2/x*y + x^2/(x + 1)
            sage: f(y)
            2*y
        """
        v = x.list()
        if self._base_morphism is not None:
            v = [self._base_morphism(a) for a in v]
        f = v[0].parent()['X'](v)
        return f(self._im_gen)

class FunctionFieldMorphism_rational(FunctionFieldMorphism):
    """
    Morphism from a rational function field to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: f = K.hom(1/x); f
        Morphism of function fields defined by x |--> 1/x
    """
    def __init__(self, parent, im_gen):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Morphism of function fields defined by x |--> 1/x
            sage: type(f)
            <class 'sage.rings.function_field.maps.FunctionFieldMorphism_rational'>
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, None)

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Morphism of function fields defined by x |--> 1/x
            sage: f(x+1)                          # indirect doctest
            (x + 1)/x
            sage: 1/x + 1
            (x + 1)/x
        """
        a = x.element()
        return a.subs({a.parent().gen():self._im_gen})
