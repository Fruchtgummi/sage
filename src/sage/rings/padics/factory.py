"""
Factories

Factory classes and constructor functions to create `p`-adic rings

AUTHORS:

- David Roe

- Julian Rueth (2012-10-24): rewrote large parts, added support for two step
  extensions and generic extensions

EXAMPLES::

    sage: K = Qp(5); K
    5-adic Field with capped relative precision 20
    sage: R = Zp(5); R
    5-adic Ring with capped relative precision 20

Unramified extensions::

    sage: K.<u> = Qq(25); K
    Unramified Extension of 5-adic Field with capped relative precision 20 in u defined by (1 + O(5^20))*x^2 + (4 + O(5^20))*x + (2 + O(5^20))
    sage: R.<u> = Zq(25); R
    Unramified Extension of 5-adic Ring with capped relative precision 20 in u defined by (1 + O(5^20))*x^2 + (4 + O(5^20))*x + (2 + O(5^20))

Eisenstein extensions::

    sage: K = Qp(5)
    sage: R.<a> = K[]
    sage: L.<a> = K.extension(a^2 + 5); L
    Eisenstein Extension of 5-adic Field with capped relative precision 20 in a defined by (1 + O(5^20))*a^2 + (5 + O(5^21))

Trivial extensions::

    sage: R.<t> = K[]
    sage: L.<t> = K.extension(t + 1); L
    Trivial extension of 5-adic Field with capped relative precision 20 in t defined by (1 + O(5^20))*t + (1 + O(5^20))

    sage: K.<u> = Qq(25)
    sage: R.<t> = K[]
    sage: L.<t> = K.extension(t + u); L
    Trivial extension of Unramified Extension of 5-adic Field with capped relative precision 20 in u defined by (1 + O(5^20))*x^2 + (4 + O(5^20))*x + (2 + O(5^20)) in t defined by (1 + O(5^20))*t + u + O(5^20)

    sage: K = Qp(5)
    sage: R.<a> = K[]
    sage: L.<a> = K.extension(a^2 + 5)
    sage: R.<t> = K[]
    sage: M.<t> = L.extension(t + a); M
    Trivial extension of Eisenstein Extension of 5-adic Field with capped relative precision 20 in a defined by (1 + O(5^20))*a^2 + (5 + O(5^21)) in t defined by (1 + O(a^40))*t + a + O(a^41)

An iterated trivial extension::

    sage: R.<s> = M[]
    sage: N.<s> = M.extension(s + t); N
    Trivial extension of Trivial extension of Eisenstein Extension of 5-adic Field with capped relative precision 20 in a defined by (1 + O(5^20))*a^2 + (5 + O(5^21)) in t defined by (1 + O(a^40))*t + a + O(a^41) in s defined by (1 + O(a^40))*s + 4*a + a^3 + O(a^41)

A totally ramified extension not defined by an Eisenstein polynomial::

    sage: K = Qp(5)
    sage: R.<a> = K[]
    sage: L.<a> = K.extension(a^2 + 5*a + 25); L

.. _padic_precision

Precision
---------

There are two notions of precision for `p`-adic elements.

The first is relative precision, which gives the number of known `p`-adic
digits::

    sage: R = Zp(5)
    sage: a = R(675); a
    2*5^2 + 5^4 + O(5^22)
    sage: a.precision_relative()
    20

The second one is absolute precision, which gives the power of the uniformizer
that this element is defined modulo::

    sage: a.precision_absolute()
    22

When creating a `p`-adic ring, one can chose from several implementations which
determine how elements keep track of their precision.

In a *capped-relative* ring, the relative precision of elements is bounded,
this is the default implementation::

    sage: R = Zp(5, 10); R
    5-adic Ring with capped relative precision 10

When an exact element is brought into such a ring, its relative precision is
truncated to the precision cap of the ring::

    sage: R(5)
    5 + O(5^11)

Consequently, an exact zero is stored to infinite absolute precision::

    sage: R(0).precision_absolute()
    +Infinity

In a *capped-absolute* ring, the absolute precision of elements is bounded::

    sage: R = Zp(5, 10, 'capped-abs')
    sage: R(5)
    5 + O(5^10)

In a *fixed-mod* ring, the absolute precision of all elements is the same. All
computations are performed module a power of the uniformizer and the elements
do not keep track of precision losses::

    sage: R = Zp(5, 10, 'fixed-mod')
    sage: R(5)
    5 + O(5^10)
    sage: R(0)
    O(5^10)

Notice how the following computation differs in these three precision models::

    sage: R = Zp(5, 10, 'fixed-mod')
    sage: a = R(5^9); a
    5^9 + O(5^10)
    sage: b = a>>9; b
    1 + O(5^10)
    sage: a + b
    1 + 5^9 + O(5^10)

    sage: R = Zp(5, 10, 'capped-abs')
    sage: a = R(5^9); a
    5^9 + O(5^10)
    sage: b = a>>9; b
    1 + O(5)
    sage: a + b
    1 + O(5)

    sage: R = Zp(5, 10, 'capped-rel')
    sage: a = R(5^9); a
    5^9 + O(5^19)
    sage: b = a>>9; b
    1 + O(5^10)
    sage: a + b
    1 + 5^9 + O(5^10)

Eventually, in a *lazy* implementation, elements will be able to increase their
precision upon request. This is currently not implemented.

.. _padic_printing

Printing
--------

Upon creation of a `p`-adic ring, several parameters control how elements
print. The following list describe the effect of different values for the
parameter ``print_mode``:

* ``'series'`` (the default):

    Elements are displayed as series in the uniformizer::

        sage: K = Qp(5, 10, print_mode='series')
        sage: K(70700)
        3*5^2 + 3*5^4 + 2*5^5 + 4*5^6 + O(5^12)
        sage: K(-70700)
        2*5^2 + 4*5^3 + 5^4 + 2*5^5 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

        sage: K.<u> = Qq(9, 10, 'capped-rel', print_mode='series')
        sage: (1+2*u)^4
        2 + (2*u + 2)*3 + (2*u + 1)*3^2 + O(3^10)
        sage: -3*(1+2*u)^4
        3 + u*3^2 + 3^3 + (2*u + 2)*3^4 + (2*u + 2)*3^5 + (2*u + 2)*3^6 + (2*u + 2)*3^7 + (2*u + 2)*3^8 + (2*u + 2)*3^9 + (2*u + 2)*3^10 + O(3^11)
        sage: ~(3*u+18)
        (u + 2)*3^-1 + 1 + 2*3 + (u + 1)*3^2 + 3^3 + 2*3^4 + (u + 1)*3^5 + 3^6 + 2*3^7 + (u + 1)*3^8 + O(3^9)

    ``print_pos`` controls whether negative integers can be used in the
    coefficients::

        sage: K = Qp(5, 10, print_mode='series', print_pos=False)
        sage: K(70700)
        -2*5^2 + 5^3 - 2*5^4 - 2*5^5 + 5^7 + O(5^12)
        sage: K(-70700)
        2*5^2 - 5^3 + 2*5^4 + 2*5^5 - 5^7 + O(5^12)

    This also works for extensions::

        sage: K.<u> = Qq(9, print_mode='series', print_pos=False)
        sage: (1+2*u)^4
        -1 - u*3 - 3^2 + (u + 1)*3^3 + O(3^20)
        sage: -3*(1+2*u)^4
        3 + u*3^2 + 3^3 + (-u - 1)*3^4 + O(3^21)

    ``names`` affects how the generator is printed. For base rings, this
    controls how the prime is printed::

        sage: K.<p> = Qp(5)
        sage: p
        p + O(p^21)

    In general, ``ram_name`` controls how the uniformizer is printed::

        sage: K.<u> = Qq(9, print_mode='series', ram_name='p')
        sage: 3*(1+2*u)^4
        2*p + (2*u + 2)*p^2 + (2*u + 1)*p^3 + O(p^21)

    ``print_max_ram_terms`` limits the number of powers of the
    uniformizer that appear::

        sage: K = Qp(5, 10, print_mode='series', print_max_ram_terms=4)
        sage: repr(R(-70700)) # have to use repr() here so the doctest does not get confused by the ...
        '2*5^2 + 4*5^3 + 5^4 + 2*5^5 + ... + O(5^12)'

        sage: K.<u> = Qq(9, print_mode='series', print_max_ram_terms=4)
        sage: repr(-3*(1+2*u)^4)
        '3 + u*3^2 + 3^3 + (2*u + 2)*3^4 + ... + O(3^21)'

    In extensions, ``print_max_unram_terms`` limits the number of terms that
    appear in a coefficient of a power of the uniformizer::

        sage: K.<u> = Qq(128, prec = 8, print_mode='series')
        sage: repr((1+u)^9)
        '(u^3 + 1) + (u^5 + u^4 + u^3 + u^2)*2 + (u^6 + u^5 + u^4 + u + 1)*2^2 + (u^5 + u^4 + u^2 + u + 1)*2^3 + (u^6 + u^5 + u^4 + u^3 + u^2 + u + 1)*2^4 + (u^5 + u^4)*2^5 + (u^6 + u^5 + u^4 + u^3 + u + 1)*2^6 + (u + 1)*2^7 + O(2^8)'
        sage: K.<u> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 3)
        sage: repr((1+u)^9)
        '(u^3 + 1) + (u^5 + u^4 + ... + u^2)*2 + (u^6 + u^5 + ... + 1)*2^2 + (u^5 + u^4 + ... + 1)*2^3 + (u^6 + u^5 + ... + 1)*2^4 + (u^5 + u^4)*2^5 + (u^6 + u^5 + ... + 1)*2^6 + (u + 1)*2^7 + O(2^8)'
        sage: K.<u> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 2)
        sage: repr((1+u)^9)
        '(u^3 + 1) + (u^5 + ... + u^2)*2 + (u^6 + ... + 1)*2^2 + (u^5 + ... + 1)*2^3 + (u^6 + ... + 1)*2^4 + (u^5 + u^4)*2^5 + (u^6 + ... + 1)*2^6 + (u + 1)*2^7 + O(2^8)'
        sage: K.<u> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 1)
        sage: repr((1+u)^9)
        '(u^3 + ...) + (u^5 + ...)*2 + (u^6 + ...)*2^2 + (u^5 + ...)*2^3 + (u^6 + ...)*2^4 + (u^5 + ...)*2^5 + (u^6 + ...)*2^6 + (u + ...)*2^7 + O(2^8)'
        sage: K.<u> = Qq(128, prec = 8, print_mode='series', print_max_unram_terms = 0)
        sage: repr((1+u)^9 - 1 - u^3)
        '(...)*2 + (...)*2^2 + (...)*2^3 + (...)*2^4 + (...)*2^5 + (...)*2^6 + (...)*2^7 + O(2^8)'

    ``print_sep``, ``print_alphabet``, and ``print_max_terse_terms`` have no effect.

* ``'val-unit'``:

    Elements are displayed as ``pi^k*u``::

        sage: K = Qp(5, print_mode='val-unit')
        sage: K(70700)
        5^2 * 2828 + O(5^22)
        sage: K(-707/5)
        5^-1 * 95367431639918 + O(5^19)

    ``print_pos`` controls whether to use negative integers::

        sage: K = Qp(5, print_mode='val-unit', print_pos=False)
        sage: K(-70700)
        5^2 * (-2828) + O(5^22)

    This also works for extensions::

        sage: K.<u> = Qq(9, 7, print_mode='val-unit', print_pos=False)
        sage: b = (1+3*u)^9 - 1; b
        3^3 * (15 - 17*u) + O(3^7)
        sage: ~b
        3^-3 * (-40 + u) + O(3)

    ``names`` affects how the generator is printed. For base rings, this
    changes how the prime is printed::

        sage: K.<pi> = Qp(5, print_mode='val-unit')
        sage: K(70700)
        pi^2 * 2828 + O(pi^22)

    In general, ``ram_name`` controls how the uniformizer is printed::

        sage: K.<u> = Qq(5, ram_name='p')
        sage: u
        p + O(p^21)

    In an extension, ``print_max_terse_terms`` controls how many terms appear
    in the unit part::

        sage: K.<u> = Qq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)
        sage: b = ~(17*(u^3-u+14)); b
        17^-1 * (22110411 + 11317400*u + 20656972*u^2 + ...) + O(17^5)
        sage: b*17*(u^3-u+14)
        1 + O(17^6)

    ``print_max_ram_terms``, ``print_sep``, ``print_max_unram_terms``, and
    ``print_alphabet`` have no effect.

* ``'terse'``

    Elements are displayed as an integer in base 10 or the quotient of an
    integer by a power of `p` (still in base 10)::

        sage: K = Qp(5, print_mode='terse')
        sage: K(70700)
        70700 + O(5^22)
        sage: K(-70700)
        2384185790944925 + O(5^22)
        sage: K(707/5^2)
        707/5^2 + O(5^18)

    ``print_pos`` controls whether to use a balanced representation::

        sage: K = Qp(5, print_mode='terse', print_pos=False)
        sage: K(-70700)
        -70700 + O(5^22)
        sage: K(-707/5)
        -707/5 + O(5^19)

    ``names`` affects how the generator is printed. For a base field, this
    changes how the prime is printed::

        sage: K.<unif> = Qp(5, print_mode='terse')
        sage: K(-707/5)
        95367431639918/unif + O(unif^19)
        sage: K(-707/5^10)
        95367431639918/unif^10 + O(unif^10)

    In an extension, ``ram_name`` controls how the uniformizer is printed::

        sage: K.<u> = Qq(125, print_mode='terse', ram_name='p')
        sage: (u - 1/5)^6
        95367431620001/p^6 + 18369/p^5*u + 1353/p^3*u^2 + O(p^14)

    ``print_max_terse_terms`` controls how many terms of the polynomial are
    shown::

        sage: K.<u> = Qq(625, print_mode='terse', print_max_terse_terms=2)
        sage: (u-1/5)^6
        106251/5^6 + 49994/5^5*u + ... + O(5^14)

    ``print_max_ram_terms``, ``print_max_unram_terms``, ``print_sep`` and
    ``print_alphabet`` have no effect.

* ``'digits'``:

    Elements are displayed as a string of base `p` digits::

        sage: K = Qp(5, print_mode='digits')
        sage: repr(K(70700)) # the repr() is necessary in these doctests to handle the leading ...
        '...4230300'
        sage: repr(K(-70700))
        '...4444444444444440214200'
        sage: repr(K(-707/5))
        '...4444444444444443413.3'
        sage: repr(K(-707/5^2))
        '...444444444444444341.33'

    Note that it's not possible to read off the precision from the
    representation in this mode.

    Naturally, this mode can not be used in unramified extensions::

        sage: K.<u> = Qq(25, print_mode='digits')
        Traceback (most recent call last):
        ...
        ValueError: digits printing mode only usable for totally ramified extensions with p at most the length of the alphabet (default 62).  Try using print_mode = 'bars' instead.

    ``print_max_ram_terms`` limits the number of digits that are printed.  Note
    that if the valuation of the element is very negative, more digits will be
    printed::

        sage: K = Qp(5, print_mode='digits', print_max_ram_terms=4)
        sage: repr(K(-70700))
        '...214200'
        sage: repr(K(-707/5^2))
        '...41.33'
        sage: repr(K(-707/5^6))
        '...?.434133'
        sage: repr(K(-707/5^6,absprec=-2))
        '...?.??4133'
        sage: repr(K(-707/5^4))
        '...?.4133'

    ``print_alphabet`` defines the symbols used to print the digits, defaults
    to the tuple of length 62 composed of ``'0'`` to ``'9'``, followed by
    ``'A'`` to ``'Z'``, followed by ``'a'`` to ``'z'``.

    The digits printing mode is only available if `p` is at most the length of
    ``print_alphabet``::

        sage: K = Qp(5, print_mode='digits', print_alphabet=('zero','one','two','three','four'))
        sage: repr(K(707))
        '...onezerothreeonetwo'

        sage: K = Qp(7, print_mode='digits', print_alphabet=('zero','one','two','three','four'))
        Traceback (most recent call last):
        ...
        ValueError: digits printing mode only usable for totally ramified extensions with p at most the length of the alphabet (default 62).  Try using print_mode = 'bars' instead.

    Restriction: you can only use the digits printing mode for
    small primes.  Namely, `p` must be less than the length of the
    alphabet tuple (default alphabet has length 62).::

    ``print_max_terse_terms``, ``print_max_unram_terms``, ``print_pos``,
    ``names`` and ``print_sep`` have no effect.

* ``'bars'``:

   Elements are displayed as a string of base `p` digits with separators::

        sage: K = Qp(5, print_mode='bars')
        sage: repr(K(70700))
        '...4|2|3|0|3|0|0'
        sage: repr(K(-70700))
        '...4|4|4|4|4|4|4|4|4|4|4|4|4|4|4|0|2|1|4|2|0|0'
        sage: repr(K(-707/5^2))
        '...4|4|4|4|4|4|4|4|4|4|4|4|4|4|4|3|4|1|.|3|3'

    This also works for extensions::

        sage: K.<u> = Qq(125)
        sage: (u+5)^6
        (4*u^2 + 3*u + 4) + (3*u^2 + 2*u)*5 + (u^2 + u + 1)*5^2 + (3*u + 2)*5^3 + (3*u^2 + u + 3)*5^4 + (2*u^2 + 3*u + 2)*5^5 + O(5^20)
        sage: K.<u> = Qq(125, print_mode='bars', prec=8)
        sage: repr((u+5)^6)
        '...[2, 3, 2]|[3, 1, 3]|[2, 3]|[1, 1, 1]|[0, 2, 3]|[4, 3, 4]'
        sage: repr((u-5)^6)
        '...[0, 4]|[1, 4]|[2, 0, 2]|[1, 4, 3]|[2, 3, 1]|[4, 4, 3]|[2, 4, 4]|[4, 3, 4]'

    Note that elements with negative valuation are shown with a
    decimal point at valuation 0.::

        sage: repr((u+1/5)^6)
        '...[3]|[4, 1, 3]|.|[1, 2, 3]|[3, 3]|[0, 0, 3]|[0, 1]|[0, 1]|[1]'
        sage: repr((u+1/5)^2)
        '...[0, 0, 1]|.|[0, 2]|[1]'

    If not enough precision is known, ``'?'`` is used instead.::

        sage: repr((u+K(1/5,relprec=3))^7)
        '...|.|?|?|?|?|[0, 1, 1]|[0, 2]|[1]'

    Again, it is not possible to read off the precision from this
    representation::

        sage: b = u + 3; repr(b)
        '...[3, 1]'
        sage: c = u + K(3, 4); repr(c)
        '...[3, 1]'
        sage: b.precision_absolute()
        8
        sage: c.precision_absolute()
        4

    ``print_pos`` controls whether the digits can be negative::

        sage: K = Qp(5, print_mode='bars',print_pos=False)
        sage: repr(K(-70700))
        '...-1|0|2|2|-1|2|0|0'

    ``print_max_ram_terms`` limits the number of digits that are printed.  Note
    that if the valuation of the element is very negative, more digits will be
    printed::

        sage: K = Qp(5, print_mode='bars', print_max_ram_terms=4)
        sage: repr(K(-70700))
        '...2|1|4|2|0|0'
        sage: repr(K(-707/5^2))
        '...4|1|.|3|3'
        sage: repr(K(-707/5^6))
        '...|.|4|3|4|1|3|3'
        sage: repr(K(-707/5^6,absprec=-2))
        '...|.|?|?|4|1|3|3'
        sage: repr(K(-707/5^4))
        '...|.|4|1|3|3'

    Note that this puts a cap on the relative precision, not the
    absolute precision::

        sage: K.<u> = Qq(125, print_mode='bars', print_max_ram_terms=3, print_pos=False)
        sage: repr((u-5)^6)
        '...[0, 0, -1]|[-2]|[-1, -2, -1]'
        sage: repr(5*(u-5)^6+50)
        '...[0, 0, -1]|[]|[-1, -2, -1]|[]'

    However, if the element has negative valuation, digits are shown
    up to the decimal point.::

        sage: repr((u-1/5)^6)
        '...|.|[-2, -1, -1]|[2, 2, 1]|[0, 0, -2]|[0, -1]|[0, -1]|[1]'

    ``print_sep`` controls the separation character::

        sage: K = Qp(5, print_mode='bars', print_sep='][')
        sage: repr(K(70700))
        '...4][2][3][0][3][0][0'

    In an extension, ``print_max_unram_terms`` controls how many terms are
    shown in each *digit*::

        sage: K.<u> = Qq(125, print_mode='bars')
        sage: b = (u+5)^6
        sage: with local_print_mode(K, {'max_unram_terms': 3}): repr(b)
        '...[2,..., 3, 2]|[3,..., 1, 3]|[2, 3]|[1,..., 1, 1]|[0,..., 2, 3]|[4,..., 3, 4]'
        sage: with local_print_mode(K, {'max_unram_terms': 2}): repr(b)
        '...[2,..., 2]|[3,..., 3]|[2, 3]|[1,..., 1]|[0,..., 3]|[4,..., 4]'
        sage: with local_print_mode(K, {'max_unram_terms': 1}): repr(b)
        '...[..., 2]|[..., 3]|[..., 3]|[..., 1]|[..., 3]|[..., 4]'
        sage: with local_print_mode(K, {'max_unram_terms':0}): repr(b-75*a)
        '...[...]|[...]|[...]|[...]|[...]|[...]|[...]|[...]|[...]|[]|[]|[]|[]|[]|[...]|[...]|[...]|[...]|[...]|[...]'

    ``names``, ``print_max_terse_terms``, and ``print_alphabet`` have no
    effect.

In general, equality of rings depends on printing options::

    sage: K = Qp(5, print_max_ram_terms=4)
    sage: L = Qp(5, print_max_ram_terms=5)
    sage: K == L
    False

However, irrelevant options are not taken into account::

    sage: K = Qp(5, print_alphabet=())
    sage: L = Qp(5, print_alphabet=('0'))
    sage: K == L
    True

:: _padic_modulus

Moduli
------

To create an extension ring, the modulus can be given in a number of ways:

* Implicitely:

    For unramified extensions, the modulus does not have to be specified::

        sage: K.<a> = Qq(27); K
        Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (2 + O(3^20))*x + (1 + O(3^20))

    In this case, the modulus is the standard lift of the generator chosen for
    `\mathbb{F}_q`::

        sage: GF(27, 'a').modulus()
        x^3 + 2*x + 1
        sage: K.modulus()
        (1 + O(3^20))*x^3 + (2 + O(3^20))*x + (1 + O(3^20))

* As a polynomial:

    The base ring can be `\mathbb{Z}`, `\mathbb{Q}`, `\mathbb{Z}_p`, `\mathbb{Q}_p`, `\mathbb{F}_p`.::

        sage: P.<x> = ZZ[]
        sage: R.<a> = Qq(27, modulus = x^3 + 2*x + 1); R.modulus()
        (1 + O(3^20))*x^3 + (2 + O(3^20))*x + (1 + O(3^20))
        sage: P.<x> = QQ[]
        sage: S.<a> = Qq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = Zp(3)[]
        sage: T.<a> = Qq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = Qp(3)[]
        sage: U.<a> = Qq(27, modulus = x^3 + 2*x + 1)
        sage: P.<x> = GF(3)[]
        sage: V.<a> = Qq(27, modulus = x^3 + 2*x + 1)

    Which form the modulus is given in has no effect on the extension
    produced::

        sage: R == S, S == T, T == U, U == V
        (True, True, True, False)

    unless the precision of the modulus differs. In the case of ``V``, the
    modulus is only given to precision 1, so the resulting field has a
    precision cap of 1::

        sage: V.precision_cap()
        1
        sage: U.precision_cap()
        20
        sage: P.<x> = Qp(3)[]
        sage: modulus = x^3 + (2 + O(3^7))*x + (1 + O(3^10))
        sage: modulus
        (1 + O(3^20))*x^3 + (2 + O(3^7))*x + (1 + O(3^10))
        sage: W.<a> = Qq(27, modulus = modulus); W.precision_cap()
        7

* As a symbolic expression:

        sage: x = var('x')
        sage: X.<a> = Qq(27, modulus = x^3 + 2*x + 1)
        sage: X.modulus()
        (1 + O(3^20))*x^3 + (2 + O(3^20))*x + (1 + O(3^20))
        sage: X == R
        True

"""
#*****************************************************************************
#       Copyright (C) 2007-2010 David Roe <roed@math.harvard.edu>
#                     2012-2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.structure.factory import UniqueFactory

class AbstractFactory(UniqueFactory):
    r"""
    Abstract base class for factories for `p`-adic rings.

    INPUT:

    - ``constructors`` -- a dictionary from strings to callables, mapping
      precision types to their respective implementation

    EXAMPLES::

        sage: from sage.rings.padics.factory import AbstractFactory
        sage: zp = AbstractFactory("zp", {"capped-rel" : ZpCR})

    """
    def __init__(self, name, constructors):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.padics.factory import AbstractFactory
            sage: isinstance(Zp, AbstractFactory)
            True
            sage: isinstance(Qp, AbstractFactory)
            True

        """
        UniqueFactory.__init__(self, name)
        self._constructors = constructors

    def _create_key_printing(self, print_mode, defaults, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms, names):
        r"""
        Shared functionality for ``create_key`` to normalize parameters related
        to printing and fill in default values.

        TESTS::

            sage: Zp._create_key_printing(None, {"print_mode":"series"}, None, None, None, None, None, None, None)
            ('series', None, None, None, None, None, None, None)

        """
        # if print_mode is a dictionary of parameters, pull out the parameters
        # if they don't contradict the parameters explicitely given
        if isinstance(print_mode, dict):
            for a,b in [("pos",print_pos), ("sep",print_sep), ("alphabet",print_alphabet), ("max_ram_terms",print_max_ram_terms), ("max_unram_terms",print_max_unram_terms), ("max_terse_terms",print_max_terse_terms)]:
                if a in print_mode:
                    if b is not None and print_mode[a] != b:
                        raise ValueError("contradicting parameters given for `%s`: `%s` and `%s`"%(a,print_mode[a],b))
                else:
                    print_mode[a] = b

            print_pos = print_mode["pos"]
            print_sep = print_mode["sep"]
            print_alphabet = print_mode["alphabet"]
            print_max_ram_terms = print_mode["max_ram_terms"]
            print_max_unram_terms = print_mode["max_unram_terms"]
            print_max_terse_terse = print_mode["max_terse_terms"]
            print_mode = print_mode["mode"]

        # fill unspecified parameters with default values
        def fill_with_default(key, value):
            if value is None:
                if key in defaults:
                    return defaults[key]
            return value

        print_mode = fill_with_default("print_mode", print_mode)
        print_pos = fill_with_default("print_pos", print_pos)
        print_sep = fill_with_default("print_sep", print_sep)
        print_alphabet = fill_with_default("print_alphabet", print_alphabet)
        if print_alphabet is not None: print_alphabet = tuple(print_alphabet)
        print_max_ram_terms = fill_with_default("print_max_ram_terms", print_max_ram_terms)
        print_max_unram_terms = fill_with_default("print_max_unram_terms", print_max_unram_terms)
        print_max_terse_terms = fill_with_default("print_max_terse_terms", print_max_terse_terms)

        # eliminate irrelevant print options
        if print_mode == "series":
            print_alphabet = None
            print_sep = None
            print_max_terse_terms = None
        elif print_mode == "val-unit":
            print_alphabet = None
            print_max_ram_terms = None
            print_sep = None
            print_max_unram_terms = None
        elif print_mode == "terse":
            print_alphabet = None
            print_max_ram_terms = None
            print_sep = None
            print_max_unram_terms = None
        elif print_mode == "digits":
            print_pos = None
            print_sep = None
            print_max_terse_terms = None
            print_max_unram_terms = None
        elif print_mode == "bars":
            print_alphabet = None
            print_max_terse_terms = None
        else:
            raise ValueError("unknown print mode `%s`"%print_mode)

        return print_mode, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms, names

    def _normalize_version(self, version, key):
        r"""
        Bring a ``key`` of an older ``version`` into the current format.

        INPUT:

        - ``version`` -- a version tuple

        - ``key`` -- a tuple of parameters

        EXAMPLES::

            sage: Qp._normalize_version((3,1), (2, 10, "capped-rel", "terse", None))
            (2, 10, 'capped-rel', 'terse', '2', None, None, None, None)

        """
        return key

    def create_object(self, version, key):
        """
        Create an object from the ``key``.

        TESTS::

            sage: Qp.create_object((3,4,2),(5, 41, 'capped-rel', 'series', '5', True, '|', (), -1))
            5-adic Field with capped relative precision 41

        """
        key = self._normalize_version(version, key)
        type, kwargs = self._decode_key(key)

        if type not in self._constructors:
            raise ValueError("unsupported precision type `%s`"%type)

        return self._constructors[type](**kwargs)

#####################################
#                                   #
#    p-adic base rings Zp and Qp    #
#                                   #
#####################################
class BaseFactory(AbstractFactory):
    r"""
    Return a `p`-adic ring.

    INPUT:

    - ``p`` -- an integer, the `p` in `\mathbb{Q}_p`

    - ``prec`` -- an integer or ``None`` (default: ``None``), the precision cap
      of the ring, see :ref:`padic_precision`

    - ``type`` -- a string (default: ``'capped-rel'``), the precision type of
      this ring, see :ref:`padic_precision`

    - ``print_mode`` -- a dictionary or one of ``'series'``, ``'val-unit'``,
      ``'terse'``, ``'digits'``, ``'bars'``, or ``None`` (default: ``None``),
      the way elements of this ring are printed, see :ref:`padic_printing`. If
      this is a dictionary, then the other printing parameters are read from it.

    - ``halt`` -- an integer or ``None`` (default: ``None``), currently
      irrelevant (to be used for lazy ring)

    - ``names`` -- a string, tuple consisting of a single string, or ``None``
      (default: ``None``), the string to print for `p`

    - ``print_pos`` -- a boolean (default: ``False``), whether to only use
      positive integers in the representations of elements, see
      :ref:`padic_printing`.

    - ``print_sep`` -- a string or ``None`` (default: ``None``), the separator
      character used when ``print_mode`` is ``'bars'``, see :ref:`padic_printing`

    - ``print_alphabet`` -- a tuple of strings or ``None`` (default: ``None``),
      the digits to use when ``print_mode`` is ``'digits'``, see
      :ref:`padic_printing`

    - ``print_max_ram_terms`` -- an integer or ``None`` (default: ``None``),
      the maximum number of terms to show in the series expansion, see
      :ref:`padic_printing`

    - ``check`` -- a boolean (default: ``True``), whether to check that ``p``
      is prime.

    EXAMPLES::

        sage: K = Qp(5); K
        5-adic Field with capped relative precision 20
        sage: R = Zp(5); R
        5-adic Ring with capped relative precision 20

    TESTS::

        sage: K = Qp(15, check=False)
        sage: a = K(999); a
        9 + 6*15 + 4*15^2 + O(15^20)

    """
    def create_key(self, p, prec = None, type = 'capped-rel', print_mode = None, halt = None, names = None, ram_name = None, print_pos = None, print_sep = None, print_alphabet = None, print_max_ram_terms = None, check = True):
        r"""
        Create a key which can be used to cache the `p`-adic ring specified by
        the parameters.

        See the classes docstring for the meaning of the parameters.

        This method essentially normalizes the parameters and fills in default
        values.

        TESTS::

            sage: Qp.create_key(5)
            (5, 20, 'capped-rel', 'series', '5', True, None, None, -1)
            sage: Zp.create_key(5)
            (5, 20, 'capped-rel', 'series', '5', True, None, None, -1)

        """
        # normalize parameters and fill in default values
        p = ZZ(p)
        if check:
            if not p.is_prime():
                raise ValueError("p must be prime")

        if prec is None:
            prec = 20
        prec = ZZ(prec)

        if halt is None:
            halt = 40
        halt = ZZ(halt)

        if isinstance(names, tuple):
            if len(names) != 1:
                raise ValueError("names must be a tuple consisting of one string")
            names = names[0]
        if ram_name is not None and names is not None:
            if ram_name != names:
                raise ValueError("ram_name and names must coincide")
        if ram_name is not None:
            names = ram_name
        if names is None:
            names = str(p)
        names = (names,)

        print_mode, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms, names = self._create_key_printing(print_mode, {"print_mode":padic_printing._printer_defaults.mode(), "print_pos":not padic_printing._printer_defaults.allow_negatives(), "print_sep":padic_printing._printer_defaults.sep(), "print_alphabet":padic_printing._printer_defaults.alphabet(), "print_max_ram_terms":padic_printing._printer_defaults.max_series_terms()}, print_pos, print_sep, print_alphabet, print_max_ram_terms, None, None, names)

        # eliminate irrelevant print options not handled by _create_key_printing
        if p == 2:
            print_pos = None

        return p, prec, type, print_mode, names[0], print_pos, print_sep, print_alphabet, print_max_ram_terms

    def _normalize_version(self, version, key):
        r"""
        Bring a ``key`` of an older ``version`` into the current format.

        INPUT:

        - ``version`` -- a version tuple

        - ``key`` -- a tuple of parameters

        EXAMPLES::

            sage: Qp._normalize_version((3,1), (2, 10, "capped-rel", "terse", None))
            (2, 10, 'capped-rel', 'terse', '2', None, None, None, None)

        """
        version = (version[0], version[1] if len(version)>1 else 0, version[2] if len(version)>2 else 0)

        if version[0] < 3 or (version[0] == 3 and version[1] < 2) or (version[0] == 3 and version[1] == 2 and version[2] < 3):
            p, prec, type, print_mode, names = key
            return self.create_key(p, prec, type, print_mode, None, names, None, None, None, None, None, False)

        if version[0] < 4 or (len(version) > 1 and version[0] == 4 and version[1] < 5) or (len(version) > 2 and version[0] == 4 and version[1] == 5 and version[2] < 3):
            # keys changed in order to reduce irrelevant duplications: e.g. two Qps with print_mode 'series' that are identical except for different 'print_alphabet' now return the same object.
            p, prec, type, print_mode, names, print_pos, print_sep, print_alphabet, print_max_ram_terms = key
            return self.create_key(p, prec, type, print_mode, None, names, None, print_pos, print_sep, print_alphabet, print_max_ram_terms, False)

        return key

    def _decode_key(self, key):
        r"""
        Split a ``key`` into the type, determining the implementation, and the
        parameters needed to construct an instance.

        INPUT:

        - ``key`` -- a tuple

        OUTPUT:

        A string determining the type of the implementation and a dictonary of
        parameters.

        EXAMPLES:

            sage: Qp._decode_key((2, 10, 'capped-rel', 'terse', '2', None, None, None, None))
            ('capped-rel', {'p': 2, 'print_mode': {'sep': None, 'alphabet': None, 'pos': None, 'mode': 'terse', 'ram_name': '2', 'max_ram_terms': None}, 'names': '2', 'prec': 10})

        """
        p, prec, type, print_mode, names, print_pos, print_sep, print_alphabet, print_max_ram_terms = key
        print_mode = { "mode":print_mode, "pos":print_pos, "sep":print_sep, "alphabet":print_alphabet, "ram_name":names, "max_ram_terms":print_max_ram_terms }
        return type, {"p":p, "prec":prec, "print_mode":print_mode, "names":names}

import padic_base_leaves
Qp = BaseFactory("Qp", { "capped-rel":padic_base_leaves.pAdicFieldCappedRelative })
Zp = BaseFactory("Zp", { "capped-rel":padic_base_leaves.pAdicRingCappedRelative,
                         "capped-abs":padic_base_leaves.pAdicRingCappedAbsolute,
                         "fixed-mod":padic_base_leaves.pAdicRingFixedMod })

def create_base_factory_with_type(base_factory, type):
    """
    Return a function which creates `p`-adic base rings.

    EXAMPLES:

        sage: from sage.rings.padics.factory import create_base_factory_with_type
        sage: zpca = create_base_factory_with_type(Zp, "capped-abs")
        sage: zpca(3)
        3-adic Ring with capped absolute precision 20

    """
    def create_ring(p, prec=None, print_mode=None, halt=None, names=None, print_pos=None, print_sep=None, print_alphabet=None, print_max_ram_terms=None, check=True):
        """
        Shortcut function to create a %s ring.

        See documentation for ``Zp`` and ``Qp`` for a description of the input
        parameters.

        EXAMPLES::

            sage: ZpCR(3, 10)
            3-adic Ring with capped relative precision 10
            sage: QpCR(3, 10)
            3-adic Field with capped relative precision 10
            sage: ZpCA(3, 10)
            3-adic Ring with capped absolute precision 10
            sage: ZpFM(3, 10)
            3-adic Ring of fixed modulus 3^10

        """%type
        return base_factory(p=p, prec=prec, print_mode=print_mode, halt=halt, names=names, print_pos=print_pos, print_sep=print_sep, print_alphabet=print_alphabet, print_max_ram_terms=print_max_ram_terms, check=check, type=type)
    return create_ring

QpCR = create_base_factory_with_type(Qp, "capped-rel")
ZpCR = create_base_factory_with_type(Zp, "capped-rel")
ZpCA = create_base_factory_with_type(Zp, "capped-abs")
ZpFM = create_base_factory_with_type(Zp, "fixed-mod")

####################################
#                                  #
#   unramified extension rings     #
#                                  #
####################################
def create_unramified_factory(base_factory):
    """
    Return a function which creates unramified extensions over a `p`-adic base
    ring.

    EXAMPLES::

        sage: from sage.rings.padics.factory import create_unramified_factory
        sage: zq = create_unramified_factory(Zp)
        sage: R.<u> = zq(9); R
        Unramified Extension of 3-adic Ring with capped relative precision 20 in u defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))

    """
    def create_unramified_ring(q, prec=None, type="capped-rel", modulus=None, names=None, print_mode=None, halt=None, ram_name=None, res_name=None, print_pos=None, print_sep=None, print_max_ram_terms=None, print_max_unram_terms=None, print_max_terse_terms=None, check=True):
        """
        Given a prime power `q = p^n`, return the unique unramified extension ring
        of degree `n`.

        INPUT:

        - ``q`` -- an integer

        - ``prec`` -- an integer or ``None`` (default: ``None``), the precision cap
          of the ring, see :ref:`padic_precision`.

        - ``type`` -- a string (default: ``'capped-rel'``), the precision type of
          this ring, see :ref:`padic_precision`

        - ``modulus`` -- a polynomial or ``None`` (default: ``None``), a polynomial
          defining an unramified extension of `\mathbb{Q}_p`, see
          :ref:`padic_modulus`.

        - ``names`` -- a string or a tuple consisting of a single string, the name
          of the generator, reducing to a generator of the residue field.

        - ``print_mode`` -- a dictionary or one of ``'series'``, ``'val-unit'``,
          ``'terse'``, ``'digits'``, ``'bars'``, or ``None`` (default: ``None``),
          the way elements of this ring are printed, see :ref:`padic_printing`. If
          this is a dictionary, then the other printing parameters are read from
          it.

        - ``halt`` -- an integer or ``None`` (default: ``None``), currently
          irrelevant (to be used for lazy ring)

        - ``ram_name`` -- a string or ``None`` (default: ``None``), controls how
          the uniformizer is printed, see :ref:`padic_printing`.

        - ``res_name`` -- a string or ``None`` (default: ``None``), the name of the
          reduction of the generator; if ``None``, ``'0'`` is appended to the name
          of the generator.

        - ``print_pos`` -- a boolean (default: ``False``), whether to only use
          positive integers in the representations of elements, see
          :ref:`padic_printing`.

        - ``print_sep`` -- a string or ``None`` (default: ``None``), the separator
          character used when ``print_mode`` is ``'bars'``, see
          :ref:`padic_printing`

        - ``print_max_ram_terms`` -- an integer or ``None`` (default: ``None``),
          the maximum number of powers of the uniformizer to print, see
          :ref:`padic_printing`

        - ``print_max_unram_terms`` -- an integer or ``None`` (default: ``None``),
          the maximum number of entries shown in a coefficient of `p`, see
          :ref:`padic_printing`

        - ``print_max_terse_terms`` -- an integer or ``None`` (default: ``None``),
          the maximum number of terms to print in a polynomial representation of an
          element, see :ref:`padic_printing`

        - ``check`` -- a boolean (default: ``True``), whether to check that ``q``
          is the power of a prime

        EXAMPLES::

            sage: R.<u> = Zq(9); R
            Unramified Extension of 3-adic Ring with capped relative precision 20 in u defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
            sage: K.<u> = Qq(9); K
            Unramified Extension of 3-adic Field with capped relative precision 20 in u defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))

        """
        q = ZZ(q)
        F = q.factor()
        if len(F) != 1 or (check and not F[0][0].is_prime()):
            raise ValueError("q must be a prime power")
        p = F[0][0]

        if prec is None:
            prec = 20
        prec = ZZ(prec)

        if halt is None:
            halt = 40
        halt = ZZ(halt)

        base = base_factory(p=p, prec=prec, type=type, print_mode=print_mode, halt=halt, names=ram_name, print_pos=print_pos, print_sep=print_sep, print_max_ram_terms=print_max_ram_terms, check=False)

        if F[0][1] == 1:
            return base

        if names is None:
            raise ValueError("generator name must be specified")
        if not isinstance(names, tuple):
            names = (names,)

        if res_name is None:
            res_name = names[0] + '0'

        if modulus is None:
            from sage.rings.finite_rings.constructor import FiniteField as GF
            modulus = PolynomialRing(base, 'x')(GF(q, res_name).modulus().change_ring(ZZ))

        return ExtensionFactory(base=base, premodulus=modulus, prec=prec, print_mode=print_mode, halt=halt, names=names, res_name=res_name, ram_name=ram_name, print_pos=print_pos, print_sep=print_sep, print_max_ram_terms=print_max_ram_terms, print_max_unram_terms=print_max_unram_terms, print_max_terse_terms=print_max_terse_terms, check=check)

    return create_unramified_ring

Qq = create_unramified_factory(Qp)
Zq = create_unramified_factory(Zp)

def create_unramified_factory_with_type(base_factory, type):
    """
    Return a function which creates an unramified `p`-adic extension ring.

    EXAMPLES:

        sage: from sage.rings.padics.factory import create_unramified_factory_with_type
        sage: zqca = create_unramified_factory_with_type(Zq, "capped-abs")
        sage: zqca(3)
        3-adic Ring with capped absolute precision 20

    """
    def create_ring(q, prec=None, modulus=None, names=None, print_mode=None, halt=None, ram_name=None, print_pos=None, print_sep=None, print_max_ram_terms=None, print_max_unram_terms=None, print_max_terse_terms=None, check=True):
        """
        Shortcut function to create an unramified %s ring.

        See documentation for ``Zq`` and ``Qq`` for a description of the input
        parameters.

        EXAMPLES::

            sage: R.<u> = ZqCR(9, 10); R
            Unramified Extension of 3-adic Ring with capped relative precision 10 in u defined by (1 + O(3^10))*x^2 + (2 + O(3^10))*x + (2 + O(3^10))
            sage: K.<u> = QqCR(9, 10); K
            Unramified Extension of 3-adic Field with capped relative precision 10 in u defined by (1 + O(3^10))*x^2 + (2 + O(3^10))*x + (2 + O(3^10))
            sage: R.<u> = ZqCA(9, 10); R
            Unramified Extension of 3-adic Ring with capped absolute precision 10 in u defined by (1 + O(3^10))*x^2 + (2 + O(3^10))*x + (2 + O(3^10))
            sage: R.<u> = ZqFM(9, 10); R
            Unramified Extension of 3-adic Ring of fixed modulus 3^10 in u defined by (1 + O(3^10))*x^2 + (2 + O(3^10))*x + (2 + O(3^10))

        """%type
        return base_factory(q=q, prec=prec, modulus=modulus, names=names, print_mode=print_mode, halt=halt, ram_name=ram_name, print_pos=print_pos, print_sep=print_sep, print_max_ram_terms=print_max_ram_terms, print_max_unram_terms=print_max_unram_terms, print_max_terse_terms=print_max_terse_terms, check=check, type=type)
    return create_ring

ZqCR = create_unramified_factory_with_type(Zq, "capped-rel")
QqCR = create_unramified_factory_with_type(Qq, "capped-rel")
ZqCA = create_unramified_factory_with_type(Zq, "capped-abs")
ZqFM = create_unramified_factory_with_type(Zq, "fixed-mod")

###########################################
#                                         #
#    Eisenstein and general extensions    #
#                                         #
###########################################
class GenericExtensionFactory(AbstractFactory):
    """
    Create an extension of ``base``.

    See :meth:`ExtensionFactory` for a description of the parameters.

    EXAMPLES::

        sage: from sage.rings.padics.factory import QpExtensionFactory
        sage: K = Qp(2)
        sage: R.<x> = K[]
        sage: QpExtensionFactory(base=K, premodulus=x^2+x+1, unram_name="u")
        Unramified Extension of 2-adic Field with capped relative precision 20 in u defined by (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20))

    """
    @staticmethod
    def krasner_check(poly):
        """
        Check whether the Eisenstein polynomial ``poly`` uniquely defines an
        extension by Krasner's lemma.

        OUTPUT:

        Raise a ValueError if ``poly`` does not uniquely determine an extension
        or if its precision is insufficient to apply the algorithm.

        ALGORITHM:

        By Krasner's lemma, given `a,b` algebraic over `K`, we have that
        `K(a)\subseteq K(b)` if `v(a-b)>v(a_i-a)` for all conjugates `a_i` of
        `a`.

        We can therefore guarantee that ``poly`` uniquely defines an extension,
        if for all `a` and `c` whose minimal polynomials `f` and `g` reduce to
        ``poly``, there is a conjugate `b` of `c` such that `a,b` satisfy
        Krasner's lemma and vice versa.

        One can show that the distance `d(f,g)=min\{v(g_i-f_i)+1/n\}` with `n`
        the degree of `f` and `g` satisfies
        `d(f,g)=\sum_{a_i}\min\{v(a-b),v(a-a_i)\}` for `b` the conjugate of `c`
        closest to `a`.

        A direct computation shows that `v(a-b)>v(a_i-a)` for all `a_i` iff
        `d(f,g)<\max_{a_i\neq a}\{v(a-a_i)\}\sum_{a_i\neq a} v(a-a-i)`. The
        right hand side of this expression can be computed from the Newton
        polygon of `f(x+a)`.

        In other words, to decide if ``poly`` uniquely defines an extension, we
        compute these `v(a-a_i)` and we compute a bound on `d(f,g)`.

        EXAMPLES::

            sage: from sage.rings.padics.factory import QpExtensionFactory

            sage: K = Qp(2,2)
            sage: R.<x> = K[]
            sage: QpExtensionFactory.krasner_check( x^2 + 2*x + 2 )
            sage: QpExtensionFactory.krasner_check( x^2 + 2 )
            Traceback (most recent call last):
            ...
            ValueError: polynomial does probably not determine a unique totally ramified extension

            sage: K = Qp(2,3)
            sage: R.<x> = K[]
            sage: QpExtensionFactory.krasner_check( x^2 + 2*x + 2 )
            sage: QpExtensionFactory.krasner_check( x^2 + 2 )

        """
        if not GenericExtensionFactory.is_eisenstein(poly):
            raise ValueError("only implemented for Eisenstein polynomials")

        K = poly.base_ring()
        n = ZZ(poly.degree())

        # a_i-a are the roots of f(x+a)
        # the valuation of a_i-a is therefore given by the slopes of the Newton
        # polygon of f(x+a)/x which can be determined from its Taylor expansion

        # to compute the Taylor expansion, we have to work in the quotient L = K[x]/(f)
        L = K.extension(poly, names="a", check=False)
        coeffs = [(poly.derivative(i)(L.gen()),L(ZZ(i).factorial())) for i in range(1,n+1)]
        vals = [a.valuation() - b.valuation() for a,b in coeffs]

        from newton_polygon import NewtonPolygon
        NP = NewtonPolygon(vals)
        assert NP == NP.principal_part()

        for i,((c,_),v) in enumerate(zip(coeffs,vals)):
            if c.is_zero() and NP[i] == v:
                raise ValueError("insufficient precision to determine whether polynomial defines a unique extension")

        # the valuations of the a-a_i
        conjugate_vals = []
        for side,slope in zip(NP.sides(),NP.slopes()):
            conjugate_vals.extend([-slope]*(side[1][0]-side[0][0]))
        assert len(conjugate_vals) == n - 1

        # now we compute a bound on d(f,g)
        dfg = min([n*poly[i].precision_absolute() + i for i in range(n)])

        if dfg <= sum(conjugate_vals) + max(conjugate_vals):
            raise ValueError("polynomial does probably not determine a unique totally ramified extension")

    @staticmethod
    def is_eisenstein(poly):
        """
        Return whether the monic irreducible polynomial ``poly`` is Eisenstein.

        A polynomial is Eisenstein if its leading coefficient is one, the
        constant term has valuation 1 and all other coefficients have positive
        valuation.

        INPUT:

        - ``poly`` -- a polynomial over a `p`-adic ring

        EXAMPLES::

            sage: R = Zp(5)
            sage: S.<x> = R[]
            sage: f = x^4 - 75*x + 15
            sage: from sage.rings.padics.factory import ZpExtensionFactory
            sage: ZpExtensionFactory.is_eisenstein(f)
            True
            sage: g = x^4 + 75
            sage: ZpExtensionFactory.is_eisenstein(g)
            False
            sage: h = x^7 + 27*x -15
            sage: ZpExtensionFactory.is_eisenstein(h)
            False
        """
        if not poly.is_monic():
            raise ValueError("poly must be monic")
        if poly[0].valuation() != 1:
            return False
        if not all([c.valuation() for c in poly.list()[:-1]]):
            return False
        return True

    @staticmethod
    def is_trivial(poly):
        """
        Return whether the monic irreducible polynomial ``poly`` defines a
        trivial extension, i.e., whether it is linear.

        EXAMPLES::

            sage: R = Zp(5)
            sage: S.<x> = R[]
            sage: f = x^4 + 14*x + 9
            sage: from sage.rings.padics.factory import ZpExtensionFactory
            sage: ZpExtensionFactory.is_trivial(f)
            False
            sage: f = x + 1
            sage: ZpExtensionFactory.is_trivial(f)
            True

        """
        if not poly.is_monic():
            raise ValueError("poly must be monic")
        return poly.degree() == 1

    @staticmethod
    def is_totally_ramified(poly):
        """
        Return whether the monic irreducible polynomial ``poly`` defines a
        totally ramified extension.

        EXAMPLES::

            sage: R = Zp(5)
            sage: S.<x> = R[]
            sage: f = x^4 + 14*x + 9
            sage: from sage.rings.padics.factory import ZpExtensionFactory
            sage: ZpExtensionFactory.is_totally_ramified(f)
            False
            sage: f = x^2 + 5*x + 25
            sage: ZpExtensionFactory.is_totally_ramified(f)

        """
        if not poly.is_monic():
            raise ValueError("poly must be monic")
        F = poly.map_coefficients(lambda c:c.residue(), poly.base_ring().residue_class_field()).factor()
        if len(F)!=1:
            raise ValueError("poly must be irreducible")

        return F[0][0].degree()==1

    @staticmethod
    def is_unramified(poly):
        """
        Return whether the monic polynomial ``poly`` is unramified.

        A monic polynomial is unramified if its reduction modulo the maximal
        ideal is irreducible.

        EXAMPLES::

            sage: R = Zp(5)
            sage: S.<x> = R[]
            sage: f = x^4 + 14*x + 9
            sage: from sage.rings.padics.factory import ZpExtensionFactory
            sage: ZpExtensionFactory.is_unramified(f)
            True
            sage: g = x^6 + 17*x + 6
            sage: ZpExtensionFactory.is_unramified(g)
            False

        TESTS:

            sage: K = QpCR(3,20)
            sage: S.<x> = K[]
            sage: L.<x> = K.extension(209/3*x^2 + 309*x + 47/9) # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: minimal polynomial is not integral after normalization

        """
        if not poly.is_monic():
            raise ValueError("poly must be monic")
        if poly[0].valuation() > 0 or poly[poly.degree()].valuation() > 0:
            return False
        if any([c.valuation()<0 for c in list(poly)]):
            return False
        if not poly.map_coefficients(lambda c:c.residue(), poly.base_ring().residue_class_field()).is_irreducible():
            return False
        return True

    def _normalize_modulus(self, base, premodulus, check):
        from sage.symbolic.expression import is_Expression
        if is_Expression(premodulus):
            if len(premodulus.variables()) != 1:
                raise ValueError("modulus must be in one variable")
            modulus = premodulus.polynomial(base)
        elif is_Polynomial(premodulus):
            if premodulus.parent().ngens() != 1:
                raise ValueError("modulus must be univariate")
            modulus = premodulus.change_ring(base)
        else:
            raise ValueError("modulus must be a univariate polynomial or symbolic expression")

        if check:
            if not modulus.is_irreducible():
                raise ValueError("modulus must be irreducible")
            if modulus.is_constant():
                raise ValueError("modulus must not be constant")
            if len(modulus.list()) > modulus.degree()+1:
                raise ValueError("modulus must not have leading zero coefficients")

        if not modulus.is_monic():
            if base.is_field():
                # TODO: modulus /= modulus.leading_coefficient() does not work. There is a bug in (padic?) polynomials: /= will coerce the rhs into the lhs parent and then do the division which puts this into the fraction field of the polynomial ring
                modulus = modulus / modulus.leading_coefficient()
            else:
                modulus *= modulus.leading_coefficient().inverse_of_unit()
                if modulus.leading_coefficient().valuation() <= min(c.valuation() for c in modulus.list()):
                    modulus = modulus.map_coefficients(lambda c:c>>modulus.leading_coefficient().valuation())
                else:
                    raise NotImplementedError("modulus with leading coefficient with high valuation")

        assert modulus.is_monic()
        return modulus

    def create_key(self, base, premodulus, prec = None, print_mode = None, halt = None, names = None, res_name = None, unram_name = None, ram_name = None, print_pos = None, print_sep = None, print_max_ram_terms = None, print_max_unram_terms = None, print_max_terse_terms = None, check=True):
        r"""
        Create a key which can be used to cache the `p`-adic ring specified by
        the parameters.

        See the classes docstring for the meaning of the parameters.

        This method essentially normalizes the parameters and fills in default
        values.

        TESTS::

            sage: from sage.rings.padics.factory import QpExtensionFactory
            sage: K = Qp(2)
            sage: R.<x> = K[]
            sage: QpExtensionFactory.create_key(base=K, premodulus=x^2+x+1, unram_name="u")
            ('u', 2-adic Field with capped relative precision 20, (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20)), (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20)), ('u', 'u0'), 20, 40, 'series', True, None, None, -1, -1, None)

        """
        # a univariate polynomial over base, the actual modulus to use for the extension
        modulus = self._normalize_modulus(base, premodulus, check)

        # decide on the extension class: unramified, Eisenstein, or general
        ext = None
        if modulus.degree() == 1:
            ext = "?"
        elif self.is_eisenstein(modulus):
            ext = "e"
        elif self.is_unramified(modulus):
            ext = "u"
        else:
            ext = "?"

        print_mode, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms, names = self._create_key_printing(print_mode, {"print_mode":base.print_mode(), "print_pos":base._printer._pos(), "print_sep":base._printer._sep(), "print_max_ram_terms":base._printer._max_ram_terms(), "print_max_unram_terms":base._printer._max_unram_terms(), "print_max_terse_terms":base._printer._max_terse_terms()}, print_pos, print_sep, None, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms, names)

        # merge ram_name, unram_name into names
        if not isinstance(names, (list, tuple)):
            names = (names,)
        names = list(names)
        if names[0] is None:
            if ext == "?":
                pass
            elif ext == "e":
                names[0] = ram_name
            elif ext == "u":
                names[0] = unram_name
            else:
                assert(False)

        if names[0] is None:
            raise ValueError("name of the generator must be specified")

        if ext=="e" and ram_name is not None and ram_name != names[0]:
            raise ValueError("name of the generator is inconsistent with ram_name")
        if ext=="u" and unram_name is not None and unram_name != names[0]:
            raise ValueError("name of the generator is inconsistent with unram_name")

        # merge res_name into names
        if ext == "u":
            if len(names) == 1:
                if res_name is None:
                    res_name = names[0]+"0"
                names.append(res_name)
            if len(names) == 2:
                if res_name is None:
                    res_name = names[1]
                if names[1] != res_name:
                    raise ValueError("name of residue extension generator inconsistent")
            if len(names) > 2:
                raise ValueError("too many generators specified")
        else:
            if len(names) != 1:
                raise ValueError("wrong number of generators specified")

        names = tuple(names)

        if halt is None:
            halt = 40

        # set default for prec - if this is a ramified extension that is not
        # Eisenstein, then this number will be multiplied with the ramification
        # index by the implementing ring
        if prec is None:
            prec = base.precision_cap()
        prec = min(prec,base.precision_cap(),min([c.precision_absolute() for c in modulus.list()]))
        if ext == "u":
            modulus = modulus.map_coefficients(lambda c:c.add_bigoh(prec))
        elif ext == "e":
            modulus = modulus.map_coefficients(lambda c:c.add_bigoh(prec + c.valuation()))
        if ext == "e":
            prec *= modulus.degree()

        return ext, base, premodulus, modulus, names, prec, halt, print_mode, print_pos, print_sep, None, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms

    def _decode_key(self, key):
        r"""
        Split a ``key`` into the type, determining the implementation, and the
        parameters needed to construct an instance.

        INPUT:

        - ``key`` -- a tuple

        OUTPUT:

        A pair of strings determining the type of the implementation and a
        dictonary of parameters.

        EXAMPLES:

            sage: K = Qp(2)
            sage: R.<x> = K[]
            sage: from sage.rings.padics.factory import QpExtensionFactory
            sage: QpExtensionFactory._decode_key(('u', K, (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20)), (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20)), ('u', 'u0'), 20, 40, 'series', True, None, None, -1, -1, None))
            (('u', 'capped-rel'), {'halt': 40, 'prec': 20, 'print_mode': {'sep': None, 'alphabet': None, 'pos': True, 'max_unram_terms': -1, 'max_terse_terms': None, 'mode': 'series', 'max_ram_terms': -1}, 'poly': (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20)), 'prepoly': (1 + O(2^20))*x^2 + (1 + O(2^20))*x + (1 + O(2^20)), 'names': ('u', 'u0')})

        """
        ext, base, premodulus, modulus, names, prec, halt, print_mode, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms = key

        precision_type = None
        if base.is_capped_relative():
            precision_type = "capped-rel"
        if base.is_capped_absolute():
            precision_type = "capped-abs"
        if base.is_fixed_mod():
            precision_type = "fixed-mod"
        assert precision_type

        print_mode = { "mode":print_mode, "pos":print_pos, "sep":print_sep, "alphabet":print_alphabet, "max_ram_terms":print_max_ram_terms, "max_unram_terms":print_max_unram_terms, "max_terse_terms":print_max_terse_terms }
        return (ext, precision_type), {"prepoly":premodulus, "poly":modulus, "names":names, "prec":prec, "halt":halt, "print_mode":print_mode}

from padic_extension_leaves import EisensteinExtensionFieldCappedRelative, UnramifiedExtensionFieldCappedRelative, GeneralExtensionFieldCappedRelative, EisensteinExtensionRingCappedRelative, EisensteinExtensionRingCappedAbsolute, EisensteinExtensionRingFixedMod, UnramifiedExtensionRingCappedRelative, UnramifiedExtensionRingCappedAbsolute, UnramifiedExtensionRingFixedMod, GeneralExtensionRingCappedRelative, GeneralExtensionRingCappedAbsolute, GeneralExtensionRingFixedMod

QpExtensionFactory = GenericExtensionFactory("Qp_ext", { ("e", "capped-rel") : EisensteinExtensionFieldCappedRelative,
                                                         ("u", "capped-rel") : UnramifiedExtensionFieldCappedRelative,
                                                         ("?", "capped-rel") : GeneralExtensionFieldCappedRelative,
                                                       })
ZpExtensionFactory = GenericExtensionFactory("Zp_ext", { ("e", "capped-rel") : EisensteinExtensionRingCappedRelative,
                                                         ("e", "capped-abs") : EisensteinExtensionRingCappedAbsolute,
                                                         ("e", "fixed-mod") : EisensteinExtensionRingFixedMod,
                                                         ("u", "capped-rel") : UnramifiedExtensionRingCappedRelative,
                                                         ("u", "capped-abs") : UnramifiedExtensionRingCappedAbsolute,
                                                         ("u", "fixed-mod") : UnramifiedExtensionRingFixedMod,
                                                         ("?", "capped-rel") : GeneralExtensionRingCappedRelative,
                                                         ("?", "capped-abs") : GeneralExtensionRingCappedAbsolute,
                                                         ("?", "fixed-mod") : GeneralExtensionRingFixedMod,
                                                       })

QpIteratedExtensionFactory = GenericExtensionFactory("Qp_ext", { ("e", "capped-rel") : GeneralExtensionFieldCappedRelative,
                                                                 ("u", "capped-rel") : GeneralExtensionFieldCappedRelative,
                                                                 ("?", "capped-rel") : GeneralExtensionFieldCappedRelative,
                                                               })
ZpIteratedExtensionFactory = GenericExtensionFactory("Zp_ext", { ("e", "capped-rel") : GeneralExtensionRingCappedRelative,
                                                                 ("e", "capped-abs") : GeneralExtensionRingCappedAbsolute,
                                                                 ("e", "fixed-mod") : GeneralExtensionRingFixedMod,
                                                                 ("u", "capped-rel") : GeneralExtensionRingCappedRelative,
                                                                 ("u", "capped-abs") : GeneralExtensionRingCappedAbsolute,
                                                                 ("u", "fixed-mod") : GeneralExtensionRingFixedMod,
                                                                 ("?", "capped-rel") : GeneralExtensionRingCappedRelative,
                                                                 ("?", "capped-abs") : GeneralExtensionRingCappedAbsolute,
                                                                 ("?", "fixed-mod") : GeneralExtensionRingFixedMod,
                                                               })

def ExtensionFactory(base, premodulus, prec = None, print_mode = None, halt = None, names = None, res_name = None, unram_name = None, ram_name = None, print_pos = None, print_sep = None, print_max_ram_terms = None, print_max_unram_terms = None, print_max_terse_terms = None, check = True):
    """
    Create an extension of ``base``.

    INPUT:

    - ``base`` -- a `p`-adic ring

    - ``premodulus`` -- a modulus for the extension, see :ref:`padic_modulus`

    - ``prec`` -- an integer or ``None`` (default: ``None``), the precision cap
      of this ring. If ``None``, the precision cap will be set to the precision
      cap of ``base`` times the ramification index; or less if ``premodulus``
      is not given to sufficient precision.

    - ``print_mode`` -- a dictionary or one of ``'series'``, ``'val-unit'``,
      ``'terse'``, ``'digits'``, ``'bars'``, or ``None`` (default: ``None``),
      the way elements of this ring are printed, see :ref:`padic_printing`. If
      this is a dictionary, then the other printing parameters are read from
      it. If ``None``, then the print mode will be the same as in ``base``.

    - ``halt`` -- an integer or ``None`` (default: ``None``), currently
      irrelevant (to be used for lazy ring)

    - ``names`` -- a string or a tuple consisting of a single string, the name
      of the generator of this extension

    - ``res_name`` -- a string or ``None`` (default: ``None``), if this
      extension has an unramified part, the name of the generator of the
      residue field extension; per default, this is ``unram_name`` with a
      ``'0'`` appended.

    - ``unram_name`` -- a string or ``None`` (default: ``None``), if this
      extension has an unramified part, the name of the generator of that part;
      if this extension is unramified, this is the same as ``names``.

    - ``ram_name`` -- a string or ``None`` (default: ``None``), if this
      extension is ramified, the name of the generator of the totally ramified
      part; if ``premodulus`` is Eisenstein, this is the same as ``names``.

    - ``print_pos`` -- a boolean or ``None`` (default: ``None``), whether to
      only use positive integers in the representations of elements, see
      :ref:`padic_printing`. If ``None``, the value of ``base`` is used.

    - ``print_sep`` -- a string or ``None`` (default: ``None``), the separator
      character used when ``print_mode`` is ``'bars'``, see
      :ref:`padic_printing`. If ``None``, the value of ``base`` is used.

    - ``print_max_ram_terms`` -- an integer or ``None`` (default: ``None``),
      the maximum number of powers of the uniformizer to print, see
      :ref:`padic_printing`. If ``None``, the value of ``base`` is used.

    - ``print_max_unram_terms`` -- an integer or ``None`` (default: ``None``),
      the maximum number of entries shown in a coefficient of the uniformizer,
      see :ref:`padic_printing`. If ``None``, the value of ``base`` is used.

    - ``print_max_terse_terms`` -- an integer or ``None`` (default: ``None``),
      the maximum number of terms to print in a polynomial representation of an
      element, see :ref:`padic_printing`. If ``None``, the value of ``base`` is
      used.

    - ``check`` --  a boolean (default: ``True``), whether to perform expensive
      checks on the input parameters

    EXAMPLES::

        TODO
    """
    args = base, premodulus, prec, print_mode, halt, names, res_name, unram_name, ram_name, print_pos, print_sep, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms, check
    if base is base.ground_ring_of_tower():
        if base.is_field():
            return QpExtensionFactory(*args)
        else:
            return ZpExtensionFactory(*args)
    else:
        if base.is_field():
            return QpIteratedExtensionFactory(*args)
        else:
            return ZpIteratedExtensionFactory(*args)

pAdicExtension = ExtensionFactory

### TODO ### everything below this line ###
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_element import is_Polynomial
from padic_base_leaves import pAdicRingCappedRelative, \
                              pAdicRingCappedAbsolute, \
                              pAdicRingFixedMod, \
                              pAdicFieldCappedRelative
from padic_valuation import pAdicValuation

import padic_printing
padic_field_cache = {}

#######################################################################################################
#
#  The Extension Factory -- creates extensions of p-adic rings and fields
#
#######################################################################################################

# class pAdicExtension_class(UniqueFactory):
#     """
#     A class for creating extensions of `p`-adic rings and fields.
# 
#     EXAMPLES::
# 
#         sage: R = Zp(5,3)
#         sage: S.<x> = ZZ[]
#         sage: W.<w> = pAdicExtension(R, x^4-15)
#         sage: W
#         Eisenstein Extension of 5-adic Ring with capped relative precision 3 in w defined by (1 + O(5^3))*x^4 + (O(5^4))*x^3 + (O(5^4))*x^2 + (O(5^4))*x + (2*5 + 4*5^2 + 4*5^3 + O(5^4))
#         sage: W.precision_cap()
#         12
# 
# ###        TODO: Stefan's example
# ###        sage: k=Qp(2,10)
# ###        sage: R.<T>=k[]
# ###        sage: F=T^4-T^2+1
# ###        sage: k1.<t>=pAdicExtension(k,F); k1
# ###        General extension of 2-adic Field with capped relative precision 10 in t defined by (1 + O(2^10))*T^4 + (1 + 2 + 2^2 + 2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^8 + 2^9 + O(2^10))*T^2 + (1 + O(2^10))
# ###
# ###        sage: k=Qp(2,10)
# ###        sage: R.<x>=k[]
# ###        sage: f=x^4+2*x^3+2*x^2-2*x+2
# ###        sage: f.is_irreducible()
# ###        True
# ###        sage: k1.<alpha>=k.extension(f)
# ###        sage: f=f.base_extend(k1)
# ###        sage: F=f.factor()
# ###        sage: g=F[1][0]
# ###        sage: k2.<beta>=k1.extension(g)
# ###        sage: g(beta)
# ###        O(pi^115)
# 
#     """
#     def create_key_and_extra_args(self, base, premodulus, prec = None, print_mode = None, halt = None, names = None, var_name = None, res_name = None, unram_name = None, ram_name = None, print_pos = None, print_sep = None, print_alphabet = None, print_max_ram_terms = None, print_max_unram_terms = None, print_max_terse_terms = None, check = True, unram = False):
#         """
#         Creates a key from input parameters for pAdicExtension.
# 
#         See the documentation for ``Qq`` for more information.
# 
#         TESTS::
# 
#             sage: R = Zp(5,3)
#             sage: S.<x> = ZZ[]
#             sage: pAdicExtension.create_key_and_extra_args(R, x^4-15,names='w')
#             (('e', 5-adic Ring with capped relative precision 3, x^4 - 15, (1 + O(5^3))*x^4 + (O(5^4))*x^3 + (O(5^4))*x^2 + (O(5^4))*x + (2*5 + 4*5^2 + 4*5^3 + O(5^4)), ('w', None, None, 'w'), 12, None, 'series', True, '|', (), -1, -1, -1), {'shift_seed': (3 + O(5^3))})
# 
#         Check that :trac:`?` has been resolved::
# 
#             sage: R = Qp(3,10)
#             sage: S.<x> = R[]
#             sage: S.<a> = R.ext(3*x^2 - 9); S
#             Eisenstein Extension of 3-adic Field with capped relative precision 10 in a defined by (1 + O(3^10))*x^2 + (O(3^11))*x + (2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + O(3^11))
# 
#         """
#         if print_mode is None:
#             print_mode = base.print_mode()
#         if print_pos is None:
#             print_pos = base._printer._pos()
#         if print_sep is None:
#             print_sep = base._printer._sep()
#         if print_alphabet is None:
#             print_alphabet = base._printer._alphabet()
#         if print_max_ram_terms is None:
#             print_max_ram_terms = base._printer._max_ram_terms()
#         if print_max_unram_terms is None:
#             print_max_unram_terms = base._printer._max_unram_terms()
#         if print_max_terse_terms is None:
#             print_max_terse_terms = base._printer._max_terse_terms()
#         from sage.symbolic.expression import is_Expression
#         if check:
#             if is_Expression(premodulus):
#                 if len(premodulus.variables()) != 1:
#                     raise ValueError, "symbolic expression must be in only one variable"
#                 modulus = premodulus.polynomial(base)
#             elif is_Polynomial(premodulus):
#                 if premodulus.parent().ngens() != 1:
#                     raise ValueError, "must use univariate polynomial"
#                 modulus = premodulus.change_ring(base)
#                 if not modulus.is_irreducible():
#                     raise ValueError, "poly must be irreducible"
#             elif isinstance(premodulus, tuple):
#                 if len(premodulus) != 2:
#                     raise ValueError, "modulus must be a polynomial or a tuple of two polynomials"
#                 if not is_Polynomial(premodulus[0]) or not is_Polynomial(premodulus[1]):
#                     raise ValueError, "entries of modulus must be polynomials"
#                 if premodulus[0].parent().ngens() != 1 or premodulus[1].parent().ngens() != 1:
#                     raise ValueError, "must use univariate polynomials"
#                 modulus = (premodulus[0].change_ring(base), premodulus[1].change_ring(premodulus[0].parent()))
#             else:
#                 raise ValueError, "modulus must be a polynomial or a tuple of two polynomials"
# 
#             if isinstance(modulus, tuple):
#                 if any([modul.degree() <= 1 for modul in modulus]):
#                     raise NotImplementedError, "degree of moduli must be at least 2"
#             elif modulus.degree() <= 0:
#                 raise NotImplementedError, "degree of modulus must be at least 1"
#             # need to add more checking here.
#             if not isinstance(modulus,tuple) and not unram: #this is not quite the condition we want for not checking these things; deal with fixed-mod sanely
#                 if not modulus.is_monic():
#                     if base.is_field():
#                         modulus = modulus / modulus.leading_coefficient()
#                     elif modulus.leading_coefficient().valuation() <= min(c.valuation() for c in modulus.list()):
#                         modulus = modulus.parent()(modulus / modulus.leading_coefficient())
#                     else:
#                         modulus = modulus / modulus.leading_coefficient()
#                         base = base.fraction_field()
#                 #Now modulus is monic
#                 if not krasner_check(modulus, prec):
#                     raise ValueError, "polynomial does not determine a unique extension.  Please specify more precision or use parameter check=False."
#             if names is None:
#                 if var_name is not None:
#                     names = var_name
#                 else:
#                     raise ValueError, "must specify name of generator of extension"
#             if not isinstance(premodulus, tuple):
#                 if isinstance(names, tuple):
#                     names = names[0]
#                 if not isinstance(names, str):
#                     names = str(names)
#         else:
#             modulus = premodulus
#         #print type(base)
#         # We now decide on the extension class: unramified, Eisenstein, two-step or general
#         if isinstance(modulus, tuple):
#             if unram_name is None:
#                 unram_name = names[0]
#             if res_name is None:
#                 res_name = unram_name + '0'
#             if ram_name is None:
#                 ram_name = names[1]
#             names = (names, res_name, unram_name, ram_name)
#             polytype = '2'
#             halt = None
#             e = modulus[1].degree()
#             if prec is None:
#                 from sage.rings.infinity import infinity
#                 prec = infinity
#             prec = min(prec, min([c.precision_absolute() for c in modulus[0].list()] + [c.precision_absolute() for d in modulus[1].list() for c in d.list()] + [base.precision_cap()]) * e)
#             shift_seed = None
#             upoly,epoly = modulus = (truncate_to_prec( modulus[0], prec), truncate_to_prec( modulus[1], prec))
#         elif unram or (modulus.degree() > 1 and is_unramified(modulus)) and base is base.ground_ring_of_tower():
#             if unram_name is None:
#                 unram_name = names
#             if res_name is None:
#                 res_name = unram_name + '0'
#             if ram_name is None:
#                 ram_name = base._printer._uniformizer_name()
#             names = (names, res_name, unram_name, ram_name)
#             polytype = 'u'
#             #if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
#             #    halt = base.halting_paramter()
#             #elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
#             #    halt = None
#             halt = None
#             if prec is None:
#                 prec = min([c.precision_absolute() for c in modulus.list()] + [base.precision_cap()])
#             else:
#                 prec = min([c.precision_absolute() for c in modulus.list()] + [base.precision_cap()] + [prec])
#             shift_seed = None
#             modulus = truncate_to_prec(modulus, prec)
#         elif modulus.degree() > 1 and is_eisenstein(modulus) and base is base.ground_ring_of_tower():
#             unram_name = None
#             res_name = None
#             if ram_name is None:
#                 ram_name = names
#             names = (names, res_name, unram_name, ram_name)
#             polytype = 'e'
#             e = modulus.degree()
#             #if halt is None and isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
#             #    halt = base.halting_paramter() * e
#             #elif not isinstance(base.ground_ring_of_tower(), (pAdicRingLazy, pAdicFieldLazy)):
#             #    halt = None
#             halt = None
#             # The precision of an eisenstein extension is governed both by the absolute precision of the polynomial,
#             # and also by the precision of polynomial with the leading term removed (for shifting).
#             # The code below is to determine the correct prec for the extension, and possibly to obtain
#             # the information needed to shift right with full precision from the premodulus.
#             if is_Expression(premodulus):
#                 # Here we assume that the output of coeffs is sorted in increasing order by exponent:
#                 coeffs = premodulus.coeffs()
#                 preseed = premodulus / coeffs[-1][0]
#                 preseed -= preseed.variables()[0]**coeffs[-1][1]
#                 preseed /= base.prime() # here we assume that the base is unramified over Qp
#                 shift_seed = -preseed.polynomial(base)
#             else: # a polynomial
#                 if not premodulus.is_monic():
#                     preseed = premodulus / premodulus.leading_coefficient()
#                 else:
#                     preseed = premodulus
#                 preseed = preseed[:preseed.degree()]
#                 if base.is_fixed_mod():
#                     shift_seed = -preseed.change_ring(base)
#                     shift_seed = shift_seed.parent()([a >> 1 for a in shift_seed.list()])
#                 else:
#                     if base.e() == 1:
#                         try:
#                             preseed *= 1/base.prime()
#                             shift_seed = -preseed.change_ring(base)
#                         except TypeError:
#                             # give up on getting more precision
#                             shift_seed = -preseed.change_ring(base)
#                             shift_seed /= base.uniformizer()
#                     else:
#                         # give up on getting more precision
#                         shift_seed = -preseed.change_ring(base)
#                         shift_seed /= base.uniformizer()
#             if prec is None:
#                 prec = min([c.precision_absolute() for c in shift_seed.list() if not c._is_exact_zero()] + [modulus.leading_coefficient().precision_absolute()] + [base.precision_cap()]) * e
#             else:
#                 prec = min([c.precision_absolute() * e for c in shift_seed.list() if not c._is_exact_zero()] + [modulus.leading_coefficient().precision_absolute() * e] + [base.precision_cap() * e] + [prec])
#             modulus = truncate_to_prec(modulus, (prec/e).ceil() + 1)
#         else:
#             if unram_name is None:
#                 unram_name = names + '_u'
#             if res_name is None:
#                 res_name = names + '0'
#             if ram_name is None:
#                 ram_name = names + '_p'
#             if var_name is None:
#                 var_name = names
#             names = (var_name, res_name, unram_name, ram_name)
# 
#             polytype = 'p'
#         if polytype == 'u' or polytype == 'e':
#             key = (polytype, base, premodulus, modulus, names, prec, halt, print_mode, print_pos, print_sep, tuple(print_alphabet), print_max_ram_terms, print_max_unram_terms, print_max_terse_terms)
#         elif polytype == '2':
#             key = (polytype, base, premodulus, upoly, epoly, names, prec, halt, print_mode, print_pos, print_sep, tuple(print_alphabet), print_max_ram_terms, print_max_unram_terms, print_max_terse_terms)
#         else:
#             polytype = 'p'
# 
#             shift_seed = None
#             key = polytype, base, modulus, names, prec, halt, print_mode, print_pos, print_sep, tuple(print_alphabet), print_max_ram_terms, print_max_unram_terms, print_max_terse_terms
# 
#         return key, {'shift_seed': shift_seed}
# 
#     def create_object(self, version, key, shift_seed):
#         """
#         Creates an object using a given key.
# 
#         See the documentation for pAdicExtension for more information.
# 
#         TESTS::
# 
#             sage: R = Zp(5,3)
#             sage: S.<x> = R[]
#             sage: pAdicExtension.create_object(version = (3,4,2), key = ('e', R, x^4 - 15, x^4 - 15, ('w', None, None, 'w'), 12, None, 'series', True, '|', (),-1,-1,-1), shift_seed = S(3 + O(5^3)))
#             Eisenstein Extension of 5-adic Ring with capped relative precision 3 in w defined by (1 + O(5^3))*x^4 + (2*5 + 4*5^2 + 4*5^3 + O(5^4))
#         """
#         polytype = key[0]
#         if polytype == 'u' or polytype == 'e':
#             polytype, base, premodulus, modulus, names, prec, halt, print_mode, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms = key
#             return ext_table[polytype, type(base.ground_ring_of_tower()).__base__](premodulus, modulus, prec, halt, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet, 'max_ram_terms': print_max_ram_terms, 'max_unram_terms': print_max_unram_terms, 'max_terse_terms': print_max_terse_terms}, shift_seed, names)
#         elif polytype == '2':
#             polytype, base, premodulus, upoly, epoly, names, prec, halt, print_mode, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms = key
#             return ext_table['2', type(base.ground_ring_of_tower()).__base__](premodulus, upoly, epoly, prec, halt, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet, 'max_ram_terms': print_max_ram_terms, 'max_unram_terms': print_max_unram_terms, 'max_terse_terms': print_max_terse_terms}, names)
#         elif polytype == 'p':
#             polytype, base, modulus, names, prec, halt, print_mode, print_pos, print_sep, print_alphabet, print_max_ram_terms, print_max_unram_terms, print_max_terse_terms = key
# 
#             from sage.categories.homset import Hom
#             # first, we reduce to the case with base being two-step/Eisenstein/unramified/no-extension over base.ground_ring_of_tower()
#             #from sage.rings.padics.eisenstein_extension_generic import EisensteinExtensionGeneric
#             if base is base.ground_ring_of_tower() or isinstance(base, UnramifiedExtensionGeneric) or isinstance(base, EisensteinExtensionGeneric):
#                 abs_base = base
#                 abs_base_primitive = abs_base.gen()
#                 to_abs_base = Hom(base,base).identity()
#                 from_abs_base = Hom(base,base).identity()
#             elif isinstance(base, GeneralExtensionGeneric) or isinstance(base, TwoStepExtensionGeneric):
#                 if isinstance(base, GeneralExtensionGeneric):
#                     abs_base = base._implementation_ring
#                     to_abs_base = base.to_implementation_ring
#                     from_abs_base = base.from_implementation_ring
#                 else:
#                     abs_base = base
#                     to_abs_base = Hom(base,base).identity()
#                     from_abs_base = Hom(base,base).identity()
#                 if isinstance(abs_base, TwoStepExtensionGeneric):
#                     ## TODO: fix the number of iterations
#                     for i in range(10):
#                         abs_base_primitive = abs_base.gen(1) + i*abs_base.gen(0)
#                         if abs_base_primitive.valuation():
#                             continue
#                         if abs_base_primitive.matrix(base=base.ground_ring_of_tower()).minpoly().degree() == abs_base.degree():
#                             break
#                     else:
#                         raise NotImplementedError("primitive element theorem failed")
#                 else:
#                     abs_base_primitive = abs_base.gen()
#             else:
#                 raise NotImplementedError,"%s"%type(base)
# 
#             modulus_over_abs_base = modulus.map_coefficients(lambda c:to_abs_base(c), abs_base)
# 
#             # find a primitive element and its modulus for the absolute
#             # extension ring over the padic base ring (this may not yet be in
#             # Eisenstein/unramified/... form)
#             if abs_base is abs_base.ground_ring_of_tower():
#                 abs_modulus = modulus_over_abs_base
#                 to_abs_ring = [abs_modulus.parent().gen()]
#                 from_abs_ring = modulus_over_abs_base.parent().gen()
#             else:
#                 # abs_base is a direct extension of a p-adic ground ring
#                 assert abs_base.base_ring() is abs_base.ground_ring_of_tower()
# 
#                 # TODO: how many attempts do we need?
#                 for i in range(100):
#                     # a number of bad things can happen if we choose our
#                     # primitive element and its minpoly unwisely:
#                     # 1) a coefficient with low absolute precision could
#                     # eliminate most information when we reduce mod the minpoly
#                     # [can we do anything about this?]
#                     # 2) a constant with high valuation makes the polynomial
#                     # almost not irreducible; in particular, powers of the
#                     # minpoly will see additional factors. However, for a
#                     # ramified extension, the constant must have non-zero
#                     # valuation.
#                     # [can we do anything about this?]
# 
#                     # guess a primitive element
#                     primitive_element = modulus_over_abs_base.parent().gen() + i*abs_base_primitive
# 
#                     formal_quotient = modulus_over_abs_base.parent().quo(modulus_over_abs_base)
#                     primitive_element = formal_quotient(primitive_element)
# 
#                     # compute the minpoly of primitive_element
#                     # what is the minpoly over an inexact ring?
#                     # 1) if there is a lift of primitive_element which is not
#                     # primitive, then the minpoly should have smaller degree
#                     # than the extension's degree
#                     # 2) if not, then for every lift of primitive_element a
#                     # lift of the minpoly must be its minpoly
#                     # 3) every lift of the polynomial must be irreducible
# 
#                     total_degree = modulus_over_abs_base.degree()*abs_base.degree()
#                     try:
#                         minpoly = _mod_minpoly(primitive_element, total_degree, abs_base.ground_ring_of_tower())
#                     except ValueError as e:
#                         # matrix equation has no solutions
#                         continue
# 
#                     # primitive_element is primitive if its minpoly has the degree of the full extension
#                     if not minpoly.is_irreducible():
#                         continue
#                     from sage.rings.polynomial.padics.factor.factoring import pfactortree
#                     # ugly hack - is_irreducible often misses factors
#                     if len(pfactortree(minpoly)) > 1:
#                         continue
#                     assert minpoly(primitive_element).is_zero(), minpoly(primitive_element)
# 
#                     # when we used change_ring(ZZ) we made a choice. This is
#                     # the minimial polynomial of some lift of primitive_element
#                     # to infinite precision but probably primitive_element has
#                     # a minpoly of lower degree? Furthermore, this polynomial
#                     # might show some special behaviour which another minimal
#                     # polynomial might not have?
#                     # TODO: How do we make sure that we got "the right one"?
# 
#                     break
#                 else:
#                     raise NotImplementedError("primitive element theorem failed")
# 
#                 abs_modulus = minpoly
# 
#                 from_abs_ring = modulus_over_abs_base.parent().gen() + i*abs_base_primitive
#                 # write a basis of abs_ring in abs_base[T]/(modulus_over_abs_base)
#                 abs_ring_basis = [formal_quotient(from_abs_ring)**i for i in range(total_degree)]
#                 abs_ring_basis = [b.lift() for b in abs_ring_basis]
#                 abs_ring_basis = [b.list() for b in abs_ring_basis]
#                 abs_ring_basis = [b[:modulus_over_abs_base.degree()] for b in abs_ring_basis]
#                 abs_ring_basis = [b + [abs_base.zero()]*(modulus_over_abs_base.degree()-len(b)) for b in abs_ring_basis]
#                 from sage.misc.flatten import flatten
#                 abs_ring_basis = [flatten([c.vector(base=base.ground_ring_of_tower()) for c in b]) for b in abs_ring_basis]
#                 from sage.matrix.constructor import matrix
#                 A = matrix(abs_ring_basis)
# 
#                 # express the generators in the basis of abs_ring
#                 to_abs_ring = []
# 
#                 if modulus_over_abs_base.degree() > 1:
#                     b = [0]*total_degree
#                     b[abs_base.degree()] = 1
#                 else:
#                     b = (-modulus_over_abs_base[0]).vector(base=base.ground_ring_of_tower())
#                 b = A.matrix_space().row_space()(b)
#                 x = A.solve_left(b)
#                 assert x*A == b
#                 x = x.change_ring(base.ground_ring_of_tower())
#                 x = x.list()
#                 while x[-1].is_zero(): x.pop()
#                 to_abs_ring.append(minpoly.parent()(x))
# 
#                 for gen in abs_base.gens():
#                     b = gen.vector(base=base.ground_ring_of_tower()) + [0]*(total_degree - abs_base.degree())
#                     b = A.matrix_space().row_space()(b)
#                     x = A.solve_left(b)
#                     assert x*A == b
#                     x = x.list()
#                     while x[-1].is_zero(): x.pop()
#                     to_abs_ring.append(minpoly.parent()(x))
# 
#             # reduce the absolute extension to an absolute extension that we
#             # can handle (two-step, Eisenstein, unramified, or trivial)
#             if abs_modulus.degree() == 1:
#                 implementation_ring = abs_base
#                 to_implementation_ring = -abs_modulus[0]
#                 from_implementation_ring = ()
#                 prec = modulus.base_ring().precision_cap()
#             else:
#                 # TODO: prec
#                 upoly, epoly, prec = split(abs_modulus, prec)
#                 if upoly.degree() >= 2 and epoly.degree() >= 2:
#                     implementation_ring = pAdicExtension(abs_base.base_ring(), (upoly, epoly), names=('u','pi'))
#                 elif upoly.degree() >= 2:
#                     implementation_ring = pAdicExtension(abs_base.base_ring(), upoly, names=('u'))
#                 else:
#                     implementation_ring = pAdicExtension(abs_base.base_ring(), epoly, names=('pi'))
#                 prec = implementation_ring.precision_cap()
# 
#                 to_implementation_ring = implementation_ring._any_root_univariate_polynomial(abs_modulus.change_ring(implementation_ring))
# 
#                 # a base of the abs_ring in the implementation_ring
#                 from itertools import product
#                 abs_basis = [to_implementation_ring**i for i in range(abs_modulus.degree())]
#                 # as a matrix over base
#                 from sage.matrix.constructor import matrix
#                 A = matrix(implementation_ring.ground_ring_of_tower(),[b.vector(base.ground_ring_of_tower()) for b in abs_basis])
# 
#                 from_implementation_ring = []
#                 # solve for the generators of the implementation_ring
#                 u_or_pi = [0]*abs_modulus.degree()
#                 u_or_pi[1]=1
#                 b = A.matrix_space().row_space()(u_or_pi)
#                 x = A.solve_left(b)
#                 assert x*A == b
#                 x = x.list()
#                 assert any(x), "The basis %s of %s permits non non-trivial solutions -- only found %s"%(abs_basis, implementation_ring, x)
#                 while x[-1].is_zero(): x.pop()
#                 from_implementation_ring.append(abs_modulus.parent()(x).change_variable_name(names[0]))
# 
#                 if upoly.degree() > 1 and epoly.degree() > 1:
#                     pi = [0]*abs_modulus.degree()
#                     pi[upoly.degree()] = 1
#                     b = A.matrix_space().row_space()(pi)
#                     x = A.solve_left(b)
#                     assert x*A == b
#                     x = x.list()
#                     while x[-1].is_zero(): x.pop()
#                     from_implementation_ring.append(abs_modulus.parent()(x).change_variable_name(names[0]))
# 
#             to_implementation_ring_base = abs_base.hom([to(to_implementation_ring) for to in to_abs_ring[1:]]) * to_abs_base
#             to_implementation_ring = to_abs_ring[0](to_implementation_ring)
#             from_implementation_ring = [ image_of_implementation_ring_generator_in_abs_ring(from_abs_ring).map_coefficients(lambda c: from_abs_base(c), base) for image_of_implementation_ring_generator_in_abs_ring in from_implementation_ring ]
# 
#             return ext_table['p', type(base.ground_ring_of_tower()).__base__](modulus, implementation_ring, to_implementation_ring, to_implementation_ring_base, from_implementation_ring, prec, halt, {'mode': print_mode, 'pos': print_pos, 'sep': print_sep, 'alphabet': print_alphabet, 'max_ram_terms': print_max_ram_terms, 'max_unram_terms': print_max_unram_terms, 'max_terse_terms': print_max_terse_terms}, names)
# 
# ExtensionFactory = pAdicExtension = pAdicExtension_class("pAdicExtension")

######################################################
# Helper functions for the Extension Factory
######################################################

#def split(poly, prec, y=None, pi=None, unram_name='u', res_name='u0'):
#    r"""
#    Given a polynomial ``poly`` and a desired precision ``prec``, computes
#    ``upoly`` and ``epoly`` so that the extension defined by ``poly`` is
#    isomorphic to the extension defined by first taking an extension by the
#    unramified polynomial ``upoly``, and then an extension by the Eisenstein
#    polynomial ``epoly``.
#
#    INPUT:
#
#        - ``poly`` -- a polynomial defined over a p-adic base ring.
#
#        - ``prec`` -- a positive integer, the coefficients of ``poly`` will
#          only be considered to that absolute precision ---- I guess that was not the intention of the original author: this should be: the resulting polynomials should support a precision of prec
#
#    OUTPUT:
#
#    A tuple ``(upoly,epoly)`` of polynomials. ``upoly`` is a univariate
#    polynomial over the p-adic base ring. ``epoly`` is a univariate polynomial
#    over a univariate polynomial ring over the p-adic base ring.
#
#    ALGORITHM:
#
#    Let `\mathcal{o}` be the ring of integers of the p-adic base ring and let
#    `\matcal{O}` be the ring of integers of the extension defined by ``poly``.
#
#    The algorithm relies on the following observation: If `y` is a unit in
#    `\matcal{O}` and `\pi` is a non-unit in `\mathcal{O}` such that
#    `\mathrm{Frac}(\mathcal{o}[y]) = \mathrm{Frac}(\mathcal{o}[\pi]) =
#    \mathrm{Frac}(\mathcal{O})` but `\mathcal{o}[\pi,y] \neq \mathcal{O}`, then
#    `\mathcal{o}[\pi,y]` is not integrally closed, and so is not regular; in
#    particular, `(\pi)` is not its maximal ideal.
#
#    Therefore not both of the following can be true:
#
#        - `p\in(\pi)\subseteq\mathcal{o}[\pi,y]`
#
#        - `G(y)\in(\pi)\subseteq\mathcal{o}[\pi,y]` where `G(y)` is a lift of
#          the minimal polynomial of `\overline{y}`
#
#    because otherwise `\mathcal{o}[y,\pi]/(\pi)` would be a field.
#
#    The idea of this algorithm is to develop `p` and `G(y)` in `(\pi)` to find
#    a new `\pi` and `y` with `\mathcal{o}[\pi,y]=\mathcal{O}`.
#
#    One can show that if `\mathcal{o}[y,\pi]=\mathcal{O}`, then
#    the valuation of `\pi` and `p` are coprime; so as soon as one has such `y`
#    and `\pi`, they can be used to find a uniformizer `\pi'` of `\mathcal{O}`.
#
#    Then `G(y)` is a minimal polynomial of the unramified extension and the
#    minimal polynomial of `\pi'` over that unramified extension will be an
#    Eisenstein polynomial.
#
#    AUTHORS:
#
#        - Julian Rueth (2012-10-28): initial version
#
#    EXAMPLES::
#
#        sage: R.<x> = Qp(13)[]
#        sage: sage.rings.padics.factory.split(x^2 - 13^2, 20)
#        Traceback (most recent call last):
#        ...
#        ValueError: poly must be irreducible
#
#    A totally ramified extension with a polynomial which is not an Eisenstein
#    polynomial::
#
#        sage: K = Qp(3,10)
#        sage: R.<x> = K[]
#        sage: f = 8/3*x^2 + 267*x + 693
#        sage: f *= 3/8
#        sage: upoly, epoly, _ = sage.rings.padics.factory.split(f, 10)
#        sage: upoly.degree(),epoly.degree()
#        (1, 2)
#        sage: L.<a> = K.extension(epoly)
#        sage: f.change_ring(L).is_irreducible() #TODO: bug in the way we use pari's factorization
#        False
#
#    An extension with an unramified and an Eisenstein part::
#
#        sage: K = Qp(3,10)
#        sage: R.<x> = K[]
#        sage: f = 53713/3*x^4 + 905*x^3 + 116728358445*x^2 + 263907*x + 81849
#        sage: f *= 3/53713
#        sage: upoly, epoly, _ = sage.rings.padics.factory.split(f,10)
#        sage: upoly.degree(), epoly.degree()
#        (2, 2)
#
#    TESTS:
#
#    This checks that ticket #6186 is still fixed:
#
#        sage: K = Qp(13)
#        sage: R.<x> = K[]
#        sage: L.<a> = K.extension(x^2 + 1)
#        Traceback (most recent call last):
#        ...
#        ValueError: poly must be irreducible
#
#    """
#    assert poly.base_ring() is poly.base_ring().ground_ring_of_tower()
#
#    from sage.rings.polynomial.padics.factor.factoring import pfactortree
#    prec = poly.parent().base_ring().precision_cap()
#
#    if not poly.leading_coefficient().is_one():
#        raise NotImplementedError
#    if any([c.valuation() < 0 for c in poly.coeffs()]):
#        raise NotImplementedError
#
#    frame = pfactortree(poly)
#    if len(frame) != 1:
#        raise ValueError("poly must be irreducible")
#    frame = frame[0]
#
#    from sage.rings.finite_rings.constructor import FiniteField
#    F = frame.F
#    assert not frame.phi_divides_Phi(), "%s is irreducible but %s divides it"%(poly, frame.phi)
#    if hasattr(frame.polygon[0],'factors'):
#        F*=frame.polygon[0].factors[0].Fplus
#    E = poly.degree() // F
#
#    upoly = FiniteField(poly.base_ring().prime()**F,names=unram_name).polynomial().change_ring(poly.base_ring()).map_coefficients(lambda c:c.lift_to_precision(prec))
#    assert upoly.degree() == F
#    if F != 1:
#        L = poly.base_ring().extension(upoly,names=unram_name)
#    else:
#        L = poly.base_ring()
#    epoly = poly.change_ring(L)
#    assert epoly.is_squarefree()
#    from padic_valuation import pAdicValuation
#    from gauss_valuation import GaussValuation
#    epoly = pAdicValuation(L).montes_factorization(epoly)
#
#    assert epoly.prod().degree() == poly.degree(), (epoly, poly)
#    assert epoly.prod() == poly
#    assert all([g.degree() == E for g,e in epoly]), "all factors in %s should have degree %s"%(epoly,E)
#    epoly = epoly[0][0]
#    assert epoly.degree() == E, "%s should have degree %s"%(epoly,E)
#
#    v = pAdicValuation(L)
#    M = epoly.parent().quo(epoly)
#    if epoly.degree() != 1:
#        is_totally_ramified, ramification_steps = v.is_totally_ramified(epoly, include_steps=True, assume_squarefree=True)
#        assert is_totally_ramified
#        assert any([v(v.phi()).denominator() == epoly.degree() for v in ramification_steps])
#        slopes = [v(v.phi()) for v in ramification_steps]
#        keys = [v.phi() for v in ramification_steps]
#        numerators = [(slope.numerator()*epoly.degree()//slope.denominator(), key) for slope,key in zip(slopes,keys)]
#        numerators.append((-epoly.degree(),L(L.prime())**-1))
#        bfs = {0:[]}
#        while not 1 in bfs:
#            for onum in bfs.keys():
#                for num, key in numerators:
#                    nnum = onum+num
#                    if nnum not in bfs or len(bfs[onum])+1 < len(bfs[nnum]):
#                        bfs[nnum] = bfs[onum]+[key]
#        uniformizer = M.one()
#        for key in bfs[1]:
#            uniformizer*=key
#        epoly = _mod_minpoly(uniformizer, E, L)
#    else:
#        slope,epoly = ZZ(1),epoly.parent().gen() - epoly.parent().base_ring().uniformizer()
#    assert epoly.degree() == E, epoly
#    assert epoly[0].valuation() == 1, "%s is not Eisenstein"%epoly
#
#    if epoly.base_ring() is poly.base_ring():
#        return upoly, epoly, upoly.base_ring().precision_cap()*epoly.degree()
#
#    uring = upoly.parent()
#    ering = uring.change_ring(uring)
#    if poly.base_ring().is_capped_relative():
#        coeffs = [[[]]*min(c.valuation(),c.parent().precision_cap())+c.list() for c in epoly.list()]
#    else:
#        coeffs = [c.list() for c in epoly.list()]
#    coeffs = [[c + [0]*(upoly.degree()-len(c)) for c in l] for l in coeffs]
#    coeffs = [[uring(d) for d in c] for c in coeffs]
#    #coeffs = [uring(ering(c)(poly.base_ring().uniformizer())).map_coefficients(lambda d:d.add_bigoh(e.precision_absolute())) for c,e in zip(coeffs,epoly.coeffs())]
#    coeffs = [uring(ering(c)(poly.base_ring().uniformizer())) for c,e in zip(coeffs,epoly.coeffs())]
#    epoly = ering(coeffs)
#
#    if upoly.degree()*epoly.degree() != poly.degree():
#        from sage.rings.padics.precision_error import PrecisionError
#        raise PrecisionError("Insufficient precision to split poly into its unramified and Eisenstein part")
#
#    return upoly, epoly, upoly.base_ring().precision_cap()*epoly.degree() #TODO: do something about prec
#
#    from sage.rings.arith import gcd,xgcd
#
#    if poly.degree() < 2:
#        raise ValueError
#
#    #TODO: do something sensible here
#    prec = poly.parent().base_ring().precision_cap()
#
#    #if prec is None:
#    #    raise ValueError("prec must not be None")
#
#    R = poly.parent()
#    K = R.base_ring()
#    if R.base_ring() is not R.base_ring().ground_ring_of_tower():
#        raise ValueError("poly must be defined over a p-adic base field")
#    if len(R.gens()) != 1:
#        raise ValueError("poly must be defined over a univariate polynomial ring")
#    F = poly.factor()
#    if len(F) != 1 or F[0][1] != 1:
#        raise ValueError("poly must be irreducible")
#
#    n = poly.degree()
#
#    if is_unramified(poly):
#        upoly, epoly = poly, R.gen()
#    if is_eisenstein(poly):
#        upoly, epoly = R.gen(), poly
#    else:
#        # since F is irreduible, L is an integral domain
#        L = R.quo(poly)
#        p = L(R.base_ring().prime())
#
#        # elements in L might have leading zero coefficients, this function strips them away.
#        normalize = lambda f: f.parent()(f.lift().list()[:f.lift().degree()+1])
#
#        # compute the valuation of f, normalized such that p has valuation n
#        def valuation(f):
#            assert f.parent() is L
#            return normalize(f).charpoly('T')[0].valuation()
#
#        # returns the minimal positive valuation that can be obtained when multiplying powers of f and g
#        def crunch_valuation(f,g):
#            valf = valuation(f)
#            assert valf>0
#            valg = valuation(g)
#            assert valg>0
#            if valf % valg == 0: # f is a power of g (up to units)
#                assert valg <= valf
#                return valg, g
#            if valg % valf == 0: # g is a power of f (up to units)
#                assert valf <= valg
#                return valf, f
#            v,s,t = xgcd(valf,valg)
#            assert v != 0
#            assert v <= valf and v <= valg
#            return v, normalize(f**s*g**t)
#
#        # compute the "current" ramification index, i.e., the valuation of p divided by the valuation of pi
#        def e(pi):
#            val = valuation(pi)
#            assert n % val == 0
#            return n//val
#
#        # compute the minpoly of the reduction of f
#        def residue_minpoly(f):
#            F = f.charpoly(poly.variable_name()).map_coefficients(lambda c:c.residue(), K.residue_field()).factor()
#            assert len(F)==1
#            return F[0][0]
#
#        # compute the "current" degree of the residue field extension, i.e., the
#        # degree of y over the residue field of the padic base ring
#        def f(y):
#            return residue_minpoly(y).degree()
#
#        # try to find an element of smaller valuation than pi using f
#        def update_pi(f,piy):
#            pi,y = piy
#            f = normalize(f)
#            assert valuation(f)>0
#            new_val, new_pi = crunch_valuation(f,pi)
#            if new_val == valuation(pi):
#                return
#            assert new_val < valuation(pi)
#            assert y.minpoly().degree() == n
#            while new_pi.minpoly().degree() != n:
#                new_pi = normalize(new_pi*y)
#            assert new_val == valuation(new_pi)
#            piy[0] = new_pi
#            raise RestartException("found new uniformizer with valuation %s"%new_val)
#
#        # try to find an element which generates a bigger residue field extension than y using f
#        def update_y(f, piy):
#            pi,y = piy
#            f = normalize(f)
#            assert valuation(f)==0
#
#            # l is the residue field of K[y]
#            G = residue_minpoly(y)
#            assert G.base_ring().is_prime_field()
#            if G.degree() == 1:
#                l = K.residue_field()
#            else:
#                from sage.rings.finite_rings.constructor import FiniteField
#                l = FiniteField(order=G.base_ring().order()**G.degree(),name="ybar",modulus=G)
#
#            # the residue field of K[f] is K[T]/(F)
#            F = residue_minpoly(f).change_ring(l).factor()
#            for FF,d in F:
#                if FF.degree() == 1:
#                    continue
#                else:
#                    if G.degree() == 1:
#                        piy[1] = f
#                        raise RestartException("extended residue field by %s"%FF)
#                    else:
#                        raise NotImplementedError("product of finite fields")
#
#        def residue_lift(f, piy):
#            pi, y = piy
#            f = normalize(f)
#            assert valuation(f)==0
#            G = residue_minpoly(y)
#            if G.degree() == 1:
#                l = K.residue_field()
#            else:
#                from sage.rings.finite_rings.constructor import FiniteField
#                l = FiniteField(order=G.base_ring().order()**G.degree(),name="ybar",modulus=G)
#
#            F = residue_minpoly(f).change_ring(l).factor()
#            # F might factor - one of the conjugate roots is our desired lift
#            assert all([h.degree() == 1 for h,e in F])
#            for h,e in F:
#                lift = (-h[0]).polynomial().map_coefficients(lambda c:K(c).lift_to_precision(prec), K)(y)
#                assert valuation(lift)==0
#                if valuation(lift-f):
#                    return lift
#            assert False
#
#        # we need a random element pi of positive valuation which generates the field extension
#        if pi is None:
#            while True:
#                pi = L.random_element()
#                if valuation(pi)>0 and pi.minpoly().degree()==n: break
#        pi = L(pi)
#        assert pi.minpoly().degree()==n
#
#        # we need a random element y of valuation zero which generates the field extension
#        if y is None:
#            while True:
#                y = L.random_element()
#                if valuation(y) % n == 0:
#                    y /= p**(valuation(y)//n)
#                if valuation(y)==0 and y.minpoly().degree()==n: break
#        y = L(y)
#        assert y.minpoly().degree()==n
#
#        pi /= p**(((valuation(pi)-1)/valuation(p)).floor())
#        piy = [pi,y]
#
#        # now we develop p and G(y) in (pi) until we have found a uniformizer and a
#        # generator of the residue field extension we want to restart this process
#        # as soon as we have found a new candidate for a uniformizer of a generator
#        # of the residue field extension; the following exception will be used to
#        # trigger this
#
#        class RestartException(Exception):pass
#
#        # develop target in (pi)
#        def develop(target, piy):
#            #print "Now developing",target
#            assert valuation(target),"target has zero valuation"
#            pi,y = piy
#            update_pi(target, piy)
#            target /= pi**(valuation(target)/valuation(pi))
#            target = normalize(target)
#            #print "target is",target,"after normalization"
#            #print "Minpoly has degree",target.minpoly().degree(),"after normalization"
#            #this is total nonsense
#            #if target.minpoly().degree() == 1: # development was successful: target is in (pi)
#            #    return None
#            if target.is_zero():
#                return None
#            update_y(target, piy)
#            target -= residue_lift(target, piy)
#            #print "Taking a way a lift of its residue we are left with",target
#            #print "Which has positive valuation",valuation(target)
#            target = normalize(target)
#            assert valuation(target) > 0, "%s does not have positive valuation, its minpoly is %s"%(target,target.minpoly())
#            if target.minpoly().degree() == 1: # development was successful: target is in (pi)
#                return None
#            return target
#
#        while e(piy[0])*f(piy[1]) != n:
#            try:
#                targets = [normalize(residue_minpoly(y).map_coefficients(lambda c:K(c).lift_to_precision(prec),K)(y)), p]
#                while True:
#                    #print "Developing",targets
#                    #print [target.minpoly().degree() for target in targets if target is not None]
#                    assert len(targets)
#                    targets = [develop(target, piy) for target in targets if target is not None]
#            except RestartException as re:
#                pass
#
#        pi,y = piy
#        upoly = residue_minpoly(y).map_coefficients(lambda c:c.lift(),K)
#
#        L = K.extension(upoly, names=(unram_name,), res_name=res_name) if upoly.degree() > 1 else K # probably we don't need these parameters if epoly is a polynomial over a polynomial ring now TODO
#
#        pim = pi.minpoly()
#        assert pim(pi).is_zero()
#        F = pim.change_ring(L).factor()
#        epoly = F[0][0]
#
#        # TODO: assert is_unram, is_Eisenstein
#
#    assert upoly.degree()*epoly.degree() == n
#
#    upoly = upoly.map_coefficients(lambda c:c.add_bigoh(prec))
#    epoly = epoly.map_coefficients(lambda c:c.add_bigoh(prec))
#
#    # turn epoly into a polynomial over a polynomial ring over the base ring
#    #print "before",epoly
#    uring = upoly.parent()
#    ering = uring.change_ring(uring)
#    coeffs = [[[]]*c.valuation()+c.list() for c in epoly.coeffs()]
#    coeffs = [[c + [0]*(upoly.degree()-len(c)) for c in l] for l in coeffs]
#    coeffs = [[uring(d) for d in c] for c in coeffs]
#    coeffs = [ering(c)(K.uniformizer()).map_coefficients(lambda d:d.add_bigoh(e.precision_absolute())) for c,e in zip(coeffs,epoly.coeffs())]
#    epoly = ering(coeffs)
#    #print "after",epoly
#
#    if upoly.degree()*epoly.degree() != n:
#        from sage.rings.padics.precision_error import PrecisionError
#        raise PrecisionError("Insufficient precision to split poly into its unramified and Eisenstein part")
#
#    return upoly, epoly, y, pi, prec #TODO: do something about prec
#
#def truncate_to_prec(poly, absprec):
#    """
#    Truncates the unused precision off of a polynomial.
#
#    EXAMPLES::
#
#        sage: R = Zp(5)
#        sage: S.<x> = R[]
#        sage: from sage.rings.padics.factory import truncate_to_prec
#        sage: f = x^4 + (3+O(5^6))*x^3 + O(5^4)
#        sage: truncate_to_prec(f, 5)
#        (1 + O(5^5))*x^4 + (3 + O(5^5))*x^3 + (O(5^5))*x^2 + (O(5^5))*x + (O(5^4))
#    """
#    R = poly.base_ring()
#    return poly.parent()([R(a, absprec=absprec) for a in poly.list()]) # Is this quite right?  We don't want flat necessarily...
#
#def krasner_check(poly, prec):
#    """
#    Returns True iff poly determines a unique isomorphism class of
#    extensions at precision prec.
#
#    Currently just returns True (thus allowing extensions that are not
#    defined to high enough precision in order to specify them up to
#    isomorphism).  This will change in the future.
#
#    EXAMPLES::
#
#        sage: from sage.rings.padics.factory import krasner_check
#        sage: krasner_check(1,2) #this is a stupid example.
#        True
#    """
#    return True #This needs to be implemented
#
