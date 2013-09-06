from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_element import Polynomial_generic_dense, Polynomial
from sage.rings.infinity import infinity
from sage.libs.all import pari_gen
from sage.structure.factorization import Factorization
import sage.rings.padics.misc

class Polynomial_padic_flat(Polynomial_generic_dense):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False, absprec=None):
        """
        TESTS:

        Check that :trac:`13620` has been fixed::

            sage: K = ZpFM(3)
            sage: R.<t> = K[]
            sage: R(R.zero())
            0

        """
        if x is None:
            Polynomial_generic_dense.__init__(self, parent, x = None, is_gen = is_gen)
            return
        R = parent.base_ring()
        if sage.rings.fraction_field_element.is_FractionFieldElement(x):
            if x.denominator() != 1:
                raise TypeError, "denominator must be 1"
            else:
                x = x.numerator()
        if isinstance(x, Polynomial):
            if x.base_ring() is R:
                x = list(x.list())
            else:
                x = [R(a) for a in x.list()]
        elif isinstance(x, list):
            if check:
                x = [R(a) for a in x]
        elif isinstance(x, dict):
            if check:
                m = infinity
                zero = R(0)
                n = max(x.keys()) if len(x) else 0
                v = [zero for _ in xrange(n+1)]
                for i, z in x.iteritems():
                    v[i] = R(z)
                    m = min(m, v[i].precision_absolute())
                x = v
            else:
                m = sage.rings.padics.misc.min(a.precision_absolute() for a in x.values())
            if not absprec is None:
                m = min(m, absprec)
            Polynomial_generic_dense.__init__(self, parent, x, absprec = m)
            return
        elif isinstance(x, pari_gen):
            x = [R(w) for w in x.list()]
        else:
            x = [R(x)]
        if absprec is None:
            absprec = infinity
        m = min([a.precision_absolute() for a in x] + [absprec])
        Polynomial_generic_dense.__init__(self, parent, x, absprec = m)

    def _mul_(self, right):
        """
        Returns the product of this Polynomial_padic_flat by right.
        """
        return self._mul_generic(right)

    def _repr(self, name=None):
        r"""
        EXAMPLES:
            sage: R.<w> = PolynomialRing(Zp(5, prec=5, type = 'capped-abs', print_mode = 'val-unit'))
            sage: f = 24 + R(4/3)*w + w^4
            sage: f._repr()
            '(1 + O(5^5))*w^4 + (1043 + O(5^5))*w + (24 + O(5^5))'
            sage: f._repr(name='z')
            '(1 + O(5^5))*z^4 + (1043 + O(5^5))*z + (24 + O(5^5))'

        AUTHOR:
            -- David Roe (2007-03-03), based on Polynomial_generic_dense._repr()
        """
        s = " "
        n = m = self.degree()
        if name is None:
            name = self.parent().variable_name()
        #atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = self.list()
        for x in reversed(coeffs):
            if x != 0:
                if n != m:
                    s += " + "
                #x = repr(x)
                x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
            n -= 1
        if s==" ":
            return "0"
        return s[1:]

    def quo_rem(self,right):
        """
        Returns the quotient and remainder of division by right

        EXAMPLES:
            sage: Kx.<x> = PolynomialRing(Zp(7))
            sage: (x^3+7*x+1).quo_rem(x-7)
            ((1 + O(7^20))*x^2 + (7 + O(7^21))*x + (7 + 7^2 + O(7^21)), (O(7^20))*x^3 + (O(7^21))*x^2 + (O(7^21))*x + (1 + 7^2 + 7^3 + O(7^20)))
        """
        return self._quo_rem_hensel(right)

    def _quo_rem_naive(self, right):
        """
        Naive quotient with remainder operating on padic polynomials as their coefficient lists
        """
        if right == 0:
            raise ZeroDivisionError, "cannot divide by a polynomial indistinguishable from 0"
        P = self.parent()
        F = list(self); G = list(right);
        fdeg = self.degree()
        gdeg = right.degree()
        Q = [0 for i in range(fdeg-gdeg+1)]
        j=1
        while fdeg >= gdeg:
            a = F[-j]
            if a!=0:
                for i in range(fdeg-gdeg,fdeg+1):
                    F[i] -= a * G[i-(fdeg-gdeg)]
                Q[fdeg-gdeg] += a
            j+=1; fdeg-=1;
        return (P(Q), P(F))
