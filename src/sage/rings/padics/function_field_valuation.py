from discrete_valuation import DiscreteValuation

from sage.rings.all import ZZ, infinity

class RationalFunctionFieldValuation(DiscreteValuation):
    def __init__(self, domain, base_valuation):
        if base_valuation.domain() is not domain._ring:
            raise ValueError

        self._base_valuation = base_valuation
        DiscreteValuation.__init__(self, domain)

    def _call_(self, x):
        return self._base_valuation(x.numerator()) - self._base_valuation(x.denominator())

    def shift(self, f, v):
        if f.parent() is not self.domain():
            raise ValueError

        if v == 0:
            return f
        elif v < 0:
            return f/self._base_valuation.element_with_valuation(-v)
        else:
            return f*self._base_valuation.element_with_valuation(v)

    def value_group(self):
        return self._base_valuation.value_group()

    def residue_field(self):
        return self._base_valuation.residue_field()

    def reduce(self, f):
        if f.parent() is not self.domain():
            raise ValueError
        if self(f) > 0:
            return self.residue_field().zero()
        if self(f) < 0:
            raise ValueError

        return self._base_valuation.reduce(f.numerator()) / self._base_valuation.reduce(f.denominator())

    def lift(self, F):
        if F.parent() is not self.residue_field():
            raise ValueError

        return self.domain()(self._base_valuation.lift(self._base_valuation.residue_ring()(F)))

class TrivialValuation(DiscreteValuation):
    def __init__(self, domain):
        DiscreteValuation.__init__(self, domain)

    def _call_(self, x):
        if x.is_zero():
            return infinity
        else:
            return ZZ.zero()

    def shift(self, x, v):
        if x.parent() is not self.domain():
            raise ValueError

        if v == 0:
            return x
        else:
            raise ValueError

    def value_group(self):
        return self._value_group(0)

    def residue_field(self):
        if self.domain().is_field():
            return self.domain()
        else:
            raise NotImplementedError("residue ring is not a field")

    def reduce(self, x):
        if x.parent() is not self.domain():
            raise ValueError

        return x

    lift = reduce
