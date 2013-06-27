from valuation import Valuation
from sage.sets.set import Set
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity

class DiscretePseudoValuation(Valuation):
    def __init__(self, domain):
        Valuation.__init__(self, domain, Set(ZZ).union(Set([infinity])))
