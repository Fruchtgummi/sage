from sage.rings.padics.factory import ZpFM
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.padics.factor.frame import Frame

def pfactortree(Phi):
    r"""
    Return a factorization of Phi

    This is accomplished by constructing a tree of Frames of approximations of factors

    INPUT::

        - ``Phi`` -- squarefree, monic padic polynomial with fixed precision coefficients

    OUTPUT::

        - list -- polynomial factors of Phi

    EXAMPLES::

    Factoring polynomials over Zp(2)[x]::

        sage: from sage.rings.polynomial.padics.factor.factoring import jorder4,pfactortree,Improve
        sage: f = ZpFM(2,24,'terse')['x']( (x^32+16)*(x^32+16+2^16*x^2)+2^34 )
        sage: len(pfactortree(f)) == 1
        True

    See the irreducibility of x^32+16 in Zp(2)[x]::

        sage: len(pfactortree(ZpFM(2)['x'](x^32+16))) == 1
        True

    Test the irreducibility of test polynomial jorder4 for Zp(3)::

        sage: len(pfactortree(jorder4(3))) == 1
        True

    AUTHORS::

        - Brian Sinclair (2012-02-22): initial version

    """
    from sage.misc.flatten import flatten

    def followsegment(next,Phi):
        """ Returns next if it corresponds to an irreducible factor of $\Phi$ and follows every branch if not """
#        print next
        if next.phi_divides_Phi():
#            print "phi_divides_Phi"
            return [next]
        # With E potentially increased, Check to see if E*F == deg(Phi) (and thus Phi is irreducible)
        if next.E * next.F * next.polygon[0].Eplus == Phi.degree():
#            print "EFEplus"
            return [next]
        # With F potentially increaded, Check to see if E*F == deg(Phi) (and thus Phi is irreducible)
        if next.E * next.F * next.polygon[0].Eplus * next.polygon[0].factors[0].Fplus == Phi.degree():
#            print "EFEplusFplus"
            return [next]
        # Check if we should begin Single Factor Lifting
        if sum([seg.length for seg in next.polygon]) == 1:
#            print "LENGTH"
            return [next]
        return [[followsegment(fact.next_frame(fact.rhoexp+1),Phi) for fact in seg.factors] for seg in next.polygon]

    if Phi.is_zero():
        raise ValueError("Factorization of zero not defined")
    if not Phi.leading_coefficient().is_one():
        raise ValueError("Phi=%s has no leading coefficient one."%Phi)

    # Construct and initialize the first frame (phi = x)
    next = Frame(Phi)
    next.first()

    # Construct the next level of the tree by following every factor of the
    # residual polynomials of every Newton polygon segment in our frame
    tree = followsegment(next, Phi)

    # tree contains the leaves of the tree of frames and each leaf corresponds
    # to an irreducible factor of Phi, so we flatten the list and start lifting
    flat = flatten(tree)

    return flat

def jorder4(p):
    r"""
    Produce a particularly complicated example of polynomials for factorization over Zp(p)

    INPUT::

        - ``p`` a prime number

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.factor.factoring import jorder4
        sage: jorder4(3)
        (1 + O(3^20))*x^24 + (24 + O(3^20))*x^23 + (276 + O(3^20))*x^22 + (2048 + O(3^20))*x^21 + (11310 + O(3^20))*x^20 + (51180 + O(3^20))*x^19 + (201652 + O(3^20))*x^18 + (709092 + O(3^20))*x^17 + (2228787 + O(3^20))*x^16 + (6232484 + O(3^20))*x^15 + (15469950 + O(3^20))*x^14 + (34143276 + O(3^20))*x^13 + (67323664 + O(3^20))*x^12 + (119268300 + O(3^20))*x^11 + (190652502 + O(3^20))*x^10 + (275456480 + O(3^20))*x^9 + (359189415 + O(3^20))*x^8 + (420635664 + O(3^20))*x^7 + (438402286 + O(3^20))*x^6 + (400618284 + O(3^20))*x^5 + (313569267 + O(3^20))*x^4 + (203945072 + O(3^20))*x^3 + (105227142 + O(3^20))*x^2 + (38341248 + O(3^20))*x + (10912597 + O(3^20))

    Input must be prime::

        sage: jorder4(4)
        Traceback (most recent call last):
        ...
        ValueError: p must be prime

    """
    K = ZpFM(p,20,print_mode='terse')
    Kx = PolynomialRing(K,names='x')
    x = Kx.gen()
    f1 = (x+1)**3+p;
    f2 = f1**2+p**2*(x+1);
    f3 = f2**2+4*p**3*f1*(x+1)**2;
    f4 = f3**2+20*p**2*f3*f2*(x+1)**2+64*p**9*f1;
    return f4;

def A(p,n,k,r,prec=20):
    """
    TESTS::

        sage: from sage.rings.polynomial.padics.factor.factoring import A
        sage: A(2,2,5,1).is_irreducible()
        True
        sage: A(3,2,5,1).is_irreducible()
        True
        sage: A(2,4,15,3).is_irreducible()
        True

    """
    if not p.is_prime():
        raise ValueError
    K = ZpFM(p,prec)
    if n <= 0:
        raise ValueError
    if k <= 0  or k >= prec:
        raise ValueError
    if r < 0 or r > (k/n).floor():
        raise ValueError
    if not n.gcd(k).is_one():
        raise ValueError

    Kx = PolynomialRing(K,names='x')
    x = Kx.gen()
    shift = K.one()
    for e in range(r):
        shift <<= 1
        shift += 1
    return (x**n + K.uniformizer_pow(k))(x + shift)

def Improve(frame,prec=None):
    r"""
    Lift a frame to a factor

    EXAMPLES::
        sage: from sage.rings.polynomial.padics.factor.factoring import Improve
        sage: from sage.rings.polynomial.padics.factor.frame import Frame
        sage: Kx.<x> = PolynomialRing(ZpFM(3,20,'terse'))
        sage: f = (x**3+x-1)*(x**2+1)
        sage: fr = Frame(f)
        sage: fr.first()
        sage: fr = fr.polygon[0].factors[0].next_frame()
        sage: fact = Improve(fr) ; fact
        (1 + O(3^20))*x + (904752403 + O(3^20))
        sage: f % fact
        0

    """
    def _reduce(poly,phi,d):
        """ returns poly mod phi and simplifies the denominator of poly """
        poly = phi.parent()(poly)._quo_rem_hensel(phi)[1]
        if d != 1:
            g = min([d] + [p.valuation() for p in poly])
            if g > 0:
                poly = poly.parent( [p >> g for p in poly] )
                d = d - g
        return poly,d
    if frame.phi_divides_Phi():
        return frame.phi

    if prec is None:
        prec = frame.Phi.base_ring().precision_cap()

    #functor, base = frame.phi.base_ring().construction()
    #functor.prec = 2*frame.O.precision_cap()
    #LiftRing = functor(base)
    LiftRing = frame.phi.base_ring()
    Lifty = PolynomialRing(LiftRing,names='y')
    Phi = Lifty(frame.Phi)
    phi = Lifty(frame.phi)
    a0,a1 = frame.elt_phi_expansion()[0:2]

    Psi = frame.find_psi(-a1.valuation())
    A0 = Psi * a0
    A1 = Psi * a1

    Psi,Psiden = Psi.polynomial(True)
    Psi = Lifty(Psi)

    C1inv = frame.polygon[0].factors[0].lift(1/(A1.residue()))
    C1inv,C1invden = C1inv.polynomial(True)
    C1inv,C1invden = _reduce(C1inv,phi,C1invden)

    A0,A0den = A0.polynomial(True)
    A0,A0den = _reduce(Lifty(A0),phi,A0den)

    C,Cden = _reduce(frame.Ox(A0*C1inv),phi,A0den+C1invden)
    phi = (phi + C)

    h = 1
    oldphi = None
    phi = phi.map_coefficients(lambda c:c.lift_to_precision(prec))
    while h < prec and phi != oldphi:
        oldphi = phi
        C1, C0 = Phi._quo_rem_hensel(phi)

        C0,C0den = _reduce((Psi*C0),phi,Psiden)
        C1,C1den = _reduce((Psi*C1),phi,Psiden)

        x1,x1den = _reduce((LiftRing(2)<<(C1den+C1invden))-C1*C1inv,phi,C1den+C1invden) ## should that really be a '2'?
        C1inv,C1invden = _reduce(C1inv*x1,phi,C1invden+x1den)

        C,Cden = _reduce((C0*C1inv),phi,C0den+C1invden)

        phi = (phi + C)
        h += 1
    return frame.Ox(phi)
