using Revise
using PuiseuxPolynomials
using Oscar

K, (t,) = puiseux_polynomial_ring(QQ, ["t"])
nu = tropical_semiring_map(K,t,max)
nu(t)
initial(t,nu)

R,(x,y) = K["x","y"]
f = x + y + t
OscarPuiseuxPolynomials.tropical_hypersurface_over_fraction_field(f,nu)
tropical_polynomial(f,nu)
tropical_hypersurface(ans)

initial(f,nu,[1,1])
has_attribute(R, :tropical_geometry_polynomial_rings_for_initial)


K, (u,v,w) = puiseux_polynomial_ring(QQ,["u","v","w"])
f = u^(1//2)*v^(2//3) + w^(1//4)
g = u^(2//3)
gh = g*h

u^(-1)

Kp, (tp1,tp2,tp3) = puiseux_polynomial_ring(QQ,["t1","t2","t3"])
a = tp2^(3//3)
b = tp1^(1//2)

Kt, (t1,t2) = puiseux_polynomial_ring(QQ, ["t1","t2"])
monomials(t1+t1^(2//4)*t2^2)
t1 * t1 - t2^2
(t1-t2)^4


puiseux_polynomial_ring_elem(Kt, (gens(Kt.underlyingPolynomialRing)[1])^2, ZZRingElem[0,0], ZZ(2))


Kt, (t1,) = puiseux_polynomial_ring(QQ, ["t1"])
valuation(t1+t1^(2//4)+t1^2) == 1//2



length(x)
x^QQ(1//2)
S = base_ring(R)
t = first(gens(S))

f = FiniteSparsePuiseuxSeries.FiniteSparsePuiseuxSeriesRingElem(R, t^2, ZZ(0), ZZ(2))
f+one(R) # Danger!!

normalize!(f)
f
x
normalize!(x)
x

hash(f)
hash(x)
hash(one(R))
hash(zero(R))

one(R)==one(R)
iszero(zero(R))


one(R)+one(R)
one(R)+zero(R)
zero(R)+zero(R)
one(R)+x

one(R)-one(R)
x - one(R)
x - zero(R)

isone(one(R))
isone(zero(R))
isone(x)

x*x
(x^2-one(R))*(x^2+one(R))

f = x^3 + x^2 - x
normalize!(f)
f

f = FiniteSparsePuiseuxSeries.FiniteSparsePuiseuxSeriesRingElem(R, t^3+t+1, ZZ(0), ZZ(2))
