

# TODO:
# Kill all standard Specs
# shortcut to get the coefficient(s) of a weil divisor
# intersect Weil with Weil on a surface
#=
julia> in_linear_system(linsys[1], D,check=false)
false
=#


add_verbosity_scope(:ellipticK3)
set_verbosity_level(:ellipticK3, 2)
k = GF(29)


# The generic fiber of the elliptic fibration
# as an elliptic curve over k(t)
kt, t = polynomial_ring(k, :t)
kP1 = fraction_field(kt)
t = gen(kP1)
E = EllipticCurve(kP1, [3*t^8 + 10*t^7 + 6*t^6 + 17*t^5 + 25*t^4 + 4*t^3 + 23*t^2 + 9*t + 14, 5*t^12 + 25*t^11 + 2*t^10 + 28*t^9 + 28*t^8 + 19*t^7 + 3*t^6 + 17*t^5 + 19*t^4 + 12*t^3 + 25*t^2 + 12*t + 6])
# A basis for the Mordell-Weil group of E
mwl_basis = [E(collect(i)) for i in [
(12*t^4 + 21*t^3 + 5*t^2 + 12*t + 18, 23*t^5 + 7*t^4 + 22*t^3 + 13*t^2),
(12*t^4 + 20*t^3 + 22*t^2 + 27*t + 18, 15*t^4 + 12*t^3 + 12*t),
(12*t^4 + 20*t^3 + 27*t^2 + 11*t + 18, -(16*t^4 + 24*t^3 + 3*t^2 + 24*t)),
(4*t^4 + 5*t^3 + 5*t^2 + 6*t + 13,  9*t^6 + 21*t^5 + 17*t^4 + 12*t^2 + 3*t + 6)]]


gram =  [
0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 0, 0, 1, -2, 1, 0, 0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0, 0, 1, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 1, 1, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 1, 0,
1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, -2, 0, 0, 2,
1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, -2, 0, 2,
1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, -2, 2,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2, 2, -2]
NS = integer_lattice(gram=matrix(ZZ, 16, 16, gram))
B = basis_matrix(NS)
fnew = matrix(QQ, 1, 16, [6, 3, -1, -2, -1, -1, -1, -2, -2, -1, 0, -1, -1, 1, -1, 0])*B
@assert inner_product(ambient_space(NS), fnew, fnew)==0
p = matrix(QQ,1,16, [5, 2, -1, -2, -1, 0, 0, -1, -1, -1, 0, -1, -1, 1, -1, 0])
@assert inner_product(ambient_space(NS), p, p)[1,1] == -2

# The section P used for the fibration hop  Fiber_new = O + P + Vertical
l = length(mwl_basis)
P = sum(ZZ(p[1,i+16-l])*mwl_basis[i] for i in 1:l)

delta = factor(discriminant(E), kt).fac
reducible_singular_fibers = [p for p in keys(delta) if delta[p]>1]



IP1 = projective_space(k, 1)
c = standard_covering(IP1)
# rename the variables on the affine charts
# to a more readable version
OO(c[1]).S = [:t]
OO(c[2]).S = [:s]


O0 = twisting_sheaf(IP1, 0)
O4 = twisting_sheaf(IP1, -4)
O6 = twisting_sheaf(IP1, -6)

bundleE = direct_sum([O0, O4, O6])

X_proj = projectivization(bundleE, var_names=["z", "x", "y"])

# remove some charts not needed to cover S
X = covered_scheme(X_proj)
OO(X[1][1]).S=[:x,:y,:t]
C = Covering([X[1][i] for i in [1,3,4,6]])
for U in patches(C)
  for V in patches(C)
    add_glueing!(C, glueings(X[1])[(U, V)])
  end
end
X = CoveredScheme(C)



# Create the singular Weierstrass model S of the elliptic K3 surface
a = a_invars(E)
U = affine_charts(X)[1]  # the standard Weierstrass chart
(x, y, t) = gens(OO(U))
@assert all(denominator(i)==1 for i in a)
a = [numerator(a)(t) for a in a]
(a1,a2,a3,a4,a6) = a
ft = y^2  + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x+a6)
I = IdealSheaf(X, U, [ft])

inc_S = Oscar.CoveredClosedEmbedding(X, I)
@test I === image_ideal(inc_S)
S = domain(inc_S)  # The ADE singular elliptic K3 surface

I_sing = Oscar.ideal_sheaf_of_singular_locus(S)

I_sing_X = radical(pushforward(inc_S)(I_sing))


# Refine the covering over the reducible singular fibers
refined_charts = AbsSpec[]
U = X[1][1]
Ising = I_sing_X(U)
if isone(Ising)
  push!(refined_charts, U)
else
  # reducible singular fibers
  disc = gens(eliminate(Ising, coordinates(U)[1:2]))[1]
  redfib = collect(keys(factor(disc).fac))
  for i in 1:length(redfib)
    r = copy(redfib)
    deleteat!(r, i)
    push!(refined_charts, PrincipalOpenSubset(U, r))
  end
end
# add the fiber at s=0 and remove all other singular fibers
V = X[1][3]
IsingV = I_sing_X(V)
if isone(IsingV)
  push!(refined_charts, V)
else
  # reducible singular fibers
  local disc = gens(eliminate(IsingV, coordinates(V)[1:2]))[1]
  (x,y,s) = coordinates(V)
  b, d = divides(disc, s)
  if b
    disc = d
  end
  redfib = collect(keys(factor(disc).fac))
  push!(refined_charts, PrincipalOpenSubset(V, redfib))
end
# no extra singularities at Z = 0
# therefore we just exclude the singularities visible here
for U in [X[1][2],X[1][4]]
  local Ising = I_sing_X(U)
  if isone(Ising)
    push!(refined_charts, U)
    continue
  end
  local (z,x,s_or_t) = coordinates(U)
  # reducible singular fibers
  local disc = gens(eliminate(Ising, [z, s_or_t]))[1]
  local redfib = collect(keys(factor(disc).fac))
  push!(refined_charts, PrincipalOpenSubset(U, redfib))
end

Cref = Covering(refined_charts)
Oscar.inherit_glueings!(Cref, C)
push!(X.coverings, Cref)

@test scheme(I_sing) === S
@test scheme(I_sing_X) === X

# The mwl_basis as ideal sheaves
UX = X[1][1]
function section_to_ideal_sheaf(X, section, U)
  (x,y,t) = coordinates(U)
  b = section
  return ideal_sheaf(X,X[1][1],[OO(UX)(i) for i in [x*denominator(b[1])(t)-numerator(b[1])(t),y*denominator(b[2])(t)-numerator(b[2])(t)]])
end
sections = [section_to_ideal_sheaf(X, b, UX) for b in mwl_basis]
PonX = section_to_ideal_sheaf(X, P, UX)

#X = CoveredScheme(C)
(x,y,t) = coordinates(UX)
fiber = IdealSheaf(X, UX, [t-1,ft]) # we want a smooth fiber
(z,x,t) = coordinates(X[1][2])
zero_section = IdealSheaf(X, X[1][2], [x,z])

# components of the singular fibers visible on S
# A3 at t=0
# A4 at s=0
# A1 at t=
# A1 at t=
# A1 at t=
(x,y,t) = coordinates(UX)
A3_0 = IdealSheaf(X, UX, [t, ft])
A1a_0 = IdealSheaf(X, UX, [t+17, ft])
A1b_0 = IdealSheaf(X, UX, [t+1, ft])
A1c_0 = IdealSheaf(X, UX, [t+20, ft])
(x,y,s) = ambient_coordinates(X[1][3])
A4_0 = IdealSheaf(X,X[1][3], OO(X[1][3]).([s,gens((modulus(OO(S[1][3]))))[1]]))


# S = K3 --> X ambient space
# all relevant divisors as ideal sheaves on X
divisors = vcat([PonX, zero_section, fiber, A3_0, A4_0, A1a_0, A1b_0,A1c_0], sections)

# initialization for the while loop
X0 = X
Y0 = S
inc_Y0 = inc_S

# blow up points (one at a time) until smooth
exceptionals = []
divisors0 = divisors
count = 0
varnames = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j,:k,:l]
projectionsX = []
projectionsY = []
while true
  global count = count+1
  @vprint :ellipticK3 1 "blowup number: $(count)\n"
  @vprint :ellipticK3 2 "computing singular locus\n"
  I_sing_Y0 = Oscar.ideal_sheaf_of_singular_locus(Y0)
  @vprint :ellipticK3 2 "decomposing singular locus\n"
  I_sing_Y0 = Oscar.maximal_associated_points(I_sing_Y0)
  @vprint :ellipticK3 1 "number of singular points: $(length(I_sing_Y0))\n"
  if length(I_sing_Y0)==0
    # stop if smooth
    break
  end
  # take the first singular point and blow it up
  I_sing_X0_1 = radical(pushforward(inc_Y0)(I_sing_Y0[1]))
  if count == 1
    cov = X0[2]
  else
    cov = Oscar.simplified_covering(X0)
  end
  prX1 = blow_up(I_sing_X0_1, covering=cov, var_name=varnames[count])
  X1 = domain(prX1)
  @vprint :ellipticK3 1 "$(X1)\n"
  E1 = exceptional_divisor(prX1)

  @vprint :ellipticK3 2 "computing strict transforms\n"
  # compute the exceptional divisors
  global exceptionals = [strict_transform(prX1, e) for e in exceptionals]
  # move the divisors coming originally from S up to the next chart
  global divisors0 = [strict_transform(prX1, e) for e in divisors0]
  push!(exceptionals, E1)

  Y1, inc_Y1, pr_Y1 = strict_transform(prX1, inc_Y0)

  push!(projectionsX, prX1)
  push!(projectionsY, pr_Y1)
  simplify!(Y1)

  # set up for the next iteration
  global Y0 = Y1
  global inc_Y0 = inc_Y1
  global X0 = X1
end
# Restrict the exceptional divisors:
# X > S <- Y0 < X0
@vprint :ellipticK3 2 "Pulling back divisors to Y0\n"
# Pulling back the divisors to Y0 may introduce multiplicities
exceptionals_res = [pullback(inc_Y0)(e) for e in exceptionals]
divisors_res = [pullback(inc_Y0)(e) for e in divisors0]

# Compute the intersection matrix of the exceptional divisors:
Ex = exceptionals_res
@vprint :ellipticK3 2 "Exceptional Cartier to Weil divisors\n"
ExWeil = weil_divisor.(Ex)   # too slow to be a test
@vprint :ellipticK3 2 "done\n"
tmp = []
ExWeil = reduce(append!, [components(i) for i in ExWeil], init= tmp)


# divisors = vcat([PonX, zero_section, fiber, A3_0, A4_0, A1a_0, A1b_0,A1c_0], sections)
NSgens = vcat(divisors_res[3], divisors_res[2], ExWeil, divisors_res[end-3:end])
NSgens = [WeilDivisor(i,ZZ, check=false) for i in NSgens]
G = zero_matrix(ZZ,16,16)
for i in 1:length(NSgens)
  for j in 1:i-1
    G[i,j] = intersect(NSgens[i],NSgens[j])
    G[j,i] = G[i,j]
  end
end
for i in 2:length(NSgens)
  G[i,i]= -2
end
@assert det(G) == -1183
@assert G == gram_matrix(NS) # apparently I have some mixups
error("")

# f: Y -- > S the total blowup
f = projectionsY[1]
for g in projectionsY[2:end]
  global f = g*f
end
(x, y, t) = coordinates(S[1][1])
fstar = pullback(f[Y0[1][1]])
# |D| = | P + O + 1 F|
linsys = prop217(Y0, E, P, 1, fstar(x), fstar(y), fstar(t), fstar(t))
kY0 = parent(linsys[1])
linsys = [kY0(fstar(t),fstar(t)^0)*i for i in linsys]
PonY, OonY = [WeilDivisor(i, ZZ,check=false) for i in divisors_res[1:2]]

error("")
A1_0onY = WeilDivisor(divisors_res[4], ZZ, check=false)
A1_1onY = WeilDivisor(ExWeil[9], ZZ, check=false)
D = PonY + OonY + A1_0onY + A1_1onY
L = linear_system(linsys, D, check=false)
@test_broken in_linear_system(gens(L)[1], D) #something is wrong
Lnew, _ = subsystem(L, A1_1onY, 1)
@assert length(gens(Lnew))==2

[order_on_divisor(g, A1_1onY) for g in gens(Lnew)]
Fnew = PonY + OonY + A1_0onY
tt = gens(Lnew)[1]//gens(Lnew)[2] # the new elliptic coordinate

(x,y,t) = coordinates(S[1][1])
Ut = hypersurface_complement(S[1][1],t)
Vt = hypersurface_complement(Y0[1][1],fstar(t))
fres = restrict(f[Y0[1][1]], Vt, Ut)
fresinv = inverse(fres)
fresinvstar = pullback(fresinv)

function map_func(phi, tt)
  nu = fresinvstar(numerator(tt))
  @assert lifted_denominator(nu)==1
  nu = lifted_numerator(nu)
  du = fresinvstar(denominator(tt))
  @assert lifted_denominator(du)==1
  du = lifted_numerator(du)
  return nu//du
end

u = map_func(fresinvstar, tt) # the new elliptic parameter

# transform to new coordinates

R = OO(X[1][1])
(x,y,t) = gens(R)
# u = (a y + b) / (cy + d)
# solve for y:
# y = (b-d*u)/ (c*u - a)
g = hom(R,R,R.([x,0,t]))
b = g(nu)
a = divexact(nu - b, y)
d = g(du)
c = divexact(du - d, y)
@assert (a*y+b) // (c*y+d) == u

Ru, (t, x, u) = polynomial_ring(base_ring(R), [:t, :x, :u])
h = hom(R, Ru, [x, 0, t])
(a,b,c,d) = h.([a,b,c,d])
Ruloc,_ = localization(Ru, c*u - a)
phi = hom(R, Ruloc, Ruloc.([x, Ruloc(b-d*u, c*u - a), t]))
g = SpecMor(Spec(Ruloc), Spec(R), phi)
U1 = closure(preimage(g,S[1][1]))

f = gens(modulus(OO(S[1][1])))[1]
fu = gens(modulus(OO(U1)))[1]


Rufrac = fraction_field(Ru)
helper = hom(Ru,Rufrac, Rufrac.([t,P[1](t),u]))
