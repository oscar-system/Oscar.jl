################################################################################
#
# Some linear systems on elliptic surfaces
#
################################################################################

@doc raw"""
    _prop217(E::EllipticCurve, P::EllipticCurvePoint, k)

Compute a basis for the linear system
``|O + P + kF|``
on the  minimal elliptic (K3) surface defined by E.
Here F is the class of a fiber O the zero section
and P any non-torsion section.

The return value is a list of pairs ``(a(t),b(t))``

```jldoctest
julia> kt,t = polynomial_ring(GF(29),:t);

julia> ktfield = fraction_field(kt);

julia> bk = [((17*t^4 + 23*t^3 + 18*t^2 + 2*t + 6, 8*t^5 + 2*t^4 + 6*t^3 + 25*t^2 + 24*t + 5 )),
             ((17*t^6 + 3*t^5 + 16*t^4 + 4*t^3 + 13*t^2 + 6*t + 5)//(t^2 + 12*t + 7), (4*t^8 + 19*t^7 + 14*t^6 + 18*t^5 + 27*t^4 + 13*t^3 + 9*t^2 + 14*t + 12)//(t^3 + 18*t^2 + 21*t + 13) ),
             ((17*t^6 + 10*t^5 + 24*t^4 + 15*t^3 + 22*t^2 + 27*t + 5)//(t^2 + 16*t + 6), (20*t^8 + 24*t^7 + 22*t^6 + 12*t^5 + 21*t^4 + 21*t^3 + 9*t^2 + 21*t + 12)//(t^3 + 24*t^2 + 18*t + 19) ),
             ((17*t^8 + 21*t^7 + 20*t^5 + 24*t^4 + 21*t^3 + 4*t^2 + 9*t + 13)//(t^4 + 17*t^3 + 12*t^2 + 28*t + 28), (23*t^11 + 25*t^10 + 8*t^9 + 7*t^8 + 28*t^7 + 16*t^6 + 7*t^5 + 23*t^4 + 9*t^3 + 27*t^2 + 13*t + 13)//(t^6 + 11*t^5 + 14*t^4 + 13*t^3 + 6*t^2 + 18*t + 12) )];

julia> E = elliptic_curve(ktfield,[3*t^8+24*t^7+22*t^6+15*t^5+28*t^4+20*t^3+16*t^2+26*t+16, 24*t^12+27*t^11+28*t^10+8*t^9+6*t^8+16*t^7+2*t^6+10*t^5+3*t^4+22*t^3+27*t^2+10*t+3]);

julia> bk = [E(collect(i)) for i in bk];

julia> Oscar._prop217(E,bk[2],2)
5-element Vector{Tuple{FqPolyRingElem, FqPolyRingElem}}:
 (t^2 + 12*t + 7, 0)
 (t^3 + 8*t + 3, 0)
 (t^4 + 23*t + 2, 0)
 (25*t + 22, 1)
 (12*t + 28, t)

julia> Oscar._prop217(E,bk[1],1)
2-element Vector{Tuple{FqPolyRingElem, FqPolyRingElem}}:
 (1, 0)
 (t, 0)
```
"""
function _prop217(E::EllipticCurve, P::EllipticCurvePoint, k)
  @req !iszero(P[3]) "P must not be torsion" # seems like we cannot check this
  xn = numerator(P[1])
  xd = denominator(P[1])
  yn = numerator(P[2])
  yd = denominator(P[2])
  OP = divexact(max(degree(xd), degree(xn) - 4), 2)
  dega = k + 2*OP
  degb = k + 2*OP - 2 - divexact(degree(xd), 2) #?
  base = base_field(E)
  Bt = base_ring(base)
  B = coefficient_ring(Bt)

  R,ab = polynomial_ring(base,vcat([Symbol(:a,i) for i in 0:dega],[Symbol(:b,i) for i in 0:degb]),cached=false)
  Rt, t1 = polynomial_ring(R,:t)
  a = reduce(+,(ab[i+1]*t1^i for i in 0:dega), init=zero(Rt))
  b = reduce(+,(ab[2+dega+j]*t1^j for j in 0:degb), init=zero(Rt))
  c = a*xn(t1) - b*yn(t1)
  r = mod(c, xd(t1))
  # setup the linear equations for coefficients of r to vanish
  # and for the degree of c to be bounded above by
  # k + 2*OP + 4 + degree(xd)
  eq1 = collect(coefficients(r))
  eq2 = [coeff(c,i) for i in (k + 2*OP + 4 + degree(xd) + 1):degree(c)]
  eqns = vcat(eq1, eq2)

  # collect the equations as a matrix
  cc = [[coeff(j, abi) for abi in ab] for j in eqns]
  mat = reduce(vcat,cc, init=elem_type(base)[])
  @assert all(is_one(denominator(x)) for x in mat)
  @assert all(is_constant(numerator(x)) for x in mat)
  mat2 = [constant_coefficient(numerator(x)) for x in mat]
  M = matrix(B, length(eqns), length(ab), mat2)
  # @assert M == matrix(base, cc) # does not work if length(eqns)==0
  K = kernel(M; side = :right)
  kerdim = ncols(K)
  result = Tuple{elem_type(Bt),elem_type(Bt)}[]
  t = gen(Bt)
  for j in 1:kerdim
    aa = reduce(+, (K[i+1,j]*t^i for i in 0:dega), init=zero(Bt))
    bb = reduce(+, (K[dega+i+2,j]*t^i for i in 0:degb), init=zero(Bt))
    push!(result, (aa, bb))
  end
  # confirm the computation
  @assert kerdim == 2*k + OP # prediced by Riemann-Roch
  for (a,b) in result
    @assert mod(a*xn - b*yn, xd) == 0
    @assert degree(a) <= k + 2*OP
    @assert degree(b) <= k + 2*OP - 2 - 1//2*degree(xd)
    @assert degree(a*xn - b*yn) <= k + 2*OP + 4 + degree(xd)
  end
  return result
end

@doc raw"""
    linear_system(X::EllipticSurface, P::EllipticCurvePoint, k::Int) -> LinearSystem

Compute the linear system ``|O + P + k F|`` on the elliptic surface ``X``.
Here ``F`` is the class of the fiber over ``[0:1]``, ``O`` the zero section
and ``P`` any section given as a point on the generic fiber.

The linear system is represented in terms of the Weierstrass coordinates.
"""
function linear_system(X::EllipticSurface, P::EllipticCurvePoint, k::Int)
  euler_characteristic(X) == 2 || error("linear system implemented only for elliptic K3s")
  #FS = function_field(weierstrass_model(X)[1])
  FS = function_field(X)
  U = weierstrass_chart_on_minimal_model(X)
  (x,y,t) = ambient_coordinates(U)

  sections = elem_type(FS)[]
  if iszero(P[3])
    append!(sections, [FS(t)^(i-k) for i in 0:k])
    append!(sections, [FS(t)^(i-k)*FS(x) for i in 0:k-4])
  else
    xn = numerator(P[1])
    xd = denominator(P[1])
    yn = numerator(P[2])
    yd = denominator(P[2])

    I = saturated_ideal(defining_ideal(U))
    IP = ideal([x*xd(t)-xn(t),y*yd(t)-yn(t)])
    @hassert :EllipticSurface 2 issubset(I, IP) || error("P does not define a point on the Weierstrasschart")

    @assert gcd(xn, xd)==1
    @assert gcd(yn, yd)==1
    ab = _prop217(generic_fiber(X), P, k)
    d = divexact(yd, xd)(t)
    den = t^k*(x*xd(t) - xn(t))
    for (a,b) in ab
      c = divexact(b*yn - a*xn, xd)
      num = a(t)*x+b(t)*d*y + c(t)
      push!(sections, FS(num//den))
    end
  end
  return sections
end

@doc raw"""
    two_neighbor_step(X::EllipticSurface, F1::Vector{QQFieldElem})

Let ``F`` be the class of a fiber of the elliptic fibration on ``X``. 
Given an isotropic nef divisor ``F_1`` with ``F_1.F = 2``,
compute the linear system ``|F_1|`` and return the corresponding generic fiber
as a double cover `C` of the projective line branched over four points.

Input:
``F_1`` is represented as a vector in the `algebraic_lattice(X)`
``X`` must be a K3 surface

Output:
A tuple `(C, (x1, y1, t1))` defined as follows.
- `C` is given by a polynomial ``y_1^2 - q(x_1)`` in ``k(t_1)[x_1,y_1]`` with ``q`` of degree ``3`` or ``4``.
- (x1,y1,t1) are expressed as rational functions in terms of the weierstrass coordinates `(x,y,t)`.
"""
function two_neighbor_step(X::EllipticSurface, F::Vector{QQFieldElem})
  E = generic_fiber(X)
  basisNS, tors, NS = algebraic_lattice(X)
  V = ambient_space(NS)
  @req inner_product(V, F, F)==0 "not an isotropic divisor"
  @req euler_characteristic(X) == 2 "not a K3 surface"
  F0 = zeros(QQ,degree(NS)); F0[1]=1

  @req inner_product(V, F, F0) == 2 "not a 2-neighbor"

  D1, D, P, l, c = horizontal_decomposition(X, F)
  u = _elliptic_parameter(X, D1, D, P, l, c)
  @assert scheme(parent(u)) === X
  pr = weierstrass_contraction(X)
  WX, _ = weierstrass_model(X)
  # The following is a cheating version of the command u = pushforward(pr)(u) (the latter has now been deprecated!)
  u = function_field(WX)(u[weierstrass_chart_on_minimal_model(X)])
  @assert scheme(parent(u)) === weierstrass_model(X)[1]

  # Helper function
  my_const(u::MPolyRingElem) = is_zero(u) ? zero(coefficient_ring(parent(u))) : first(coefficients(u))

  # transform to a quartic y'^2 = q(x)
  if iszero(P[3])  #  P = O
    eqn1, phi1 = _elliptic_parameter_conversion(X, u, case=:case1)
    eqn2, phi2 = _normalize_hyperelliptic_curve(eqn1)
#   function phi_func(x)
#     y = phi1(x)
#     n = numerator(y)
#     d = denominator(y)
#     return phi2(n)//phi2(d)
#   end
#   phi = MapFromFunc(domain(phi1), codomain(phi2), phi_func)
#   # TODO: Verify that the construction below also works and replace by that, eventually.
#   phi_alt = compose(phi1, extend_domain_to_fraction_field(phi2))
#   @assert phi.(gens(domain(phi))) == phi_alt.(gens(domain(phi)))
   phi = compose(phi1, extend_domain_to_fraction_field(phi2))
  elseif iszero(2*P) # P is a 2-torsion section
    eqn1, phi1 = _elliptic_parameter_conversion(X, u, case=:case3)
    #eqn1, phi1 = _conversion_case_3(X, u)
    (x2, y2) = gens(parent(eqn1))

    # Make sure the coefficient of y² is one (or a square) so that 
    # completing the square works. 
    c = my_const(coeff(eqn1, [x2, y2], [0, 2]))::AbstractAlgebra.Generic.FracFieldElem
    eqn1 = inv(unit(factor(c)))*eqn1

    eqn2, phi2 = _normalize_hyperelliptic_curve(eqn1)
    phi = compose(phi1, extend_domain_to_fraction_field(phi2))
  else  # P has infinite order
    eqn1, phi1 = _elliptic_parameter_conversion(X, u, case=:case2)
    #eqn1, phi1 = _conversion_case_2(X, u)
    (x2, y2) = gens(parent(eqn1))
    
    # Make sure the coefficient of y² is one (or a square) so that 
    # completing the square works. 
    c = my_const(coeff(eqn1, [x2, y2], [0, 2]))::AbstractAlgebra.Generic.FracFieldElem
    eqn1 = inv(unit(factor(c)))*eqn1

    eqn2, phi2 = _normalize_hyperelliptic_curve(eqn1)
    phi = compose(phi1, extend_domain_to_fraction_field(phi2))
  end

  return eqn2, phi
end

@doc raw"""
    horizontal_decomposition(X::EllipticSurface, L::Vector{QQFieldElem}) -> AbsWeilDivisor, EllipticCurvePoint

Given a divisor ``L`` as a vector in the `algebraic_lattice(X)`
find a linearly equivalent divisor ``(n-1) O + P + V = D ~ L`` where
``O`` is the zero section, ``P`` is any section and ``V`` is vertical.

Return a tuple `(D1, D, P, l, c)` where `D` and `P` are as above and
``D <= D1 = (n-1)O + P + n_1F_1 + ... n_k F_k`` with ``l = n_1 + ... n_k`` minimal
and the `F_i` are some other fibers.
The rational function `c=c(t)` has divisor of zeros and poles``
(c) = lF - n_0F_1 + ... n_k F_k``
"""
function horizontal_decomposition(X::EllipticSurface, F::Vector{QQFieldElem})
  E = generic_fiber(X)
  basisNS, tors, NS = algebraic_lattice(X)
  V = ambient_space(NS)
  @req F in algebraic_lattice(X)[3] "not in the algebraic lattice"
  @req inner_product(V, F, F)==0 "not an isotropic divisor"
  @req euler_characteristic(X) == 2 "not a K3 surface"
  # how to give an ample divisor automagically in general?
  # @req is_nef(X, F) "F is not nef"
  l = F[1]
  rk_triv = nrows(trivial_lattice(X)[2])
  n = rank(NS)
  @assert degree(NS) == rank(NS)
  p, P = _vertical_part(X,F)
  D = section(X, P)
  F2 = F - p
  @vprint :EllipticSurface 4 "F2 = $(F2)\n"
  D = D + ZZ(F2[2])*zero_section(X)
  D1 = D
  F2 = ZZ.(F2); F2[2] = 0
  l = F2[1] # number of fibers that we need
  # find the fiber components meeting O necessary
  F3 = F2
  (_,_,t) = ambient_coordinates(weierstrass_chart_on_minimal_model(X))
  c = t^0
  for (pt, rt, fiber, comp, gram) in reducible_fibers(X)
    Fib0 = comp[1]
    f0 = zeros(QQFieldElem, length(basisNS))
    for i in 1:length(basisNS)
      if !isone(components(Fib0)[1]+components(basisNS[i])[1])
        if length(comp)==2 && 2<i<=rk_triv
          f0[i] = 2
        else
          f0[i] = 1
        end
      end
    end
    f0 = ZZ.(f0 * inv(gram_matrix(ambient_space(NS))))
    @assert inner_product(ambient_space(NS), f0,f0) == -2
    nonzero = [i for i in 3:rk_triv if f0[i]!=0]
    if pt[2]==0 # at infinity
      t0 = t
    else
      t0 = t//(t*pt[2]-pt[1])
    end
    while any(F3[i]<0 for i in nonzero)
      F3 = F3 - f0
      D = D + Fib0
      D1 = D1 + fiber
      c = c*t0
    end
  end
  pt, _ = fiber(X)
  if pt[2]==0 # at infinity
    t0 = t
  else
    t0 = t//(t*pt[2]-pt[1])
  end
  c = c*(t0//1)^ZZ(F3[1])
  D = D + F3[1]*basisNS[1]
  D1 = D1 + F3[1]*basisNS[1]
  F4 = copy(F3); F4[1]=0
  @assert all(F4[i]>=0 for i in 1:length(basisNS))
  D = D + sum(ZZ(F4[i])*basisNS[i] for i in 1:length(basisNS))
  @assert D<=D1
  return D1, D, P, Int(l), c
end
  
# internal method used for two neighbor steps
# in the horizontal_decomposition
function _vertical_part(X::EllipticSurface, v::QQMatrix)
  @req nrows(v)==1 "not a row vector"
  _,tors, NS = algebraic_lattice(X)
  E = generic_fiber(X)
  @req ncols(v)==degree(NS) "vector of wrong size $(ncols(v))"
  @req v in NS "not an element of the lattice"
  mwl_rank = length(X.MWL)
  rk_triv = rank(NS)-mwl_rank
  n = rank(NS)
  P = sum([ZZ(v[1,i])*X.MWL[i-rk_triv] for i in (rk_triv+1):n], init = E([0,1,0]))
  p = zero_matrix(QQ, 1, rank(NS)) # the section part
  p[1,end-mwl_rank+1:end] = v[1,end-mwl_rank+1:end]
  p[1,2] = 1 - sum(p)  # assert p.F = 1 by adding a multiple of the zero section O
  
  # P meets exactly one fiber component per fiber 
  # and that one must be simple, it can be the one meeting O or not
  # assert this by adding fiber components under the additional condition that p stays in the algebraic lattice
  simples = []
  E = identity_matrix(QQ,rank(NS))
  z = zero_matrix(QQ,1, rank(NS))
  r = 2
  for fiber in _trivial_lattice(X)[3]
    fiber_type = fiber[2]
    fiber_rk = fiber_type[2]
    h = highest_root(fiber_type...)
    simple_indices = [r+i for i in 1:ncols(h) if isone(h[1,i])]
    simple_or_zero = [E[i:i,:] for i in simple_indices]
    push!(simple_or_zero, z)
    push!(simples,simple_or_zero)
    r += fiber_rk
  end
  G = gram_matrix(ambient_space(NS))
  pG = (p*G)[1:1,3:r]
  T = Tuple(simples)
  GF = G[3:r,3:r]
  candidates = QQMatrix[]
  for s in Base.Iterators.ProductIterator(T)
    g = sum(s)[:,3:r]
    y = (g - pG)
    xx = solve(GF,y;side=:left)
    x = zero_matrix(QQ,1, rank(NS))
    x[:,3:r] = xx
    if x in NS
      push!(candidates,p+x)
    end 
  end
  @assert length(candidates)>0
  # Select the candidate congruent to v modulo Triv 
  mwg = _mordell_weil_group(X)
  vmwg = mwg(vec(collect(v)))
  candidates2 = [mwg(vec(collect(x))) for x in candidates]
  i = findfirst(==(vmwg), candidates2)
  t = mwg(vec(collect(p - candidates[i])))
  mwl_tors_gens = [mwg(vec(collect(i[2]))) for i in tors]
  ag = abelian_group(zeros(ZZ,length(tors)))
  mwlAb = abelian_group(mwg)
  phi = hom(ag, mwlAb, mwlAb.(mwl_tors_gens))
  a = preimage(phi, mwlAb(t))
  for i in 1:ngens(ag)
    P += a[i]*rational_point(tors[i][1])
  end 
  
  p = candidates[i]
  k = (p*G*transpose(p))[1,1]
  # assert p^2 = -2
  p[1,1] = -k/2-1
  V = ambient_space(NS)
  @hassert :EllipticSurface 1 inner_product(V, p, p)[1,1]== -2
  @hassert :EllipticSurface 1 mwg(vec(collect(p))) == mwg(vec(collect(p)))
  @hassert :EllipticSurface 3 basis_representation(X,section(X,P))==vec(collect(p))
  return p, P
end
  
function _vertical_part(X::EllipticSurface, v::Vector{QQFieldElem}) 
  vv = matrix(QQ,1,length(v),v)
  p, P = _vertical_part(X,vv)
  pp = vec(collect(p))
  return pp, P
end

@doc raw"""
    elliptic_parameter(X::EllipticSurface, F::Vector{QQFieldElem}) -> LinearSystem

Return the elliptic parameter ``u`` of the divisor class `F`. 

The input `F` must be given with respect to the basis of
`algebraic_lattice(X)` and be an isotropic nef divisor. 
This method assumes that $X$ is a K3 surface.
"""
function elliptic_parameter(X::EllipticSurface, F::Vector{QQFieldElem})
  D1, D, P, l, c = horizontal_decomposition(X, F)
  return _elliptic_parameter(X, D1, D, P, l, c)
end

@doc raw"""
    _elliptic_parameter(X::EllipticSurface, D::AbsWeilDivisor, l, c)

Compute the linear system of ``D = (n-1) O + P + V``.
where V is vertical and `l` is the coefficient of the fiber class.
Assumes `D` nef and `D^2=0`.
Typically ``D`` is the output of `horizontal_decomposition`.
"""
function _elliptic_parameter(X::EllipticSurface, D1::AbsWeilDivisor, D::AbsWeilDivisor, P::EllipticCurvePoint, l::Int, c)
  S, piS = weierstrass_model(X);
  piX = weierstrass_contraction(X)
  c = function_field(X)(c)
  L = [i*c for i in linear_system(X, P, l)];
  LonX = linear_system(L, D1, check=false);

  LsubF, Tmat = subsystem(LonX, D);
  LsubFonS = [sum(Tmat[i,j]*L[j] for j in 1:ncols(Tmat)) for i in 1:nrows(Tmat)]

  @assert length(LsubFonS)==2
  u2 = LsubFonS[2]//LsubFonS[1]
  return u2
end


########################################################################
# Internal functionality for Weierstrass transformation 
########################################################################

@doc raw"""
    _normalize_hyperelliptic_curve(g::MPolyRingElem, parent=nothing)

Transform ``a(x)y^2 + b(x)y - h(x)`` in ``K(t)[x,y]`` to ``y'^2 - h(x')``
"""
function _normalize_hyperelliptic_curve(g::MPolyRingElem; parent::Union{MPolyRing, Nothing}=parent(g))
  R = Oscar.parent(g)
  @assert ngens(R) == 2 "polynomial must be bivariate"
  F = fraction_field(R)
  kt = coefficient_ring(R)
  (x, y) = gens(R)

  # Prepare the output ring
  if parent===nothing
    R1, (x1, y1) = R, gens(R)
  else
    R1 = parent
    @assert coefficient_ring(R1) == coefficient_ring(R) "coefficient ring of output is incompatible with input"
    (x1, y1) = gens(R1)
  end

  # Get the coefficients of g as a univariate polynomial in y
  ktx, X = polynomial_ring(kt, :X, cached=false)
  ktxy, Y = polynomial_ring(ktx, :y, cached=false)

  # Maps to transform to univariate polynomials in y
  split_map_R = hom(R, ktxy, [ktxy(X), Y])
  split_map_R1 = hom(R1, ktxy, [ktxy(X), Y])
  G = split_map_R(g)
  @assert degree(G) == 2 "polynomial must be of degree 2 in its second variable"

  #complete the square
  h, b, a = collect(coefficients(G))
  h = -h
  u = unit(factor(a))
  a = inv(u)*a
  b = inv(u)*b
  success, sqa = is_square_with_sqrt(a)
  @assert success "leading coefficient as univariate polynomial in the second variable must be a square"

  F1 = fraction_field(R1)
  psi = hom(R1, F, F.([x, (2*evaluate(a, x)*y + evaluate(b, x))//(2*evaluate(sqa, x))]))
  conv = MapFromFunc(ktx, R1, f->evaluate(f, x1))
  (a1, b1, sqa1) = conv.([a, b, sqa])
  phi = hom(R, F1, F1.([x1, (2*sqa1*y1-b1)//(2*a1)]))
  phiF = MapFromFunc(F, F1, x-> phi(numerator(x))//phi(denominator(x)))
  # the inverse map if wanted
  # psiF = MapFromFunc(F1, F, x-> psi(numerator(x))//psi(denominator(x)))
  # @assert all(phiF(psiF(F1(i)))==i for i in gens(R1))

  # absorb squares into y1
  g1 = numerator(phi(g))
  G1 = split_map_R1(g1)
  ff = factor(first(coefficients(G1)))
  c = prod([p^div(i, 2) for (p, i) in ff], init=one(ktx))
  #d = sqrt(my_coeff(g1, y1, 2))
  d = last(coefficients(split_map_R1(g1)))
  success, d = is_square_with_sqrt(d)
  @assert success "leading coefficient must be a square"

  phi1 = hom(R1, F1, [F1(x1), F1(evaluate(c, x1), evaluate(d, x1))*y1])
  phiF1 = MapFromFunc(F1, F1, x-> phi1(numerator(x))//phi1(denominator(x)))
  phi2 = compose(phi, phiF1)
  g2 = numerator(phi1(g1))
  #c = my_coeff(g2, y1, 2)
  c = last(coefficients(split_map_R1(g2)))
  g2 = divexact(g2, evaluate(c, x1))
  @assert degree(g2, gen(parent, 1)) <= 4 "degree in the first variable is too high"
  @assert degree(g2, gen(parent, 1)) >= 3 "degree in the first variable is too low"
  return g2, phi2
end

@doc raw"""
    transform_to_weierstrass(g::MPolyRingElem, x::MPolyRingElem, y::MPolyRingElem, P::Vector{<:RingElem})

Transform a bivariate polynomial `g` of the form `y^2 - Q(x)` with `Q(x)` of degree ``≤ 4``
to Weierstrass form. This returns a pair `(f, trans)` where `trans` is an endomorphism of the 
`fraction_field` of `parent(g)` and `f` is the transform. The input `P` must be a rational point 
on the curve defined by `g`, i.e. `g(P) == 0`.
"""
function transform_to_weierstrass(g::MPolyRingElem, x::MPolyRingElem, y::MPolyRingElem, P::Vector{<:RingElem})
  R = parent(g)
  F = fraction_field(R)
  @assert ngens(R) == 2 "input polynomial must be bivariate"
  @assert x in gens(R) "second argument must be a variable of the parent of the first"
  @assert y in gens(R) "third argument must be a variable of the parent of the first"
  # In case of variables in the wrong order, switch and transform the result.
  if x == R[2] && y == R[1]
    switch = hom(R, R, reverse(gens(R)))
    g_trans, trans = transform_to_weierstrass(switch(g), y, x, reverse(P))
    new_trans = MapFromFunc(F, F, f->begin
                                switch_num = switch(numerator(f))
                                switch_den = switch(denominator(f))
                                interm_res = trans(F(switch_num))//trans(F(switch_den))
                                num = numerator(interm_res)
                                den = denominator(interm_res)
                                switch(num)//switch(den)
                            end
                           )
    return switch(g_trans), new_trans
  end

  g = inv(coeff(g,[0,2]))*g # normalise g
  kk = coefficient_ring(R)
  kkx, X = polynomial_ring(kk, :x, cached=false)
  kkxy, Y = polynomial_ring(kkx, :y, cached=false)

  imgs = [kkxy(X), Y]
  split_map = hom(R, kkxy, imgs)

  G = split_map(g)
  @assert degree(G) == 2 "input polynomial must be of degree 2 in y"
  @assert all(h->degree(h)<=4, coefficients(G)) "input polynomial must be of degree <= 4 in x"
  @assert iszero(coefficients(G)[1]) "coefficient of linear term in y must be zero"
  @assert isone(coefficients(G)[2]) "leading coefficient in y must be one"
  
  if length(P) == 3 && isone(P[3])
      P = P[1:2]
  end 
      

  if length(P) == 2
    @assert iszero(evaluate(g, P)) "point does not lie on the hypersurface"
    (px, py) = P
  else 
    px = P[1]
  end
  #    assert g.subs({x:px,y:py})==0
  gx = -evaluate(g, [X + px, zero(X)])
  coeff_gx = collect(coefficients(gx))
  A = coeff(gx, 4)
  B = coeff(gx, 3)
  C = coeff(gx, 2)
  D = coeff(gx, 1)
  E = coeff(gx, 0)
  #E, D, C, B, A = coeff_gx
  if length(P)==3
    @req all(h->degree(h)<=3, coefficients(G)) "infinity (0:1:0) is not a point of this hypersurface"
    # y^2 = B*x^3+C*x^2+C*x+D
    x1 = F(inv(B)*x)
    y1 = F(inv(B)*y)
    trans = MapFromFunc(F, F, f->evaluate(numerator(f), [x1, y1])//evaluate(denominator(f), [x1, y1]))
    f_trans = B^2*trans(F(g))
    result = numerator(B^2*f_trans)
    return result, trans
  elseif !iszero(E)
    b = py
    a4, a3, a2, a1, a0 = A,B,C,D,E
    A = b
    B = a1//(2*b)
    C = (4*a2*b^2-a1^2)//(8*b^3)
    D = -2*b

    x1 = x//y
    y1 = (A*y^2+B*x*y+C*x^2+D*x^3)//y^2
    x1 = x1+px

    # TODO: The following are needed for the inverse. To be added eventually.
    # x2 = (y-(A+B*x+C*x^2))//(D*x^2)
    # y2 = x2//x
    # x2 = evaluate(x2, [x-px, y])
    # y2 = evaluate(y2, [x-px, y])

    # @assert x == evaluate(x1, [x2, y2])
    # @assert y == evaluate(y1, [x2, y2])
  else
    # TODO compute the inverse transformation (x2,y2)
    x1 = 1//x
    y1 = y//x^2
    g1 = numerator(evaluate(g, [x1, y1]))
    c = coeff(g1, [x], [3])
    x1 = evaluate(x1, [-x//c, y//c])
    y1 = evaluate(y1, [-x//c, y//c])
    x1 = x1+px
    #@assert x == evaluate(x1, [x2, y2])
    #@assert y == evaluate(y1, [x2, y2])
  end
  @assert F === parent(x1) "something is wrong with caching of fraction fields"
  # TODO: eventually add the inverse.
  trans = MapFromFunc(F, F, f->evaluate(numerator(f), [x1, y1])//evaluate(denominator(f), [x1, y1]))
  f_trans = trans(F(g))
  fac = [a[1] for a in factor(numerator(f_trans)) if isone(a[2]) && _is_in_weierstrass_form(a[1])]
  isone(length(fac)) || error("transform to weierstrass form did not succeed")

  # normalize the output
  result = first(fac)
  result = inv(first(coefficients(coeff(result, gens(parent(result)), [3, 0]))))*result

  return result, trans
end

function _is_in_weierstrass_form(f::MPolyRingElem)
  R = parent(f)
  @req ngens(R) == 2 "polynomial must be bivariate"
  # Helper function
  my_const(u::MPolyRingElem) = is_zero(u) ? zero(coefficient_ring(parent(u))) : first(coefficients(u))

  (x, y) = gens(R)
  f = -inv(my_const(coeff(f, [x, y], [0, 2]))) * f
  isone(-coeff(f, [x, y], [0, 2])) || return false
  isone(coeff(f, [x, y], [3, 0])) || return false
  
  a6 = coeff(f, [x,y], [0,0])
  a4 = coeff(f, [x,y], [1,0])
  a2 = coeff(f, [x,y], [2,0])
  a3 = -coeff(f, [x,y], [0,1])
  a1 = -coeff(f, [x,y], [1,1])
  a_invars = [my_const(i) for i in [a1,a2,a3,a4,a6]]
  (a1,a2,a3,a4,a6) = a_invars
  return f == (-(y^2 + a1*x*y + a3*y) + (x^3 + a2*x^2 + a4*x + a6))
end

########################################################################
# The three conversions from Section 39.1 in 
#   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface" 
# pp. 44--45.
########################################################################

function _elliptic_parameter_conversion(X::EllipticSurface, u::VarietyFunctionFieldElem; 
    case::Symbol=:case1, names=[:x, :y, :t]
  )
  @req variety(parent(u)) === weierstrass_model(X)[1] "function field element must live on the weierstrass model of the first argument"
  @req length(names) == 3 "need 3 variable names x, y, t"
  U = weierstrass_chart(X)
  R = ambient_coordinate_ring(U)
  x, y, t = gens(R)
  loc_eqn = first(gens(modulus(OO(U))))
  E = generic_fiber(X)::EllipticCurve
  f = equation(E)
  kk = base_ring(X)
  kkt_frac_XY = parent(f)::MPolyRing
  (xx, yy) = gens(kkt_frac_XY)
  kkt_frac = coefficient_ring(kkt_frac_XY)::AbstractAlgebra.Generic.FracField
  kkt = base_ring(kkt_frac)::PolyRing
  T = first(gens(kkt))

# kk = base_ring(U)
# kkt, T = polynomial_ring(kk, :T, cached=false)
# kkt_frac = fraction_field(kkt)
# kkt_frac_XY, (xx, yy) = polynomial_ring(kkt_frac, [:X, :Y], cached=false)
  R_to_kkt_frac_XY = hom(R, kkt_frac_XY, [xx, yy, kkt_frac_XY(T)])

  f_loc = first(gens(modulus(OO(U))))
  @assert f == R_to_kkt_frac_XY(f_loc) && _is_in_weierstrass_form(f) "local equation is not in Weierstrass form"
  a = a_invariants(E)

  u_loc = u[U]::AbstractAlgebra.Generic.FracFieldElem # the representative on the Weierstrass chart

  # Set up the ambient_coordinate_ring of the new Weierstrass-chart
  kkt2, t2 = polynomial_ring(kk, names[3], cached=false)
  kkt2_frac = fraction_field(kkt2)
  S, (x2, y2) = polynomial_ring(kkt2_frac, names[1:2], cached=false)
  FS = fraction_field(S)

  # Helper function
  my_const(u::MPolyRingElem) = is_zero(u) ? zero(coefficient_ring(parent(u))) : first(coefficients(u))

  # We verify the assumptions made on p. 44 of
  #   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface"
  # for the first case considered there.
  @assert all(x->isone(denominator(x)), a) "local equation does not have the correct form"
  a = numerator.(a)
  @assert iszero(a[1]) "local equation does not have the correct form"
  @assert degree(a[2]) <= 4 "local equation does not have the correct form"
  @assert iszero(a[3]) "local equation does not have the correct form"
  @assert degree(a[4]) <= 8 "local equation does not have the correct form"
  @assert degree(a[5]) <= 12 "local equation does not have the correct form" # This is really a₆ in the notation of the paper, a₅ does not exist.
  # reduce fraction
  u_frac = R_to_kkt_frac_XY(numerator(u_loc))//R_to_kkt_frac_XY(denominator(u_loc))
  u_num = numerator(u_frac)
  u_den = denominator(u_frac)
  if case == :case1
    # D = 2O
    u_poly = u_num*inv(u_den) # Will throw if the latter is not a unit
    # Extract a(t) and b(t) as in the notation of the paper
    a_t = my_const(coeff(u_poly, [xx, yy], [0, 0]))
    b_t = my_const(coeff(u_poly, [xx, yy], [1, 0]))

    a_t = evaluate(a_t, x2)
    b_t = evaluate(b_t, x2)
    phi = hom(R, FS, FS.([(t2 - a_t)//b_t, y2, x2]))
    f_trans = phi(f_loc)
    return numerator(f_trans), phi
  elseif case == :old
    # D = O + P
    @assert degree(u_num, 2) == 1 && degree(u_num, 1) <= 1 "numerator does not have the correct degree"
    @assert degree(u_den, 1) == 1 && degree(u_den, 2) == 0 "denominator does not have the correct degree"

    # We expect a form as on p. 44, l. -4
    denom_unit = my_const(coeff(u_den, [xx, yy], [1, 0]))
    x0 = -inv(denom_unit)*my_const(coeff(u_den, [xx, yy], [0, 0]))
    b_t = inv(denom_unit)*my_const(coeff(u_num, [xx, yy], [0, 1]))
    u_num = u_num - denom_unit * b_t * yy
    a_t = inv(denom_unit)*my_const(coeff(u_num, [xx, yy], [1, 0]))
    u_num = u_num - denom_unit * a_t * (xx - x0)
    @assert is_constant(u_num) "numerator is not in the correct form"
    y0 = my_const(coeff(u_num, [xx, yy], [0, 0])) * inv(denom_unit * b_t)

    @assert a_t + b_t*(yy + y0)//(xx - x0) == u_frac "decomposition failed"
    # We have 
    #
    #   y ↦ (u - a_t) * (x - x₀) / b_t - y₀ = (t₂ - a_t(x₂)) * (y₂ - x₀(x₂)) / b_t(x₂) - y₀(x₂)
    #   x ↦ y₂
    #   t ↦ x₂
    phi = hom(R, FS, FS.([y2, (t2 - evaluate(a_t, x2)) * (y2 - evaluate(x0, x2)) // evaluate(b_t, x2) - evaluate(y0, x2), x2]))
    f_trans = phi(f_loc)
    eqn1 = numerator(f_trans)
    # According to 
    #   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface" 
    # p. 45, l. 1 we expect the following cancelation to be possible:
    divisor_num = evaluate(numerator(x0), x2)
    divisor_den = evaluate(denominator(x0), x2)
    divisor = divisor_den * y2 - divisor_num
    success, eqn1 = divides(eqn1, divisor) # This division must only be possible in the ring K(x2)[y2].
                                           # Hence, multiplying by the denominator `divisor_den` is 
                                           # merely an educated guess.
    @assert success "division failed"
    return eqn1, phi
  elseif case == :case3
    # D = O + T

    @assert u_den == xx "elliptic parameter was not brought to the correct form"
    @assert degree(u_num, 1) <= 1 && degree(u_num, 2) <= 1 "numerator does not have the correct degrees"
    a_t = my_const(coeff(u_num, [xx, yy], [1, 0]))
    b_t = my_const(coeff(u_num, [xx, yy], [0, 1]))

    # New Weierstrass equation is of the form 
    #
    #   x^2 = h(t, u)
    #
    # so y₂ = x, x₂ = t, and t₂ = u.
    #
    # We have u = a_t + b_t * y/x ⇒ y = (u - a_t) * x / b_t = (t₂ - a_t(x₂)) * y₂ / b_t(x₂)
    phi = hom(R, FS, FS.([y2, (t2 - evaluate(a_t, x2)) * y2 // evaluate(b_t, x2), x2]))
    f_trans = phi(f_loc)
    eqn1 = numerator(f_trans)
    # According to 
    #   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface" 
    # p. 45, l. 15 we expect the following cancelation to be possible:
    success, eqn1 = divides(eqn1, y2)
    @assert success "equation did not come out in the anticipated form"
    return eqn1, phi
  elseif case == :case2
    # D = O + P
    @assert degree(u_num, 2) == 1 && degree(u_num, 1) <= 1 "numerator does not have the correct degree"
    @assert degree(u_den, 1) == 1 && degree(u_den, 2) <= 1 "denominator does not have the correct degree"

    # u = (ax + by + c)/(a'x + b'y + c')
    an = my_const(coeff(u_num, [xx, yy], [1, 0]))
    bn = my_const(coeff(u_num, [xx, yy], [0, 1]))
    cn = my_const(coeff(u_num, [xx, yy], [0, 0]))

    ad = my_const(coeff(u_den, [xx, yy], [1, 0]))
    bd = my_const(coeff(u_den, [xx, yy], [0, 1]))
    cd = my_const(coeff(u_den, [xx, yy], [0, 0]))

    @assert (an*xx+bn*yy+cn)//(ad*xx+bd*yy+cd) == u_frac "decomposition failed"


    v = solve(matrix(parent(an), 2, 2, [-an, bn,-ad, bd]), matrix(parent(an), 2, 1, [cn, cd]); side=:right)
    x0 = v[1,1]
    y0 = v[2,1]
    @assert evaluate(f_loc,[x0,y0,gen(parent(x0))])==0

    ad = evaluate(ad,x2)
    an = evaluate(an,x2)
    bd = evaluate(bd,x2)
    bn = evaluate(bn,x2)
    cn = evaluate(cn,x2)
    cd = evaluate(cd,x2)
    #x0 = evaluate(x0,x2)
    #y0 = evaluate(y0,x2)


    imgy = -FS(((ad*t2 - an )*y2 + (cd*t2 -cn)) //(bd*t2 -bn))

    # We have
    #
    #   y ↦ -((ad u - an )x + (cd u -cn)) // (bd*u -bn)
    #   x ↦ y₂
    #   t ↦ x₂
    phi = hom(R, FS, [y2, imgy, x2])
    f_trans = phi(f_loc)
    eqn1 = numerator(f_trans)
    # According to
    #   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface"
    # p. 45, l. 1 we expect the following cancellation to be possible:
    divisor_num = evaluate(numerator(x0), x2)
    divisor_den = evaluate(denominator(x0), x2)
    divisor = divisor_den * y2 - divisor_num
    success, eqn1 = divides(eqn1, divisor) # This division must only be possible in the ring K(x2)[y2].
                                           # Hence, multiplying by the denominator `divisor_den` is
                                           # merely an educated guess.
    @assert success "division failed"
    return eqn1, phi
  else
    error("case not recognized")
  end
end


# Given a bivariate polynomial over a univariate function field, 
# normalize the associated elliptic curve so that the usual constructor 
# for elliptic surfaces digests it, and then return it, together with the 
# transformation on the algebraic side. 
#
# The transformation is a morphism from the fraction field of the 
# parent of g to the fraction field of the `ambient_coordinate_ring` 
# of the `weierstrass_chart` of the resulting surface.
function _elliptic_surface_with_trafo(g::MPolyRingElem{<:AbstractAlgebra.Generic.FracFieldElem}; minimize::Bool=true)
  x, y = gens(parent(g))
  E = elliptic_curve(g, x, y)
  kkt = base_field(E)
  kk = coefficient_ring(base_ring(kkt))

  FFt, t = rational_function_field(kk, :t)

  # The following three commands won't work unless we convert to a rational_function_field
  EE = base_change(x->evaluate(x, t), E)
  if minimize
    EE = tates_algorithm_global(EE)
    EE, _ = short_weierstrass_model(EE)
    EE, _ = integral_model(EE)
  end
  
  # ...and back.
  E2 = base_change(x->evaluate(x, gen(kkt)), EE)

  @assert is_isomorphic(E, E2)
  a, b, _ = rational_maps(isomorphism(E2, E))

  eq_E = equation(E)
  eq_E2 = equation(E2)

  h = evaluate(eq_E, [a, b])
  @assert divides(h, eq_E2)[1]

  cod = parent(a)::MPolyRing

  #phi = hom(R, cod, cod.([a, b]))
  #Phi = extend_domain_to_fraction_field(phi)
  
  result = elliptic_surface(E2, 2)
  W = weierstrass_chart(result)
  R = ambient_coordinate_ring(W)
  FR = fraction_field(R)

  help_map = hom(cod, FR, t->evaluate(t, FR(R[3])), FR.([R[1], R[2]]))
  A = help_map(a)
  B = help_map(b)

  res_map = hom(parent(g), FR, t->evaluate(t, FR(R[3])), [A, B])
  return result, extend_domain_to_fraction_field(res_map)
end



# TODO: Instead return an abelian group A and two maps. 
# algebraic_lattice -> A 
# A -> MWL= E(k(t))
function _mordell_weil_group(X)
  N = algebraic_lattice(X)[3]
  V = ambient_space(N)
  t = nrows(trivial_lattice(X)[2])
  Triv = lattice(V, identity_matrix(QQ,dim(V))[1:t,:])
  return torsion_quadratic_module(N, Triv;modulus=1, modulus_qf=1, check=false)
end 

