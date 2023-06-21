# TODO:
# shortcut to get the coefficient(s) of a weil divisor
# einsetzen von VarietyFunctionFieldElem in univariate polynome
# in_linear_system hakt noch total
#=
julia> in_linear_system(linsys[1], D,check=false)
false
=#


"""
    prop217(E::EllCrv, P::EllCrvPt, k)

Compute a basis for the linear system
|O + P + kF|
on the  minimal elliptic (K3) surface defined by E.
Here F is the class of a fiber O the zero section
and P any non-torsion section.

```jldoctest
julia> kt,t = polynomial_ring(GF(29),:t);

julia> ktfield = fraction_field(kt);

julia> bk = [((17*t^4 + 23*t^3 + 18*t^2 + 2*t + 6, 8*t^5 + 2*t^4 + 6*t^3 + 25*t^2 + 24*t + 5 )),
             ((17*t^6 + 3*t^5 + 16*t^4 + 4*t^3 + 13*t^2 + 6*t + 5)//(t^2 + 12*t + 7), (4*t^8 + 19*t^7 + 14*t^6 + 18*t^5 + 27*t^4 + 13*t^3 + 9*t^2 + 14*t + 12)//(t^3 + 18*t^2 + 21*t + 13) ),
             ((17*t^6 + 10*t^5 + 24*t^4 + 15*t^3 + 22*t^2 + 27*t + 5)//(t^2 + 16*t + 6), (20*t^8 + 24*t^7 + 22*t^6 + 12*t^5 + 21*t^4 + 21*t^3 + 9*t^2 + 21*t + 12)//(t^3 + 24*t^2 + 18*t + 19) ),
             ((17*t^8 + 21*t^7 + 20*t^5 + 24*t^4 + 21*t^3 + 4*t^2 + 9*t + 13)//(t^4 + 17*t^3 + 12*t^2 + 28*t + 28), (23*t^11 + 25*t^10 + 8*t^9 + 7*t^8 + 28*t^7 + 16*t^6 + 7*t^5 + 23*t^4 + 9*t^3 + 27*t^2 + 13*t + 13)//(t^6 + 11*t^5 + 14*t^4 + 13*t^3 + 6*t^2 + 18*t + 12) )];

julia> E = EllipticCurve(ktfield,[3*t^8+24*t^7+22*t^6+15*t^5+28*t^4+20*t^3+16*t^2+26*t+16, 24*t^12+27*t^11+28*t^10+8*t^9+6*t^8+16*t^7+2*t^6+10*t^5+3*t^4+22*t^3+27*t^2+10*t+3]);

julia> bk = [E(collect(i)) for i in bk];

julia> prop217(E,bk[2],2)
5-element Vector{Any}:
 (t^2 + 12*t + 7, 0)
 (t^3 + 8*t + 3, 0)
 (t^4 + 23*t + 2, 0)
 (25*t + 22, 1)
 (12*t + 28, t)

julia> prop217(E,bk[1],1)
2-element Vector{Any}:
 (1, 0)
 (t, 0)
```
"""
function prop217(E::EllCrv, P::EllCrvPt, k)
  @req !iszero(P[3]) "P must not be torsion" # seems like we cannot check this
  xn = numerator(P[1])
  xd = denominator(P[1])
  yn = numerator(P[2])
  yd = denominator(P[2])
  OP = divexact(max(degree(xd), degree(xn) - 4), 2)
  dega = k + 2*OP
  degb = k + 2*OP - 2 - divexact(degree(xd), 2) #?
  base = base_ring(X)

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
  M = matrix(base, length(eqns), length(ab), reduce(vcat,cc, init=elem_type(base)[]))
  # @assert M == matrix(base, cc) # does not work if length(eqns)==0
  kerdim, K = kernel(M)
  result = []
  Bt = base_ring(base_field(E))
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

function prop217(X::AbsCoveredScheme, Weierstrasschart::AbsSpec, E::EllCrv, P::EllCrvPt, k)
  FX = function_field(X)
  U = Weierstrasschart
  (x,y,t) = ambient_coordinates(Weierstrasschart)

  sections = elem_type(FX)[]
  xn = numerator(P[1])
  xd = denominator(P[1])
  yn = numerator(P[2])
  yd = denominator(P[2])

  I = Oscar.ambient_closure_ideal(U)
  IP = ideal([x*xd(t)-xn(t),y*yd(t)-yn(t)])
  issubset(I, IP) || error("P does not define a point on the Weierstrasschart")

  @assert gcd(xn, xd)==1
  @assert gcd(yn, yd)==1
  ab = prop217(E, P, k)
  d = divexact(yd, xd)(t)
  den = t^k(x*xd(t) - xn(t))
  for (a,b) in ab
    c = divexact(b*yn - a*xn, xd)
    num = a(t)*x+b(t)*d*y + c(t)
    push!(sections, FX(num//den))
  end
  return sections
end

function Oscar.saturation(I::IdealSheaf, J::IdealSheaf)
  X = scheme(I)
  K = IdDict{AbsSpec,Ideal}()
  for U in X[1]
    K[U] = saturation(I(U),J(U))
  end
  return IdealSheaf(X, K, check=false)
end

function (f::AbstractAlgebra.Generic.Frac)(t)
  return numerator(f)(t)//denominator(f)(t)
end

function my_coeff(g::MPolyRingElem, x, deg)
  R = parent(g)
  @req parent(x)=== R "parent missmatch"
  i = findfirst(==(x), gens(R))
  c = MPolyBuildCtx(R)
  for (co, mon) in coefficients_and_exponents(g)
    if mon[i] == deg
      mon[i] = 0
      push_term!(c, co, mon)
    end
  end
  return finish(c)
end

function my_degree(g::MPolyRingElem, x)
  R = parent(g)
  @req R === parent(x) "g and x must have the same parent"
  i = findfirst(==(x), gens(R))
  return maximum(c[i] for c in exponents(g))
end

"""
Transform
a(x)y^2 + b(x) y = h(x)
to y'^2 = h(x')
"""
function normalize_quartic(g::MPolyRingElem, parent_out=nothing)
  R = parent(g)
  F = fraction_field(R)
  kt = base_ring(R)
  (x, y) = gens(R)

  #complete the square
  a = my_coeff(g, y, 2)
  b = my_coeff(g, y, 1)
  u = unit(factor(a))
  a = inv(u)*a
  b = inv(u)*b
  sqa = sqrt(a)
  # inverse map
  if parent_out isa Nothing
    R1, (x1,y1) = polynomial_ring(kt, [:x, :y])
  else
    R1 = parent_out
    @assert base_ring(R1) == base_ring(R)
    (x1, y1) = gens(R1)
  end
  F1 = fraction_field(R1)
  psi = hom(R1, F, F.([x, (2*a*y + b)//(2*sqa)]))
  conv = hom(R, R1, [x1, 0])
  (a1,b1,sqa1) = conv.((a,b,sqa))
  phi = hom(R, F1, F1.([x1, (2*sqa1*y1-b1)//(2*a1)]))
  phiF = map_from_func(x-> phi(numerator(x))//phi(denominator(x)), F, F1)
  psiF = map_from_func(x-> psi(numerator(x))//psi(denominator(x)), F1, F)
  @assert all(phiF(psiF(F1(i)))==i for i in gens(R1))

  # absorb squares into y1
  g1 = numerator(phi(g))
  ff = factor(hom(R1,R1,[x1,0])(g1))
  c = prod([p^divexact(i,2) for (p,i) in ff if mod(i,2)==0],init=R1(1))
  d = sqrt(my_coeff(g1, y1, 2))

  phi1 = hom(R1, F1, [x1, (c//d)*y1])
  phiF1 = map_from_func(x-> phi1(numerator(x))//phi1(denominator(x)), F1, F1)
  phi2 = compose(phi, phiF1)
  g2 = numerator(phi1(g1))
  c = my_coeff(g2, y1, 2)
  g2 = divexact(g2, c)
  return g2, phi2
end


function transform_to_weierstrass(g::MPolyElem, x::MPolyElem, y::MPolyElem, P::Vector{T};
    return_inverse::Bool=false
  ) where {T<:FieldElem}
#    r"""
#    y^2 - quartic(x)
#    """
    @assert my_degree(g, y)== 2
    @assert my_degree(g, x)== 4
    @assert my_coeff(g, y, 1) == 0
    @assert my_coeff(g, y, 2) == 1
    R = parent(g)
    kk = coefficient_ring(R)
    S, X = polynomial_ring(kk, "X", cached=false)
    length(P) == 2 || error("need precisely two point coordinates")
    (px, py) = P
#    assert g.subs({x:px,y:py})==0
    @assert iszero(evaluate(g, P)) "point does not lie on the hypersurface"
    gx = -evaluate(g, [X + px, zero(X)])
    coeff_gx = collect(coefficients(gx))
    E, D, C, B, A = coeff_gx
    if !iszero(E)
      b = py
      a4, a3, a2, a1, a0 = A,B,C,D,E
      A = b
      B = a1//(2*b)
      C = (4*a2*b^2-a1^2)//(8*b^3)
      D = -2*b

      x1 = x//y
      y1 = (A*y^2+B*x*y+C*x^2+D*x^3)//y^2
      x1 = x1+px

      x2 = (y-(A+B*x+C*x^2))//(D*x^2)
      y2 = x2//x
      x2 = evaluate(x2, [x-px, y])
      y2 = evaluate(y2, [x-px, y])

      @assert x == evaluate(x1, [x2, y2])
      @assert y == evaluate(y1, [x2, y2])
    else
      error("formulas to be confirmed before use")
      # TODO compute the inverse transformation (x2,y2)
      x1 = 1//x
      y1 = y//x^2
      g1 = numerator(evaluate(g, [x1, y1]))
      c = coeff(g1, [x], [3])
      x1 = evaluate(x1, [-x//c, y//c])
      y1 = evaluate(y1, [-x//c, y//c])
      x1 = x1+px
      @assert x == evaluate(x1, [x2, y2])
      @assert y == evaluate(y1, [x2, y2])
    end
    if return_inverse
      return x1, y1, x2, y2
    else
      return x1, y1
    end
end

function elliptic_curve(f::MPolyRingElem, x,y)
  # asserts
  k = coefficient_ring(f)
  c = coeff(f, [x,y], [3,0])
  c = c(k(0),k(0))
  @assert parent(c)===k
  f = inv(c)*f
  @assert coeff(f, [x,y], [0,2]) == -1
  a6 = coeff(f, [x,y], [0,0])
  a4 = coeff(f, [x,y], [1,0])
  a2 = coeff(f, [x,y], [2,0])
  a3 = -coeff(f, [x,y], [0,1])
  a1 = -coeff(f, [x,y], [1,1])
  a_invars = [i(k(0),k(0)) for i in (a1,a2,a3,a4,a6)]
  (a1,a2,a3,a4,a6) = a_invars
  @assert f  == (-(y^2 + a1*x*y + a3*y) + (x^3 + a2*x^2 + a4*x + a6))
  E = EllipticCurve(k, [a1,a2,a3,a4,a6])
  return E
end

function is_isomorphic_with_map(G1::Graph, G2::Graph)
  f12 = Polymake.graph.find_node_permutation(G1.pm_graph, G2.pm_graph)
  if isnothing(f12)
    return false, Vector{Int}()
  end
  return true, Polymake.to_one_based_indexing(f12)
end

function graph(G::MatElem)
  n = nrows(G)
  g = Graph{Undirected}(n)
  for i in 1:n
    # small hack to single out the fiber
    if G[i,i]==0
      add_edge!(g,i,i)
    end
    for j in 1:i-1
      if G[i,j] == 1
        add_edge!(g,i,j)
      end
    end
  end
  return g
end

function is_isomorphic_with_permutation(A1::MatElem, A2::MatElem)
  b, T = is_isomorphic_with_map(graph(A1),graph(A2))
  @assert b || A1[T] == A2
  return b, T
end

function extend_domain(phi::Map, K::Field)
  F = codomain(phi)
  @assert base_ring(K) == domain(phi)
  return map_from_func(x-> phi(numerator(x))*inv(phi(denominator(x))), K, F)
end

function hom(dom::FracField, cod::FracField, imgs::Vector)
  imgs = cod.(imgs)
  phi = hom(base_ring(dom), cod, imgs)
  return extend_domain(phi, cod)
end



