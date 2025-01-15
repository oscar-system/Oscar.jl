function _parse_symbol(symbol::String, I::UnitRange)
  return ["$symbol[$i]" for i in I]
end
function _parse_symbol(symbol::String, n::Int, I::UnitRange)
  return [symbol*"[$n, $i]" for i in I]
end

@doc raw"""
      present_finite_extension_ring(F::Oscar.AffAlgHom)

Given a finite homomorphism `F` $:$ `A` $\rightarrow$ `B`  of algebras of type `<: Union{MPolyRing, MPolyQuoRing}` over a field, return a presentation

$A^r \rightarrow A^s\rightarrow B \rightarrow 0$

of `B` as an `A`-module.

More precisely, return a tuple `(gs, PM, sect)`, say, where
- `gs` is a vector of polynomials representing generators for `B` as an `A`-module,
- `PM` is an `r` $\times$ `s`-matrix of polynomials defining the map $A^r \rightarrow A^s$, and
- `sect` is a function which gives rise to a section of the augmentation map $ A^s\rightarrow B$.

!!! note
    The finiteness condition on `F` is checked by the function.

!!! note
    The function is implemented so that the last element of `gs` is `one(B)`.

# Examples
```jldoctest
julia> RA, (h,) = polynomial_ring(QQ, [:h]);

julia> A, _ = quo(RA, ideal(RA, [h^9]));

julia> RB, (k, l) = polynomial_ring(QQ, [:k, :l]);

julia> B, _ = quo(RB, ideal(RB, [k^3, l^3]));

julia> F = hom(A, B, [k+l])
Ring homomorphism
  from quotient of multivariate polynomial ring by ideal (h^9)
  to quotient of multivariate polynomial ring by ideal (k^3, l^3)
defined by
  h -> k + l

julia> gs, PM, sect = present_finite_extension_ring(F);

julia> gs
3-element Vector{QQMPolyRingElem}:
 l^2
 l
 1

julia> PM
3×3 Matrix{QQMPolyRingElem}:
 h^3     0       0
 -3*h^2  h^3     0
 3*h     -3*h^2  h^3

julia> sect(k*l)
3-element Vector{QQMPolyRingElem}:
 -1
 h
 0

```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [z^2-y^2*(y+1)]);

julia> A, _ = quo(R, I);

julia> B, (s,t) =  polynomial_ring(QQ, [:s, :t]);

julia> F = hom(A,B, [s, t^2-1, t*(t^2-1)])
Ring homomorphism
  from quotient of multivariate polynomial ring by ideal (-y^3 - y^2 + z^2)
  to multivariate polynomial ring in 2 variables over QQ
defined by
  x -> s
  y -> t^2 - 1
  z -> t^3 - t

julia> gs, PM, sect = present_finite_extension_ring(F);

julia> gs
2-element Vector{QQMPolyRingElem}:
 t
 1

julia> PM
2×2 Matrix{QQMPolyRingElem}:
 y   -z
 -z  y^2 + y

julia> sect(t)
2-element Vector{QQMPolyRingElem}:
 1
 0

julia> sect(one(B))
2-element Vector{QQMPolyRingElem}:
 0
 1

julia> sect(s)
2-element Vector{QQMPolyRingElem}:
 0
 x

```

```jldoctest
julia> A, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c]);

julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x*y]);

julia> B, _ = quo(R, I);

julia> (x, y, z) = gens(B);

julia> F = hom(A, B, [x^2+z, y^2-1, z^3])
Ring homomorphism
  from multivariate polynomial ring in 3 variables over QQ
  to quotient of multivariate polynomial ring by ideal (x*y)
defined by
  a -> x^2 + z
  b -> y^2 - 1
  c -> z^3

julia> gs, PM, sect = present_finite_extension_ring(F);

julia> gs
2-element Vector{QQMPolyRingElem}:
 y
 1

julia> PM
2×2 Matrix{QQMPolyRingElem}:
 a^3 - c  0
 0        a^3*b + a^3 - b*c - c

julia> sect(y)
2-element Vector{QQMPolyRingElem}:
 1
 0

julia> sect(one(B))
2-element Vector{QQMPolyRingElem}:
 0
 1

```
"""
function  present_finite_extension_ring(F::Oscar.AffAlgHom)
  A, B = F.domain, F.codomain
  a, b = ngens(A), ngens(B)
  
  if A isa MPolyQuoRing
    AR = base_ring(A)
  else
    AR = A
  end
  if B isa MPolyQuoRing
    BR = base_ring(B)
    M = [F(gens(A)[i]).f for i = 1:a]
  else
    BR = B
    M = [F(gens(A)[i]) for i = 1:a]
  end

  @assert base_ring(AR) == base_ring(BR)

  I = ideal(BR, isdefined(B, :I) ? vcat(gens(B.I), M) : M)
  C, _ = quo(BR, I)
  gs = monomial_basis(C) # monomials whose residue classes form a K-basis of B
  @assert gs[end] == 1 # the last one should always be 1
  g = length(gs)

  R, _ = tensor_product(BR, AR, use_product_ordering = true)
  ba = gens(R)
  ARtoR = hom(AR, R, ba[b+1:end], check = false)
  BRtoR = hom(BR, R, ba[1:b], check = false)
  RtoAR = hom(R, AR, vcat(repeat([AR()], b), gens(AR)))
  gs_lift = [BRtoR(g) for g in gs]
  
  # compute the ideal J of the graph of F
  Rels = [ba[b+i]-BRtoR(m) for (i,m) in enumerate(M)]
  if isdefined(A, :I) for g in gens(A.I) push!(Rels, ARtoR(g)) end end
  if isdefined(B, :I) for g in gens(B.I) push!(Rels, BRtoR(g)) end end
  J = ideal(R, Rels) # the ideal of the graph of F
  V = groebner_basis(J)

  sect = x -> (y = reduce(BRtoR(x), gens(V));
	      ans = elem_type(AR)[];
	      for i in 1:g
	        q = div(y, gs_lift[i])
	        push!(ans, RtoAR(q))
	        y -= q * gs_lift[i]
	      end; ans)

  FM = free_module(R, g)
  gB = elem_type(FM)[FM(push!([j == i ? R(1) : R() for j in 1:g-1], -gs_lift[i])) for i in 1:g-1]
  gJ = elem_type(FM)[FM([j==i ? x : R() for j in 1:g]) for x in gens(V) for i in 1:g]
  U  = vcat(gB, gJ)
  S, _ = sub(FM, U)
  P = groebner_basis(S, ordering = default_ordering(R)*lex(FM))
  Rw, _ = grade(R, vcat(repeat([1], b), repeat([0], a)))
  RtoRw = hom(R, Rw, gens(Rw))
  inA = x -> x == zero(Rw) ?  true : (degree(Int, leading_term(RtoRw(x)))) <= 0
  PM = vcat([(RtoAR.(transpose(Vector(P[i])))) for i in 1:ngens(P) if all(inA, Vector(P[i]))]...)
  return gs, PM, sect
end

######################################
@doc raw"""
      blowup(i::AbstractVarietyMap; symbol::String="e")

Given an inclusion `i`$ : $ `X` $\rightarrow$ `Y`, say, return the blowup of `Y` along `X`.

More precisely, return a tuple `(Bl, E, j)`, say, where
- `Bl`, an abstract variety, is the blowup,
- `E`, an abstract variety, is the exceptional divisor, and
- `j`, a map of abstract varieties, is the inclusion of `E` into `Bl`.

!!! note
    The resulting maps `Bl` $\rightarrow$ `Y` and `E` $\rightarrow$ `X` are obtained entering `structure_map(Bl)` and `structure_map(E)`, respectively.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P5 = abstract_projective_space(5, symbol = "H")
AbstractVariety of dim 5

julia> h = gens(P2)[1]
h

julia> H = gens(P5)[1]
H

julia> i = map(P2, P5, [2*h])
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 5

julia> Bl, E, j = blowup(i)
(AbstractVariety of dim 5, AbstractVariety of dim 4, AbstractVarietyMap from AbstractVariety of dim 4 to AbstractVariety of dim 5)

julia> e, HBl = gens(chow_ring(Bl))
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 e
 H

julia> integral((6*HBl-2*e)^5)
3264

```
"""
function blowup(i::AbstractVarietyMap; symbol::String = "e")
  # use construction via generators and relations as in Eisenbud-Harris, Section 13.6
  #  E ---j---> Bl
  #  |          |
  #  g          f
  #  |          |
  #  Z ---i---> X

  SED = symbol # SED = "e"
  Z, X = i.domain, i.codomain
  AZ, RZ = Z.ring, base_ring(Z.ring)
  AX, RX = X.ring, base_ring(X.ring)

  # we write N for the normal bundle N_{ZX}

  N = -i.T # normal bundle
  rN = rank(N) # codimension of Z in X
  rN <= 0 && error("not a proper subvariety")

  # we construct E as the projective bundle PN of N, see [E-H, p 472]

  E = abstract_projective_bundle(N) 
  AE, RE = E.ring, base_ring(E.ring)
  g = E.struct_map
  ζ = g.O1  # the first Chern class of O_PN(1)
  ### Q =  E.bundles[2] # the universal quotient bundle on PN
  ### ctopQ = top_chern_class(Q)

  # we set up the generators of ABl

  # we first represent AZ as an AX-module

  gs, PM, sect = present_finite_extension_ring(i.pullback)
  ngs = length(gs)
  # note: gs is a vector of polynomials representing generators for AZ as an AX-module;
  #       write Z(gs[i]) for the class in AZ represented by gs[i];
  #       if e[i] = jₓgˣ(Z(gs[i])), then e[end] is the class of the exceptional divisor in ABl;
  #       here, present_finite_extension_ring needs to return the generators accordingly;
  #       that is, Z(gs[end]) must by the generator 1_{AZ}; note that ζ = -jˣ(e[end]).

  # by E-H, Prop. 13.12, ABl is generated by jₓ(AE) and  fˣ(AX) as a K-algebra;
  # equivalently, as an AX-algebra, ABl is generated by the e[i].
  syms = vcat(push!(_parse_symbol(SED, 1:ngs-1), SED), string.(gens(RX))) # set e = e[end]
  degs = [degree(Int, Z(gs[i])) + 1 for i in 1:ngs]
  degsRX = [Int(degree(gens(RX)[i])[1])  for i = 1:ngens(RX)]
  RBl = graded_polynomial_ring(X.base, syms, vcat(degs, degsRX))[1]
  
  ev, xv = gens(RBl)[1:ngs], gens(RBl)[ngs+1:end]
  RXtoRBl = hom(RX, RBl, xv) # fˣ on polynomial ring level
  jₓgˣ = x -> sum(ev .* RXtoRBl.(sect(x.f))) # AZ --> RBl

  # we set up the relations of ABl

  Rels = elem_type(RBl)[]

  # 1) relations from AX

  if isdefined(AX, :I)
    for r in RXtoRBl.(gens(AX.I)) push!(Rels, r) end
  end

  # 2) relations for AZ as an AX-module

  for r in RXtoRBl.(PM)*ev push!(Rels, r) end

  # 3) relation from AE: ∑ gˣcₖ(N) ⋅ ζᵈ⁻ᵏ = 0, see EH Thm. 9.6

  cN = total_chern_class(N)[0:rN] # cN[k] = c_k ₋₁(N)
  push!(Rels, sum([jₓgˣ(cN[k+1]) * (-ev[end])^(rN-k) for k in 0:rN]))

  # 4) jₓx ⋅ jₓy = -jₓ(x⋅y⋅ζ) # rule for multiplication, EH, Prop. 13.12
  # recall that e[i] = jₓgˣ(Z(gs[i]))

  for j in 1:ngs-1, k in j:ngs-1
    push!(Rels, ev[j] * ev[k] + jₓgˣ(Z(gs[j] * Z(gs[k]))) * (-ev[end]))
  end

  # 5) relation as in the proof of [E-H], Theorem 13.14:
  # RXtoRBliₓx = jₓ(gˣx ⋅ ctop(Q)) where Q is the tautological quotient bundle on E
  # we have ctop(Q) = ∑ gˣcₖ₋₁(N) ⋅ ζᵈ⁻ᵏ, EH p. 477

  for j in 1:ngs
    lhs = RXtoRBl(i.pushforward(Z(gs[j])).f) # this is the crucial step where iₓ is needed
    rhs = sum([jₓgˣ(Z(gs[j]) * cN[k]) * (-ev[end])^(rN-k) for k in 1:rN])
    push!(Rels, lhs - rhs)
  end
 
  Rels = minimal_generating_set(ideal(RBl, Rels)) ### TODO: make this an option?
  ABl, _ = quo(RBl, Rels)
  Bl = abstract_variety(X.dim, ABl)

  # Bl and g being constructed, we set up the morphisms f and j

  # pushforward of f and more

  RBltoRX = hom(RBl, RX, vcat(repeat([RX()], ngs), gens(RX)))
  fₓ = x -> (xf = simplify(x).f;
	     X(RBltoRX(xf));)
  fₓ = map_from_func(fₓ, ABl, AX)
  f = AbstractVarietyMap(Bl, X, Bl.(xv), fₓ)
  Bl.struct_map = f
  if isdefined(X, :point) Bl.point = f.pullback(X.point) end

  # pullback of j
  
  jˣ = vcat([-ζ * g.pullback(Z(xi)) for xi in gs], [g.pullback(i.pullback(f)) for f in gens(AX)])
  
  # pushforward of j: write as a polynomial in ζ, and compute term by term

  REtoRZ = hom(RE, RZ, pushfirst!(gens(RZ), RZ()))
  jₓ = x -> (xf = simplify(x).f;
            ans = RBl();
	    for k in rN-1:-1:0
	      q = div(xf, ζ.f^k)
	      ans += jₓgˣ(Z(REtoRZ(q))) * (-ev[end])^k
	      xf -= q * ζ.f^k
	    end;
	    Bl(ans))
  jₓ = map_from_func(jₓ, AE, AX)
  j = AbstractVarietyMap(E, Bl, jˣ, jₓ)

  # the normal bundle of E in Bl is O(-1)

  j.T = -E.bundles[1]

  # finally, compute the tangent bundle of Bl

  # 0 → Bl.T → RXtoRBl(X.T) → jₓ(Q) → 0 where Q is the tautological quotient bundle

  f.T = -pushforward(j, E.bundles[2])
  Bl.T = pullback(f, X.T) + f.T

# chern(Bl.T) can be readily computed from its Chern character, but the following is faster
  α = sum(sum((binomial(ZZ(rN-j), ZZ(k)) - binomial(ZZ(rN-j), ZZ(k+1))) * ζ^k for k in 0:rN-j) * g.pullback(chern_class(N, j)) for j in 0:rN)
  Bl.T.chern = simplify(f.pullback(total_chern_class(X.T)) + j.pushforward(g.pullback(total_chern_class(Z.T)) * α))
  set_attribute!(E, :projections => [j, g])
  set_attribute!(Bl, :exceptional_divisor => E)
  set_attribute!(Bl, :description => "Blowup of $X with center $Z")
  if get_attribute(Z, :alg) == true && get_attribute(X, :alg) == true
    set_attribute!(Bl, :alg => true)
  end
  return Bl, E, j
end


@doc raw"""
    function blowup_points(X::AbstractVariety, n::Int; symbol::String = "e")

Return the blowup of `X` at `n` points.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> Bl = blowup_points(P2, 1)
AbstractVariety of dim 2

julia> chow_ring(Bl)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    e -> [1]
    h -> [1]
  by ideal (e*h, e^2 + h^2)

```
"""
function blowup_points(X::AbstractVariety, n::Int; symbol::String = "e")
  SED = symbol # SED = "e"
  if n == 1
    symbs = [SED]
  else
    symbs = _parse_symbol(SED, 1:n)
  end
  Bl = X
  P = abstract_point(base = X.base)
  for i in 1:n
    Bl = blowup(map(P, Bl, [zero(P.ring) for j = 1:i]), symbol=symbs[i])[1]
  end
  set_attribute!(Bl, :description => "Blowup of $X at $n points")
  Bl.struct_map = map(Bl, X)
  if get_attribute(X, :alg) == true
    set_attribute!(Bl, :alg => true)
  end
  return Bl
end







