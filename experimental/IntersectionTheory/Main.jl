###############################################################################
#
# AbsBundle
#
@doc Markdown.doc"""
    bundle(X::AbsVariety, ch)
    bundle(X::AbsVariety, r, c)
Construct a bundle on $X$ by specifying its Chern character, or its rank and
total Chern class.
"""
bundle(X::V, ch::MPolyDecRingOrQuoElem) where V <: AbsVarietyT = AbsBundle(X, ch)
bundle(X::V, r::RingElement, c::MPolyDecRingOrQuoElem) where V <: AbsVarietyT = AbsBundle(X, r, c)

==(F::AbsBundle, G::AbsBundle) = ch(F) == ch(G)

@doc Markdown.doc"""
    ch(F::AbsBundle)
Return the Chern character."""
ch(F::AbsBundle) = (
  if !isdefined(F, :ch) F.ch = F.rank + _logg(F.chern) end;
  F.ch)

@doc Markdown.doc"""
    chern(F::AbsBundle)
    chern(F::TnBundle)
Compute the total Chern class.
"""
chern(F::AbsBundle) = (
  if !isdefined(F, :chern) F.chern = _expp(F.ch) end;
  F.chern)

@doc Markdown.doc"""
    chern(k::Int, F::AbsBundle)
    chern(k::Int, F::TnBundle)
Compute the $k$-th Chern class.
"""
chern(k::Int, F::AbsBundle) = (
  isdefined(F, :chern) && return chern(F)[k];
  _expp(F.ch, truncate=k)[k])

@doc Markdown.doc"""
    ctop(F::AbsBundle)
    ctop(F::TnBundle)
Compute the top Chern class.
"""
ctop(F::AbsBundle) = chern(F.rank, F)

@doc Markdown.doc"""
    segre(F::AbsBundle)
Compute the total Segre class."""
segre(F::AbsBundle) = inv(chern(F))

@doc Markdown.doc"""
    segre(k::Int, F::AbsBundle)
Compute the $k$-th Segre class."""
segre(k::Int, F::AbsBundle) = segre(F)[k]

@doc Markdown.doc"""
    todd(F::AbsBundle)
Compute the Todd class."""
todd(F::AbsBundle) = _todd(ch(F))

@doc Markdown.doc"""
    pontryagin(F::AbsBundle)
Compute the total Pontryagin class."""
function pontryagin(F::AbsBundle)
  n = F.parent.dim
  x = chern(F) * chern(dual(F))
  comps = x[0:n]
  sum([(-1)^i*comps[2i+1] for i in 0:n÷2])
end

@doc Markdown.doc"""
    pontryagin(k::Int, F::AbsBundle)
Compute the $k$-th Pontryagin class."""
pontryagin(k::Int, F::AbsBundle) = pontryagin(F)[2k]

@doc Markdown.doc"""
    chi(F::AbsBundle)
    chi(F::AbsBundle, G::AbsBundle)
Compute the holomorphic Euler characteristic $\chi(F)$, or the Euler pairing
$\chi(F,G)$.
"""
chi(F::AbsBundle) = integral(ch(F) * todd(F.parent)) # Hirzebruch-Riemann-Roch
chi(F::AbsBundle, G::AbsBundle) = begin
  F, G = _coerce(F, G)
  integral(ch(dual(F)) * ch(G) * todd(F.parent))
end

###############################################################################
#
# AbsVarietyHom
#
@doc Markdown.doc"""
    hom(X::AbsVariety, Y::AbsVariety, fˣ::Vector)
    hom(X::AbsVariety, Y::AbsVariety, fˣ::Vector, fₓ)

Construct a variety morphism from $X$ to $Y$, by specifying the pullbacks of
the generators of the Chow ring of $Y$. The pushforward can be automatically
computed in certain cases.

In case of an inclusion $i:X\hookrightarrow Y$ where the class of $X$ is not
present in the Chow ring of $Y$, use the argument `inclusion=true`.
A copy of $Y$ will be created, with extra classes added so that one can
pushforward classes on $X$.
"""
function hom(X::AbsVariety, Y::AbsVariety, fˣ::Vector, fₓ=nothing; inclusion::Bool=false, symbol::String="x")
  AbsVarietyHom(X, Y, fˣ, fₓ)
  # !inclusion && return AbsVarietyHom(X, Y, fˣ, fₓ)
  # _inclusion(AbsVarietyHom(X, Y, fˣ), symbol=symbol)
end

@doc Markdown.doc"""
    dim(f::AbsVarietyHom)
Return the relative dimension."""
dim(f::AbsVarietyHom) = f.dim

@doc Markdown.doc"""
    tangent_bundle(f::AbsVarietyHom)
Return the relative tangent bundle."""
tangent_bundle(f::AbsVarietyHom) = f.T

@doc Markdown.doc"""
    cotangent_bundle(f::AbsVarietyHom)
Return the relative cotangent bundle."""
cotangent_bundle(f::AbsVarietyHom) = dual(f.T)

@doc Markdown.doc"""
    todd(f::AbsVarietyHom)
Compute the Todd class of the relative tangent bundle."""
todd(f::AbsVarietyHom) = todd(f.T)

@doc Markdown.doc"""
    pullback(f::AbsVarietyHom, x::MPolyDecRingElem)
    pullback(f::AbsVarietyHom, F::AbsBundle)
Compute the pullback of a Chow ring element $x$ or a bundle $F$ by a morphism $f$.
"""
pullback(f::AbsVarietyHom, x::MPolyDecRingOrQuoElem) = f.pullback(x)
pullback(f::AbsVarietyHom, F::AbsBundle) = AbsBundle(f.domain, f.pullback(ch(F)))

@doc Markdown.doc"""
    pushforward(f::AbsVarietyHom, x::MPolyDecRingElem)
    pushforward(f::AbsVarietyHom, F::AbsBundle)
Compute the pushforward of a Chow ring element $x$ or a bundle $F$ by a
morphism $f$. For abstract bundles, the pushforward is derived, e.g., for a
bundle $F$ it is understood as the alternating sum of all direct images.
"""
pushforward(f::AbsVarietyHom, x::MPolyDecRingOrQuoElem) = f.pushforward(x)
pushforward(f::AbsVarietyHom, F::AbsBundle) = AbsBundle(f.codomain, f.pushforward(ch(F) * todd(f))) # Grothendieck-Hirzebruch-Riemann-Roch

function identity_hom(X::V) where V <: AbsVarietyT
  AbsVarietyHom(X, X, gens(X.ring), map_from_func(identity, X.ring, X.ring))
end

@doc Markdown.doc"""
    *(f::AbsVarietyHom, g::AbsVarietyHom)
Construct the composition morphism $g\circ f: X\to Z$ for $f: X\to Y$ and $g:Y\to Z$.
"""
function *(f::AbsVarietyHom, g::AbsVarietyHom)
  X, Y = f.domain, f.codomain
  @assert g.domain == Y
  Z = g.codomain
  gofₓ = nothing
  if isdefined(f, :pushforward) && isdefined(g, :pushforward)
    gofₓ = map_from_func(g.pushforward ∘ f.pushforward, X.ring, Z.ring)
  end
  gof = AbsVarietyHom(X, Z, g.pullback * f.pullback, gofₓ)
  return gof
end

###############################################################################
#
# AbsVariety
#
# generic variety with some classes in given degrees
@doc Markdown.doc"""
    variety(n::Int, symbols::Vector{String}, degs::Vector{Int})

Construct a generic variety of dimension $n$ with some classes in given degrees.

Return the variety and the list of classes.
"""
function variety(n::Int, symbols::Vector{String}, degs::Vector{Int}; base::Ring=QQ)
  @assert length(symbols) > 0
  R, x = grade(PolynomialRing(base, symbols)[1], degs)
  return AbsVariety(n, R), x
end

# generic variety with some bundles in given ranks
@doc Markdown.doc"""
    variety(n::Int, bundles::Vector{Pair{Int, T}}) where T

Construct a generic variety of dimension $n$ with some bundles of given ranks.

Return the variety and the list of bundles.
"""
function variety(n::Int, bundles::Vector{Pair{Int, T}}; base::Ring=QQ) where T
  symbols = vcat([_parse_symbol(s,1:r) for (r,s) in bundles]...)
  degs = vcat([collect(1:r) for (r,s) in bundles]...)
  X = variety(n, symbols, degs, base=base)[1]
  i = 1
  X.bundles = AbsBundle[]
  for (r,s) in bundles
    push!(X.bundles, AbsBundle(X, r, 1 + sum(gens(X.ring)[i:i+r-1])))
    i += r
  end
  return X, X.bundles
end

# generic variety with tangent bundle
@doc Markdown.doc"""
    variety(n::Int)

Construct a generic variety of dimension $n$ and define its tangent bundle.

Return the variety.
"""
function variety(n::Int; base::Ring=QQ)
  n == 0 && return point()
  X, (T,) = variety(n, [n=>"c"], base=base)
  X.T = T
  return X
end

(X::AbsVariety)(f::RingElement) = X.ring(f)
gens(X::AbsVariety) = gens(X.ring)

@doc Markdown.doc"""
    OO(X::AbsVariety)
    OO(X::TnVariety)
Return the trivial bundle $\mathcal O_X$ on $X$.
"""
OO(X::AbsVariety) = AbsBundle(X, X(1))

@doc Markdown.doc"""
    OO(X::AbsVariety, n)
    OO(X::AbsVariety, D)
Return the line bundle $\mathcal O_X(n)$ on $X$ if $X$ has been given a
polarization, or a line bundle $\mathcal O_X(D)$ with first Chern class $D$.
"""
OO(X::AbsVariety, n::RingElement) = AbsBundle(X, 1, 1+n*X.O1)
OO(X::AbsVariety, D::MPolyDecRingElem) = AbsBundle(X, 1, 1+D[1])

@doc Markdown.doc"""
    degree(X::AbsVariety)
Compute the degree of $X$ with respect to its polarization."""
degree(X::AbsVariety) = integral(X.O1^X.dim)

@doc Markdown.doc"""
    tangent_bundle(X::AbsVariety)
    tangent_bundle(X::TnVariety)
Return the tangent bundle of a variety $X$. Same as `X.T`.
"""
tangent_bundle(X::AbsVariety) = X.T

@doc Markdown.doc"""
    cotangent_bundle(X::AbsVariety)
    cotangent_bundle(X::TnVariety)
Return the cotangent bundle of a variety $X$.
"""
cotangent_bundle(X::AbsVariety) = dual(X.T)

@doc Markdown.doc"""
    canonical_class(X::AbsVariety)
Return the canonical class of a variety $X$."""
canonical_class(X::AbsVariety) = -chern(1, X.T)

@doc Markdown.doc"""
    canonical_bundle(X::AbsVariety)
Return the canonical bundle of a variety $X$."""
canonical_bundle(X::AbsVariety) = det(cotangent_bundle(X))

@doc Markdown.doc"""
    bundles(X::AbsVariety)
    bundles(X::TnVariety)
Return the tautological bundles of a variety $X$. Same as `X.bundles`.
"""
bundles(X::AbsVariety) = X.bundles

@doc Markdown.doc"""
    chern(X::AbsVariety)
    chern(X::TnVariety)
Compute the total Chern class of the tangent bundle of $X$.
"""
chern(X::AbsVariety) = chern(X.T)

@doc Markdown.doc"""
    chern(k::Int, X::AbsVariety)
    chern(k::Int, X::TnVariety)
Compute the $k$-th Chern class of the tangent bundle of $X$.
"""
chern(k::Int, X::AbsVariety) = chern(k, X.T)

@doc Markdown.doc"""
    euler(X::AbsVariety)
    euler(X::TnVariety)
Compute the Euler number of a variety $X$.
"""
euler(X::AbsVariety) = integral(chern(X.T))

@doc Markdown.doc"""
    todd(X::AbsVariety)
Compute the Todd class of the tangent bundle of $X$."""
todd(X::AbsVariety) = todd(X.T)

@doc Markdown.doc"""
    pontryagin(X::AbsVariety)
Compute the total Pontryagin class of the tangent bundle of $X$."""
pontryagin(X::AbsVariety) = pontryagin(X.T)

@doc Markdown.doc"""
    pontryagin(k::Int, X::AbsVariety)
Compute the $k$-th Pontryagin class of the tangent bundle of $X$."""
pontryagin(k::Int, X::AbsVariety) = pontryagin(k, X.T)

chi(p::Int, X::AbsVariety) = chi(exterior_power(p, dual(X.T))) # generalized Todd genus

function todd_polynomial(n::Int)
  X = variety(n)
  R, z = X.ring["z"]
  sum(chi(p, X) * (z-1)^p for p in 0:n)
end

@doc Markdown.doc"""
    chern_number(X::AbsVariety, λ::Int...)
    chern_number(X::AbsVariety, λ::Vector{Int})
    chern_number(X::AbsVariety, λ::Partition)
Compute the Chern number $c_\lambda (X):=\int_X c_{\lambda_1}(X)\cdots
c_{\lambda_k}(X)$, where $\lambda:=(\lambda_1,\dots,\lambda_k)$ is a partition
of the dimension of $X$.
"""
chern_number(X::AbsVariety, λ::Int...) = chern_number(X, collect(λ))
chern_number(X::AbsVariety, λ::Partition) = chern_number(X, collect(λ))
function chern_number(X::AbsVariety, λ::Vector{Int})
  @assert sum(λ) == X.dim
  c = chern(X)[1:X.dim]
  integral(prod([c[i] for i in λ]))
end

@doc Markdown.doc"""
    chern_numbers(X::AbsVariety)
Compute all the Chern numbers of $X$ as a list of pairs $\lambda\Rightarrow
c_\lambda(X)$.
"""
function chern_numbers(X::AbsVariety)
  c = chern(X)[1:X.dim]
  [λ => integral(prod([c[i] for i in λ])) for λ in partitions(X.dim)]
end

for g in [:a_hat_genus, :l_genus]
  @eval function $g(k::Int, X::AbsVariety)
    R = X.ring
    k == 0 && return R(1)
    p = pontryagin(X.T)[1:2k]
    R isa MPolyDecRing && return R($g(k).f([p[2i].f for i in 1:k]...))
    R isa MPolyQuoRing && return R(base_ring(R)($g(k).f([p[2i].f.f for i in 1:k]...)))
  end
  @eval function $g(X::AbsVariety)
    !iseven(X.dim) && error("the variety is not of even dimension")
    integral($g(X.dim÷2, X))
  end
end

@doc Markdown.doc"""
    a_hat_genus(k::Int, X::AbsVariety)
Compute the $k$-th $\hat A$ genus of a variety $X$."""
a_hat_genus(k::Int, X::AbsVariety)

@doc Markdown.doc"""
    l_genus(k::Int, X::AbsVariety)
Compute the $k$-th L genus of a variety $X$."""
l_genus(k::Int, X::AbsVariety)

@doc Markdown.doc"""
    a_hat_genus(X::AbsVariety)
Compute the top $\hat A$ genus of a variety $X$ of even dimension."""
a_hat_genus(X::AbsVariety)

@doc Markdown.doc"""
    l_genus(X::AbsVariety)
Compute the top L genus of a variety $X$ of even dimension."""
l_genus(X::AbsVariety)

@doc Markdown.doc"""
    signature(X::AbsVariety)
Compute the signature of a variety $X$ of even dimension."""
signature(X::AbsVariety) = l_genus(X) # Hirzebruch signature theorem

@doc Markdown.doc"""
    hilbert_polynomial(F::AbsBundle)
    hilbert_polynomial(X::AbsVariety)
Compute the Hilbert polynomial of a bundle $F$ or the Hilbert polynomial of $X$
itself, with respect to the polarization $\mathcal O_X(1)$ on $X$.
"""
function hilbert_polynomial(F::AbsBundle)
  !isdefined(F.parent, :O1) && error("no polarization is specified for the variety")
  X, O1 = F.parent, F.parent.O1
  # extend the coefficient ring to QQ(t)
  # TODO should we use FunctionField here?
  Qt, t = PolynomialRing(QQ, "t")
  @assert X.ring isa MPolyQuoRing
  R = parent(change_base_ring(Qt, base_ring(X.ring).R()))
  GR = grade(R, gradings(base_ring(X.ring)))[1]
  toR = x -> GR(change_base_ring(Qt, x, parent=R))
  I = ideal(toR.(gens(X.ring.I)))
  R_ = quo(GR, I)[1]
  set_attribute!(R_, :variety_dim => X.dim)
  ch_O_t = 1 + _logg(1 + t * R_(toR(O1.f)))
  ch_F = R_(toR(ch(F).f))
  td = R_(toR(todd(X).f))
  pt = R_(toR(X.point.f))
  hilb = constant_coefficient(div(simplify(ch_F * ch_O_t * td).f, simplify(pt).f))
  return hilb
end
hilbert_polynomial(X::AbsVariety) = hilbert_polynomial(OO(X))

# find canonically defined morphism from X to Y
function _hom(X::AbsVariety, Y::AbsVariety)
  X == Y && return identity_hom(X)
  # first handle the case where X is a (fibered) product
  projs = get_attribute(X, :projections)
  if projs !== nothing
    for p in projs
      p.codomain == Y && return p
    end
  else
    # follow the chain of structure maps to see if we can arrive at Y
    homs = AbsVarietyHom[]
    while isdefined(X, :struct_map) && X != Y
      push!(homs, X.struct_map)
      X = X.struct_map.codomain
    end
    X == Y && return reduce(*, homs)
  end
  error("no canonical homomorphism between the given varieties")
end

# morphisms for points are convenient, but are not desired when doing coercion
@doc Markdown.doc"""
    hom(X::AbsVariety, Y::AbsVariety)
Return a canonicallly defined morphism from $X$ to $Y$."""
function hom(X::AbsVariety, Y::AbsVariety)
  get_attribute(Y, :point) !== nothing && return hom(X, Y, [X(0)]) # Y is a point
  get_attribute(X, :point) !== nothing && return hom(X, Y, repeat([X(0)], length(gens(Y.ring)))) # X is a point
  _hom(X, Y)
end
→(X::AbsVariety, Y::AbsVariety) = hom(X, Y)

# product variety
@doc Markdown.doc"""
    *(X::AbsVariety, Y::AbsVariety)
Construct the product variety $X\times Y$. If both $X$ and $Y$ have a
polarization, $X\times Y$ will be endowed with the polarization of the Segre
embedding.
"""
function *(X::AbsVariety, Y::AbsVariety)
  prod_cache = get_attribute(X, :prod_cache)
  prod_cache !== nothing && Y in keys(prod_cache) && return prod_cache[Y]
  if prod_cache === nothing
    prod_cache = Dict{AbsVariety, AbsVariety}()
    set_attribute!(X, :prod_cache => prod_cache)
  end
  @assert X.base == Y.base
  base = X.base
  A, B = X.ring, Y.ring
  symsA, symsB = string.(gens(A)), string.(gens(B))
  a = length(symsA)
  R, x = grade(PolynomialRing(base, vcat(symsA, symsB))[1], vcat(gradings(A), gradings(B)))
  # TODO: fails with check = true
  AtoR = hom(A, R, x[1:a], check = false)
  BtoR = hom(B, R, x[a+1:end], check = false)
  IA = ideal(A isa MPolyQuoRing ? AtoR.(A.(gens(A.I))) : [R()])
  IB = ideal(B isa MPolyQuoRing ? BtoR.(B.(gens(B.I))) : [R()])
  AXY, _ = quo(R, IA + IB)
  XY = AbsVariety(X.dim+Y.dim, AXY)
  if isdefined(X, :point) && isdefined(Y, :point)
    XY.point = XY(AtoR(X.point) * BtoR(Y.point))
  end
  p = AbsVarietyHom(XY, X, XY.(x[1:a]))
  q = AbsVarietyHom(XY, Y, XY.(x[a+1:end]))
  if isdefined(X, :T) && isdefined(Y, :T)
    XY.T = pullback(p, X.T) + pullback(q, Y.T)
  end
  if isdefined(X, :O1) && isdefined(Y, :O1) # Segre embedding
    XY.O1 = p.pullback(X.O1) + q.pullback(Y.O1)
  end
  if get_attribute(X, :alg) == true && get_attribute(Y, :alg) == true
    set_attribute!(XY, :alg => true)
  end
  set_attribute!(XY, :projections => [p, q])
  set_attribute!(XY, :description => "Product of $X and $Y")
  prod_cache[Y] = XY
  return XY
end

@doc Markdown.doc"""
    graph(f::AbsVarietyHom)
Given a morphism $f: X\to Y$, construct $i:\Gamma_f\to X\times Y$, the
inclusion of the graph into the product.
"""
function graph(f::AbsVarietyHom)
  X, Y = f.domain, f.codomain
  hom(X, X * Y, vcat(gens(X), f.pullback.image))
end

###############################################################################
#
# Operators on AbsBundle
#
function adams(k::Int, x::MPolyDecRingOrQuoElem)
  R = parent(x)
  n = get_attribute(R, :variety_dim)
  comps = x[0:n]
  sum([ZZ(k)^i*comps[i+1] for i in 0:n])
end

@doc Markdown.doc"""
    dual(F::AbsBundle)
    dual(F::TnBundle)
Return the dual bundle.
"""
function dual(F::AbsBundle)
  Fdual = AbsBundle(F.parent, adams(-1, ch(F)))
  if isdefined(F, :chern)
    Fdual.chern = adams(-1, chern(F))
  end
  return Fdual
end
+(n::RingElement, F::AbsBundle) = AbsBundle(F.parent, n + ch(F))
*(n::RingElement, F::AbsBundle) = AbsBundle(F.parent, n * ch(F))
+(F::AbsBundle, n::RingElement) = n + F
*(F::AbsBundle, n::RingElement) = n * F
-(F::AbsBundle) = AbsBundle(F.parent, -ch(F))
^(F::AbsBundle, n::Int) = AbsBundle(F.parent, ch(F)^n)

@doc Markdown.doc"""
    det(F::AbsBundle)
    det(F::TnBundle)
Return the determinant bundle.
"""
det(F::AbsBundle) = AbsBundle(F.parent, 1, 1 + chern(1, F))
function _coerce(F::AbsBundle, G::AbsBundle)
  X, Y = F.parent, G.parent
  X == Y && return F, G
  try
    return F, pullback(_hom(X, Y), G)
  catch
    try
      return pullback(_hom(Y, X), F), G
    catch
      error("the sheaves are not on compatible varieties")
    end
  end
end

for O in [:(+), :(-), :(*)]
  @eval ($O)(F::AbsBundle, G::AbsBundle) = (
    (F, G) = _coerce(F, G);
    AbsBundle(F.parent, $O(ch(F), ch(G))))
end
hom(F::AbsBundle, G::AbsBundle) = dual(F) * G

@doc Markdown.doc"""
    exterior_power(k::Int, F::AbsBundle)
    exterior_power(k::Int, F::TnBundle)
Return the $k$-th exterior power.
"""
function exterior_power(k::Int, F::AbsBundle)
  AbsBundle(F.parent, _wedge(k, ch(F))[end])
end

function exterior_power(F::AbsBundle)
  AbsBundle(F.parent, sum([(-1)^(i-1) * w for (i, w) in enumerate(_wedge(F.rank, ch(F)))]))
end

@doc Markdown.doc"""
    symmetric_power(k, F::AbsBundle)
    symmetric_power(k::Int, F::TnBundle)
Return the $k$-th symmetric power. For an `AbsBundle`, $k$ can contain parameters.
"""
function symmetric_power(k::Int, F::AbsBundle)
  AbsBundle(F.parent, _sym(k, ch(F))[end])
end

function symmetric_power(k::RingElement, F::AbsBundle)
  X = F.parent
  PF = proj(dual(F))
  p = PF.struct_map
  AbsBundle(X, p.pushforward(sum((ch(OO(PF, k)) * todd(p))[0:PF.dim])))
end

@doc Markdown.doc"""
    schur_functor(λ::Vector{Int}, F::AbsBundle)
    schur_functor(λ::Partition, F::AbsBundle)
Return the result of the Schur functor $\mathbf S^\lambda$.
"""
function schur_functor(λ::Vector{Int}, F::AbsBundle) schur_functor(Partition(λ), F) end
function schur_functor(λ::Partition, F::AbsBundle)
  λ = conj(λ)
  X = F.parent
  w = _wedge(sum(λ), ch(F))
  S, ei = PolynomialRing(QQ, ["e$i" for i in 1:length(w)])
  e = i -> i < 0 ? S() : ei[i+1]
  M = [e(λ[i]-i+j) for i in 1:length(λ), j in 1:length(λ)]
  sch = det(matrix(S, M)) # Jacobi-Trudi
  R = X.ring
  if R isa MPolyQuoRing
    StoX = hom(S, R.R.R, [wi.f.f for wi in w])
    return AbsBundle(X, X(R.R(StoX(sch))))
  else
    StoX = hom(S, R.R, [wi.f for wi in w])
    return AbsBundle(X, X(StoX(sch)))
  end
end
function giambelli(λ::Vector{Int}, F::AbsBundle)
  R = F.parent.ring
  M = [chern(λ[i]-i+j, F).f for i in 1:length(λ), j in 1:length(λ)]
  R(det(matrix(R, M)))
end

###############################################################################
#
# Various computations
#
@doc Markdown.doc"""
    basis(X::AbsVariety)
Return an additive basis of the Chow ring of $X$, grouped by increasing
degree (i.e., increasing codimension)."""
function basis(X::AbsVariety)
  # it is important for this to be cached!
  return get_attribute!(X, :basis) do
    R = X.ring
    try_trim = "Try use `trim!`."
    !(R isa MPolyQuoRing) && error("the ring has no ideal. "*try_trim)
    dim(R.I) > 0 && error("the ideal is not 0-dimensional. "*try_trim)
    b = Oscar._kbase(R)
    ans = [MPolyQuoRingElem[] for i in 0:X.dim]
    for bi in b
      push!(ans[_total_degree(bi)+1], R(bi))
    end
    return ans
  end
end

@doc Markdown.doc"""
    basis(k::Int, X::AbsVariety)
Return an additive basis of the Chow ring of $X$ in codimension $k$."""
basis(k::Int, X::AbsVariety) = basis(X)[k+1]

@doc Markdown.doc"""
    betti(X::AbsVariety)
Return the Betti numbers of the Chow ring of $X$. Note that these are not
necessarily equal to the usual Betti numbers, i.e., the dimensions of
(co)homologies."""
betti(X::AbsVariety) = length.(basis(X))

@doc Markdown.doc"""
    integral(x::MPolyDecRingElem)

Compute the integral of a Chow ring element.

If the variety $X$ has a (unique) point class `X.point`, the integral will be a
number (an `fmpq` or a function field element). Otherwise the 0-dimensional
part of $x$ is returned.
"""
function integral(x::MPolyDecRingOrQuoElem)
  X = get_attribute(parent(x), :variety)
  if isdefined(X, :point) && length(basis(X.dim, X)) == 1
    return constant_coefficient(div(simplify(x).f, simplify(X.point).f))
  else
    return x[X.dim]
  end
end

@doc Markdown.doc"""
    intersection_matrix(a::Vector)
    intersection_matrix(a::Vector, b::Vector)
    intersection_matrix(X::AbsVariety)
Compute the intersection matrix among entries of a vector $a$ of Chow ring
elements, or between two vectors $a$ and $b$. For a variety $X$, this computes
the intersection matrix of the additive basis given by `basis(X)`.
"""
function intersection_matrix(X::AbsVariety) intersection_matrix(vcat(basis(X)...)) end
function intersection_matrix(a::Vector{}, b=nothing)
  if b === nothing b = a end
  matrix([integral(ai*bj) for ai in a, bj in b])
end

@doc Markdown.doc"""
    dual_basis(k::Int, X::AbsVariety)
Compute the dual basis of the additive basis in codimension $k$ given by
`basis(k, X)` (the returned elements are therefore in codimension
$\dim X-k$)."""
function dual_basis(k::Int, X::AbsVariety)
  d = get_attribute!(X, :dual_basis) do
    d = Dict{Int, Vector{elem_type(X.ring)}}()
    return d
  end
  if !(k in keys(d))
    B = basis(X)
    b_k = B[k+1]
    b_comp = B[X.dim-k+1]
    M = Matrix(inv(intersection_matrix(b_comp, b_k)))
    d[k] = M * b_comp
    d[X.dim-k] = transpose(M) * b_k
  end
  return d[k]
end

@doc Markdown.doc"""
    dual_basis(X::AbsVariety)
Compute the dual basis with respect to the additive basis given by `basis(X)`,
grouped by decreasing degree (i.e., decreasing codimension)."""
dual_basis(X::AbsVariety) = [dual_basis(k, X) for k in 0:X.dim]

# the parameter for truncation is usually the dimension, but can also be set
# manually, which is used when computing particular Chern classes (without
# computing the total Chern class)
function _expp(x::MPolyDecRingOrQuoElem; truncate::Int=-1)
  R = parent(x)
  n = truncate < 0 ? get_attribute(R, :variety_dim) : truncate
  comps = x[0:n]
  p = [(-1)^i * factorial(ZZ(i)) * comps[i+1] for i in 0:n]
  e = repeat([R(0)], n+1)
  e[1] = R(1)
  for i in 1:n
    e[i+1] = QQ(-1, i) * sum(p[j+1] * e[i-j+1] for j in 1:i)
  end
  simplify(sum(e))
end

function _logg(x::MPolyDecRingOrQuoElem)
  R = parent(x)
  n = get_attribute(R, :variety_dim)
  n == 0 && return R()
  e = x[1:n]
  p = pushfirst!(repeat([R()], n-1), -e[1])
  for i in 1:n-1
    p[i+1] = -ZZ(i+1)*e[i+1] - sum(e[j] * p[i-j+1] for j in 1:i)
  end
  simplify(sum((-1)^i//factorial(ZZ(i))*p[i] for i in 1:n))
end

# returns all the wedge from 0 to k
function _wedge(k::Int, x::MPolyDecRingOrQuoElem)
  R = parent(x)
  k == 0 && return [R(1)]
  n = get_attribute(R, :variety_dim)
  wedge = repeat([R(0)], k+1)
  wedge[1] = R(1)
  wedge[2] = x
  for j in 2:k
    wedge[j+1] = 1//ZZ(j) * sum(sum((-1)^(j-i+1) * wedge[i+1] * adams(j-i, x) for i in 0:j-1)[0:n])
  end
  wedge
end

# returns all the sym from 0 to k
function _sym(k::Int, x::MPolyDecRingOrQuoElem)
  R = parent(x)
  k == 0 && return [R(1)]
  n = get_attribute(R, :variety_dim)
  r = min(k, Int(ZZ(QQ(constant_coefficient(x.f)))))
  wedge = _wedge(r, x)
  sym = repeat([R(0)], k+1)
  sym[1] = R(1)
  sym[2] = x
  for j in 2:k
    sym[j+1] = sum(sum((-1)^(i+1) * wedge[i+1] * sym[j-i+1] for i in 1:min(j,r))[0:n])
  end
  sym
end

function _genus(x::MPolyDecRingOrQuoElem, taylor::Vector{})
  R = parent(x)
  iszero(x) && return R(1)
  n = get_attribute(R, :variety_dim)
  R, (t,) = grade(PolynomialRing(QQ, ["t"])[1])
  set_attribute!(R, :variety_dim, n)
  lg = _logg(R(sum(taylor[i+1] * t^i for i in 0:n)))
  comps = lg[1:n]
  lg = [leading_coefficient(comps[i].f) for i in 1:n]
  comps = x[1:n]
  _expp(sum(factorial(ZZ(i)) * lg[i] * comps[i] for i in 1:n))
end

function _todd(x::MPolyDecRingOrQuoElem)
  n = get_attribute(parent(x), :variety_dim)
  # the Taylor series of t/(1-exp(-t))
  taylor = [(-1)^i//factorial(ZZ(i))*bernoulli(i) for i in 0:n]
  _genus(x, taylor)
end

function _l_genus(x::MPolyDecRingOrQuoElem)
  n = get_attribute(parent(x), :variety_dim)
  # the Taylor series of sqrt(t)/tanh(sqrt(t))
  taylor = [ZZ(2)^2i//factorial(ZZ(2i))*bernoulli(2i) for i in 0:n]
  _genus(x, taylor)
end

function _a_hat_genus(x::MPolyDecRingOrQuoElem)
  n = get_attribute(parent(x), :variety_dim)
  # the Taylor series of (sqrt(t)/2)/sinh(sqrt(t)/2)
  R, t = PowerSeriesRing(QQ, 2n+1, "t")
  s = divexact(t, exp(QQ(1//2)*t)-exp(-QQ(1//2)*t))
  taylor = [coeff(s, 2i) for i in 0:n]
  _genus(x, taylor)
end

for (g,s) in [:a_hat_genus=>"p", :l_genus=>"p", :todd=>"c"]
  _g = Symbol("_", g)
  @eval function $g(n::Int)
    n == 0 && return QQ(1)
    R, p = grade(PolynomialRing(QQ, _parse_symbol($s, 1:n))[1], collect(1:n))
    set_attribute!(R, :variety_dim, n)
    $_g(_logg(R(1+sum(p))))[n]
  end
end

@doc Markdown.doc"""
    todd(n::Int)
Compute the (generic) $n$-th Todd genus."""
todd(n::Int)

@doc Markdown.doc"""
    l_genus(n::Int)
Compute the (generic) $n$-th L genus."""
l_genus(n::Int)

@doc Markdown.doc"""
    a_hat_genus(n::Int)
Compute the (generic) $n$-th $\hat A$ genus."""
a_hat_genus(n::Int)

@doc Markdown.doc"""
    section_zero_locus(F::AbsBundle)

Construct the zero locus of a general section of a bundle $F$.

Use the argument `class=true` to only compute the class of the zero locus (same
as `ctop(F)`).
"""
function section_zero_locus(F::AbsBundle; class::Bool=false)
  X = parent(F)
  R = X.ring
  cZ = ctop(F)
  # return only the class of Z in the chow ring of X
  class && return cZ
  if R isa MPolyQuoRing
    I = quotient(modulus(R), ideal(base_ring(R), [cZ.f]))
    AZ = quo(base_ring(R), I)[1]
  else
    AZ = R
  end
  Z = AbsVariety(X.dim - F.rank, AZ)
  if isdefined(X, :point)
    ps = basis(Z.dim, Z) # the 0-cycles
    @assert length(ps) == 1 # make sure that the 0-cycle is unique
    p = ps[1]
    degp = integral(R(p.f) * cZ) # compute the degree of iₓp
    Z.point = Z(inv(degp) * p.f)
  end
  if isdefined(X, :T)
    Z.T = AbsBundle(Z, Z((ch(X.T) - ch(F)).f))
  end
  if isdefined(X, :O1)
    Z.O1 = Z(X.O1.f)
  end
  iₓ = x -> x.f * cZ
  iₓ = map_from_func(iₓ, Z.ring, X.ring)
  @assert R isa MPolyQuoRing
  i = AbsVarietyHom(Z, X, Z.(gens(base_ring(R))), iₓ)
  i.T = pullback(i, -F)
  Z.struct_map = i
  set_attribute!(Z, :description, "Zero locus of a section of $F")
  return Z
end

@doc Markdown.doc"""
    complete_intersection(X::AbsVariety, degs::Int...)
    complete_intersection(X::AbsVariety, degs::Vector{Int})
Construct the complete intersection in $X$ of general hypersurfaces with
degrees $d_1,\dots,d_k$.
"""
complete_intersection(X::AbsVariety, degs::Int...) = complete_intersection(X, collect(degs))
complete_intersection(X::AbsVariety, degs::Vector{Int}) = (
  Y = section_zero_locus(sum(OO(X, d) for d in degs));
  set_attribute!(Y, :description => "Complete intersection of degree $(tuple(degs...)) in $X");
  Y)

@doc Markdown.doc"""
    degeneracy_locus(k::Int, F::AbsBundle, G::AbsBundle)

Construct the $k$-degeneracy locus for a general bundle map from $F$ to $G$.

Use the argument `class=true` to only compute the class of the degeneracy locus.
"""
function degeneracy_locus(k::Int, F::AbsBundle, G::AbsBundle; class::Bool=false)
  F, G = _coerce(F, G)
  m, n = rank(F), rank(G)
  @assert k < min(m,n)
  if class
    # return only the class of D in the chow ring of X
    if (m-k)*(n-k) <= F.parent.dim # Porteous' formula
      return ch(schur_functor(repeat([m-k], n-k), G-F))[(m-k)*(n-k)]
    else # expected dimension is negative
      return F.parent.ring(0)
    end
  end
  Gr = (m-k == 1) ? proj(F) : flag(m-k, F)
  S = Gr.bundles[1]
  D = section_zero_locus(dual(S) * G)
  D.struct_map = D → F.parent # skip the flag variety
  set_attribute!(D, :description, "Degeneracy locus of rank $k from $F to $G")
  return D
end

###############################################################################
@doc Markdown.doc"""
    point()
Construct a point as an abstract variety."""
function point(; base::Ring=QQ)
  R, (p,) = grade(PolynomialRing(base, ["p"])[1])
  I = ideal([p])
  pt = AbsVariety(0, quo(R, I)[1])
  pt.point = pt(1)
  pt.T = AbsBundle(pt, pt(0))
  pt.O1 = pt(0)
  set_attribute!(pt, :description, "Point")
  set_attribute!(pt, :point, true)
  return pt
end

@doc Markdown.doc"""
    proj(n::Int)
Construct an abstract projective space of dimension $n$, parametrizing
1-dimensional *subspaces* of a vector space of dimension $n+1$."""
function proj(n::Int; base::Ring=QQ, symbol::String="h")
  R, (h,) = grade(PolynomialRing(base, [symbol])[1])
  I = ideal([h^(n+1)])
  AP = quo(R, I)[1]
  P = AbsVariety(n, AP)
  h = P(h)
  P.point = h^n
  P.O1 = h
  chTP = P(n)
  for i in 1:n chTP += ZZ(n+1)//factorial(ZZ(i))*h^i end
  P.T = AbsBundle(P, chTP)
  P.T.chern = (1+h)^(n+1)
  S = AbsBundle(P, 1, 1-h)
  Q = OO(P)*(n+1) - S
  P.bundles = [S, Q]
  P.struct_map = hom(P, point(base=base), [P(1)])
  set_attribute!(P, :description => "Projective space of dim $n")
  set_attribute!(P, :grassmannian => :absolute)
  set_attribute!(P, :alg => true)
  return P
end

@doc Markdown.doc"""
    proj(F::AbsBundle)
Construct the projectivization of a bundle $F$, parametrizing 1-dimensional
*subspaces*."""
function proj(F::AbsBundle; symbol::String="h")
  X, r = F.parent, F.rank
  !(r isa Int) && error("expect rank to be an integer")
  R = X.ring
  syms = vcat([symbol], string.(gens(R)))
  # FIXME add product ordering
  # ord = ordering_dp(1) * R.R.ord
  # construct the ring
  w = vcat([1], gradings(R))
  R1, (h,) = grade(PolynomialRing(X.base, syms)[1], w)
  # TODO: why does this fail with check = true
  pback = hom(R, R1, gens(R1)[2:end], check = false)
  pfwd = hom(R1, R, pushfirst!(gens(R), R()))
  # construct the ideal
  rels = [sum(pback(chern(i, F)) * h^(r-i) for i in 0:r)]
  if R isa MPolyQuoRing rels = vcat(pback.(R.(gens(R.I))), rels) end
  APF = quo(R1, ideal(rels))[1]
  h = APF(h)
  # construct the variety
  PF = AbsVariety(X.dim+r-1, APF)
  pₓ = x -> X(pfwd(div(simplify(x).f, simplify(PF(h^(r-1))).f)))
  pₓ = map_from_func(pₓ, PF.ring, X.ring)
  p = AbsVarietyHom(PF, X, PF.(gens(R1)[2:end]), pₓ)
  if isdefined(X, :point) PF.point = p.pullback(X.point) * h^(r-1) end
  p.O1 = PF(h)
  PF.O1 = PF(h)
  S = AbsBundle(PF, 1, 1-h)
  Q = pullback(p, F) - S
  p.T = dual(S)*Q
  if isdefined(X, :T) PF.T = pullback(p, X.T) + p.T end
  PF.bundles = [S, Q]
  PF.struct_map = p
  set_attribute!(PF, :description => "Projectivization of $F")
  set_attribute!(PF, :grassmannian => :relative)
  return PF
end

@doc Markdown.doc"""
    grassmannian(k::Int, n::Int; bott::Bool=false)

Construct an abstract Grassmannian $\mathrm{Gr}(k, n)$, parametrizing
$k$-dimensional subspaces of an $n$-dimensional vector space.

Use the argument `bott=true` to construct the Grassmannian as a `TnVariety` for
computing integrals using Bott's formula.
"""
function grassmannian(k::Int, n::Int; bott::Bool=false, weights=:int, base::Ring=QQ, symbol::String="c")
  # combine the interface for AbsVariety and TnVariety versions
  bott && return tn_grassmannian(k, n, weights=weights)
  abs_grassmannian(k, n, base=base, symbol=symbol)
end

function abs_grassmannian(k::Int, n::Int; base::Ring=QQ, symbol::String="c")
  @assert k < n
  d = k*(n-k)
  R, c = grade(PolynomialRing(base, _parse_symbol(symbol, 1:k))[1], collect(1:k))
  inv_c = sum((-sum(c))^i for i in 1:n) # this is c(Q) since c(S)⋅c(Q) = 1
  # Q is of rank n-k: the vanishing of Chern classes in higher degrees provides all the relations for the Chow ring
  AGr = quo(R, ideal(inv_c[n-k+1:n]))[1]
  c = gens(AGr)
  Gr = AbsVariety(d, AGr)
  Gr.O1 = Gr(-c[1])
  S = AbsBundle(Gr, k, 1 + sum(c))
  Q = OO(Gr)*n - S
  Q.chern = 1 + Gr(sum(inv_c[1:n-k]))
  Gr.point = Gr((-1)^d*c[end]^(n-k))
  Gr.T = dual(S) * Q
  Gr.bundles = [S, Q]
  Gr.struct_map = hom(Gr, point(base=base), [Gr(1)])
  set_attribute!(Gr, :description => "Grassmannian Gr($k, $n)")
  set_attribute!(Gr, :grassmannian => :absolute)
  set_attribute!(Gr, :alg => true)
  return Gr
end

@doc Markdown.doc"""
    flag(dims::Int...; bott::Bool=false)
    flag(dims::Vector{Int}; bott::Bool=false)

Construct an abstract flag variety $\mathrm{Fl}(d_1,\dots,d_k)$, parametrizing
flags of subspaces $V_{d_1}\subset V_{d_2}\subset\cdots\subset V_{d_k}=V$.

Use the argument `bott=true` to construct the flag variety as a `TnVariety` for
computing integrals using Bott's formula.
"""
function flag(dims::Int...; bott::Bool=false, weights=:int, base::Ring=QQ, symbol::String="c")
  # combine the interface for AbsVariety and TnVariety versions
  bott && return tn_flag(collect(dims), weights=weights)
  abs_flag(collect(dims), base=base, symbol=symbol)
end

function flag(dims::Vector{Int}; bott::Bool=false, weights=:int, base::Ring=QQ, symbol::String="c")
  bott && return tn_flag(dims, weights=weights)
  abs_flag(dims, base=base, symbol=symbol)
end

function abs_flag(dims::Vector{Int}; base::Ring=QQ, symbol::String="c")
  n, l = dims[end], length(dims)
  ranks = pushfirst!([dims[i+1]-dims[i] for i in 1:l-1], dims[1])
  @assert all(r->r>0, ranks)
  d = sum(ranks[i] * sum(dims[end]-dims[i]) for i in 1:l-1)
  syms = vcat([_parse_symbol(symbol, i, 1:r) for (i,r) in enumerate(ranks)]...)
  # FIXME ordering
  # ord = prod(ordering_dp(r) for r in ranks)
  R = grade(PolynomialRing(base, syms)[1], vcat([collect(1:r) for r in ranks]...))[1]
  c = pushfirst!([1+sum(gens(R)[dims[i]+1:dims[i+1]]) for i in 1:l-1], 1+sum(gens(R)[1:dims[1]]))
  gi = prod(c)[0:n]
  # XXX cannot mod using graded ring element
  Rx, x = R.R["x"]
  g = sum(gi[i+1].f * x^(n-i) for i in 0:n)
  q = mod(x^n, g)
  rels = [R(coeff(q, i)) for i in 0:n-1]
  AFl = quo(R, ideal(rels))[1]
  c = AFl.(c)
  Fl = AbsVariety(d, AFl)
  Fl.bundles = [AbsBundle(Fl, r, ci) for (r,ci) in zip(ranks, c)]
  Fl.O1 = simplify(sum((i-1)*chern(1, Fl.bundles[i]) for i in 1:l))
  Fl.point = prod(ctop(E)^dims[i] for (i,E) in enumerate(Fl.bundles[2:end]))
  Fl.T = sum(dual(Fl.bundles[i]) * sum([Fl.bundles[j] for j in i+1:l]) for i in 1:l-1)
  Fl.struct_map = hom(Fl, point(base=base), [Fl(1)])
  set_attribute!(Fl, :description => "Flag variety Flag$(tuple(dims...))")
  if l == 2 set_attribute!(Fl, :grassmannian => :absolute) end
  set_attribute!(Fl, :alg => true)
  # if all(r->r==1, ranks)
  #   set_attribute!(Fl, :weyl_group => WeylGroup("A$(n-1)"))
  #   set_attribute!(Fl, :roots => -[c[i] - c[i+1] for i in 1:n-1])
  # end
  return Fl
end

@doc Markdown.doc"""
    flag(d::Int, F::AbsBundle)
    flag(dims::Vector{Int}, F::AbsBundle)
Construct the relative flag variety of a bundle $F$, parametrizing
flags of subspaces $V_{d_1}\subset V_{d_2}\subset\cdots\subset V_{d_k}$. The
last dimension (i.e., the rank of $F$) can be omitted.
"""
function flag(d::Int, F::AbsBundle; symbol::String="c") flag([d], F, symbol=symbol) end
function flag(dims::Vector{Int}, F::AbsBundle; symbol::String="c")
  X, n = F.parent, F.rank
  !(n isa Int) && error("expect rank to be an integer")
  # compute the ranks and relative dim
  l = length(dims)
  ranks = pushfirst!([dims[i+1]-dims[i] for i in 1:l-1], dims[1])
  @assert all(r->r>0, ranks) && dims[end] <= n
  if dims[end] < n # the last dim can be omitted
    dims = vcat(dims, [n])
    push!(ranks, n-dims[l])
    l += 1
  end
  d = sum(ranks[i] * sum(dims[end]-dims[i]) for i in 1:l-1)
  # FIXME ordering
  # construct the ring
  R = X.ring
  syms = vcat([_parse_symbol(symbol, i, 1:r) for (i,r) in enumerate(ranks)]..., string.(gens(R.R)))
  # ord = prod(ordering_dp(r) for r in ranks) * ordering(X.ring.R.R)
  w = vcat([collect(1:r) for r in ranks]..., gradings(R))
  R1 = grade(PolynomialRing(X.base, syms)[1], w)[1]
  pback = hom(R, R1, gens(R1)[n+1:end])
  pfwd = hom(R1, R, vcat(repeat([R()], n), gens(R)))
  
  # compute the relations
  c = pushfirst!([1+sum(gens(R1)[dims[i]+1:dims[i+1]]) for i in 1:l-1], 1+sum(gens(R1)[1:dims[1]]))
  # XXX cannot mod using graded ring element
  Rx, x = R1.R["x"]
  fi = pback(chern(F))[0:n]
  f = sum(fi[i+1].f * x^(n-i) for i in 0:n)
  gi = prod(c)[0:n]
  g = sum(gi[i+1].f * x^(n-i) for i in 0:n)
  rels = [R1(coeff(mod(f, g), i)) for i in 0:n-1]
  if R isa MPolyQuoRing rels = vcat(pback.(R.(gens(R.I))), rels) end
  AFl = quo(R1, ideal(rels))[1]
  c = AFl.(c)
  Fl = AbsVariety(X.dim + d, AFl)
  
  Fl.bundles = [AbsBundle(Fl, r, ci) for (r,ci) in zip(ranks, c)]
  section = prod(ctop(E)^sum(dims[i]) for (i, E) in enumerate(Fl.bundles[2:end]))
  if isdefined(X, :point)
    Fl.point = pback(X.point) * section
  end
  pˣ = Fl.(gens(R1)[n+1:end])
  pₓ = x -> (@warn("possibly wrong ans"); X(pfwd(div(simplify(x).f, simplify(section).f))))
  pₓ = map_from_func(pₓ, Fl.ring, X.ring)
  p = AbsVarietyHom(Fl, X, pˣ, pₓ)
  p.O1 = simplify(sum((i-1)*chern(1, Fl.bundles[i]) for i in 1:l))
  Fl.O1 = p.O1
  p.T = sum(dual(Fl.bundles[i]) * sum([Fl.bundles[j] for j in i+1:l]) for i in 1:l-1)
  if isdefined(X, :T)
    Fl.T = pullback(p, X.T) + p.T
  end
  Fl.struct_map = p
  set_attribute!(Fl, :description => "Relative flag variety Flag$(tuple(dims...)) for $F")
  set_attribute!(Fl, :section => section)
  if l == 2 set_attribute!(Fl, :grassmannian => :relative) end
  return Fl
end

@doc Markdown.doc"""
    schubert_class(G::AbsVariety, λ::Int...)
    schubert_class(G::AbsVariety, λ::Vector{Int})
    schubert_class(G::AbsVariety, λ::Partition)
Return the Schubert class $\sigma_\lambda$ on a (relative) Grassmannian $G$.
"""
function schubert_class(G::AbsVariety, λ::Int...) schubert_class(G, collect(λ)) end
function schubert_class(G::AbsVariety, λ::Partition) schubert_class(G, collect(λ)) end
function schubert_class(G::AbsVariety, λ::Vector{Int})
  get_attribute(G, :grassmannian) === nothing && error("the variety is not a Grassmannian")
  (length(λ) > rank(G.bundles[1]) || sort(λ, rev=true) != λ) && error("the Schubert input is not well-formed")
  giambelli(λ, G.bundles[2])
end

@doc Markdown.doc"""
    schubert_classes(m::Int, G::AbsVariety)
Return all the Schubert classes in codimension $m$ on a (relative) Grassmannian $G$."""
function schubert_classes(m::Int, G::AbsVariety)
  get_attribute(G, :grassmannian) === nothing && error("the variety is not a Grassmannian")
  S, Q = G.bundles
  [schubert_class(G, λ) for λ in partitions(m, rank(S), rank(Q))]
end
