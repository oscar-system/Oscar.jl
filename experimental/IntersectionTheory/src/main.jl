
###############################################################################
#
# AbstractBundle
#
@doc raw"""
    abstract_bundle(X::AbstractVariety, ch::Union{MPolyDecRingElem, MPolyQuoRingElem})
    abstract_bundle(X::AbstractVariety, r::RingElement, c::Union{MPolyDecRingElem, MPolyQuoRingElem})

Return an abstract vector bundle on `X` by specifying its Chern character `ch`.
Equivalently, specify its rank `r` and total Chern class `c`.

# Examples

We show two ways of constructing the Horrocks-Mumford bundle `F` [HM73](@cite).
First, we create `F` as the cohomology bundle of its Beilinson monad
(see Equation (2.1) in [HM73](@cite)):

$0 \rightarrow \mathcal O_{\mathbb P^4}^5(2)\rightarrow \bigwedge^2\mathrm{T}^*_{\mathbb P^4}(5)
\rightarrow \mathcal O_{\mathbb P^4}^5(3)\rightarrow 0.$

Then, we show the constructor above at work.

```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> A = 5*OO(P4, 2)
AbstractBundle of rank 5 on AbstractVariety of dim 4

julia> B = 2*exterior_power(cotangent_bundle(P4), 2)*OO(P4, 5)
AbstractBundle of rank 12 on AbstractVariety of dim 4

julia> C = 5*OO(P4, 3)
AbstractBundle of rank 5 on AbstractVariety of dim 4

julia> F = B - A - C
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> total_chern_class(F)
10*h^2 + 5*h + 1

julia> h = gens(P4)[1]
h

julia> F == abstract_bundle(P4, 2, 10*h^2 + 5*h + 1)
true

```
"""
abstract_bundle(X::AbstractVariety, ch::MPolyDecRingOrQuoElem) = AbstractBundle(X, ch)
abstract_bundle(X::AbstractVariety, r::RingElement, c::MPolyDecRingOrQuoElem) = AbstractBundle(
  X, r, c
)

#######################################################
@doc raw"""
    ==(F::AbstractBundle, G::AbstractBundle)

Return `true` if `F` is equal to `G`, and `false` otherwise.
This means that the underlying varieties of `F` and `G` are the same,
and that the Chern characters of `F` and `G` are equal.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> 3*OO(P2, 1) - OO(P2) == tangent_bundle(P2) # Euler sequence
true

```
"""
==(F::AbstractBundle, G::AbstractBundle) =
  parent(F) == parent(G) && chern_character(F) == chern_character(G)

function Base.hash(F::AbstractBundle, h::UInt)
  return hash(chern_character(F), h)
end

@doc raw"""
    chern_character(F::AbstractBundle)

Return the Chern character of `F`.

!!! note "Chern classes vs. Chern character"
    The **Chern classes** $\operatorname{c}_k(F)\in A^k(X)$ are the coefficients of the total
    Chern class $\operatorname{c}(F)=1+\operatorname{c}_1(F)+\operatorname{c}_2(F)+\cdots$, which is
    *multiplicative*: $\operatorname{c}(F\oplus G)=\operatorname{c}(F)\cdot\operatorname{c}(G)$.

    The **Chern character** $\operatorname{ch}(F)\in A^*(X)\otimes\mathbb Q$ is
    instead *additive* and *multiplicative*:
    $\operatorname{ch}(F\oplus G)=\operatorname{ch}(F)+\operatorname{ch}(G)$ and
    $\operatorname{ch}(F\otimes G)=\operatorname{ch}(F)\cdot\operatorname{ch}(G)$.
    It is related to the Chern classes by the Newton identity
    $\operatorname{ch}(F)=\operatorname{rk}(F)+\sum_{k\ge1}\frac{1}{k!}s_k$
    where $s_k$ is the $k$-th Newton polynomial in the Chern roots.
    In particular $\operatorname{ch}_0=\operatorname{rk}$,
    $\operatorname{ch}_1=\operatorname{c}_1$,
    $\operatorname{ch}_2=\frac12(\operatorname{c}_1^2-2\operatorname{c}_2)$.

    Because of these ring-homomorphism properties, the Chern character — not the total
    Chern class — is used as the internal representation of a bundle in this package
    (see `AbstractBundle`).

# Examples
```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 6

julia> chern_character(Q)
-1//2*c[1]^2 + 1//6*c[1]*c[2] - 1//24*c[1]*c[3] - c[1] + c[2] - 1//3*c[3] + 2

```
"""
function chern_character(F::AbstractBundle)
  if !isdefined(F, :ch)
    F.ch = rank(F) + _logg(total_chern_class(F))
  end
  return F.ch
end

@doc raw"""
    total_chern_class(F::AbstractBundle)

Return the total Chern class of `F`.

# Examples
```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 6

julia> total_chern_class(Q)
c[1]^2 - c[1] - c[2] + 1

julia> chern_class(Q, 1)
-c[1]

julia> chern_class(Q, 2)
c[1]^2 - c[2]

```
"""
function total_chern_class(F::AbstractBundle)
  if !isdefined(F, :chern)
    F.chern = _expp(chern_character(F))
  end
  return F.chern
end

@doc raw"""
    chern_class(F::AbstractBundle, k::Int)

Return the `k`-th Chern class of `F`.

# Examples
```jldoctest
julia> T, (d,) = polynomial_ring(QQ, [:d])
(Multivariate polynomial ring in 1 variable over QQ, QQMPolyRingElem[d])

julia> QT = fraction_field(T)
Fraction field
  of multivariate polynomial ring in 1 variable over QQ

julia> P4 = abstract_projective_space(4, base = QT)
AbstractVariety of dim 4

julia> h = gens(P4)[1]
h

julia> F = abstract_bundle(P4, 2, 10*h^2 + 5*h + 1) # Horrocks-Mumford bundle
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> chern_class(F*OO(P4, d), 1)
(2*d + 5)*h

julia> chern_class(F*OO(P4, d), 2)
(d^2 + 5*d + 10)*h^2

julia> chern_class(F*OO(P4, -3), 1)
-h

julia> chern_class(F*OO(P4, -3), 2)
4*h^2

```
"""
function chern_class(F::AbstractBundle, k::Int)
  if isdefined(F, :chern)
    return total_chern_class(F)[k]
  end
  # we don't compute `total_chern_class` here, in case it is an expensive computation
  # whether this is actually a useful choice is not supported by evidence at this point, so it may be changed in the future
  return _expp(chern_character(F), truncate=k)[k]
end

@doc raw"""
    top_chern_class(F::AbstractBundle)

Return the top Chern class of `F`.

# Examples
```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> h = gens(P4)[1]
h

julia> F = abstract_bundle(P4, 2, 10*h^2 + 5*h + 1)
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> top_chern_class(F)
10*h^2

```
"""
top_chern_class(F::AbstractBundle) = chern_class(F, rank(F))

@doc raw"""
    total_segre_class(F::AbstractBundle)

Return the total Segre class of `F`.

# Examples
```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 6

julia> C = total_chern_class(Q)
c[1]^2 - c[1] - c[2] + 1

julia> S = total_segre_class(Q)
c[1] + c[2] + c[3] + 1

julia> C*S
1

```
"""
total_segre_class(F::AbstractBundle) = inv(total_chern_class(F))

@doc raw"""
    segre_class(F::AbstractBundle, k::Int)

Return the `k`-th Segre class of `F`.

# Examples
```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 6

julia> segre_class(Q, 0)
1

julia> segre_class(Q, 1)
c[1]

julia> segre_class(Q, 2)
c[2]

julia> segre_class(Q, 3)
c[3]

```
"""
segre_class(F::AbstractBundle, k::Int) = total_segre_class(F)[k]

@doc raw"""
    todd_class(F::AbstractBundle)

Return the Todd class of `F`.

# Examples
```jldoctest
julia> P = abstract_projective_space(4, symbol = "H"); # Hartshorne, p. 433

julia> F = exterior_power(cotangent_bundle(P), 3)*OO(P, 3);

julia> G = OO(P, 1) + 4*OO(P);

julia> Z = degeneracy_locus(F, G, 3) # rational surface in P4
AbstractVariety of dim 2

julia> TZ = tangent_bundle(Z);

julia> K = canonical_class(Z)
z - H

julia> chern_class(TZ, 1) == -K
true

julia> tc = todd_class(TZ)
-1//2*z + 1//8*H^2 + 1//2*H + 1

julia> tc == 1 - 1//2*K + 1//12*(K^2 + chern_class(TZ, 2))
true

```
"""
todd_class(F::AbstractBundle) = _todd_class(chern_character(F))

@doc raw"""
    total_pontryagin_class(F::AbstractBundle)

Return the total Pontryagin class of `F`.

# Examples
```jldoctest
julia> S = complete_intersection(abstract_projective_space(3), 4);

julia> total_pontryagin_class(tangent_bundle(S))
-12*h^2 + 1
```
"""
function total_pontryagin_class(F::AbstractBundle)
  n = dim(parent(F))
  x = total_chern_class(F) * total_chern_class(dual(F))
  comps = x[0:n]
  sum([(-1)^i*comps[2i + 1] for i in 0:div(n, 2)])
end

@doc raw"""
    pontryagin_class(F::AbstractBundle, k::Int)

Return the `k`-th Pontryagin class of `F`.

# Examples
```jldoctest
julia> S = complete_intersection(abstract_projective_space(3), 4);

julia> pontryagin_class(tangent_bundle(S), 1)
-12*h^2
```
"""
pontryagin_class(F::AbstractBundle, k::Int) = total_pontryagin_class(F)[2k]

@doc raw"""
    euler_characteristic(F::AbstractBundle)

Return the Euler characteristic $\chi(F)$ computed via the Hirzebruch-Riemann-Roch formula.

# Examples
```jldoctest
julia> P = abstract_projective_space(4, symbol = "H"); # Hartshorne, p. 433

julia> F = exterior_power(cotangent_bundle(P), 3)*OO(P, 3);

julia> G = OO(P, 1) + 4*OO(P);

julia> Z = degeneracy_locus(F, G, 3) # rational surface in P4
AbstractVariety of dim 2

julia> TZ = tangent_bundle(Z);

julia> K = canonical_class(Z)
z - H

julia> H = polarization(Z)
H

julia> ec = euler_characteristic(OO(Z, H))
4

julia> ec == integral(1//2*H*(H - K) + 1//12*(K^2 + chern_class(TZ, 2)))
true

```
"""
Oscar.euler_characteristic(F::AbstractBundle) = integral(
  chern_character(F) * todd_class(parent(F))
) # Hirzebruch-Riemann-Roch

@doc raw"""
    euler_pairing(F::AbstractBundle, G::AbstractBundle)

Return the Euler pairing $\chi(F, G) = \int \operatorname{ch}(F^\vee) \cdot \operatorname{ch}(G) \cdot \operatorname{td}(X)$.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> euler_pairing(OO(P2, 1), OO(P2, 2))
3

```
"""
euler_pairing(F::AbstractBundle, G::AbstractBundle) = begin
  F, G = _coerce(F, G)
  euler_characteristic(dual(F) * G)
end

###############################################################################
#
# AbstractVarietyMap
#

### utility functions

function _trim(R::MPolyDecRingOrQuo)
  d = get_attribute(R, :abstract_variety_dim)
  isnothing(d) && error("ring is not the Chow ring of an abstract variety")

  if isdefined(R, :I)
    gI = gens(R.I) # will consist of the generators of the modulus I of the desired Chow ring
    S = base_ring(R)
  else
    gI = MPolyRingElem[] # will consist of the generators of the modulus I of the desired Chow ring
    S = R
  end

  w = weights(Int, S)
  if !isdefined(R, :I) || krull_dim(R.I) > 0
    # add powers of variables to gI so that I = <gI> becomes 0-dimensional
    for (i, x) in enumerate(gens(R))
      push!(gI, x^Int(ceil((d+1)//w[i])))
    end
  end
  # only keep elements of degree > dim in gI
  b = monomial_basis(quo(S, ideal(S, gI))[1])
  gI = vcat(gI, filter(x -> degree(Int, R(x)) > d, b))
  result = quo(R, gI)[1]
  return result
end

@doc raw"""
    trim!(X::AbstractVariety)

Return `X` with its Chow ring modified as follows:
- Mod out suitable powers of the generators of `base_ring(R)` so that the resulting new ring is zero-dimensional (in the example below, these are the powers  `h^3`, `c1^3`, and `c2^2`).

- Then, if `K = base(X)`, the new ring has a finite `K`-basis. Now, additionally mod out the given generators of each class in the basis which should be zero due to dimension.

!!! note
    For the construction of a new abstract variety, it is occasionally convenient to start from the underlying graded polynomial ring of the Chow ring, possibly with some relations already given, and adding further relations step by step. The `trim!` function offers one way of doing this.

!!! warning
    As this function changes an entry in the collection of data making up `X`, it should only be used by experts who are able to foresee the potential impact of such a change.

# Examples

```jldoctest
julia> RS, _ = graded_polynomial_ring(QQ, ["h", "c1", "c2"], [1, 1, 2]);

julia> S = abstract_variety(2, RS) # no relations yet
AbstractVariety of dim 2

julia> trim!(S);

julia> chow_ring(S)
Quotient
  of multivariate polynomial ring in 3 variables over QQ graded by
    h -> [1]
    c1 -> [1]
    c2 -> [2]
  by ideal with 14 generators

julia> basis(S)
3-element Vector{Vector{MPolyQuoRingElem}}:
 [1]
 [c1, h]
 [c2, c1^2, h*c1, h^2]

```
"""
function trim!(X::AbstractVariety)
  R = _trim(chow_ring(X))
  set_attribute!(R, :abstract_variety, "AbstractVariety of dim $(dim(X))")
  set_attribute!(R, :abstract_variety_dim, dim(X))
  X.ring = R
  return X
end

### constructor

@doc raw"""
    map(X::AbstractVariety, Y::AbstractVariety, f_pullback::Vector, f_pushforward = nothing; inclusion::Bool = false, symbol::String = "x")

Return an abstract variety map $f:X \rightarrow Y$ by specifying the pullbacks of the generators of
the Chow ring $\mathrm{N}^*(Y)_{\mathbb Q}.$ If needed, also specify the pushforward map $\mathrm{N}^*(X)_{\mathbb Q} \rightarrow \mathrm{N}^*(Y)_{\mathbb Q}.$

!!! note
    The pullback is relatively easy to describe since it is functorial, while the pushforward is usually more complicated.
    In some cases, the pushforward of an element $x \in \mathrm{N}^*(X)_{\mathbb Q}$ can be automatically computed via pullback.
    Specifically, the projection formula tells us that

    $f_\ast(f^\ast(y)\cdot x) = f_\ast(x)\cdot y \;\text{ for all }\; x\in \mathrm{N}^*(X)_{\mathbb Q}, y\in \mathrm{N}^*(Y)_{\mathbb Q}.$

    Since we are working with numerical equivalence, to determine $f_\ast(x)$, it suffices to know all its intersection numbers.
    These can be readily computed via pullback, provided that all classes in $\mathrm{N}^*(Y)_{\mathbb Q}$ are known (or at least those
    classes having non-zero intersection numbers with $f_\ast(x)$).
    This is the case in the following situations:
    - When $Y$ is a point or a curve;
    - when all classes in $\mathrm{N}^*(Y)_{\mathbb Q}$ are known;
    - when `:alg` is passed as the fourth argument. This can be done when we are certain that the computed pushforward
      is correct, even though not all classes in $\mathrm{N}^*(X)_{\mathbb Q}$ are known.
    In the other cases, if no pushforward map has been specified, a warning will be given when trying to do pushforward.

!!! note
    In the case of an inclusion `X` $\hookrightarrow$ `Y` where the class of `X` in $\mathrm{N}^*(Y)_{\mathbb Q}$
    is not known, use the argument `extend_inclusion = true`. Then,
    a modified version of `Y` will be created, with extra classes added so that one can
    pushforward all classes on `X`. See the subsection [Example: Cubic fourfolds](@ref)
    of the documentation for an example.

# Examples

```jldoctest
julia> P2xP2 = abstract_projective_space(2, symbol = "k")*abstract_projective_space(2, symbol = "l")
AbstractVariety of dim 4

julia> k, l = gens(P2xP2)
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 k
 l

julia> P8 = abstract_projective_space(8)
AbstractVariety of dim 8

julia> h = gens(P8)[1]
h

julia> Se = map(P2xP2, P8, [k+l]) # Segre embedding
AbstractVarietyMap from AbstractVariety of dim 4 to AbstractVariety of dim 8

julia> pullback(Se, h)
k + l

julia> pushforward(Se, k+l)
6*h^5

julia> pushforward(Se, pullback(Se, h)*(k+l)) == h*pushforward(Se, k+l)
true

```
"""
function map(
  X::AbstractVariety,
  Y::AbstractVariety,
  f_pullback::Vector,
  f_pushforward=nothing;
  inclusion::Bool=false,
  symbol::String="x",
)
  #AbstractVarietyMap(X, Y, f_pullback, f_pushforward)
  !inclusion && return AbstractVarietyMap(X, Y, f_pullback, f_pushforward)
  extend_inclusion(AbstractVarietyMap(X, Y, f_pullback), symbol=symbol)
end

@doc raw"""
    dim(f::AbstractVarietyMap)

Return the relative dimension of `f`, that is, return `dim(domain(f)) - dim(codomain(f))`.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P5 = abstract_projective_space(5, symbol = "H")
AbstractVariety of dim 5

julia> h = gens(P2)[1]
h

julia> i = map(P2, P5, [2*h])
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 5

julia> dim(i)
-3

```
"""
dim(f::AbstractVarietyMap) = f.dim

@doc raw"""
    tangent_bundle(f::AbstractVarietyMap)

If `domain(f)` and `codomain(f)` are given tangent bundles, return the relative tangent bundle of `f`.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> T = tangent_bundle(P2)
AbstractBundle of rank 2 on AbstractVariety of dim 2

julia> PT = projective_bundle(T)
AbstractVariety of dim 3

julia> pr = structure_map(PT)
AbstractVarietyMap from AbstractVariety of dim 3 to AbstractVariety of dim 2

julia> PBT = pullback(pr, T)
AbstractBundle of rank 2 on AbstractVariety of dim 3

julia> tangent_bundle(PT) - PBT == tangent_bundle(pr)
true

julia> PBT*OO(PT, 1) - OO(PT) == tangent_bundle(pr) # relative Euler sequence
true

```
"""
tangent_bundle(f::AbstractVarietyMap) = f.T

@doc raw"""
    cotangent_bundle(f::AbstractVarietyMap)

If `domain(f)` and `codomain(f)` are given tangent bundles, return the relative cotangent bundle of `f`.
"""
cotangent_bundle(f::AbstractVarietyMap) = dual(tangent_bundle(f))

@doc raw"""
    todd_class(f::AbstractVarietyMap)

Return the Todd class of the relative tangent bundle of `f`.

# Examples

The fibre variable `z` is the relative polarization of the projective bundle:

```jldoctest
julia> P2 = abstract_projective_space(2);

julia> PF = projective_bundle(tautological_bundles(P2)[2]);

julia> z = polarization(PF)  # relative O(1) class on the fibre
z

julia> todd_class(structure_map(PF))
z - 1//4*h^2 + 1//2*h + 1

```
"""
todd_class(f::AbstractVarietyMap) = todd_class(tangent_bundle(f))

@doc raw"""
    pullback(f::AbstractVarietyMap, y::MPolyDecRingOrQuoElem)

Return the pullback of `y` via `f`.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P5 = abstract_projective_space(5, symbol = "H")
AbstractVariety of dim 5

julia> h = gens(P2)[1]
h

julia> i = map(P2, P5, [2*h])
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 5

julia> H = gens(P5)[1]
H

julia> pullback(i, H)
2*h

```
"""
pullback(f::AbstractVarietyMap, y::MPolyDecRingOrQuoElem) = f.pullback(y)

@doc raw"""
    pullback_morphism(f::AbstractVarietyMap)

Return the underlying ring homomorphism (an `AffAlgHom`) representing the
pullback on Chow rings.
"""
pullback_morphism(f::AbstractVarietyMap) = f.pullback

@doc raw"""
    pushforward_morphism(f::AbstractVarietyMap)

Return the underlying map representing the pushforward on Chow rings.
"""
pushforward_morphism(f::AbstractVarietyMap) = f.pushforward

@doc raw"""
    polarization(f::AbstractVarietyMap)

Return the relative polarization class (the first Chern class of
$\mathcal O(1)$) associated to `f`.
"""
polarization(f::AbstractVarietyMap) = f.O1

@doc raw"""
    pullback(f::AbstractVarietyMap, F::AbstractBundle)

Return the pullback of `F` via `f`.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P5 = abstract_projective_space(5, symbol = "H")
AbstractVariety of dim 5

julia> h = gens(P2)[1]
h

julia> i = map(P2, P5, [2*h])
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 5

julia> E = pullback(i, OO(P5,1))
AbstractBundle of rank 1 on AbstractVariety of dim 2

julia> total_chern_class(E)
2*h + 1

```
"""
pullback(f::AbstractVarietyMap, F::AbstractBundle) = AbstractBundle(
  domain(f), pullback(f, chern_character(F))
)

@doc raw"""
    pushforward(f::AbstractVarietyMap, x::MPolyDecRingOrQuoElem)

Return the pushforward of `x` via `f`.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P5 = abstract_projective_space(5, symbol = "H")
AbstractVariety of dim 5

julia> h = gens(P2)[1]
h

julia> i = map(P2, P5, [2*h])
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 5

julia> pushforward(i, h)
2*H^4
```
"""
pushforward(f::AbstractVarietyMap, x::MPolyDecRingOrQuoElem) = f.pushforward(x)

@doc raw"""
    pushforward(f::AbstractVarietyMap, F::AbstractBundle)

Return the pushforward of `F` via `f`, that is, return the alternating sum of all direct images of `F` via `f`.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P5 = abstract_projective_space(5, symbol = "H")
AbstractVariety of dim 5

julia> h = gens(P2)[1]
h

julia> f = map(P2, P5, [2*h])
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 5

julia> F = OO(P2,1)
AbstractBundle of rank 1 on AbstractVariety of dim 2

julia> c = pushforward(f, chern_character(F) * todd_class(f))
7*H^5 - 7*H^4 + 4*H^3

julia> chern_character(pushforward(f, F)) == c # Grothendieck-Riemann-Roch
true

```
"""
pushforward(f::AbstractVarietyMap, F::AbstractBundle) = AbstractBundle(
  codomain(f), pushforward(f, chern_character(F) * todd_class(f))
) # Grothendieck-Hirzebruch-Riemann-Roch

@doc raw"""
    identity_map(X::AbstractVariety)

Return the identity map on `X`.
"""
function identity_map(X::AbstractVariety)
  R = chow_ring(X)
  AbstractVarietyMap(X, X, gens(R), MapFromFunc(R, R, identity))
end

@doc raw"""
    compose(f::AbstractVarietyMap, g::AbstractVarietyMap)

Given abstract variety maps `f` (from `X` to `Y`) and `g` (from `Y` to `Z`), return their composition.
"""
function compose(f::AbstractVarietyMap, g::AbstractVarietyMap)
  X, Y = domain(f), codomain(f)
  @assert domain(g) == Y

  Z = codomain(g)

  gof_pushforward = nothing
  if isdefined(f, :pushforward) && isdefined(g, :pushforward)
    gof_pushforward = MapFromFunc(
      chow_ring(X), chow_ring(Z), x -> pushforward_morphism(g)(pushforward_morphism(f)(x))
    )
  end

  gof = AbstractVarietyMap(
    X, Z, pullback_morphism(g) * pullback_morphism(f), gof_pushforward
  )

  return gof
end

###############################################################################
#
# AbstractVariety
#
# generic abstract variety with some classes in given degrees
@doc raw"""
    abstract_variety(n::Int, symbols::Vector{String}, degs::Vector{Int}; base::Ring=QQ)

Construct a generic abstract variety of dimension $n$ with some classes in given degrees.

Return the abstract variety and the list of classes.
"""
function abstract_variety(n::Int, symbols::Vector{String}, degs::Vector{Int}; base::Ring=QQ)
  @assert length(symbols) > 0
  R, x = graded_polynomial_ring(base, symbols, degs)
  return AbstractVariety(n, R), x
end

# generic abstract variety with some bundles in given ranks
@doc raw"""
    abstract_variety(n::Int, bundles::Vector{Pair{Int, T}}) where T

Construct a generic abstract variety of dimension $n$ with some bundles of given ranks.

Return the abstract variety and the list of bundles.
"""
function abstract_variety(n::Int, bundles::Vector{Pair{Int,T}}; base::Ring=QQ) where {T}
  symbols = reduce(vcat, [_parse_symbol(s, 1:r) for (r, s) in bundles])
  degs = reduce(vcat, [1:r for (r, _) in bundles])
  X = abstract_variety(n, symbols, degs, base=base)[1]
  i = 1
  X.bundles = AbstractBundle[]
  # we create a bundle for each rank
  # the Chern character of the bundle is set to be the sum of the generators of the corresponding degrees
  for (r, _) in bundles
    push!(X.bundles, AbstractBundle(X, r, 1 + sum(gens(X)[i:(i + r - 1)])))
    i += r
  end
  return X, tautological_bundles(X)
end

# generic abstract variety with tangent bundle
@doc raw"""
    abstract_variety(n::Int)

Construct a generic abstract variety of dimension $n$ and define its tangent bundle.

Return the abstract variety.

# Examples
```jldoctest
julia> X = abstract_variety(3)
AbstractVariety of dim 3

julia> chern_character(tangent_bundle(X))
1//6*c[1]^3 + 1//2*c[1]^2 - 1//2*c[1]*c[2] + c[1] - c[2] + 1//2*c[3] + 3
```
"""
function abstract_variety(n::Int; base::Ring=QQ)
  n == 0 && return abstract_point()
  X, (T,) = abstract_variety(n, [n=>"c"], base=base)
  X.T = T
  return X
end

# abstract variety with dimension and Chow ring
@doc raw"""
    abstract_variety(n::Int, A::Union{MPolyDecRing, MPolyQuoRing{<:MPolyDecRingElem}})

Return an abstract variety by specifying its dimension `n` and its Chow ring `A`.

!!! note
    We allow (graded) polynomial rings here since for the construction of a new abstract variety,
    the expert user may find it useful to start from the underlying graded polynomial ring of the Chow ring,
    and add its defining relations step by step.

!!! note
    In addition to the dimension and the Chow ring, further data making up an abstract variety can be set.
    See the corresponding setter functions in the section [Some Particular Constructions](@ref) of the documentation.

# Examples
```jldoctest
julia> R, (h,) = graded_polynomial_ring(QQ, [:h])
(Graded multivariate polynomial ring in 1 variable over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[h])

julia> A, _ = quo(R, ideal(R, [h^3]));

julia> P2 = abstract_variety(2, A)
AbstractVariety of dim 2
```

!!! note
    The example above shows one way of setting up a version of the abstract projective plane. A more convenient way is to use the built-in command `abstract_projective_space` which implements additional data such as the tangent bundle:

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> chow_ring(P2)
Quotient
  of multivariate polynomial ring in 1 variable over QQ graded by
    h -> [1]
  by ideal (h^3)

julia> TP2 = tangent_bundle(P2)
AbstractBundle of rank 2 on AbstractVariety of dim 2

julia> chern_character(TP2)
3//2*h^2 + 3*h + 2

```

!!! note
    Typically, the coefficient ring of a Chow ring in OSCAR will be $\mathbb Q$. In order to
    allow the use of parameters, the coefficient ring may be extended, say to a function field of type $\mathbb Q(t_1, \dots, t_r).$
    The example below shows how to do this when using built-in constructors such as `abstract_projective_space`:

```jldoctest
julia> T, (t, ) = polynomial_ring(QQ, [:t])
(Multivariate polynomial ring in 1 variable over QQ, QQMPolyRingElem[t])

julia> QT = fraction_field(T)
Fraction field
  of multivariate polynomial ring in 1 variable over QQ

julia> P2 = abstract_projective_space(2, base = QT)
AbstractVariety of dim 2

julia> chow_ring(P2)
Quotient
  of multivariate polynomial ring in 1 variable over QT graded by
    h -> [1]
  by ideal (h^3)

julia> TP2 = tangent_bundle(P2)
AbstractBundle of rank 2 on AbstractVariety of dim 2

julia> chern_character(TP2*OO(P2, t))
(t^2 + 3*t + 3//2)*h^2 + (2*t + 3)*h + 2

```
"""
function abstract_variety(n::Int, R::MPolyDecRingOrQuo)
  return AbstractVariety(n::Int, R::MPolyDecRingOrQuo)
end

###############################################
### getter functions
###############################################

@doc raw"""
    chow_ring(X::AbstractVariety)

Return the Chow ring of `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P3 = abstract_projective_space(3, symbol = "H")
AbstractVariety of dim 3

julia> chow_ring(P2*P3)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    h -> [1]
    H -> [1]
  by ideal (h^3, H^4)

```
"""
chow_ring(X::AbstractVariety) = X.ring

@doc raw"""
    (X::AbstractVariety)(f::RingElement)

Coerce the ring element `f` into the Chow ring of `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> h = gens(P2)[1]
h

julia> P2(h^2 + 1)
h^2 + 1

```
"""
(X::AbstractVariety)(f::RingElement) = chow_ring(X)(f)

@doc raw"""
    gens(X::AbstractVariety)

Return the generators of the Chow ring of `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> gens(P2)
1-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 h

```
"""
gens(X::AbstractVariety) = gens(chow_ring(X))

@doc raw"""
    base(X::AbstractVariety)

Return the coefficient ring of the polynomial ring underlying the Chow ring of `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> chow_ring(P2)
Quotient
  of multivariate polynomial ring in 1 variable over QQ graded by
    h -> [1]
  by ideal (h^3)

julia> base(P2)
Rational field

```
"""
base(X::AbstractVariety) = X.base

@doc raw"""
    point_class(X::AbstractVariety)

If `X` has been given a point class, return that class.

!!! note
    A *point class* is a top-degree element of the Chow ring of `X` which integrates to 1.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P3 = abstract_projective_space(3, symbol = "H")
AbstractVariety of dim 3

julia> p = point_class(P2*P3)
h^2*H^3

julia> degree(p)
[5]

julia> integral(p)
1

```
"""
point_class(X::AbstractVariety) = X.point

@doc raw"""
    polarization(X::AbstractVariety)

If `X` has been given a polarization, return that polarization.

!!! note
    To implement a polarization $\mathcal O_X(1)$ means to implement its first Chern class.
    For Grassmannians, this is the polarization of the Plücker embedding.
    For the product of two abstract varieties with given polarizations, it is the polarization of the Segre embedding.

# Examples

```jldoctest
julia> G = abstract_grassmannian(2, 5)
AbstractVariety of dim 6

julia> D = polarization(G)
-c[1]

julia> degree(G) == integral(D^dim(G)) == 5
true

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 3 on AbstractVariety of dim 6

julia> OO(G, D) == det(Q)
true

```

```jldoctest
julia> P1s = abstract_projective_space(1, symbol = "s")
AbstractVariety of dim 1

julia> P1t = abstract_projective_space(1, symbol = "t")
AbstractVariety of dim 1

julia> P = P1s*P1t
AbstractVariety of dim 2

julia> D = polarization(P)
s + t

julia> degree(P) == integral(D^dim(P)) == 2
true

```
"""
polarization(X::AbstractVariety) = X.O1

@doc raw"""
    trivial_line_bundle(X::AbstractVariety)

Return the trivial line bundle $\mathcal O_X$ on `X`.

Alternatively, use `OO` instead of `trivial_line_bundle`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> OX = trivial_line_bundle(P2)
AbstractBundle of rank 1 on AbstractVariety of dim 2

julia> chern_character(OX)
1

```
"""
trivial_line_bundle(X::AbstractVariety) = AbstractBundle(X, X(1))
(OO)(X::AbstractVariety) = trivial_line_bundle(X)

@doc raw"""
    tautological_bundles(X::AbstractVariety)

If `X` has been given tautological bundles, return these bundles.

!!! note
    The tautological bundles depend on the construction of `X`. Here are the
    conventions for the standard constructors:
    - `abstract_projective_space(n)`: returns `[O(-1), Q]`, the tautological
      line bundle and the tautological quotient bundle of rank `n`.
    - `abstract_grassmannian(k, n)`: returns `[S, Q]`, the universal subbundle
      of rank `k` and the universal quotient bundle of rank `n - k`.
    - `abstract_flag_variety(d₁, d₂, ..., dₛ, n)`: returns the subquotient
      bundles of ranks `d₁, d₂ - d₁, ..., n - dₛ`.
    - `projective_bundle(F)`: returns `[O(-1), Q]` on ``\mathbb P(F)``.
    - `flag_bundle(F, d₁, d₂, ...)`: returns the subquotient bundles of `F`.
    - `zero_locus_section(s)`: inherits the tautological bundles from the ambient
      variety, restricted to the zero locus.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> TB = tautological_bundles(P2)
2-element Vector{AbstractBundle}:
 AbstractBundle of rank 1 on AbstractVariety of dim 2
 AbstractBundle of rank 2 on AbstractVariety of dim 2

julia> TB[1] == OO(P2, -1)
true

julia> TB[2] == tangent_bundle(P2)*OO(P2, -1)
true

```

```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> tautological_bundles(G)
2-element Vector{AbstractBundle}:
 AbstractBundle of rank 3 on AbstractVariety of dim 6
 AbstractBundle of rank 2 on AbstractVariety of dim 6

```
"""
tautological_bundles(X::AbstractVariety) = X.bundles

@doc raw"""
    structure_map(X::AbstractVariety)

If `X` has been given a structure map, return that map.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> structure_map(P2)
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 0

julia> T = tangent_bundle(P2)
AbstractBundle of rank 2 on AbstractVariety of dim 2

julia> E = projective_bundle(T, symbol = "H")
AbstractVariety of dim 3

julia> chow_ring(E)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    H -> [1]
    h -> [1]
  by ideal (h^3, H^2 + 3*H*h + 3*h^2)

julia> structure_map(E)
AbstractVarietyMap from AbstractVariety of dim 3 to AbstractVariety of dim 2

```
"""
structure_map(X::AbstractVariety) = X.structure_map

@doc raw"""
    line_bundle(X::AbstractVariety, n::RingElement)

If `X` has been given a polarization $\mathcal O_X(1)$, return the line bundle $\mathcal O_X(n)$ on `X`.

    line_bundle(X::AbstractVariety, D::Union{MPolyDecRingElem, MPolyQuoRingElem})

Given an element `D` of the Chow ring of `X`, return the line bundle $\mathcal O_X(D)$ with first Chern class $D[1]$.
 Here, $D[1]$ is the degree-1 part of `D` (geometrically, this is the codimension 1 part of $D$).

Alternatively, use `OO` instead of `line_bundle`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> tautological_bundles(P2)[1] == OO(P2, -1)
true

julia> h = gens(P2)[1]
h

julia> OO(P2, h) == OO(P2, 1)
true

```
"""
line_bundle(X::AbstractVariety, n::RingElement) = AbstractBundle(X, 1, 1+n*polarization(X))
line_bundle(X::AbstractVariety, D::Union{MPolyDecRingElem,MPolyQuoRingElem}) = AbstractBundle(
  X, 1, 1+D[1]
)

(OO)(X::AbstractVariety, n::RingElement) = line_bundle(X, n)
OO(X::AbstractVariety, D::Union{MPolyDecRingElem,MPolyQuoRingElem}) = line_bundle(X, D)

@doc raw"""
    degree(X::AbstractVariety)

If `X` has been given a polarization, return the corresponding degree of `X`.

# Examples
```jldoctest
julia> G = abstract_grassmannian(2,5)
AbstractVariety of dim 6

julia> degree(G)
5

```
"""
degree(X::AbstractVariety) = integral(polarization(X)^dim(X))

@doc raw"""
    tangent_bundle(X::AbstractVariety)

If `X` has been given a tangent bundle, return that bundle.

# Examples
```jldoctest
julia> G = abstract_grassmannian(2,5)
AbstractVariety of dim 6

julia> TG = tangent_bundle(G)
AbstractBundle of rank 6 on AbstractVariety of dim 6

julia> chern_character(TG)
-5//6*c[1]^3 + 5//24*c[1]^2*c[2] + 3//2*c[1]^2 + 5//2*c[1]*c[2] - 5*c[1] + 7//72*c[2]^3 - 25//24*c[2]^2 - c[2] + 6

```
"""
tangent_bundle(X::AbstractVariety) = X.T

@doc raw"""
    cotangent_bundle(X::AbstractVariety)

If `X` has been given a tangent bundle, return the dual bundle.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> cotangent_bundle(P2) == dual(tangent_bundle(P2))
true

```
"""
cotangent_bundle(X::AbstractVariety) = dual(tangent_bundle(X))

@doc raw"""
    canonical_class(X::AbstractVariety)

If `X` has been given a tangent bundle, return the canonical class of `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> canonical_class(P2) == chern_class(cotangent_bundle(P2), 1)
true

julia> canonical_class(P2)
-3*h

```
"""
canonical_class(X::AbstractVariety) = -chern_class(tangent_bundle(X), 1)

@doc raw"""
    canonical_bundle(X::AbstractVariety)

If `X` has been given a tangent bundle, return the canonical bundle on `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> canonical_bundle(P2) == det(cotangent_bundle(P2))
true

```
"""
canonical_bundle(X::AbstractVariety) = det(cotangent_bundle(X))

@doc raw"""
    total_chern_class(X::AbstractVariety)
    total_chern_class(X::TnVariety)

Return the total Chern class of the tangent bundle of `X`.
"""
total_chern_class(X::AbstractVariety) = total_chern_class(tangent_bundle(X))

@doc raw"""
    chern_class(X::AbstractVariety, k::Int)
    chern_class(X::TnVariety, k::Int)

Return the `k`-th Chern class of the tangent bundle of `X`.
"""
chern_class(X::AbstractVariety, k::Int) = chern_class(tangent_bundle(X), k)

@doc raw"""
    euler_number(X::AbstractVariety)

Return the Euler number of `X`.

!!! note
    Recall that in our geometric interpretation, we think of `X` as a smooth projective complex variety. The returned number is then the Euler number of `X` considered as a compact complex manifold. Note that this number coincides with the topological Euler characteristic of the manifold.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> euler_number(P2)
3

julia> euler_number(P2) == integral(total_chern_class(tangent_bundle(P2)))
true

```
"""
euler_number(X::AbstractVariety) = integral(total_chern_class(tangent_bundle(X)))

@doc raw"""
    todd_class(X::AbstractVariety)

Compute the Todd class of the tangent bundle of `X`.

# Examples

```jldoctest
julia> P2 = abstract_projective_space(2);

julia> todd_class(P2)
h^2 + 3//2*h + 1

```
"""
todd_class(X::AbstractVariety) = todd_class(tangent_bundle(X))

@doc raw"""
    total_pontryagin_class(X::AbstractVariety)

Compute the total Pontryagin class of the tangent bundle of `X`.

# Examples
```jldoctest
julia> S = complete_intersection(abstract_projective_space(3), 4);

julia> total_pontryagin_class(S)
-12*h^2 + 1

```
"""
total_pontryagin_class(X::AbstractVariety) = total_pontryagin_class(tangent_bundle(X))

@doc raw"""
    pontryagin_class(X::AbstractVariety, k::Int)

Compute the `k`-th Pontryagin class of the tangent bundle of `X`.

# Examples
```jldoctest
julia> S = complete_intersection(abstract_projective_space(3), 4);

julia> pontryagin_class(S, 1)
-12*h^2

```
"""
pontryagin_class(X::AbstractVariety, k::Int) = pontryagin_class(tangent_bundle(X), k)

@doc raw"""
    chi(p::Int, X::AbstractVariety)

Return the generalized Todd genus $\chi_p(X) = \chi(\Omega^p_X)$, the holomorphic
Euler characteristic of the `p`-th exterior power of the cotangent bundle of `X`.

For $p = 0$ this is the arithmetic genus $\chi(\mathcal{O}_X)$.

# Examples
```jldoctest
julia> P3 = abstract_projective_space(3)
AbstractVariety of dim 3

julia> chi(0, P3)
1

julia> chi(1, P3)
-1

```
"""
chi(p::Int, X::AbstractVariety) = euler_characteristic(
  exterior_power(dual(tangent_bundle(X)), p)
) # generalized Todd genus

@doc raw"""
    todd_polynomial(n::Int)

Return the Todd polynomial $T_n(y) = \sum_{p=0}^{n} \chi_p \cdot (y-1)^p$ of a
generic variety of dimension `n`.

!!! warning "Experimental"
    This function requires a generic abstract variety with symbolic Chern classes.
    It may not work correctly in all cases due to the symbolic integral computation.
"""
function todd_polynomial(n::Int)
  X = abstract_variety(n)
  R, z = chow_ring(X)[:z]
  sum(chi(p, X) * (z-1)^p for p in 0:n)
end

@doc raw"""
    chern_number(X::AbstractVariety, λ::Int...)
    chern_number(X::AbstractVariety, λ::Vector{Int})
    chern_number(X::AbstractVariety, λ::Partition)

Compute the Chern number $\operatorname{c}_\lambda (X):=\int_X \operatorname{c}_{\lambda_1}(X)\cdots
\operatorname{c}_{\lambda_k}(X)$, where $\lambda:=(\lambda_1,\dots,\lambda_k)$ is a partition
of the dimension of `X`.

# Examples

A quartic surface in $\mathbb P^3$ (a K3 surface) and a quintic threefold in $\mathbb P^4$:

```jldoctest
julia> K3 = zero_locus_section(OO(abstract_projective_space(3), 4));

julia> chern_number(K3, 2)
24

julia> chern_number(K3, 1, 1)
0

julia> Q = zero_locus_section(OO(abstract_projective_space(4), 5));

julia> chern_number(Q, 3)
-200

julia> chern_number(Q, 2, 1)
0

```
"""
chern_number(X::AbstractVariety, lambda::Int...) = chern_number(X, collect(lambda))
chern_number(X::AbstractVariety, lambda::Partition) = chern_number(X, Vector(lambda))

function chern_number(X::AbstractVariety, lambda::Vector{Int})
  @assert sum(lambda) == dim(X)
  c = total_chern_class(X)[1:dim(X)]
  integral(prod([c[i] for i in lambda]))
end

@doc raw"""
    chern_numbers(X::AbstractVariety)

Compute all the Chern numbers of `X` as a list of pairs $\lambda\Rightarrow
\operatorname{c}_\lambda(X)$.

# Examples

A K3 surface (quartic in $\mathbb P^3$):

```jldoctest
julia> K3 = complete_intersection(abstract_projective_space(3), 4);

julia> chern_numbers(K3)
2-element Vector{Pair{Partition{Int64}, QQFieldElem}}:
    [2] => 24
 [1, 1] => 0

```

The Fano variety of lines on a cubic fourfold, a hyperkähler fourfold:

```jldoctest
julia> G = abstract_grassmannian(2, 6);

julia> F = zero_locus_section(symmetric_power(dual(tautological_bundles(G)[1]), 3));

julia> chern_numbers(F)
5-element Vector{Pair{Partition{Int64}, QQFieldElem}}:
          [4] => 324
       [3, 1] => 0
       [2, 2] => 828
    [2, 1, 1] => 0
 [1, 1, 1, 1] => 0

```
"""
function chern_numbers(X::AbstractVariety)
  c = total_chern_class(X)[1:dim(X)]
  [lambda => integral(prod([c[i] for i in lambda])) for lambda in partitions(dim(X))]
end

for g in [:a_hat_genus, :l_genus]
  @eval function $g(k::Int, X::AbstractVariety)
    R = chow_ring(X)
    k == 0 && return R(1)
    p = total_pontryagin_class(tangent_bundle(X))[1:2k]
    R isa MPolyDecRing && return R($g(k).f([p[2i].f for i in 1:k]...))
    R isa MPolyQuoRing && return R(base_ring(R)($g(k).f([p[2i].f.f for i in 1:k]...)))
  end

  @eval function $g(X::AbstractVariety)
    !iseven(dim(X)) && error("the abstract variety is not of even dimension")
    integral($g(div(dim(X), 2), X))
  end
end

@doc raw"""
    a_hat_genus(k::Int, X::AbstractVariety)

Compute the `k`-th $\hat A$ genus of the abstract variety `X`.

# Examples
```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> a_hat_genus(1, P4)
-5//24*h^2

```
"""
a_hat_genus(k::Int, X::AbstractVariety)

@doc raw"""
    l_genus(k::Int, X::AbstractVariety)

Compute the `k`-th L genus of the abstract variety `X`.

# Examples
```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> l_genus(1, P4)
5//3*h^2

```
"""
l_genus(k::Int, X::AbstractVariety)

@doc raw"""
    a_hat_genus(X::AbstractVariety)

Compute the top $\hat A$ genus of the abstract variety `X` of even dimension.

# Examples
```jldoctest
julia> a_hat_genus(abstract_projective_space(4))
3//128

```
"""
a_hat_genus(X::AbstractVariety)

@doc raw"""
    l_genus(X::AbstractVariety)

Compute the top L genus of the abstract variety `X` of even dimension.

# Examples
```jldoctest
julia> l_genus(abstract_projective_space(4))
1

```
"""
l_genus(X::AbstractVariety)

@doc raw"""
    signature(X::AbstractVariety)

Compute the signature of the abstract variety `X` by evaluating the Hirzebruch
L-genus in the Chow ring.

For smooth compact complex varieties of even complex dimension $2m$, this agrees
with the topological signature of the intersection form on middle cohomology,
i.e. the difference $b_+ - b_-$ on $H^{2m}(X, \mathbb Q)$.

# Examples
```jldoctest
julia> signature(complete_intersection(abstract_projective_space(3), 4))
-16

```
"""
signature(X::AbstractVariety) = l_genus(X) # Hirzebruch signature theorem

@doc raw"""
    hilbert_polynomial(F::AbstractBundle)

If `F` is an abstract vector bundle on an abstract variety with a given polarization,
then return the corresponding Hilbert polynomial of `F`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> hilbert_polynomial(OO(P2))
1//2*t^2 + 3//2*t + 1

julia> euler_characteristic(OO(P2))
1

julia> euler_characteristic(OO(P2, 1))
3

julia> euler_characteristic(OO(P2, 2))
6

julia> euler_characteristic(OO(P2, 3))
10

```
"""
function hilbert_polynomial(F::AbstractBundle)
  X = parent(F)
  !isdefined(X, :O1) && error("no polarization is specified for the abstract variety")
  O1 = polarization(X)

  # extend the coefficient ring to QQ[t] to allow for the parameter t in the Hilbert polynomial
  Qt, t = polynomial_ring(QQ, :t)
  R_chow = chow_ring(X)
  @assert R_chow isa MPolyQuoRing
  R = parent(change_base_ring(Qt, base_ring(R_chow).R()))
  GR = grade(R, gradings(base_ring(R_chow)))[1]
  toR = x -> GR(change_base_ring(Qt, x, parent=R))
  I = ideal(toR.(gens(R_chow.I)))
  R_ = quo(GR, I)[1]
  set_attribute!(R_, :abstract_variety_dim => dim(X))

  # translate the bundles to new coefficient ring
  ch_O_t = 1 + _logg(1 + t * R_(toR(O1.f)))
  ch_F = R_(toR(chern_character(F).f))
  td = R_(toR(todd_class(X).f))
  pt = R_(toR(point_class(X).f))

  constant_coefficient(div(simplify(ch_F * ch_O_t * td).f, simplify(pt).f))
end

@doc raw"""
    hilbert_polynomial(X::AbstractVariety)

If `X` has been given a polarization, return the corresponding Hilbert polynomial of `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> hilbert_polynomial(P2) == hilbert_polynomial(OO(P2))
true

julia> hilbert_polynomial(P2)
1//2*t^2 + 3//2*t + 1

```
"""
hilbert_polynomial(X::AbstractVariety) = hilbert_polynomial(trivial_line_bundle(X))

# find canonically defined morphism from X to Y
function _map(X::AbstractVariety, Y::AbstractVariety)
  X == Y && return identity_map(X)
  # first handle the case where X is a (fibered) product
  projs = get_attribute(X, :projections)
  if !isnothing(projs)
    for p in projs
      codomain(p) == Y && return p
    end
  else
    # follow the chain of structure maps to see if we can arrive at Y
    homs = AbstractVarietyMap[]
    while isdefined(X, :structure_map) && X != Y
      push!(homs, structure_map(X))
      X = codomain(structure_map(X))
    end
    X == Y && return reduce(compose, homs)
  end
  error("no canonical homomorphism between the given varieties")
end

# morphisms for points are convenient, but are not desired when doing coercion
@doc raw"""
    map(X::AbstractVariety, Y::AbstractVariety)

Return a canonically defined map from `X` to `Y`.
"""
function map(X::AbstractVariety, Y::AbstractVariety)
  !isnothing(get_attribute(Y, :point)) && return map(X, Y, [X(0)]) # Y is a point
  !isnothing(get_attribute(X, :point)) && return map(X, Y, repeat([X(0)], length(gens(Y)))) # X is a point
  _map(X, Y)
end

# product abstract variety
@doc raw"""
    product(X::AbstractVariety, Y::AbstractVariety)

Return the product $X\times Y$. Alternatively, use `*`.

!!! note
    If both `X` and `Y` have been given a polarization, $X\times Y$ will be endowed with the polarization corresponding to the Segre embedding.

```jldoctest
julia> P2 = abstract_projective_space(2);

julia> P3 = abstract_projective_space(3, symbol = "H");

julia> P2xP3 = P2*P3
AbstractVariety of dim 5

julia> chow_ring(P2xP3)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    h -> [1]
    H -> [1]
  by ideal (h^3, H^4)

```
"""
function product(X::AbstractVariety, Y::AbstractVariety)
  # caching products for when we deal with maps
  product_cache = get_attribute(X, :prod_cache)
  !isnothing(product_cache) && Y in keys(product_cache) && return product_cache[Y]
  if isnothing(product_cache)
    product_cache = Dict{AbstractVariety,AbstractVariety}()
    set_attribute!(X, :prod_cache => product_cache)
  end

  @assert base(X) == base(Y)
  b = base(X)
  A, B = chow_ring(X), chow_ring(Y)
  R, x, y = graded_polynomial_ring(
    b, symbols(A), symbols(B); weights=vcat(gradings(A), gradings(B))
  )
  # we must bypass the check because R is not yet the appropriate quotient
  AtoR = Oscar.hom(A, R, x, check=false)
  BtoR = Oscar.hom(B, R, y, check=false)
  IA = ideal(A isa MPolyQuoRing ? AtoR.(A.(gens(A.I))) : [R()])
  IB = ideal(B isa MPolyQuoRing ? BtoR.(B.(gens(B.I))) : [R()])
  AXY, _ = quo(R, IA + IB)
  XY = AbstractVariety(dim(X)+dim(Y), AXY)
  if isdefined(X, :point) && isdefined(Y, :point)
    XY.point = XY(AtoR(point_class(X)) * BtoR(point_class(Y)))
  end
  p = AbstractVarietyMap(XY, X, XY.(x))
  q = AbstractVarietyMap(XY, Y, XY.(y))
  if isdefined(X, :T) && isdefined(Y, :T)
    XY.T = pullback(p, tangent_bundle(X)) + pullback(q, tangent_bundle(Y))
  end
  if isdefined(X, :O1) && isdefined(Y, :O1) # Segre embedding
    XY.O1 = pullback(p, polarization(X)) + pullback(q, polarization(Y))
  end
  if get_attribute(X, :alg) == true && get_attribute(Y, :alg) == true
    set_attribute!(XY, :alg => true)
  end
  set_attribute!(XY, :projections => [p, q])
  set_attribute!(XY, :description => "Product of $X and $Y")
  product_cache[Y] = XY
  return XY
end
*(X::AbstractVariety, Y::AbstractVariety) = product(X, Y)

@doc raw"""
    graph(f::AbstractVarietyMap)

Given a morphism $f: X\to Y$, construct $i:\Gamma_f\to X\times Y$, the
inclusion of the graph into the product.
"""
function graph(f::AbstractVarietyMap)
  X, Y = domain(f), codomain(f)
  map(X, X * Y, vcat(gens(X), pullback_morphism(f).image))
end

###############################################################################
#
# Operators on AbstractBundle
#
@doc raw"""
    adams(k::RingElement, x::MPolyDecRingOrQuoElem)

Apply the $k$-th Adams operation to the Chern character or total Chern class `x`.
The Adams operation acts on degree-$i$ components by multiplication with $k^i$.

The parameter `k` can be an integer or a symbolic ring element.
"""
function adams(k::RingElement, x::MPolyDecRingOrQuoElem)
  R = parent(x)
  n = get_attribute(R, :abstract_variety_dim)::Int
  comps = x[0:n]
  sum([ZZ(k)^i*comps[i + 1] for i in 0:n])
end

adams(k::Int, x::MPolyDecRingOrQuoElem) = adams(ZZ(k), x)

@doc raw"""
    adams(k::RingElement, F::AbstractBundle)

Return the abstract bundle obtained by applying the $k$-th Adams operation to `F`.
The Adams operation $\psi^k$ acts on line bundles by $\psi^k(L) = L^{\otimes k}$
and extends to all bundles by additivity.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> T = tangent_bundle(P2)
AbstractBundle of rank 2 on AbstractVariety of dim 2

julia> adams(2, T) == abstract_bundle(P2, adams(2, chern_character(T)))
true

julia> adams(-1, T) == dual(T)
true

```
"""
function adams(k::RingElement, F::AbstractBundle)
  AbstractBundle(parent(F), adams(k, chern_character(F)))
end

adams(k::Int, F::AbstractBundle) = adams(ZZ(k), F)

@doc raw"""
    cannibalistic(k::RingElement, x::MPolyDecRingOrQuoElem)

Compute the cannibalistic class $\theta^k(x)$ of a total Chern class `x`.
The cannibalistic class is defined by $\theta^k(E) = \psi^k(E) / E$ in the
K-theory ring, expressed in terms of Chern classes.

More precisely, if `x` is the total Chern class of a bundle $E$ of rank $r$,
then `cannibalistic(k, x)` returns the total Chern class of the virtual bundle
whose Chern character is $\psi^k(\operatorname{ch}(E)) / \operatorname{ch}(E)$.
"""
function cannibalistic(k::RingElement, x::MPolyDecRingOrQuoElem)
  R = parent(x)
  n = get_attribute(R, :abstract_variety_dim)::Int
  # Convert Chern class to Chern character, apply Adams, and divide
  ch = _logg(x)
  adams_ch = adams(k, ch)
  # The cannibalistic class is exp(adams_ch) / exp(ch) = exp(adams_ch - ch)
  # But actually θ^k is defined via: if ch = Σ e^{xᵢ} then ψ^k(ch) = Σ e^{k·xᵢ}
  # and θ^k = ψ^k / id in K-theory, i.e. Π (1-e^{k·xᵢ})/(1-e^{xᵢ}) for the cannibalistic class acting on (1-E)
  # Here we implement it simply: return the Chern class of the bundle with Chern character adams(k, ch)
  _expp(adams_ch)
end

cannibalistic(k::Int, x::MPolyDecRingOrQuoElem) = cannibalistic(ZZ(k), x)

@doc raw"""
    cannibalistic(k::RingElement, F::AbstractBundle)

Return the cannibalistic class $\theta^k(F)$.
"""
function cannibalistic(k::RingElement, F::AbstractBundle)
  cannibalistic(k, total_chern_class(F))
end

cannibalistic(k::Int, F::AbstractBundle) = cannibalistic(ZZ(k), F)

@doc raw"""
    dual(F::AbstractBundle)

Return the dual bundle of `F`.

# Examples
```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> h = gens(P4)[1]
h

julia> F = abstract_bundle(P4, 2, 10*h^2 + 5*h + 1) # Horrocks-Mumford bundle
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> c1 = chern_class(F, 1)
5*h

julia> Fd = dual(F)
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> chern_class(Fd, 1)
-5*h

julia> F == Fd*OO(P4, 5) # self-duality up to twist
true

```
"""
function dual(F::AbstractBundle)
  Fdual = AbstractBundle(parent(F), adams(-1, chern_character(F)))
  if isdefined(F, :chern)
    Fdual.chern = adams(-1, total_chern_class(F))
  end
  return Fdual
end

@doc raw"""
    -(F::AbstractBundle)
    *(n::Integer, F::AbstractBundle)
    *(F::AbstractBundle, n::Integer)
    ^(F::AbstractBundle, n::Integer)
    +(F::AbstractBundle, G::AbstractBundle)
    -(F::AbstractBundle, G::AbstractBundle)
    *(F::AbstractBundle, G::AbstractBundle)

Compute the negation, scalar product, n-fold tensor product, sum, difference, and tensor product of abstract bundles, respectively.

# Examples
```jldoctest
julia> P3 = abstract_projective_space(3)
AbstractVariety of dim 3

julia> 4*OO(P3, 1) - OO(P3) == tangent_bundle(P3) # Euler sequence
true

```
"""
-(F::AbstractBundle) = AbstractBundle(parent(F), -chern_character(F))
*(n::RingElement, F::AbstractBundle) = AbstractBundle(parent(F), n * chern_character(F))
*(F::AbstractBundle, n::RingElement) = n * F
^(F::AbstractBundle, n::Integer) = AbstractBundle(parent(F), chern_character(F)^n)
for O in [:(+), :(-), :(*)]
  @eval ($O)(F::AbstractBundle, G::AbstractBundle) = (
    (F, G)=_coerce(F, G);
    AbstractBundle(parent(F), $O(chern_character(F), chern_character(G))))
end

@doc raw"""
    det(F::AbstractBundle)

Return the determinant bundle of `F`.

# Examples
```jldoctest
julia> P3 = abstract_projective_space(3)
AbstractVariety of dim 3

julia> T = tangent_bundle(P3)
AbstractBundle of rank 3 on AbstractVariety of dim 3

julia> chern_class(T, 1)
4*h

julia> det(T) == OO(P3, 4)
true

```
"""
det(F::AbstractBundle) = AbstractBundle(parent(F), 1, 1 + chern_class(F, 1))

function _coerce(F::AbstractBundle, G::AbstractBundle)
  X, Y = parent(F), parent(G)
  X == Y && return F, G
  try
    return F, pullback(_map(X, Y), G)
  catch
    try
      return pullback(_map(Y, X), F), G
    catch
      error("the sheaves are not on compatible varieties")
    end
  end
end

hom(F::AbstractBundle, G::AbstractBundle) = dual(F) * G

@doc raw"""
    exterior_power(F::AbstractBundle, k::Int)

Return the `k`-th exterior power of `F`.
"""
function exterior_power(F::AbstractBundle, k::Int)
  AbstractBundle(parent(F), _wedge(k, chern_character(F))[end])
end

@doc raw"""
    exterior_power(F::AbstractBundle)

Return the total exterior power $\lambda_{-1}(F) = \sum_{k=0}^{r} (-1)^k \bigwedge^k F$,
where $r$ is the rank of `F`.
"""
function exterior_power(F::AbstractBundle)
  AbstractBundle(
    parent(F),
    sum([(-1)^(i-1) * w for (i, w) in enumerate(_wedge(rank(F), chern_character(F)))]),
  )
end

@doc raw"""
    symmetric_power(F::AbstractBundle, k::Int)
    symmetric_power(F::AbstractBundle, k::RingElement)

Return the `k`-th symmetric power of `F`. Here, `k` can contain parameters.
"""
function symmetric_power(F::AbstractBundle, k::Int)
  AbstractBundle(parent(F), _sym(k, chern_character(F))[end])
end

function symmetric_power(F::AbstractBundle, k::RingElement)
  X = parent(F)
  PF = projective_bundle(dual(F))
  p = structure_map(PF)
  AbstractBundle(
    X, pushforward(p, sum((chern_character(line_bundle(PF, k)) * todd_class(p))[0:dim(PF)]))
  )
end

@doc raw"""
    schur_functor(F::AbstractBundle, λ::Vector{Int})
    schur_functor(F::AbstractBundle, λ::Partition)

Return the result of the Schur functor $\mathbf S^\lambda$.
"""
function schur_functor(F::AbstractBundle, lambda::Vector{Int})
  schur_functor(F, partition(lambda))
end
function schur_functor(F::AbstractBundle, lambda::Partition)
  lambda = conjugate(lambda)
  X = parent(F)
  w = _wedge(sum(lambda), chern_character(F))
  S, ei = polynomial_ring(QQ, "e#" => 1:length(w))
  e = i -> i < 0 ? S() : ei[i + 1]
  M = [e(lambda[i]-i+j) for i in 1:length(lambda), j in 1:length(lambda)]
  sch = det(matrix(S, M)) # Jacobi-Trudi
  R = chow_ring(X)
  if R isa MPolyQuoRing
    # StoX = Oscar.hom(S, R.R.R, [wi.f.f for wi in w])
    StoX = Oscar.hom(S, base_ring(R).R, [wi.f.f for wi in w])
    # return AbstractBundle(X, X(R.R(StoX(sch))))
    return AbstractBundle(X, X(base_ring(R)(StoX(sch))))
  else
    StoX = Oscar.hom(S, R.R, [wi.f for wi in w])
    return AbstractBundle(X, X(StoX(sch)))
  end
end

function giambelli(F::AbstractBundle, lambda::Vector{Int})
  R = chow_ring(parent(F))
  M = [chern_class(F, lambda[i]-i+j).f for i in 1:length(lambda), j in 1:length(lambda)]
  R(det(matrix(R, M)))
end

###############################################################################
#
# Various computations
#
@doc raw"""
    basis(X::AbstractVariety)
    basis(X::AbstractVariety, k::Int)

If `K = base(X)`, return a `K`-basis of the Chow ring of `X` (return the elements of degree `k` in that basis).

!!! note
    The basis elements are ordered by increasing degree (geometrically, by increasing codimension).

# Examples
```jldoctest
julia> G = abstract_grassmannian(2,4)
AbstractVariety of dim 4

julia> chow_ring(G)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    c[1] -> [1]
    c[2] -> [2]
  by ideal (-c[1]^3 + 2*c[1]*c[2], c[1]^4 - 3*c[1]^2*c[2] + c[2]^2)

julia> basis(G)
5-element Vector{Vector{MPolyQuoRingElem}}:
 [1]
 [c[1]]
 [c[2], c[1]^2]
 [c[1]*c[2]]
 [c[2]^2]

julia> basis(G, 2)
2-element Vector{MPolyQuoRingElem}:
 c[2]
 c[1]^2

```
"""
@attr Vector{Vector{MPolyQuoRingElem}} function basis(X::AbstractVariety)
  # it is important for this to be cached!
  R = chow_ring(X)
  try_trim = "Try use `trim!`."
  !(R isa MPolyQuoRing) && error("the ring has no ideal. "*try_trim)
  krull_dim(R.I) > 0 && error("the ideal is not 0-dimensional. "*try_trim)
  b = Oscar._kbase(R)
  ans = [MPolyQuoRingElem[] for i in 0:dim(X)]
  for bi in b
    push!(ans[_total_degree(bi) + 1], R(bi))
  end
  return ans
end

@doc raw"""
    basis(X::AbstractVariety, k::Int)

If `K = base(X)`, return the elements of degree `k` in a `K`-basis of the Chow ring of `X`.
"""
basis(X::AbstractVariety, k::Int) = basis(X)[k + 1]

@doc raw"""
    betti_numbers(X::AbstractVariety)

Return the Betti numbers of the Chow ring of `X`.

!!! note
    The Betti number of `X` in a given degree is the number of elements of `basis(X)` in that degree.

!!! note
    These are the Betti numbers of the *numerical* Chow ring, not the topological Betti numbers.
    For example, for a quartic surface in $\mathbb P^3$, the numerical Chow ring sees only the
    algebraic classes, so the Betti numbers are `[1, 1, 1]`, while the topological Betti numbers
    would be `[1, 0, 22, 0, 1]` (since `h^{1,1} = 20` but only 1 algebraic class is detected
    by the Chow ring).

# Examples
```jldoctest
julia> P2xP2 = abstract_projective_space(2, symbol = "k")*abstract_projective_space(2, symbol = "l")
AbstractVariety of dim 4

julia> betti_numbers(P2xP2)
5-element Vector{Int64}:
 1
 2
 3
 2
 1

julia> basis(P2xP2)
5-element Vector{Vector{MPolyQuoRingElem}}:
 [1]
 [l, k]
 [l^2, k*l, k^2]
 [k*l^2, k^2*l]
 [k^2*l^2]

```

```jldoctest
julia> betti_numbers(complete_intersection(abstract_projective_space(3), 4))
3-element Vector{Int64}:
 1
 1
 1

```
"""
betti_numbers(X::AbstractVariety) = length.(basis(X))

@doc raw"""
    integral(c:::Union{MPolyDecRingElem, MPolyQuoRingElem})

Given an element `c` of the Chow ring of an abstract variety, say, `X`, return the integral of `c`.

!!! note
    If `X` has been given a (unique) point class,
    then the integral will be an element of the coefficient ring of the Chow ring.
    That is, typically, in the applications we discuss here, it will be a rational number (the degree of the 0-dimensional part
    of `c`) or an element of a function field of type $\mathbb Q(t_1, \dots, t_r)$.  If one of the conditions is not fulfilled, the 0-dimensional
    part of `c` is returned.

# Examples

```jldoctest
julia> G = abstract_grassmannian(2, 5)
AbstractVariety of dim 6

julia> p = point_class(G)
c[2]^3

julia> integral(p)
1

```

```jldoctest
julia> T, (t, ) = polynomial_ring(QQ, [:t])
(Multivariate polynomial ring in 1 variable over QQ, QQMPolyRingElem[t])

julia> QT = fraction_field(T)
Fraction field
  of multivariate polynomial ring in 1 variable over QQ

julia> P3 = abstract_projective_space(3, base = QT)
AbstractVariety of dim 3

julia> h = gens(P3)[1]
h

julia> integral(t^2*h^3+t*h)
t^2

```

# Lines on a General Cubic Hypersurface in $\mathbb P^3$

```jldoctest
julia> G = abstract_grassmannian(2, 4)
AbstractVariety of dim 4

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> E = symmetric_power(Q, 3)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> integral(top_chern_class(E))
27

```

# Lines on a General Complete Intersection Calabi-Yau Threefold of Type (2,2,2,2)

```jldoctest
julia> G = abstract_grassmannian(2, 4+4)
AbstractVariety of dim 12

julia> S = tautological_bundles(G)[1]
AbstractBundle of rank 2 on AbstractVariety of dim 12

julia> E = symmetric_power(S, 2)
AbstractBundle of rank 3 on AbstractVariety of dim 12

julia> integral(top_chern_class(E)^4)
512

```
"""
function integral(x::MPolyDecRingOrQuoElem)
  X = get_attribute(parent(x), :abstract_variety)
  if isdefined(X, :point) && length(basis(X, dim(X))) == 1
    return constant_coefficient(div(simplify(x).f, simplify(point_class(X)).f))
  else
    return x[dim(X)]
  end
end

@doc raw"""
    intersection_matrix(X::AbstractVariety)

If `b = vcat(basis(X)...)`, return `matrix([integral(bi*bj) for bi in b, bj in b])`.

    intersection_matrix(a::Vector, b::Vector)

Return `matrix([integral(ai*bj) for ai in a, bj in b])`.

    intersection_matrix(a::Vector)

As above, with `b = a`.

# Examples
```jldoctest
julia> G = abstract_grassmannian(2,4)
AbstractVariety of dim 4

julia> basis(G)
5-element Vector{Vector{MPolyQuoRingElem}}:
 [1]
 [c[1]]
 [c[2], c[1]^2]
 [c[1]*c[2]]
 [c[2]^2]

julia> b = vcat(basis(G)...)
6-element Vector{MPolyQuoRingElem}:
 1
 c[1]
 c[2]
 c[1]^2
 c[1]*c[2]
 c[2]^2

julia> intersection_matrix(G)
[0   0   0   0   0   1]
[0   0   0   0   1   0]
[0   0   1   1   0   0]
[0   0   1   2   0   0]
[0   1   0   0   0   0]
[1   0   0   0   0   0]

julia> integral(b[4]*b[4])
2

```
"""
intersection_matrix(X::AbstractVariety) = intersection_matrix(reduce(vcat, basis(X)))

function intersection_matrix(a::Vector{}, b=nothing)
  if isnothing(b)
    b = a
  end
  matrix([integral(ai*bj) for ai in a, bj in b])
end

@doc raw"""
    dual_basis(X::AbstractVariety)

If `X` has been given a point class, return a `K`-basis of the Chow ring of `X` which is dual to `basis(X)` with respect to the `K`-bilinear form defined by `intersection_matrix(X)`. Here, `K = base(X)`.

    dual_basis(X::AbstractVariety, k::Int)

If `X` has been given a point class, return the elements of degree `k` in the dual basis.

!!! note
    The basis elements are ordered by decreasing degree (geometrically, by decreasing codimension).

# Examples
```jldoctest
julia> G = abstract_grassmannian(2,4)
AbstractVariety of dim 4

julia> b = basis(G)
5-element Vector{Vector{MPolyQuoRingElem}}:
 [1]
 [c[1]]
 [c[2], c[1]^2]
 [c[1]*c[2]]
 [c[2]^2]

julia> intersection_matrix(G)
[0   0   0   0   0   1]
[0   0   0   0   1   0]
[0   0   1   1   0   0]
[0   0   1   2   0   0]
[0   1   0   0   0   0]
[1   0   0   0   0   0]

julia> bd = dual_basis(G)
5-element Vector{Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 [c[2]^2]
 [c[1]*c[2]]
 [-c[1]^2 + 2*c[2], c[1]^2 - c[2]]
 [c[1]]
 [1]

julia> integral(b[3][2]*b[3][2])
2

julia> integral(b[3][2]*bd[3][2])
1

```
"""
dual_basis(X::AbstractVariety) = [dual_basis(X, k) for k in 0:dim(X)]

function dual_basis(X::AbstractVariety, k::Int)
  isdefined(X, :point) || error("point class not defined")

  R = chow_ring(X)
  T = Dict{Int,Vector{elem_type(R)}}
  d = get_attribute!(X, :dual_basis) do
    T()
  end::T
  if !(k in keys(d))
    B = basis(X)
    b_k = B[k + 1]
    b_comp = B[dim(X) - k + 1]
    M = Matrix(inv(intersection_matrix(b_comp, b_k)))
    d[k] = M * b_comp
    d[dim(X) - k] = transpose(M) * b_k
  end
  return d[k]
end

# the parameter for truncation is usually the dimension, but can also be set
# manually, which is used when computing particular Chern classes (without
# computing the total Chern class)
function _expp(x::MPolyDecRingOrQuoElem; truncate::Int=-1)
  R = parent(x)
  n = truncate < 0 ? get_attribute(R, :abstract_variety_dim)::Int : truncate
  comps = x[0:n]
  p = [(-1)^i * factorial(ZZ(i)) * comps[i + 1] for i in 0:n]
  e = repeat([R(0)], n+1)
  e[1] = R(1)
  for i in 1:n
    e[i + 1] = QQ(-1, i) * sum(p[j + 1] * e[i - j + 1] for j in 1:i)
  end
  simplify(sum(e))
end

function _logg(x::MPolyDecRingOrQuoElem)
  R = parent(x)
  n = get_attribute(R, :abstract_variety_dim)::Int
  n == 0 && return R()
  e = x[1:n]
  p = pushfirst!(repeat([R()], n-1), -e[1])
  for i in 1:(n - 1)
    p[i + 1] = -ZZ(i+1)*e[i + 1] - sum(e[j] * p[i - j + 1] for j in 1:i)
  end
  simplify(sum((-1)^i//factorial(ZZ(i))*p[i] for i in 1:n))
end

# returns all the wedge from 0 to k
function _wedge(k::Int, x::MPolyDecRingOrQuoElem)
  R = parent(x)
  k == 0 && return [R(1)]
  n = get_attribute(R, :abstract_variety_dim)::Int
  wedge = repeat([R(0)], k+1)
  wedge[1] = R(1)
  wedge[2] = x
  for j in 2:k
    wedge[j + 1] =
      1//ZZ(j) *
      sum(sum((-1)^(j-i+1) * wedge[i + 1] * adams(j-i, x) for i in 0:(j - 1))[0:n])
  end
  wedge
end

# returns all the sym from 0 to k
function _sym(k::Int, x::MPolyDecRingOrQuoElem)
  R = parent(x)
  k == 0 && return [R(1)]
  n = get_attribute(R, :abstract_variety_dim)::Int
  r = min(k, Int(ZZ(QQ(constant_coefficient(x.f)))))
  wedge = _wedge(r, x)
  sym = repeat([R(0)], k+1)
  sym[1] = R(1)
  sym[2] = x
  for j in 2:k
    sym[j + 1] = sum(
      sum((-1)^(i+1) * wedge[i + 1] * sym[j - i + 1] for i in 1:min(j, r))[0:n]
    )
  end
  sym
end

function _genus(x::MPolyDecRingOrQuoElem, taylor::Vector{})
  R = parent(x)
  iszero(x) && return R(1)
  n = get_attribute(R, :abstract_variety_dim)
  R, (t,) = graded_polynomial_ring(QQ, [:t])
  set_attribute!(R, :abstract_variety_dim, n)
  lg = _logg(R(sum(taylor[i + 1] * t^i for i in 0:n)))
  comps = lg[1:n]
  lg = [
    if iszero(comps[i].f)
      zero(coefficient_ring(comps[i].f))
    else
      leading_coefficient(comps[i].f)
    end for i in 1:n
  ]
  comps = x[1:n]
  _expp(sum(factorial(ZZ(i)) * lg[i] * comps[i] for i in 1:n))
end

function _todd_class(x::MPolyDecRingOrQuoElem)
  n = get_attribute(parent(x), :abstract_variety_dim)::Int
  # the Taylor series of t/(1-exp(-t))
  taylor = [(-1)^i//factorial(ZZ(i))*bernoulli(i) for i in 0:n]
  _genus(x, taylor)
end

function _l_genus(x::MPolyDecRingOrQuoElem)
  n = get_attribute(parent(x), :abstract_variety_dim)::Int
  # the Taylor series of sqrt(t)/tanh(sqrt(t))
  taylor = [ZZ(2)^2i//factorial(ZZ(2i))*bernoulli(2i) for i in 0:n]
  _genus(x, taylor)
end

function _a_hat_genus(x::MPolyDecRingOrQuoElem)
  n = get_attribute(parent(x), :abstract_variety_dim)::Int
  # the Taylor series of (sqrt(t)/2)/sinh(sqrt(t)/2)
  R, t = power_series_ring(QQ, 2n+1, "t")
  s = divexact(t, exp(QQ(1//2)*t)-exp(-QQ(1//2)*t))
  taylor = [coeff(s, 2i) for i in 0:n]
  _genus(x, taylor)
end

for (g, s) in [:a_hat_genus=>"p", :l_genus=>"p", :todd_class=>"c"]
  _g = Symbol("_", g)
  @eval function $g(n::Int)
    n == 0 && return QQ(1)
    R, p = graded_polynomial_ring(QQ, $s => 1:n; weights=1:n)
    set_attribute!(R, :abstract_variety_dim, n)
    $_g(_logg(R(1+sum(p))))[n]
  end
end

@doc raw"""
    todd_class(n::Int)

Compute the (generic) $n$-th Todd genus as a polynomial in the Chern classes
$\operatorname{c}_1, \ldots, \operatorname{c}_n$.

# Examples

```jldoctest
julia> todd_class(1)
1//2*c[1]

julia> todd_class(2)
1//12*c[1]^2 + 1//12*c[2]

julia> todd_class(3)
1//24*c[1]*c[2]

julia> todd_class(4)
-1//720*c[1]^4 + 1//180*c[1]^2*c[2] + 1//720*c[1]*c[3] + 1//240*c[2]^2 - 1//720*c[4]

```
"""
todd_class(n::Int)

@doc raw"""
    l_genus(n::Int)

Compute the (generic) $n$-th L genus as a polynomial in the Pontryagin classes
$p_1, \ldots, p_n$. The first L genus $L_1 = p_1/3$ recovers the Hirzebruch
signature formula in dimension 4.

# Examples

```jldoctest
julia> l_genus(1)
1//3*p[1]

julia> l_genus(2)
-1//45*p[1]^2 + 7//45*p[2]

julia> l_genus(3)
2//945*p[1]^3 - 13//945*p[1]*p[2] + 62//945*p[3]

```
"""
l_genus(n::Int)

@doc raw"""
    a_hat_genus(n::Int)

Compute the (generic) $n$-th $\hat A$ genus as a polynomial in the Pontryagin
classes $p_1, \ldots, p_n$. The first $\hat A$ genus $\hat A_1 = -p_1/24$
appears in the Atiyah-Singer index theorem (Dirac operator).

# Examples

```jldoctest
julia> a_hat_genus(1)
-1//24*p[1]

julia> a_hat_genus(2)
7//5760*p[1]^2 - 1//1440*p[2]

julia> a_hat_genus(3)
-31//967680*p[1]^3 + 11//241920*p[1]*p[2] - 1//60480*p[3]

```
"""
a_hat_genus(n::Int)

###############################################
### setter functions
###############################################

@doc raw"""
    set_point_class(X::AbstractVariety, p::MPolyDecRingOrQuoElem)

Set the point class of `X` to `p`, an element of the Chow ring of `X` representing
the class of a point.

Once set, `integral` can be used to compute intersection numbers.

!!! note
    The point class can only be set once.
"""
function set_point_class(X::AbstractVariety, p::MPolyDecRingOrQuoElem)
  if isdefined(X, :point)
    error("Point class already set")
  end
  @assert parent(p) == chow_ring(X)
  X.point = p
end

@doc raw"""
    set_tangent_bundle(X::AbstractVariety, t::AbstractBundle)

Set the tangent bundle of `X` to `t`.

Once set, characteristic classes (Chern classes, Todd class, etc.) of `X`
can be computed.

!!! note
    The tangent bundle can only be set once.
"""
function set_tangent_bundle(X::AbstractVariety, t::AbstractBundle)
  if isdefined(X, :T)
    error("Tangent bundle already set")
  end
  @assert parent(t) == X
  X.T = t
end

@doc raw"""
    set_polarization(X::AbstractVariety, o1::MPolyDecRingOrQuoElem)

Set the polarization of `X` to `o1`, an element of degree 1 in the Chow ring.
This represents the first Chern class of an ample line bundle.

Once set, Hilbert polynomials can be computed.

!!! note
    The polarization can only be set once.
"""
function set_polarization(X::AbstractVariety, o1::MPolyDecRingOrQuoElem)
  if isdefined(X, :O1)
    error("Polarization already set")
  end
  @assert parent(o1) == chow_ring(X)
  X.O1 = o1
end

@doc raw"""
    set_tautological_bundles(X::AbstractVariety, vb::Vector{<:AbstractBundle})

Set the tautological bundles of `X` to `vb`.

Tautological bundles are the distinguished bundles that arise naturally from
the construction of `X` (e.g., the universal sub- and quotient bundles on a
Grassmannian).

!!! note
    The tautological bundles can only be set once.
"""
function set_tautological_bundles(
  X::AbstractVariety, vb::Vector{AbstractBundle{T}}
) where {T}
  if isdefined(X, :bundles)
    error("Tautological bundles already set")
  end
  for b in vb
    @assert parent(b) == X
  end
  X.bundles = vb
end

@doc raw"""
    set_structure_map(X::AbstractVariety, f::AbstractVarietyMap)

Set the structure map of `X` to `f`. The structure map is a morphism from `X`
to some base variety, used for push-forward (`integral`) and pull-back computations.

!!! note
    The structure map can only be set once.
"""
function set_structure_map(X::AbstractVariety, f::AbstractVarietyMap)
  if isdefined(X, :structure_map)
    error("Structure map already set")
  end
  @assert domain(f) == X
  X.structure_map = f
end

###############################################

@doc raw"""
    zero_locus_section(F::AbstractBundle; class::Bool = false)

Return the zero locus of a general section of `F`.

Use the argument `class = true` to only compute the class of the zero locus (same
as `top_chern_class(F)`).

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> C = zero_locus_section(OO(P2, 3)) # a plane cubic curve
AbstractVariety of dim 1

julia> dim(C)
1

julia> degree(C)
3

```

```jldoctest
julia> P3 = abstract_projective_space(3)
AbstractVariety of dim 3

julia> C = zero_locus_section(OO(P3, 2) + OO(P3, 2)) # a complete intersection
AbstractVariety of dim 1

julia> dim(C)
1

julia> degree(C)
4

```

```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> h = gens(P4)[1]
h

julia> F = abstract_bundle(P4, 2, 10*h^2 + 5*h + 1) # Horrocks-Mumford bundle
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> A = zero_locus_section(F) # abelian surface
AbstractVariety of dim 2

julia> dim(A)
2

julia> degree(A)
10

```
"""
function zero_locus_section(F::AbstractBundle; class::Bool=false)
  X = parent(F)
  R = chow_ring(X)
  cZ = top_chern_class(F)

  # return only the class of Z in the Chow ring of X
  class && return cZ

  if R isa MPolyQuoRing
    I = quotient(modulus(R), ideal(base_ring(R), [cZ.f]))
    AZ = quo(base_ring(R), I)[1]
  else
    AZ = R
  end
  Z = AbstractVariety(dim(X) - rank(F), AZ)

  if isdefined(X, :point)
    ps = basis(Z, dim(Z)) # the 0-cycles
    @assert length(ps) == 1 # make sure that the 0-cycle is unique
    p = ps[1]
    degp = integral(R(p.f) * cZ) # compute the degree of i_x p
    Z.point = Z(inv(degp) * p.f)
  end

  if isdefined(X, :T)
    Z.T = AbstractBundle(Z, Z((chern_character(tangent_bundle(X)) - chern_character(F)).f))
  end

  if isdefined(X, :O1)
    Z.O1 = Z(polarization(X).f)
  end

  i_x = x -> x.f * cZ
  i_x = MapFromFunc(chow_ring(Z), chow_ring(X), i_x)

  @assert R isa MPolyQuoRing
  i = AbstractVarietyMap(Z, X, Z.(gens(base_ring(R))), i_x)
  i.T = pullback(i, -F)
  Z.structure_map = i
  set_attribute!(Z, :description, "Zero locus of a section of $F")

  return Z
end

@doc raw"""
    complete_intersection(X::AbstractVariety, degs::Int...)
    complete_intersection(X::AbstractVariety, degs::Vector{Int})

If `X` has been given a polarization, return the complete intersection in
`X` of general hypersurfaces with the degrees given by `degs`.

# Examples
```jldoctest
julia> P3 = abstract_projective_space(3)
AbstractVariety of dim 3

julia> CI = complete_intersection(P3, 2, 2)
AbstractVariety of dim 1

julia> dim(CI)
1

julia> degree(CI)
4

julia> chow_ring(CI)
Quotient
  of multivariate polynomial ring in 1 variable over QQ graded by
    h -> [1]
  by ideal (h^2)

```
"""
complete_intersection(X::AbstractVariety, degs::Int...) = complete_intersection(
  X, collect(degs)
)
function complete_intersection(X::AbstractVariety, degs::Vector{Int})
  Y = zero_locus_section(sum(line_bundle(X, d) for d in degs))
  set_attribute!(
    Y, :description => "Complete intersection of degree $(tuple(degs...)) in $X"
  )
  Y
end

@doc raw"""
    degeneracy_locus(F::AbstractBundle, G::AbstractBundle, k::Int; class::Bool=false)

Return the `k`-th degeneracy locus of a general map from `F` to `G`.

The `k`-th degeneracy locus $\mathrm{D}_k(\varphi)$ of a map $\varphi\colon F \to G$
is the locus where $\operatorname{rank}(\varphi) \leq k$.

Use the argument `class = true` to only compute the class of the degeneracy locus
(Porteous' formula). This corresponds to `degeneracyLocus2` in Schubert2 (Macaulay2).
Without `class = true`, the actual variety is constructed, corresponding to
`degeneracyLocus` in Schubert2.

# Examples
```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> F = 3*OO(P4, -1)
AbstractBundle of rank 3 on AbstractVariety of dim 4

julia> G = cotangent_bundle(P4)*OO(P4,1)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> CZ = degeneracy_locus(F, G, 2, class = true) # only class of degeneracy locus
4*h^2

julia> CZ == chern_class(G-F, 2) # Porteous' formula
true

julia> Z = degeneracy_locus(F, G, 2) # Veronese surface in P4
AbstractVariety of dim 2

julia> degree(Z)
4

julia> chow_ring(Z)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    z -> [1]
    h -> [1]
  by ideal (2*z - 3*h, h^3)

```

```jldoctest
julia> P = abstract_projective_space(4, symbol = "H")
AbstractVariety of dim 4

julia> F = exterior_power(cotangent_bundle(P), 3)*OO(P, 3)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> G = OO(P, 1) + 4*OO(P)
AbstractBundle of rank 5 on AbstractVariety of dim 4

julia> Z = degeneracy_locus(F, G, 3) # rational surface in P4
AbstractVariety of dim 2

julia> K = canonical_class(Z)
z - H

julia> integral(K^2)
-7

julia> H = polarization(Z)
H

julia> integral(H^2) # degree of surface
8

julia> A = K+H
z

julia> integral(A^2) # degree of first adjoint surface which is a del Pezzo surface in P5
5

```
"""
function degeneracy_locus(F::AbstractBundle, G::AbstractBundle, k::Int; class::Bool=false)
  F, G = _coerce(F, G)
  m, n = rank(F), rank(G)
  @assert 0 <= k < min(m, n)

  if class
    # return only the class of D in the Chow ring of X
    if (m-k)*(n-k) <= dim(parent(F)) # Porteous' formula
      return chern_character(schur_functor(G-F, repeat([m-k], n-k)))[(m - k) * (n - k)]
    else # expected dimension is negative
      return chow_ring(parent(F))(0)
    end
  end

  if k == 0
    # k=0: the degeneracy locus where the map is zero.
    # Use zero_locus_section of dual(F) * G directly on the base variety.
    D = zero_locus_section(dual(F) * G)
    if isdefined(parent(F), :O1)
      D.O1 = pullback(D.structure_map, polarization(parent(F)))
    end
    set_attribute!(D, :description, "Degeneracy locus of rank $k from $F to $G")
    return D
  end

  Gr = (m-k == 1) ? projective_bundle(F) : flag_bundle(F, m-k)
  S = tautological_bundles(Gr)[1]
  D = zero_locus_section(dual(S) * G)
  D.structure_map = map(D, parent(F)) # skip the flag abstract variety

  if isdefined(parent(F), :O1)
    D.O1 = pullback(D.structure_map, polarization(parent(F)))
  end
  set_attribute!(D, :description, "Degeneracy locus of rank $k from $F to $G")
  return D
end

@doc raw"""
    kernel_bundle(F::AbstractBundle, G::AbstractBundle, k::Int)

Given a general map $\varphi\colon F \to G$, construct the `k`-th degeneracy
locus $\mathrm{D}_k$ where $\operatorname{rank}(\varphi) \leq k$, and return a pair
`(D, K)` where `D` is the degeneracy locus and `K` is the kernel bundle
$\ker(\varphi|_{\mathrm{D}_k})$ on `D`.

This corresponds to `kernelBundle` in Schubert2 (Macaulay2).

# Examples

The universal kernel on the degeneracy locus of a map between bundles on $\mathbb{P}^4$:

```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> F = 3*OO(P4, -1)
AbstractBundle of rank 3 on AbstractVariety of dim 4

julia> G = cotangent_bundle(P4)*OO(P4,1)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> D, K = kernel_bundle(F, G, 2);

julia> dim(D)
2

julia> rank(K)
1

```
"""
function kernel_bundle(F::AbstractBundle, G::AbstractBundle, k::Int)
  F, G = _coerce(F, G)
  m, n = rank(F), rank(G)
  @assert 0 < k < min(m, n) "k must satisfy 0 < k < min(rank(F), rank(G))"

  # Build the Grassmannian Gr(m-k, F) and the zero locus inside it
  Gr = (m-k == 1) ? projective_bundle(F) : flag_bundle(F, m-k)
  S = tautological_bundles(Gr)[1]   # tautological subbundle of rank m-k on Gr
  D = zero_locus_section(dual(S) * G)

  # Pull back the tautological subbundle to D (before resetting structure map)
  D_to_Gr = D.structure_map
  K = pullback(D_to_Gr, S)

  # Reset the structure map to go directly to X
  D.structure_map = map(D, parent(F))
  if isdefined(parent(F), :O1)
    D.O1 = pullback(D.structure_map, polarization(parent(F)))
  end
  set_attribute!(D, :description, "Degeneracy locus of rank $k from $F to $G")
  return D, K
end

###############################################################################
@doc raw"""
    abstract_point(; base::Ring = QQ)

Return an abstract variety consisting of a point.

# Examples
```jldoctest
julia> p = abstract_point()
AbstractVariety of dim 0

julia> chow_ring(p)
Quotient
  of multivariate polynomial ring in 1 variable over QQ graded by
    p -> [1]
  by ideal (p)
```
"""
function abstract_point(; base::Ring=QQ)
  R, (p,) = graded_polynomial_ring(base, [:p])
  I = ideal([p])
  pt = AbstractVariety(0, quo(R, I)[1])
  pt.point = pt(1)
  pt.T = AbstractBundle(pt, pt(0))
  pt.O1 = pt(0)
  set_attribute!(pt, :description, "Point")
  set_attribute!(pt, :point, true)
  return pt
end

@doc raw"""
    abstract_projective_space(n::Int; base::Ring = QQ, symbol::String = "h")

Return the abstract projective space of lines in an `n+1`-dimensional vector space.

!!! note
    The string `symbol` specifies how to print the generators of the Chow ring.

# Examples
```jldoctest
julia> P3 = abstract_projective_space(3)
AbstractVariety of dim 3

julia> chow_ring(P3)
Quotient
  of multivariate polynomial ring in 1 variable over QQ graded by
    h -> [1]
  by ideal (h^4)
```

```jldoctest
julia> T, (s,t) = polynomial_ring(QQ, [:s, :t])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[s, t])

julia> QT = fraction_field(T)
Fraction field
  of multivariate polynomial ring in 2 variables over QQ

julia> P3 = abstract_projective_space(3, base = QT)
AbstractVariety of dim 3

julia> chow_ring(P3)
Quotient
  of multivariate polynomial ring in 1 variable over QT graded by
    h -> [1]
  by ideal (h^4)

julia> TB = tangent_bundle(P3)
AbstractBundle of rank 3 on AbstractVariety of dim 3

julia> CTB = cotangent_bundle(P3)
AbstractBundle of rank 3 on AbstractVariety of dim 3

julia> chern_character((s*TB)*(t*CTB))
-4*s*t*h^2 + 9*s*t

```
"""
function abstract_projective_space(n::Int; base::Ring=QQ, symbol::String="h")
  R, (h,) = graded_polynomial_ring(base, [symbol])
  I = ideal([h^(n+1)])
  AP = quo(R, I)[1]
  P = AbstractVariety(n, AP)
  h = P(h)
  P.point = h^n
  P.O1 = h
  chTP = P(n)
  for i in 1:n
    chTP += ZZ(n+1)//factorial(ZZ(i))*h^i
  end
  P.T = AbstractBundle(P, chTP)
  P.T.chern = (1+h)^(n+1)
  S = AbstractBundle(P, 1, 1-h)
  Q = trivial_line_bundle(P)*(n+1) - S
  P.bundles = [S, Q]
  P.structure_map = map(P, abstract_point(base=base), [P(0)])
  set_attribute!(P, :description => "Projective space of dim $n")
  set_attribute!(P, :grassmannian => :absolute)
  set_attribute!(P, :alg => true)
  return P
end

@doc raw"""
    projective_bundle(F::AbstractBundle; symbol::String = "z")

Return the projective bundle of 1-dimensional subspaces in the fibers of `F`.

!!! note
    The string `symbol` specifies how to print the first Chern class of the line bundle $\mathcal O_{\mathbb P(\mathcal F)}(1)$.

!!! note
    For generators and relations of Chow rings of projective bundles see [EH16](@cite).
    For the pushforward of cycle classes from a projective bundle to its base variety
    see [EH16](@cite) and [DP17](@cite).

# Examples
```jldoctest
julia> P4 = abstract_projective_space(4)
AbstractVariety of dim 4

julia> T = tangent_bundle(P4)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> PT = projective_bundle(T)
AbstractVariety of dim 7

julia> R = chow_ring(PT)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    z -> [1]
    h -> [1]
  by ideal (h^5, z^4 + 5*z^3*h + 10*z^2*h^2 + 10*z*h^3 + 5*h^4)

julia> z, h = gens(R)
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 z
 h

julia> polarization(PT)
z

julia> p = structure_map(PT)
AbstractVarietyMap from AbstractVariety of dim 7 to AbstractVariety of dim 4

julia> pushforward(p, h)
0

julia> pushforward(p, z*h)
0

julia> pushforward(p, z^2*h)
0

julia> pushforward(p, z^3*h)
h

julia> s = total_segre_class(T)
70*h^4 - 35*h^3 + 15*h^2 - 5*h + 1

julia> sum(pushforward(p, (z^(i+rank(T)-1))) for i in 0:dim(P4)) == s
true

```

*Number of conics on the general quintic hypersurface in $\mathbb P^4$:*

```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> USBd = dual(tautological_bundles(G)[1])
AbstractBundle of rank 3 on AbstractVariety of dim 6

julia> F = symmetric_power(USBd, 2)
AbstractBundle of rank 6 on AbstractVariety of dim 6

julia> PF = projective_bundle(F)
AbstractVariety of dim 11

julia> A = symmetric_power(USBd, 5) - symmetric_power(USBd, 3)*OO(PF, -1)
AbstractBundle of rank 11 on AbstractVariety of dim 11

julia> integral(top_chern_class(A))
609250

```
"""
function projective_bundle(F::AbstractBundle; symbol::String="z")
  X, r = parent(F), rank(F)
  !(r isa Int) && error("expect rank to be an integer")
  R = chow_ring(X)

  # construct the ring

  w = vcat([1], gradings(R))
  R1, (z,), imgs_in_R1 = graded_polynomial_ring(base(X), [symbol], symbols(R); weights=w)
  if R isa MPolyQuoRing
    PR = base_ring(R)
  else
    PR = R
  end
  pback = Oscar.hom(PR, R1, imgs_in_R1)
  pfwd = Oscar.hom(R1, R, pushfirst!(gens(R), R()))

  # construct the relations

  rels = [sum(pback(chern_class(F, i).f) * z^(r-i) for i in 0:r)]
  if R isa MPolyQuoRing
    rels = vcat(pback.(gens(R.I)), rels)
  end
  APF = quo(R1, ideal(rels))[1]
  z = APF(z)

  # construct the abstract variety

  PF = AbstractVariety(dim(X)+r-1, APF)
  p_pushforward = x -> X(pfwd(div(simplify(x).f, simplify(PF(z^(r-1))).f)))
  p_pushforward = MapFromFunc(chow_ring(PF), chow_ring(X), p_pushforward)
  p = AbstractVarietyMap(PF, X, PF.(imgs_in_R1), p_pushforward)
  if isdefined(X, :point)
    PF.point = pullback(p, point_class(X)) * z^(r-1)
  end
  p.O1 = PF(z)
  PF.O1 = PF(z)
  S = AbstractBundle(PF, 1, 1-z)
  Q = pullback(p, F) - S
  p.T = dual(S)*Q
  isdefined(X, :T) && (PF.T = pullback(p, tangent_bundle(X)) + tangent_bundle(p))
  PF.bundles = [S, Q]
  PF.structure_map = p
  set_attribute!(PF, :description => "Projectivization of $F")
  set_attribute!(PF, :grassmannian => :relative)
  get_attribute(X, :alg) == true && set_attribute!(PF, :alg => true)
  return PF
end

@doc raw"""
    abstract_hirzebruch_surface(n::Int)

Return the `n`-th Hirzebruch surface.

Recall that the `n`-th Hirzebruch surface is the projective bundle associated to
the bundle $\mathcal O_{\mathbb P_1} \oplus O_{\mathbb P_1(-n)}$.

# Examples

```jldoctest
julia> H2 =  abstract_hirzebruch_surface(2)
AbstractVariety of dim 2

julia> chow_ring(H2)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    z -> [1]
    h -> [1]
  by ideal (h^2, z^2 - 2*z*h)

```
"""
function abstract_hirzebruch_surface(n::Int)
  P1 = abstract_projective_space(1)
  E = OO(P1) + OO(P1, -n)
  return projective_bundle(E)
end

###############################################################################

@doc raw"""
    abstract_curve(g::IntegerUnion; base::Ring = QQ)

Return an abstract smooth curve of genus `g`.

The Chow ring is $\mathbb{Q}[p]/(p^2)$ where $p$ is the point class.
The tangent bundle has total Chern class $1 + (2-2g) \cdot p$, reflecting
the Euler number $\chi_{\text{top}} = 2 - 2g$.

# Examples
```jldoctest
julia> C = abstract_curve(0)
AbstractVariety of dim 1

julia> euler_number(C)
2

julia> euler_characteristic(OO(C))
1

julia> C2 = abstract_curve(2)
AbstractVariety of dim 1

julia> euler_number(C2)
-2

```
"""
function abstract_curve(g::IntegerUnion; base::Ring=QQ)
  R, (p,) = graded_polynomial_ring(base, ["p"], [1])
  I = ideal(R, [p^2])
  C = AbstractVariety(1, quo(R, I)[1])
  p = gens(chow_ring(C))[1]
  C.point = p
  C.O1 = chow_ring(C)(0)
  # Tangent bundle: rank 1, total Chern class 1 + (2-2g)·p
  C.T = AbstractBundle(C, 1, 1 + (2 - 2*g)*p)
  C.structure_map = map(C, abstract_point(base=base), [C(0)])
  set_attribute!(C, :description, "Curve of genus $g")
  set_attribute!(C, :alg, true)
  return C
end

###############################################################################

@doc raw"""
    abstract_K3_surface(g::IntegerUnion; base::Ring = QQ)

Return an abstract K3 surface polarized by a line bundle $L$ of genus $g$.

The *genus* $g$ of a polarized K3 surface $(S, L)$ is related to the degree of $L$ by
the formula $L^2 = 2g - 2$. For example, a quartic surface in $\mathbb{P}^3$ has $g = 3$.

The Chow ring is generated by the divisor class $L$ and the point class $p$ in degrees
1 and 2, with the single relation $L^2 = (2g-2) \cdot p$.

The tangent bundle has total Chern class $1 + 24p$ (trivial $c_1$, $c_2 = 24$), reflecting
$\chi_{\text{top}}(S) = 24$ for any K3 surface.

# Examples
```jldoctest
julia> S = abstract_K3_surface(3)
AbstractVariety of dim 2

julia> L = polarization(S)
L

julia> integral(L^2)
4

julia> euler_number(S)
24

julia> euler_characteristic(OO(S))
2

julia> euler_characteristic(OO(S, L))
4

```
"""
function abstract_K3_surface(g::IntegerUnion; base::Ring=QQ)
  R, (L, p) = graded_polynomial_ring(base, ["L", "p"], [1, 2])
  I = ideal(R, [L^2 - (2*g - 2)*p, L*p, p^2])
  S = AbstractVariety(2, quo(R, I)[1])
  L, p = gens(chow_ring(S))
  S.point = p
  S.O1 = L
  # Tangent bundle of a K3 surface: c₁ = 0 (trivial canonical class), c₂ = 24 (Euler number).
  # Total Chern class: c(T_S) = 1 + c₁ + c₂ = 1 + 24p.
  S.T = AbstractBundle(S, 2, 1 + 24*p)
  set_attribute!(S, :description, "K3 surface of genus $g")
  set_attribute!(S, :alg, true)
  return S
end

@doc raw"""
    abstract_quadric(n::Int; base::Ring = QQ)

Return an abstract smooth quadric hypersurface of dimension $n$.

For odd $n$, the quadric is constructed as a degree-2 complete intersection in $\mathbb{P}^{n+1}$.

For even $n = 2m$, the Chow ring has an extra generator in the middle degree:
the two families of $m$-dimensional linear subspaces on the quadric give classes
$l_1$ and $l_2$, with $h^m = l_1 + l_2$.

# Examples
```jldoctest
julia> Q1 = abstract_quadric(1)
AbstractVariety of dim 1

julia> euler_number(Q1)
2

julia> Q2 = abstract_quadric(2)
AbstractVariety of dim 2

julia> euler_number(Q2)
4

julia> betti_numbers(Q2)
3-element Vector{Int64}:
 1
 2
 1

julia> Q3 = abstract_quadric(3)
AbstractVariety of dim 3

julia> euler_number(Q3)
4

```
"""
function abstract_quadric(n::Int; base::Ring=QQ)
  @assert n >= 1 "Dimension must be at least 1"
  if isodd(n)
    # odd-dimensional quadric: complete intersection of degree 2 in P^{n+1}
    Pn1 = abstract_projective_space(n + 1, base=base)
    Q = complete_intersection(Pn1, 2)
    set_attribute!(Q, :description, "Quadric hypersurface of dim $n")
    return Q
  else
    m = div(n, 2)
    R, (h, l1, l2) = graded_polynomial_ring(base, ["h", "l1", "l2"], [1, m, m])
    p = (1 // ZZ(2)) * h^n
    # Relations depend on parity of m
    if isodd(m)
      rels = [h^m - (l1 + l2), h * (l1 - l2), l1^2, l2^2, l1*l2 - p]
    else
      rels = [h^m - (l1 + l2), h * (l1 - l2), l1^2 - p, l2^2 - p, l1*l2]
    end
    # Add dimension-trimming relations
    rels = vcat(rels, [h^(n + 1)])
    I = ideal(R, rels)
    Q = AbstractVariety(n, quo(R, I)[1])
    h, l1, l2 = gens(chow_ring(Q))
    p = (1 // ZZ(2)) * h^n
    Q.O1 = h
    Q.point = p
    # Tangent bundle: pull back from P^{n+1} and subtract the normal bundle O(2)
    Pn1 = abstract_projective_space(n + 1, base=base)
    Q.structure_map = map(Q, Pn1, [h])
    Q.T = pullback(Q.structure_map, tangent_bundle(Pn1)) - OO(Q, 2)
    set_attribute!(Q, :description, "Quadric hypersurface of dim $n")
    set_attribute!(Q, :alg, true)
    return Q
  end
end

@doc raw"""
    abstract_cayley_plane(; base::Ring = QQ)

Return the abstract Cayley plane $\mathbf{OP}^2$, also known as the
$E_6$-Grassmannian.

This is a 16-dimensional smooth projective variety. Its Chow ring is generated by
classes `H, σ₄', σ₈, p` in degrees `1, 4, 8, 16` with explicit relations
from [IM05](@cite).

# Examples
```jldoctest
julia> X = abstract_cayley_plane()
AbstractVariety of dim 16

julia> euler_number(X)
27

julia> betti_numbers(X)
17-element Vector{Int64}:
 1
 1
 1
 1
 2
 2
 2
 2
 3
 2
 2
 2
 2
 1
 1
 1
 1

```
"""
function abstract_cayley_plane(; base::Ring=QQ)
  R, (H, s4p, s8, p) = graded_polynomial_ring(base, ["H", "s4p", "s8", "p"], [1, 4, 8, 16])
  s4pp = H^4 - s4p
  s6p = 2*s4p*H^2 - s4pp*H^2
  s6pp = s4pp*H^2 - s4p*H^2
  s8p = (s4p*H^4 - s8)*3 - (s4pp*H^4 - s8)*2
  s8pp = (s4pp*H^4 - s8)*3 - (s4p*H^4 - s8)*4
  s12p = s4p*s8
  s12pp = s4pp*s8
  rels = [s8^2 - p, s8p^2 - p, s8pp^2 - p, s8*s8p, s8p*s8pp, s8*s8pp,
    s4p*s12p - p, s4pp*s12pp - p, s4pp*s12p, s4p*s12pp,
    s4p^2 - s8 - s8p - s8pp, s4pp^2 - s8 - 2*s8p - 2*s8pp,
    s4p*s4pp - 2*s8p - s8pp,
    6*H*s4p*s8 - H*s4p^3,
    s8pp*H^4 - 2*s12p - 3*s12pp,
    s8p*H^4 - 3*s12p - 4*s12pp,
    2*H*s8 - 5*H*s4p^2 + 2*H^5*s4p,
    H^17, p^2]
  I = ideal(R, rels)
  X = AbstractVariety(16, quo(R, I)[1])
  H, s4p, s8, p = gens(chow_ring(X))
  s4pp = H^4 - s4p
  s6p = 2*s4p*H^2 - s4pp*H^2
  s6pp = s4pp*H^2 - s4p*H^2
  s8p = (s4p*H^4 - s8)*3 - (s4pp*H^4 - s8)*2
  s8pp = (s4pp*H^4 - s8)*3 - (s4p*H^4 - s8)*4
  X.O1 = H
  X.point = p
  # Normal bundle of embedding in P^26: rank 10, with explicit total Chern class
  N_chern =
    1 + 15*H + 102*H^2 + 414*H^3 +
    1107*s4p + 1113*s4pp +
    2025*s4p*H + 2079*s4pp*H +
    5292*s6p + 8034*s6pp +
    4698*s6p*H + 7218*s6pp*H +
    2751*s8 + 9786*s8p + 7032*s8pp +
    963*s8*H + 3438*s8p*H + 2466*s8pp*H +
    153*s8*H^2 + 549*s8p*H^2 + 387*s8pp*H^2
  N = AbstractBundle(X, 10, N_chern)
  P26 = abstract_projective_space(26, base=base)
  X.structure_map = map(X, P26, [H])
  X.T = pullback(X.structure_map, tangent_bundle(P26)) - N
  set_attribute!(X, :description, "Cayley plane OP^2")
  set_attribute!(X, :alg, true)
  return X
end

@doc raw"""
    abstract_cayley_grassmannian(; base::Ring = QQ)

Return the abstract Cayley Grassmannian $\mathbf{CG}$.

This is an 8-dimensional smooth projective variety, constructed as the
zero locus of a section of $\bigwedge^3 S^*$ on $\mathrm{Gr}(4, 7)$,
where $S$ is the tautological subbundle.

# Examples
```jldoctest
julia> CG = abstract_cayley_grassmannian()
AbstractVariety of dim 8

julia> euler_number(CG)
15

```
"""
function abstract_cayley_grassmannian(; base::Ring=QQ)
  G = abstract_grassmannian(4, 7, base=base)
  S = tautological_bundles(G)[1]
  CG = zero_locus_section(exterior_power(dual(S), 3))
  set_attribute!(CG, :description, "Cayley Grassmannian CG")
  return CG
end

@doc raw"""
    abstract_grassmannian(k::Int, n::Int; base::Ring = QQ, symbol::String = "c")

Return the abstract Grassmannian $\mathrm{Gr}(k, n)$ of `k`-dimensional subspaces of an
`n`-dimensional vector space.

!!! note
    The string `symbol` specifies how to print the generators of the Chow ring.

!!! tip
    There are several useful ways of representing the Chow ring of $\mathrm{Gr}(k, n)$ in terms
    of generators and relations. The function `abstract_grassmannian` yields a minimal
    set of generators, the Chern classes of the universal subbundle on $\mathrm{Gr}(k, n)$, with
    relations as described, for example, in [EH16](@cite). Another way is to handle the Grassmannian as a
    special case of a flag variety: Entering `abstract_flag_variety(k, n)`, the returned generators
    will be the Chern classes of the tautological subquotient bundles, with relations as described,
    for example, in [HK-MW24](@cite). Depending on whether $n-k$ is small (large), working with the
    latter (former) recipe may be more performant.

# Examples
```jldoctest
julia> G = abstract_grassmannian(2,4)
AbstractVariety of dim 4

julia> CR = chow_ring(G)
Quotient
  of multivariate polynomial ring in 2 variables over QQ graded by
    c[1] -> [1]
    c[2] -> [2]
  by ideal (-c[1]^3 + 2*c[1]*c[2], c[1]^4 - 3*c[1]^2*c[2] + c[2]^2)

julia> S = tautological_bundles(G)[1]
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> V = [chern_class(S, i) for i in 1:2]
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 c[1]
 c[2]

julia> is_regular_sequence(gens(modulus(CR)))
true

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 4

julia> tangent_bundle(G) == dual(S)*Q
true

```

```jldoctest
julia> G = abstract_flag_variety(2,4)
AbstractVariety of dim 4

julia> CR = chow_ring(G)
Quotient
  of multivariate polynomial ring in 4 variables over QQ graded by
    c[1, 1] -> [1]
    c[1, 2] -> [2]
    c[2, 1] -> [1]
    c[2, 2] -> [2]
  by ideal with 4 generators

julia> modulus(CR)
Ideal generated by
  -c[1, 2]*c[2, 2]
  -c[1, 1]*c[2, 2] - c[1, 2]*c[2, 1]
  -c[1, 1]*c[2, 1] - c[1, 2] - c[2, 2]
  -c[1, 1] - c[2, 1]

```
"""
function abstract_grassmannian(k::Int, n::Int; base::Ring=QQ, symbol::String="c")
  @assert k < n
  d = k*(n-k)
  R, c = graded_polynomial_ring(base, symbol => 1:k; weights=1:k)

  l = [-sum(c)]
  for i in 1:(n - 1)
    push!(l, l[1]*l[end])
  end
  inv_c = sum(l) # this is c(Q) since c(S) * c(Q) = 1

  # Q is of rank n-k: the vanishing of Chern classes in higher degrees provides all the relations for the Chow ring

  AGr = quo(R, ideal(inv_c[(n - k + 1):n]))[1]
  c = gens(AGr)
  Gr = AbstractVariety(d, AGr)
  Gr.O1 = Gr(-c[1])
  S = AbstractBundle(Gr, k, 1 + sum(c))
  Q = trivial_line_bundle(Gr)*n - S
  Q.chern = 1 + Gr(sum(inv_c[1:(n - k)]))
  Gr.point = Gr((-1)^d*c[end]^(n-k))
  Gr.T = dual(S) * Q
  Gr.bundles = [S, Q]
  Gr.structure_map = map(Gr, abstract_point(base=base), [Gr(0)])
  set_attribute!(Gr, :description => "Grassmannian Gr($k, $n)")
  set_attribute!(Gr, :grassmannian => :absolute)
  set_attribute!(Gr, :alg => true)
  return Gr
end

@doc raw"""
    abstract_flag_variety(dims::Int...; base::Ring = QQ, symbol::String = "c")
    abstract_flag_variety(dims::Vector{Int}; base::Ring = QQ, symbol::String = "c")

Given integers, say, $d_1, \dots, d_{k}, n$ with $0 < d_1 < \dots < d_{k} < n$ or a vector of such integers,
return the abstract flag variety $\mathrm{F}(d_1, \dots, d_{k}; n)$ of nested sequences of subspaces of
dimensions $d_1, \dots, d_{k}$ of an $n$-dimensional vector space.

!!! note
    The string `symbol` specifies how to print the generators of the Chow ring.

!!! note
    The returned generators are the Chern classes of the tautological subquotient bundles, with relations as described,
    for example, in [HK-MW24](@cite).

# Examples
```jldoctest
julia> F = abstract_flag_variety(1,3,4)
AbstractVariety of dim 5

julia> chow_ring(F)
Quotient
  of multivariate polynomial ring in 4 variables over QQ graded by
    c[1, 1] -> [1]
    c[2, 1] -> [1]
    c[2, 2] -> [2]
    c[3, 1] -> [1]
  by ideal with 4 generators

julia> modulus(chow_ring(F))
Ideal generated by
  -c[1, 1]*c[2, 2]*c[3, 1]
  -c[1, 1]*c[2, 1]*c[3, 1] - c[1, 1]*c[2, 2] - c[2, 2]*c[3, 1]
  -c[1, 1]*c[2, 1] - c[1, 1]*c[3, 1] - c[2, 1]*c[3, 1] - c[2, 2]
  -c[1, 1] - c[2, 1] - c[3, 1]

```
"""
function abstract_flag_variety(
  dims::Int...; bott::Bool=false, weights=:int, base::Ring=QQ, symbol::String="c"
)
  abs_flag(collect(dims), base=base, symbol=symbol)
end

function abstract_flag_variety(
  dims::Vector{Int}; bott::Bool=false, weights=:int, base::Ring=QQ, symbol::String="c"
)
  abs_flag(dims; base=base, symbol=symbol)
end

function abs_flag(dims::Vector{Int}; base::Ring=QQ, symbol::String="c")
  n, l = dims[end], length(dims)
  ranks = pushfirst!([dims[i + 1]-dims[i] for i in 1:(l - 1)], dims[1])
  @assert all(>(0), ranks)
  d = sum(ranks[i] * sum(dims[end]-dims[i]) for i in 1:(l - 1))
  syms = reduce(vcat, [_parse_symbol(symbol, i, 1:r) for (i, r) in enumerate(ranks)])
  # FIXME ordering
  # ord = prod(ordering_dp(r) for r in ranks)
  R = graded_polynomial_ring(base, syms, reduce(vcat, [1:r for r in ranks]))[1]
  c = pushfirst!(
    [1+sum(gens(R)[(dims[i] + 1):dims[i + 1]]) for i in 1:(l - 1)],
    1+sum(gens(R)[1:dims[1]]),
  )
  gi = prod(c)[0:n]
  # XXX cannot mod using graded ring element
  Rx, x = R.R[:x]
  g = sum(gi[i + 1].f * x^(n-i) for i in 0:n)
  q = mod(x^n, g)
  rels = [R(coeff(q, i)) for i in 0:(n - 1)]
  AFl = quo(R, ideal(rels))[1]
  c = AFl.(c)
  Fl = AbstractVariety(d, AFl)
  Fl.bundles = [AbstractBundle(Fl, r, ci) for (r, ci) in zip(ranks, c)]
  Fl.O1 = simplify(sum((i-1)*chern_class(tautological_bundles(Fl)[i], 1) for i in 1:l))
  Fl.point = prod(
    top_chern_class(E)^dims[i] for (i, E) in enumerate(tautological_bundles(Fl)[2:end])
  )
  Fl.T = sum(
    dual(tautological_bundles(Fl)[i]) *
    sum([tautological_bundles(Fl)[j] for j in (i + 1):l]) for i in 1:(l - 1)
  )
  Fl.structure_map = map(Fl, abstract_point(; base=base), [Fl(0)])
  set_attribute!(Fl, :description => "Flag abstract variety Flag$(tuple(dims...))")
  if l == 2
    set_attribute!(Fl, :grassmannian => :absolute)
  end
  set_attribute!(Fl, :alg => true)
  # if all(r->r==1, ranks)
  #   set_attribute!(Fl, :weyl_group => WeylGroup("A$(n-1)"))
  #   set_attribute!(Fl, :roots => -[c[i] - c[i+1] for i in 1:n-1])
  # end
  return Fl
end

@doc raw"""
    flag_bundle(F::AbstractBundle, dims::Int...; symbol::String = "c")
    flag_bundle(F::AbstractBundle, dims::Vector{Int}; symbol::String = "c")

Given integers, say, $d_1, \dots, d_{k}, n$ with $0 < d_1 < \dots < d_{k} < n$ or a vector of such integers,
and given an abstract bundle $F$ of rank $n$, return the abstract flag bundle of nested sequences
of subspaces of dimensions $d_1, \dots, d_{k}$ in the fibers of $F$.

!!! note
    Entering the number $n$ can be omitted since this number can be recovered as the rank of $F$.

!!! note
    Generators and relations for the Chow ring of a flag bundle are described in [Gro58](@cite).
    For the pushforward of cycle classes from a flag bundle to its base variety see [GSS22](@cite).

# Examples
```jldoctest
julia> P = abstract_projective_space(4)
AbstractVariety of dim 4

julia> F = exterior_power(cotangent_bundle(P),  3)*OO(P,3)
AbstractBundle of rank 4 on AbstractVariety of dim 4

julia> FB = flag_bundle(F, 1, 3)
AbstractVariety of dim 9

julia> CR = chow_ring(FB)
Quotient
  of multivariate polynomial ring in 5 variables over QQ graded by
    c[1, 1] -> [1]
    c[2, 1] -> [1]
    c[2, 2] -> [2]
    c[3, 1] -> [1]
    h -> [1]
  by ideal with 5 generators

julia> modulus(CR)
Ideal generated by
  h^5
  -c[1, 1]*c[2, 2]*c[3, 1] + h^4
  -c[1, 1]*c[2, 1]*c[3, 1] - c[1, 1]*c[2, 2] - c[2, 2]*c[3, 1] - 2*h^3
  -c[1, 1]*c[2, 1] - c[1, 1]*c[3, 1] - c[2, 1]*c[3, 1] - c[2, 2] + 4*h^2
  -c[1, 1] - c[2, 1] - c[3, 1] - 3*h

julia> tautological_bundles(FB)
3-element Vector{AbstractBundle}:
 AbstractBundle of rank 1 on AbstractVariety of dim 9
 AbstractBundle of rank 2 on AbstractVariety of dim 9
 AbstractBundle of rank 1 on AbstractVariety of dim 9

```
"""
function flag_bundle(F::AbstractBundle, dims::Int...; symbol::String="c")
  flag_bundle(F, collect(dims), symbol=symbol)
end

function flag_bundle(F::AbstractBundle, dims::Vector{Int}; symbol::String="c")
  X, n = parent(F), rank(F)
  !(n isa Int) && error("expect rank to be an integer")

  # compute the ranks of successive subquotients and the relative dimension

  if dims[end] < n # the last dim can be omitted
    dims = vcat(dims, [n])
  end
  l = length(dims)
  ranks = pushfirst!([dims[i + 1]-dims[i] for i in 1:(l - 1)], dims[1])
  @assert all(>(0), ranks) && dims[end] <= n
  d = sum(ranks[i] * sum(dims[end]-dims[i]) for i in 1:(l - 1))

  # construct the ring

  R = chow_ring(X)
  if R isa MPolyQuoRing
    PR = base_ring(R)
  else
    PR = R
  end
  syms = reduce(vcat, [_parse_symbol(symbol, i, 1:r) for (i, r) in enumerate(ranks)])
  w = reduce(vcat, [1:r for r in ranks])
  append!(w, gradings(R))
  R1, gens_for_rels_R1, imgs_in_R1 = graded_polynomial_ring(
    base(X), syms, symbols(PR); weights=w
  )
  pback = Oscar.hom(PR, R1, imgs_in_R1)
  pfwd = Oscar.hom(R1, R, vcat(repeat([R()], n), gens(R)))

  # compute the relations

  c = [1+sum(gens_for_rels_R1[(dims[i] + 1):dims[i + 1]]) for i in 1:(l - 1)]
  pushfirst!(c, 1+sum(gens_for_rels_R1[1:dims[1]]))
  Rx, x = R1[:x]
  fi = pback(total_chern_class(F).f)[0:n]
  f = sum(fi[i + 1].f * x^(n-i) for i in 0:n)
  gi = prod(c)[0:n]
  g = sum(gi[i + 1].f * x^(n-i) for i in 0:n)
  rels = [R1(coeff(mod(f, g), i)) for i in 0:(n - 1)]
  if R isa MPolyQuoRing
    rels = vcat(pback.(gens(R.I)), rels)
  end
  AFl = quo(R1, ideal(rels))[1]
  c = AFl.(c)

  # construct the abstract variety

  Fl = AbstractVariety(dim(X) + d, AFl)
  Fl.bundles = [AbstractBundle(Fl, r, ci) for (r, ci) in zip(ranks, c)]
  section = prod(
    top_chern_class(E)^sum(dims[i]) for (i, E) in enumerate(tautological_bundles(Fl)[2:end])
  )### TODO skip when fglm is used
  if isdefined(X, :point)
    Fl.point = pback(point_class(X).f) * section ### TODO better when fglm is used (take section from _find_sect as second output
  end
  p_pullback = Fl.(imgs_in_R1)

  ###p_pushforward = x -> (@warn("possibly wrong ans"); X(pfwd(div(simplify(x).f, simplify(section).f))))

  ##################################
  ##alternative approach
  ##TODO: To be discussed. Use fglm in _find_sect as soon as fglm is fixed.

  ME = _matrix_of_exponents(ranks)
  nl = size(ME)[1]
  wcs = reduce(vcat, [1:r for r in ranks])
  Rcs, _ = graded_polynomial_ring(base(X), syms; weights=wcs)
  RcstoR1 = Oscar.hom(Rcs, R1, gens_for_rels_R1)
  gs = [monomial(Rcs, [Int(ME[i, j]) for j in 1:n]) for i in nl:-1:1]
  gs = [RcstoR1(gs[i]) for i in 1:length(gs)]
  ds = [degree(Int, gs[i]) for i in 1:1:nl]
  dm = argmax(ds)
  ### TODO find and return section
  fm = Oscar.hom(R, AFl, p_pullback)
  p_pushforward = x -> X(_find_sect(fm, gs)(simplify(x).f)[dm])
  ##################################

  p_pushforward = MapFromFunc(chow_ring(Fl), chow_ring(X), p_pushforward)
  p = AbstractVarietyMap(Fl, X, p_pullback, p_pushforward)
  p.O1 = simplify(sum((i-1)*chern_class(tautological_bundles(Fl)[i], 1) for i in 1:l))
  Fl.O1 = p.O1
  p.T = sum(
    dual(tautological_bundles(Fl)[i]) *
    sum([tautological_bundles(Fl)[j] for j in (i + 1):l]) for i in 1:(l - 1)
  )
  if isdefined(X, :T)
    Fl.T = pullback(p, tangent_bundle(X)) + tangent_bundle(p)
  end
  Fl.structure_map = p
  set_attribute!(
    Fl, :description => "Relative flag abstract variety Flag$(tuple(dims...)) for $F"
  )
  set_attribute!(Fl, :section => section)
  get_attribute(X, :alg) == true && set_attribute!(Fl, :alg => true)

  l == 2 && set_attribute!(Fl, :grassmannian => :relative)
  return Fl
end

@doc raw"""
    section_class(f::AbstractVarietyMap)

Return the section class of the flag bundle map `f`. The section class is a
cycle class on the flag variety of the correct codimension such that its
push-forward to the base is 1.

In particular, `pushforward(f, section_class(f)) == 1` in the Chow ring of the
base. Integrating that class then depends on the base dimension; for example,
on `\mathbb{P}^3` it is `0`, while multiplying by the point class gives `1`.

This corresponds to `sectionClass` in Schubert2 (Macaulay2), where the method
is defined for abstract variety maps in general. In this implementation, the
class is currently provided for maps whose source carries a precomputed section
class attribute (notably structure maps of flag bundles).

# Examples
```jldoctest
julia> X = abstract_projective_space(3);

julia> Fl = flag_bundle(tangent_bundle(X), 1, 2);

julia> p = structure_map(Fl);

julia> sc = section_class(p);

julia> pushforward(p, sc)
1

julia> integral(pushforward(p, sc))
0

julia> integral(pushforward(p, sc) * point_class(X))
1

```

```jldoctest
julia> pt = abstract_point();

julia> F = 4*trivial_line_bundle(pt);

julia> G = flag_bundle(F, 2);

julia> f = structure_map(G);

julia> sc = section_class(f);

julia> pushforward(f, sc)
1

julia> integral(sc)
1

```
"""
function section_class(f::AbstractVarietyMap)
  get_attribute(domain(f), :section)
end

###########################################################
#helper functions above
###########################################################

function _matrix_of_exponents(ni::Vector{Int}) # [GSS22](@cite) Thm 2.15
  n = sum(ni)
  r = length(ni)
  A = zero_matrix(ZZ, n, n)
  b = zero_matrix(ZZ, n, 1)

  C = zero_matrix(ZZ, r, n)
  d = zero_matrix(ZZ, r, 1)
  s = 0
  for i in 1:r
    for j in 1:ni[i]
      C[i, j + s] = 1
    end
    s += ni[i]
    d[i] = n-sum(ni[i:r])
  end

  C = vcat(-C, identity_matrix(ZZ, n))
  d = vcat(-d, zero_matrix(ZZ, n, 1))

  return solve_mixed(A, b, C, d)
end

function _find_sect(F::Oscar.AffAlgHom, gs::Vector) # see function present_finite_extension_ring
  A, B = F.domain, F.codomain
  a, b = ngens(A), ngens(B)

  if A isa MPolyQuoRing
    AR = base_ring(A)
  else
    AR = A
  end
  if B isa MPolyQuoRing
    BR = base_ring(B)
    M = [F(gens(A)[i]).f for i in 1:a]
  else
    BR = B
    M = [F(gens(A)[i]) for i in 1:a]
  end

  @assert base_ring(AR) == base_ring(BR)

  I = ideal(BR, isdefined(B, :I) ? vcat(gens(B.I), M) : M)
  C, _ = quo(BR, I)
  @assert gs[end] == 1 # the last one should always be 1
  g = length(gs)

  R, _ = tensor_product(BR, AR; use_product_ordering=true)
  ba = gens(R)
  ARtoR = Oscar.hom(AR, R, ba[(b + 1):end])
  BRtoR = Oscar.hom(BR, R, ba[1:b])
  RtoAR = Oscar.hom(R, AR, vcat(repeat([AR()], b), gens(AR)))
  gs_lift = [BRtoR(g) for g in gs]

  # compute the ideal J of the graph of F

  Rels = [ba[b + i]-BRtoR(m) for (i, m) in enumerate(M)]
  if isdefined(A, :I)
    for g in gens(A.I)
      push!(Rels, ARtoR(g))
    end
  end
  if isdefined(B, :I)
    for g in gens(B.I)
      push!(Rels, BRtoR(g))
    end
  end
  J = ideal(R, Rels) # the ideal of the graph of F
  V = groebner_basis(J) ###TODO replace by fglm below as soon as fglm is fixed
  ###W = vcat(weights(Int, BR), weights(Int, AR)) ####
  ###V = fglm(J, start_ordering = wdegrevlex(R, W), destination_ordering = default_ordering(R))

  sect =
    x -> (y=reduce(BRtoR(x), gens(V); complete_reduction=true);
      ans=elem_type(AR)[];
      for i in 1:g
        q = div(y, gs_lift[i])
        push!(ans, RtoAR(q))
        y -= q * gs_lift[i]
      end; ans)
  return sect
end

###############################################################################
#
#
function _parse_symbol(symbol::String, I::AbstractUnitRange)
  return [string(symbol, "[", i, "]") for i in I]
end

function _parse_symbol(symbol::String, n::Int, I::AbstractUnitRange)
  return [string(symbol, "[", n, ", ", i, "]") for i in I]
end
