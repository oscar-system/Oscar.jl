#module MPolyModule

##############################################################################
#
# could/ should be in AbstractAlgebra
#
# some sugar to make creation of strange rings easier
# possibly lacks the passing of the ordering...
##############################################################################

#TODO: reduce = divrem in Nemo. Should be faster - if we have the correct basis

# To make a list of variables safe for Singular to use
# Allowable schemes should be:
# (:x), (:x, :y), (:x, :y, :z)
# (:x1, :x2, ...)
function _variables_for_singular(n::Int)
  n > 3 && return [Symbol("x[$i]") for i in 1:n]
  return [ :x, :y, :z ][1:n]
end
_variables_for_singular(S::Vector{Symbol}) = _variables_for_singular(length(S))


######################################################################
# pretty printing for iJulia notebooks..
#

function Base.show(io::IO, mime::IJuliaMime, R::MPolyRing)
  io = IOContext(io, :compact => true)
  print(io, "\$")
  math_html(io, R)
  print(io, "\$")
end

function math_html(io::IO, R::MPolyRing)
  print(io, "\\text{Multivariate Polynomial Ring in $(nvars(R)) variables:} ")
  math_html(io, gens(R))
  print(io, "\\text{ over }")
  math_html(io, base_ring(R))
end

###################################################

using .Orderings


##############################################################################
#
# workhorse: IdealGens
# ideals are (mostly) generated on the Nemo side, but structural computations
# are in Singular. To avoid permanent conversion, the list of generators = sideal
# is captured in IdealGens: for Ocsar this is Vector{MPoly}
#                                 Singular      sideal
#
#TODO/ to think
#  default in Nemo is     :lex
#             Singular is :degrevlex -> better for std
#by default, use different orders???
#make IdealGens use different orders for both? make the type depend on it?
#
#for std: abstraction to allow Christian to be used
#
#type for orderings, use this...
#in general: all algos here needs revision: do they benefit from gb or not?

@doc raw"""
    default_ordering(R::MPolyRing)

Return the monomial ordering that is used for computations with ideals in `R`
if no other ordering is specified -- either directly by the user or by
requirements of a specific algorithm.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> default_ordering(R)
degrevlex([x, y, z])

julia> F = free_module(R, 2)
Free module of rank 2 over R

julia> default_ordering(F)
degrevlex([x, y, z])*lex([gen(1), gen(2)])

julia> S, _ = grade(R, [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> default_ordering(S)
wdegrevlex([x, y, z], [1, 2, 3])
```
"""
@attr MonomialOrdering{T} function default_ordering(R::T) where {T<:MPolyRing}
  return degrevlex(R)
end

# Only for internal use
function set_default_ordering!(R::MPolyRing, o::MonomialOrdering)
  @assert R === base_ring(o)
  set_attribute!(R, :default_ordering, o)
  return nothing
end

@doc raw"""
    with_ordering(f, R::MPolyRing, o::MonomialOrdering)

Use the monomial ordering `o` for computations in `R` during the execution of
`f`.
This may be used with `do` block syntax, see the example.

This functionality is meant for advanced users. In general it should not be
necessary to explicitly set a monomial ordering.
Further, there is no guarantee that `o` is actually used. For example, if an
algorithm requires an elimination ordering, `o` might be ignored.

# Example
```jldoctest withordering
julia> R, (x, y, z) = QQ[:x, :y, :z];

julia> f = x + y^2;

julia> I = ideal(R, [y^2 - z, x - z^2]);

julia> normal_form(f, I) # this uses degrevlex
x + z

julia> with_ordering(R, lex(R)) do
           # this uses lex
           normal_form(f, I)
       end
z^2 + z
```
Notice that in this small example we could have achieved the same by using the
keyword argument `ordering`:
```jldoctest withordering
julia> normal_form(f, I, ordering = lex(R))
z^2 + z
```
"""
function with_ordering(f, R::MPolyRing, o::MonomialOrdering)
  old = default_ordering(R)
  set_default_ordering!(R, o)
  x = try
    f()
  finally
    set_default_ordering!(R, old)
  end
  return x
end

mutable struct BiPolyArray{S}
  Ox::NCRing #Oscar Poly Ring or Algebra
  O::Vector{S}
  Sx::NCRing # Singular Poly Ring or Algebra, poss. with different ordering
  f #= isomorphism Ox -> Sx =#
  S::Singular.sideal

  function BiPolyArray(Ox::T) where {T <: NCRing}
    if T <: MPolyQuoRing
      return new{elem_type(base_ring_type(T))}(Ox)
    end
    return new{elem_type(T)}(Ox)
  end

  function BiPolyArray(O::Vector{<: NCRingElem})
    return BiPolyArray(parent(O[1]), O)
  end

  function BiPolyArray(Ox::NCRing, O::Vector{<: NCRingElem})
    r = BiPolyArray(Ox)
    r.O = O
    return r
  end

  function BiPolyArray(Ox::T, S::Singular.sideal) where {T <: NCRing}
    r = BiPolyArray(Ox)
    r.Sx = base_ring(S)
    r.S = S
    return r
  end
end

mutable struct IdealGens{S}
  gensBiPolyArray::BiPolyArray{S}
  isGB::Bool
  isReduced::Bool
  ord::Orderings.MonomialOrdering
  keep_ordering::Bool

  # internal constructor
  function IdealGens(B::BiPolyArray{S}, isGB::Bool, isReduced::Bool) where S
    return new{S}(B, isGB, isReduced)
  end

  function IdealGens(O::Vector{T}; keep_ordering::Bool = true) where {T <: NCRingElem}
    return IdealGens(parent(O[1]), O; keep_ordering = keep_ordering)
  end

  function IdealGens(O::Vector{T}, ordering::Orderings.MonomialOrdering; keep_ordering::Bool = true, isGB::Bool = false) where {T <: NCRingElem}
    return IdealGens(parent(O[1]), O, ordering; keep_ordering = keep_ordering, isGB)
  end

  function IdealGens(Ox::NCRing, O::Vector{T}, ordering::Orderings.MonomialOrdering; keep_ordering::Bool = true, isGB::Bool = false, isReduced::Bool = false) where {T <: NCRingElem}
    r = IdealGens(BiPolyArray(Ox, O), isGB, isReduced)
    r.ord = ordering
    r.keep_ordering = keep_ordering
    return r
  end

  function IdealGens(Ox::NCRing, O::Vector{T}; keep_ordering::Bool = true) where {T <: NCRingElem}
    r = IdealGens(BiPolyArray(Ox, O), false, false)
    r.keep_ordering = keep_ordering
    return r
  end

  function IdealGens(Ox::T, S::Singular.sideal, isReduced::Bool = false) where {T <: NCRing}
    r = IdealGens(BiPolyArray(Ox, S), S.isGB, isReduced)
    if T <: MPolyQuoRing
      r.ord = Ox.ordering
    else
      r.ord = monomial_ordering(Ox, Singular.ordering(base_ring(S)))
    end
    r.keep_ordering = true
    return r
  end
end

function show(io::IO, ::MIME"text/plain", I::IdealGens)
  io = pretty(io)
  if I.isGB
    if is_global(I.ord)
      print(io, LowercaseOff(), "Gröbner basis")
    else
      print(io, "Standard basis")
    end
    print(io, " with elements", Indent())
    for (i, g) in enumerate(gens(I))
        print(io, "\n", i, ": ", OscarPair(g, I.ord))
    end
    print(io, Dedent())
    print(io, "\nwith respect to the ordering")
    print(io, Indent(), "\n", I.ord, Dedent())
  else
    print(io, "Ideal generating system with elements")
    print(io, Indent())
    for (i,g) in enumerate(gens(I))
      print(io, "\n", i, ": ", g)
    end
    print(io, Dedent())
    if isdefined(I, :ord)
      print(io, "\nwith associated ordering")
      print(io, Indent(), "\n", I.ord, Dedent())
    end
  end
end

function show(io::IO, I::IdealGens)
  if I.isGB
    io = pretty(io)
    if is_global(I.ord)
      print(io, LowercaseOff(), "Gröbner basis")
    else
      print(io, "Standard basis")
    end
    if !is_terse(io)
      print(io, " with $(ItemQuantity(length(I), "element"))")
      print(io, " w.r.t. ", I.ord)
    end
  else
    print(io, "Ideal generating system")
    if !is_terse(io)
      print(io, " with $(ItemQuantity(length(I), "element"))")
      print(io, " with associated ordering ", I.ord)
    end
  end
end


function Base.length(A::BiPolyArray)
  if isdefined(A, :S)
    return Singular.ngens(A.S)
  else
    return length(A.O)
  end
end

function Base.iterate(A::BiPolyArray, s::Int = 1)
  if s > length(A)
    return nothing
  end
  return oscar_generators(A)[s], s+1
end

Base.eltype(::Type{<:BiPolyArray{S}}) where S = S

function gen(A::IdealGens, i::Int)
  return oscar_generators(A)[i]
end

Base.getindex(A::IdealGens, i::Int) = gen(A, i)


function Base.length(A::IdealGens)
  return length(A.gensBiPolyArray)
end

function base_ring(A::IdealGens)
  return A.gensBiPolyArray.Ox
end

function Base.iterate(A::IdealGens, s::Int = 1)
  if s > length(A)
    return nothing
  end
  return gen(A, s), s+1
end

Base.eltype(::Type{IdealGens{S}}) where S = S

function Base.keys(A::IdealGens)
  return 1:length(A)
end

function gens(I::IdealGens)
  return collect(I)
end

function elements(I::IdealGens)
  return collect(I)
end

function ordering(G::IdealGens)
  isdefined(G, :ord) ? G.ord : error("The ideal generating system does not have an associated ordering")
end

# for internal use only
function set_ordering!(G::IdealGens, monord::MonomialOrdering)
  @assert base_ring(G) == monord.R "Base rings of IdealGens and MonomialOrdering are inconsistent"
  isdefined(G, :ord) && G.ord !== monord && error("Monomial ordering is already set to a different value")
  G.ord = monord
end

function singular_generators(B::IdealGens, monorder::MonomialOrdering=default_ordering(base_ring(B)))
  if !isdefined(B.gensBiPolyArray, :S)
    g = iso_oscar_singular_poly_ring(base_ring(B); keep_ordering = B.keep_ordering)
    B.gensBiPolyArray.Sx = codomain(g)
    B.gensBiPolyArray.f = g
    B.gensBiPolyArray.S = Singular.Ideal(B.gensBiPolyArray.Sx, elem_type(B.gensBiPolyArray.Sx)[g(x) for x in oscar_generators(B)])
  end
  if B.isGB && B.ord == monomial_ordering(base_ring(B), internal_ordering(B.gensBiPolyArray.Sx))
    B.gensBiPolyArray.S.isGB = true
  end

  # in case of quotient rings, monomial ordering is ignored so far in singular_poly_ring
  isa(base_ring(B), MPolyQuoRing) && return B.gensBiPolyArray.S
  isdefined(B, :ord) && B.ord == monorder && monomial_ordering(base_ring(B), Singular.ordering(base_ring(B.gensBiPolyArray.S))) == B.ord && return B.gensBiPolyArray.S
  g = iso_oscar_singular_poly_ring(base_ring(B), monorder)
  SR = codomain(g)
  f = Singular.AlgebraHomomorphism(B.gensBiPolyArray.Sx, SR, gens(SR))
  S = Singular.map_ideal(f, B.gensBiPolyArray.S)
  if isdefined(B, :ord) && B.ord == monorder
    S.isGB = B.isGB
  end
  return S
end

@doc raw"""
    set_ordering(I::IdealGens, monord::MonomialOrdering)

Return an ideal generating system with an associated monomial ordering.

# Examples
```jldoctest
julia> R, (x0, x1, x2) = polynomial_ring(QQ, [:x0,:x1,:x2]);

julia> I = ideal([x0*x1, x2]);

julia> g = generating_system(I);

julia> set_ordering(g, degrevlex(gens(R)))
Ideal generating system with elements
  1: x0*x1
  2: x2
with associated ordering
  degrevlex([x0, x1, x2])
```
"""
function set_ordering(G::IdealGens, monord::MonomialOrdering)
  @assert base_ring(G) == monord.R "Base rings of IdealGens and MonomialOrdering are inconsistent"
  H = IdealGens(base_ring(G), gens(G))
  if isdefined(G, :ord) && G.ord != monord
      H.isGB = false
      H.isReduced = false
  end
  H.keep_ordering = H.keep_ordering
  H.ord = monord
  return H
end


function Base.:(==)(G1::IdealGens, G2::IdealGens)
  @assert !G1.isGB || isdefined(G1, :ord)  "Gröbner basis must have an ordering"
  @assert !G2.isGB || isdefined(G2, :ord)  "Gröbner basis must have an ordering"
  if isdefined(G1, :ord) != isdefined(G2, :ord)
      return false
  end
  if G1.isGB != G2.isGB
      return false
  end
  if isdefined(G1, :ord) && G1.ord != G2.ord
      return false
  end
  return G1.gensBiPolyArray == G2.gensBiPolyArray
end

function Base.hash(G::IdealGens, h::UInt)
  h = hash(isdefined(G, :ord), h)
  h = hash(G.isGB, h)
  if isdefined(G, :ord)
    h = hash(G.ord, h)
  end
  h = hash(G.gensBiPolyArray, h)
  return h
end


##############################################################################
#
# Conversion to and from Singular: in particular, some Rings are
# special as they exist natively in Singular and thus should be used
#
##############################################################################
#
# Needs convert(Target(Ring), elem)
# Ring(s.th.)
#
# Singular's polynomial rings are not recursive:
# 1. singular_poly_ring(R::Ring) tries to create a Singular.PolyRing (with
#    elements of type Singular.spoly) isomorphic to R
# 2. singular_coeff_ring(R::Ring) tries to create a ring isomorphic to R that is
#    acceptable to Singular.jl as 'coefficients'
#
# a real native Singular polynomial ring with Singular's native QQ as coefficients:
#  singular_poly_ring(QQ[t]) => Singular.PolyRing{Singular.n_Q}
#
# Singular's native Fp(5):
#  singular_coeff_ring(GF(5)) => Singular.N_ZpField
#
# Singular wrapper of the Oscar type QQPolyRingElem:
#  singular_coeff_ring(QQ[t]) => Singular.N_Ring{QQPolyRingElem}
#
# even more wrappings of the immutable Oscar type FpFieldElem:
#  singular_coeff_ring(GF(ZZRingElem(5))) => Singular.N_Field{Singular.FieldElemWrapper{FpField, FpFieldElem}}

for T in [:MPolyRing, :(AbstractAlgebra.Generic.MPolyRing)]
@eval function (Ox::$T)(f::Singular.spoly)
  O = base_ring(Ox)
  Sx = parent(f)
  @assert ngens(Sx) == ngens(Ox)
  g = MPolyBuildCtx(Ox)
  for (c, e) = Base.Iterators.zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    push_term!(g, O(c), e)
  end
  return finish(g)
end
end

#Note: Singular crashes if it gets Nemo.ZZ instead of Singular.ZZ ((Coeffs(17)) instead of (ZZ))
singular_coeff_ring(::ZZRing) = Singular.Integers()
singular_coeff_ring(::QQField) = Singular.Rationals()

# if the characteristic overflows an Int, Singular doesn't support it anyways
singular_coeff_ring(F::fpField) = Singular.Fp(Int(characteristic(F)))

function singular_coeff_ring(F::Union{zzModRing, ZZModRing})
  return Singular.residue_ring(Singular.Integers(), BigInt(modulus(F)))[1]
end

singular_poly_ring(R::Singular.PolyRing; keep_ordering::Bool = true) = R

# Note: Several Singular functions crash if they get the catch-all
# Singular.CoefficientRing(F) instead of the native Singular equivalent as
# conversions to/from factory are not implemented.
singular_coeff_ring(R::Union{AbsSimpleNumField, fqPolyRepField, FqField}) = codomain(iso_oscar_singular_coeff_ring(R))

function (K::FqField)(a::Singular.n_algExt)
  SK = parent(a)
  SF = parent(Singular.modulus(SK))
  SFa = SF(a)
  numSa = Singular.n_transExt_to_spoly(numerator(SFa))
  denSa = first(AbstractAlgebra.coefficients(Singular.n_transExt_to_spoly(denominator(SFa))))
  @assert isone(denSa)
  res = zero(K)
  Ka = gen(K)
  for (c, e) in zip(AbstractAlgebra.coefficients(numSa), AbstractAlgebra.exponent_vectors(numSa))
    res += K(Int(c))*Ka^e[1]
  end
  return res
end

function (SF::Singular.N_AlgExtField)(a::FqFieldElem)
  F = parent(a)
  Sa = gen(SF)
  res = SF(lift(ZZ, coeff(a, 0)))
  var = one(SF)
  for i in 1:degree(F)-1
    var = mul!(var, Sa)
    res = addmul!(res, SF(lift(ZZ, coeff(a, i))), var)
  end
  return res
end

function (F::FqField)(a::Singular.n_Zp)
  return F(Int(a))
end

function (SF::Singular.N_ZpField)(a::FqFieldElem)
   return SF(lift(ZZ, a))
end

#### TODO stuff to move to singular.jl
function (F::Singular.N_FField)(a::Union{zzModRingElem, fpFieldElem})
  return F(a.data)
end

function (K::fqPolyRepField)(a::Singular.n_algExt)
  SK = parent(a)
  SF = parent(Singular.modulus(SK))
  SFa = SF(a)
  numSa = Singular.n_transExt_to_spoly(numerator(SFa))
  denSa = first(AbstractAlgebra.coefficients(Singular.n_transExt_to_spoly(denominator(SFa))))
  @assert isone(denSa)
  res = zero(K)
  Ka = gen(K)
  for (c, e) in zip(AbstractAlgebra.coefficients(numSa), AbstractAlgebra.exponent_vectors(numSa))
    res += K(Int(c))*Ka^e[1]
  end
  return res
end

function (SF::Singular.N_AlgExtField)(a::fqPolyRepFieldElem)
  F = parent(a)
  Sa = gen(SF)
  res = SF(coeff(a, 0))
  var = one(SF)
  for i in 1:degree(F)-1
    var = mul!(var, Sa)
    res = addmul!(res, SF(coeff(a, i)), var)
  end
  return res
end
#### end stuff to move to singular.jl

function singular_poly_ring(Rx::MPolyRing{T}; keep_ordering::Bool = false) where {T <: RingElem}
  return _create_singular_poly_ring(singular_coeff_ring(base_ring(Rx)), Rx; keep_ordering)
end

function singular_poly_ring(Rx::MPolyRing{T}, ord::Union{Symbol, Singular.sordering, MonomialOrdering}) where {T <: RingElem}
  return _create_singular_poly_ring(singular_coeff_ring(base_ring(Rx)), Rx, ord)
end

#catch all for generic nemo rings
function Oscar.singular_coeff_ring(F::AbstractAlgebra.Ring)
  return Singular.CoefficientRing(F)
end

#??? needs to coerce into b? assert parent?
function (b::AbstractAlgebra.Ring)(a::Singular.n_unknown)
  Singular.libSingular.julia(Singular.libSingular.cast_number_to_void(a.ptr))::elem_type(b)
end

##############################################################################
#
# Multivariate ideals - also used for the decorated stuff
#
##############################################################################

@doc raw"""
    mutable struct MPolyIdeal{S} <: Ideal{S}

Ideal in a multivariate polynomial ring R with elements of type `S`.

Fields:
  * `gens::IdealGens{S}`, a bi-list of generators of the ideal. This is not supposed to be altered ever after assignment of the ideal;
  * `gb::IdealGens{S}`, a field used for caching of Groebner basis computations;
  * `dim::Int`, a field used for caching the dimension of the ideal.
"""
@attributes mutable struct MPolyIdeal{S} <: Ideal{S}
  gens::IdealGens{S}
  gb::Dict{MonomialOrdering, IdealGens{S}}
  dim::Union{Int, NegInf, Nothing}

  function MPolyIdeal(R::Ring, g::Vector{T}) where {T <: MPolyRingElem}
    return MPolyIdeal(IdealGens(R, g, keep_ordering = false))
  end

  function MPolyIdeal(Ox::T, s::Singular.sideal) where {T <: MPolyRing}
    r = MPolyIdeal(IdealGens(Ox, s))
    #=if s.isGB
      ord = monomial_ordering(Ox, Singular.ordering(base_ring(s)))
      r.ord = ord
      r.isGB = true
      r.gb[ord] = r.gens
    end=#
    return r
  end

  function MPolyIdeal(B::IdealGens{T}) where T
    if length(B) >= 1
      R = base_ring(B)
      if is_graded(R)
        @req all(is_homogeneous, oscar_generators(B)) "The generators of an ideal in a graded ring must be homogeneous"
      end
    end
    r = new{T}()
    r.gens = B
    r.dim = nothing
    r.gb = Dict()
    if B.isGB
      r.gb[B.ord] = B
    end
    return r
  end
end

function ideal(g::Vector{Any})
  return ideal(typeof(g[1])[x for x = g])
end

function ideal(Rx::MPolyRing, s::Singular.sideal)
  return MPolyIdeal(Rx, s)
end

function singular_generators(I::MPolyIdeal, monorder::MonomialOrdering=default_ordering(base_ring(I)))
  return singular_generators(generating_system(I), monorder)
end

function singular_polynomial_ring(I::MPolyIdeal, monorder::MonomialOrdering=default_ordering(base_ring(I)))
  return (singular_generators(I, monorder)).base_ring
end

function singular_polynomial_ring(G::IdealGens, monorder::MonomialOrdering=G.ord)
  return (singular_generators(G, monorder)).base_ring
end

oscar_generators(I::MPolyIdeal) = oscar_generators(I.gens)

oscar_generators(IG::IdealGens) = oscar_generators(IG.gensBiPolyArray)

function oscar_generators(B::BiPolyArray)
  if !isdefined(B, :O)
    if B.Ox isa MPolyQuoRing
      R = oscar_origin_ring(B.Ox)
    else
      R = B.Ox
    end
    B.O = [R(x) for x in gens(B.S)]
  end
  return B.O
end

function map_entries(R, M::Singular.smatrix)
  s = nrows(M), ncols(M)
  S = parent(R(zero(base_ring(M))))
  return matrix(S, s[1], s[2], elem_type(S)[R(M[i,j]) for i=1:s[1] for j=1:s[2]])
end

@doc raw"""
    generating_system(I::MPolyIdeal)

Return the system of generators of `I`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, [:x,:y]);

julia> I = ideal([x*(x+1), x^2-y^2+(x-2)*y]);

julia> generating_system(I)
Ideal generating system with elements
  1: x^2 + x
  2: x^2 + x*y - y^2 - 2*y
```
"""
function generating_system(I::MPolyIdeal)
  return I.gens
end

function syzygy_module(a::Vector{MPolyRingElem})
  #only graded modules exist
  error("not implemented yet")
end

function (F::Generic.FreeModule)(s::Singular.svector)
  pv = Tuple{Int, elem_type(base_ring(F))}[]
  pos = Int[]
  values = []
  Rx = base_ring(F)
  R = base_ring(Rx)
  for (i, e, c) = s
    f = findfirst(==(i), pos)
    if f === nothing
      push!(values, MPolyBuildCtx(Rx))
      f = length(values)
      push!(pos, i)
    end
    push_term!(values[f], R(c), e)
  end
  pv = Tuple{Int, elem_type(Rx)}[(pos[i], finish(values[i])) for i=1:length(pos)]
  e = zero(F)
  for (k,v) = pv
    e += v*gen(F, k)
  end
  return e
end

@doc raw"""
    jacobian_matrix(f::MPolyRingElem)

Given a polynomial $f$, return the Jacobian matrix ``J_f=(\partial_{x_1}f,...,\partial_{x_n}f)^T`` of $f$.
"""
function jacobian_matrix(f::MPolyRingElem)
  R = parent(f)
  n = nvars(R)
  return matrix(R, n, 1, [derivative(f, i) for i=1:n])
end

@doc raw"""
    jacobian_ideal(f::MPolyRingElem)

Given a polynomial $f$, return the Jacobian ideal of $f$.
"""
function jacobian_ideal(f::MPolyRingElem)
  R = parent(f)
  n = nvars(R)
  return ideal(R, [derivative(f, i) for i=1:n])
end

@doc raw"""
    jacobian_matrix([R::MPolyRing,] g::Vector{<:MPolyRingElem})

Given an array ``g=[f_1,...,f_m]`` of polynomials over the same base ring `R`,
return the Jacobian matrix ``J=(\partial_{x_i}f_j)_{i,j}`` of ``g``.
"""
function jacobian_matrix(g::Vector{<:MPolyRingElem})
  @req length(g) > 0 "specify the common parent as first argument"
  R = parent(g[1])
  return jacobian_matrix(R, g)
end

function jacobian_matrix(R::MPolyRing, g::Vector{<:MPolyRingElem})
  n = nvars(R)
  @req all(x->parent(x) === R, g) "polynomials must be elements of R"
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

##########################################
#
# Singular library related functions
#
##########################################

##########################
#
# basic maps
#
##########################
function im_func(f::MPolyRingElem, S::MPolyRing, i::Vector{Int})
  O = base_ring(S)
  g = MPolyBuildCtx(S)
  for (c, e) = Base.Iterators.zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    f = zeros(Int, nvars(S))
    for j=1:length(e)
      if i[j] == 0
        e[j] != 0 && error("illegal map: var $(j) is used")
      else
        f[i[j]] = e[j]
      end
    end
    push_term!(g, O(c), f)
  end
  return finish(g)
end


abstract type OscarMap <: SetMap end

mutable struct MPolyHom_vars{T1, T2}  <: Map{T1, T2, Hecke.HeckeMap, MPolyHom_vars}
  header::Hecke.MapHeader
  i::Vector{Int}

  function MPolyHom_vars{T1, T2}(R::T1, S::T2, i::Vector{Int}) where {T1 <: MPolyRing, T2 <: MPolyRing}
    p = sortperm(i)
    j = Int[]
    for h = 1:length(p)
      if i[p[h]] != 0
        j = p[h:length(p)]
        break
      end
    end
    header = MapHeader{T1, T2}(R, S, x -> im_func(x, S, i), y-> im_func(y, R, j))
    return new(header, i)
  end

  function MPolyHom_vars{T1, T2}(R::T1, S::T2; type::Symbol = :none) where {T1 <: MPolyRing, T2 <: MPolyRing}

    if type == :names
      i = Int[]
      for h = symbols(R)
        push!(i, findfirst(==(h), symbols(S)))
      end
      return MPolyHom_vars{T1, T2}(R, S, i)
    end
    error("type not supported")
  end
end

(f::MPolyHom_vars)(g::MPolyRingElem) = image(f, g)

function Hecke.hom(R::MPolyRing, S::MPolyRing, i::Vector{Int})
  return MPolyHom_vars{typeof(R), typeof(S)}(R, S, i)
end

function _lift(S::Singular.sideal, T::Singular.sideal)
  R = base_ring(S)
  @assert base_ring(T) == R
  c, r = Singular.libSingular.id_Lift(S.ptr, T.ptr, R.ptr)
  M = Singular.Module(R, c)

  if Singular.ngens(M) == 0 || iszero(M[1])
    error("elem not in module")
  end
  return M
end

#TODO: return a matrix??
@doc raw"""
    coordinates(a::Vector{<:MPolyRingElem}, b::Vector{<:MPolyRingElem})

Tries to write the entries of `b` as linear combinations of `a`.
"""
function coordinates(a::Vector{<:MPolyRingElem}, b::Vector{<:MPolyRingElem})
  ia = ideal(a)
  ib = ideal(b)
  c = _lift(singular_generators(ia), singular_generators(ib))
  F = free_module(parent(a[1]), length(a))
  return [F(c[x]) for x = 1:Singular.ngens(c)]
end

function coordinates(a::Vector{<:MPolyRingElem}, b::MPolyRingElem)
  return coordinates(a, [b])[1]
end

# used in expressify for mpolys and modules wrt to ordering
function _push_monomial_expr!(prod::Expr, x::AbstractVector, e::AbstractVector)
  for i in 1:length(e)
    ei = e[i]
    if iszero(ei)
    elseif isone(ei)
      push!(prod.args, x[i])
    else
      push!(prod.args, Expr(:call, :^, x[i], ei))
    end
  end
end

# expressify wrt to an arbitrary permutation
function expressify(a::OscarPair{<:MPolyRingElem, Vector{Int}}; context = nothing)
  f = a.first
  x = symbols(parent(f))
  s = Expr(:call, :+)
  for j in a.second
    prod = Expr(:call, :*)
    c = coeff(f, j)
    if !isone(c)
      push!(prod.args, expressify(c, context = context))
    end
    _push_monomial_expr!(prod, x, exponent_vector(f, j))
    push!(s.args, prod)
  end
  return s
end
@enable_all_show_via_expressify OscarPair{<:MPolyRingElem, Vector{Int}}

# expressify wrt to an ordering
function expressify(a::OscarPair{<:MPolyRingElem, <:MonomialOrdering}; context = nothing)
  perm = permutation_of_terms(a.first, a.second)
  return expressify(OscarPair(a.first, perm); context = context)
end
@enable_all_show_via_expressify OscarPair{<:MPolyRingElem, <:MonomialOrdering}

@doc raw"""
    GeneralPermutedIterator{S, T, P}

Iterator for the `S` of `elem` in the order `perm`.
"""
struct GeneralPermutedIterator{S, T, P}
  elem::T
  perm::P
end

function GeneralPermutedIterator{S}(e::T, p::P) where {S, T, P}
  return GeneralPermutedIterator{S,T,P}(e, p)
end

function Base.length(a::GeneralPermutedIterator{S, T, <:AbstractVector}) where {S, T}
  return length(a.perm)
end

function expressify(a::GeneralPermutedIterator{S}; context = nothing) where S
  msg = replace(String(S), "_" => " ") * " iterator of "
  return Expr(:sequence, Expr(:text, msg),
                     expressify(OscarPair(a.elem, a.perm); context = context))
end

@enable_all_show_via_expressify GeneralPermutedIterator


@doc raw"""
    coefficients(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return an iterator for the coefficients of `f` with respect to the order `ordering`.
"""
function coefficients(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:coefficients}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:coefficients, <:MPolyRingElem}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  return coeff(a.elem, a.perm[state]), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:coefficients, <:MPolyRingElem{C}}}) where C
  return C
end

@doc raw"""
    coefficients_and_exponents(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return an iterator whose elements are tuples of coefficients of `f` and exponent
vectors (as `Vector{Int}`) with respect to the order `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> coefficients_and_exponents((1 + x + 2*y)^2; ordering = neglex(R))
coefficients and exponents iterator of 1 + 4*y + 4*y^2 + 2*x + 4*x*y + x^2

julia> collect(ans)
6-element Vector{Tuple{QQFieldElem, Vector{Int64}}}:
 (1, [0, 0])
 (4, [0, 1])
 (4, [0, 2])
 (2, [1, 0])
 (4, [1, 1])
 (1, [2, 0])
```
"""
function coefficients_and_exponents(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:coefficients_and_exponents}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:coefficients_and_exponents, <:MPolyRingElem}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  return (coeff(a.elem, a.perm[state]), exponent_vector(a.elem, a.perm[state])), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:coefficients_and_exponents, <:MPolyRingElem{C}}}) where C
  return Tuple{C, Vector{Int}}
end

@doc raw"""
    exponents(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return an iterator for the exponent vectors (as `Vector{Int}`) of `f` with
respect to the order `ordering`.
"""
function exponents(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:exponents}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:exponents, <:MPolyRingElem}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  return exponent_vector(a.elem, a.perm[state]), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:exponents, <:MPolyRingElem}})
  return Vector{Int}
end

@doc raw"""
    monomials(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return an iterator for the monomials of `f` with respect to the order `ordering`.
"""
function monomials(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:monomials}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:monomials, <:MPolyRingElem}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  return monomial(a.elem, a.perm[state]), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:monomials, T}}) where T <: MPolyRingElem
  return T
end

@doc raw"""
    terms(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return an iterator for the terms of `f` with respect to the order `ordering`.
"""
function terms(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:terms}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:terms, <:MPolyRingElem}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  return term(a.elem, a.perm[state]), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:terms, T}}) where T <: MPolyRingElem
  return T
end

@doc raw"""
    leading_coefficient(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return the leading coefficient of `f` with respect to the order `ordering`.
"""
function leading_coefficient(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return coeff(f, index_of_leading_term(f, ordering))
end

@doc raw"""
    leading_exponent(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return the leading exponent vector (as `Vector{Int}`) of `f` with
respect to the order `ordering`.
"""
function leading_exponent(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return AbstractAlgebra.exponent_vector(f, index_of_leading_term(f, ordering))
end

@doc raw"""
    leading_coefficient_and_exponent(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return the leading coefficient of `f` with respect to the order `ordering`.
"""
function leading_coefficient_and_exponent(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  i = index_of_leading_term(f, ordering)
  return (coeff(f, i), AbstractAlgebra.exponent_vector(f, i))
end

@doc raw"""
    leading_monomial(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return the leading monomial of `f` with respect to the order `ordering`.
"""
function leading_monomial(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return monomial(f, index_of_leading_term(f, ordering))
end

@doc raw"""
    leading_term(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return the leading term of `f` with respect to the order `ordering`.
"""
function leading_term(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  return term(f, index_of_leading_term(f, ordering))
end

function _delete_index(f::MPolyRingElem, i::Int)
  z = MPolyBuildCtx(parent(f))
  for (c, e) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    i -= 1
    if i != 0
       push_term!(z, c, e)
    end
  end
  return finish(z)
end

@doc raw"""
    tail(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))

Return the tail of `f` with respect to the order `ordering`.
"""
function tail(f::MPolyRingElem; ordering::MonomialOrdering = default_ordering(parent(f)))
  # f - leading_term(f) is too easy and might have problems with inexact
  i = index_of_leading_term(f, ordering)
  return _delete_index(f, i)
end

##############################################################################
#
##############################################################################

################################################################################

@doc raw"""
    divrem(a::Vector{T}, b::Vector{T}) where T <: MPolyRingElem{S} where S <: RingElem

Return an array of tuples (qi, ri) consisting of an array of polynomials qi, one
for each polynomial in b, and a polynomial ri such that
a[i] = sum_i b[i]*qi + ri.
"""
function divrem(a::Vector{T}, b::Vector{T}) where T <: MPolyRingElem{S} where S <: RingElem
  return [divrem(x, b) for x in a]
end

################################################################################

function Base.:*(f::MPolyRingElem, I::MPolyIdeal)
  R = base_ring(I)
  R == parent(f) || error("base rings do not match")
  return ideal(R, f.*gens(I))
end

*(I::MPolyIdeal, f::MPolyRingElem) = f*I

################################################################################

function _is_integral_domain(R::MPolyRing)
  return true
end

@doc raw"""
    weighted_degree(f::MPolyRingElem, w::Vector{Int})

Given a multivariate polynomial `f` and a weight vector `w`
return the total degree of `f` with respect to the weights `w`.
"""
function weighted_degree(f::MPolyRingElem, w::Vector{Int})
  x = gens(parent(f))
  n = length(x)
  n == length(w) || error("weight vector does not have the correct length")
  vals = [sum([degree(m, j)*w[j] for j in 1:n]) for m in monomials(f)]
  return maximum(vals)
end

# assure compatibility with generic code for MPolyQuos:
lift(f::MPolyRingElem) = f

### Hessian matrix
function hessian_matrix(f::MPolyRingElem)
  R = parent(f)
  n = nvars(R)
  df = jacobian_matrix(f)
  result = zero_matrix(R, n, n)
  for i in 1:n
    for j in i:n
      result[i, j] = result[j, i] = derivative(df[i], j)
    end
  end
  return result
end

hessian(f::MPolyRingElem) = det(hessian_matrix(f))
