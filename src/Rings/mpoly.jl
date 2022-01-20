#module MPolyModule

export PolynomialRing, total_degree, degree,  MPolyIdeal, MPolyElem, ideal, coordinates,
       jacobi_matrix, jacobi_ideal,  normalize, divrem, isprimary, isprime

##############################################################################
#
# could/ should be in AbstractAlgebra
#
# some sugar to make creation of strange rings easier
# possibly lacks the passing of the ordering...
##############################################################################

#TODO: reduce = divrem in Nemo. Should be faster - if we have the correct basis

#allows
# PolynomialRing(QQ, :a=>1:3, "b"=>1:3, "c=>1:5:10)
# -> QQx, [a1, a2, a3], [b1 ,b2, b3], ....

function PolynomialRing(R::AbstractAlgebra.Ring, v::Vararg{<:Pair{<:Union{String, Symbol}, <:Any}, N}; cached::Bool = false, ordering::Symbol = :lex) where {M, N, S}
  str = _make_strings(v)
  strings = vcat(str...)
  Rx, c = PolynomialRing(R, strings, cached = cached, ordering = ordering)
  # Now we need to collect the variables
  # We do it recursively to make it type stable
  Rx, _collect_variables(c, v)...
end

# To print [1, 2, 3] or (1, 2, 3) as "1, 2, 3"
function _print_comma_list(i)
  s = IOBuffer()
  print(s, i[1])
  for j in 2:length(i)
    print(s, ", ", i[j])
  end
  return String(take!(s))
end

# To turn "x", 'x' or :x, (1, 2, 3) into x[1, 2, 3]

_make_variable(a, i) = _make_variable(String(a), i)

function _make_variable(a::String, i)
  ii = _print_comma_list(i)
  if occursin('#', a)
    aa = replace(a, '#' => "$ii")
  else
    if Hecke.inNotebook()
      aa = "$(a)_{$ii}"
    else
      aa = "$a[$ii]"
    end
  end
  return aa
end

# Type stable recursive function to create strings from "a" => 1:2 or
# "a" => (1:3, 1:3)
function _make_strings(v::Pair{<:Union{String, Symbol}, <: Any})
  lv = last(v)
  if lv isa Tuple
    p = Iterators.product(lv...)
  else
    p = lv
  end
  res = String[]
  a = first(v)
  for i in p
    push!(res, _make_variable(a, i))
  end
  return res
end

function _make_strings(v)
  s = _make_strings(first(v))
  if length(v) == 1
    return (s, )
  end
  return tuple(s, _make_strings(Base.tail(v))...)
end

# Type stable recursive function that given a vector of
# variables (or any polynomials) and v = "a" => 1:2 or "a" =>
# (1:2, 1:3) extracts the first variables into an
# n-dimensional array with the given dimensions.
# For example, _collect_variables([x1, x2, x3, x4, x5], "a" => (1:2, 1:2))
# returns [x1 x3; x2 x4], 5
function _collect_variables(c::Vector, v::Pair, start = 1)
  lv = last(v)
  if lv isa Tuple
    res = Array{eltype(c)}(undef, map(length, lv))
  else
    res = Array{eltype(c)}(undef, length(lv))
  end
  for i in eachindex(res)
    res[i] = c[start]
    start += 1
  end
  return res, start
end

function _collect_variables(c, v, start = 1)
  s, next = _collect_variables(c, first(v), start)
  if length(v) == 1
    return (s, )
  end
  return tuple(s, _collect_variables(c, Base.tail(v), next)...)
end

function Base.getindex(R::MPolyRing, i::Int)
  i == 0 && return zero(R)
  return gen(R, i)
end

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

module Orderings

using Oscar, Markdown
import Oscar: Ring, MPolyRing, MPolyElem, weights, IntegerUnion
import Base.==
export anti_diagonal, lex, degrevlex, deglex, weights, MonomialOrdering, singular

abstract type AbsOrdering end
"""
Ring-free monomial ordering: just the indices of the variables are given.
`T` can be a `UnitRange` to make Singular happy or any `Array` if the
  variables are not consequtive
"""
mutable struct GenOrdering{T} <: AbsOrdering
  vars::T
  ord::Symbol
  wgt::fmpz_mat
  function GenOrdering(u::T, s::Symbol) where {T <: AbstractVector{Int}}
    r = new{typeof(u)}()
    r.vars = u
    r.ord = s
    return r
  end
  function GenOrdering(u::T, m::fmpz_mat; ord::Symbol = :weight) where {T <: AbstractVector{Int}}
    r = new{typeof(u)}()
    @assert ncols(m) == length(u)
    r.vars = u
    r.ord = ord
    r.wgt = m
    return r
  end
end

"""
The product of `a` and `b` (`vcat` of the the matrices)
"""
mutable struct ProdOrdering <: AbsOrdering
  a::AbsOrdering
  b::AbsOrdering
end

Base.:*(a::AbsOrdering, b::AbsOrdering) = ProdOrdering(a, b)

#not really user facing
function ordering(a::AbstractVector{Int}, s::Union{Symbol, fmpz_mat})
  i = minimum(a)
  I = maximum(a)
  if I-i+1 == length(a) #testif variables are consecutive or not.
    return GenOrdering(i:I, s)
  end
  return GenOrdering(collect(a), s)
end

#not really user facing
function ordering(a::AbstractVector{Int}, s::Symbol, w::fmpz_mat)
  i = minimum(a)
  I = maximum(a)
  if I-i+1 == length(a)
    return GenOrdering(i:I, w, ord = s)
  end
  return GenOrdering(collect(a), w, ord = s)
end


#not really user facing, flattens a product of product orderings into an array 
function flat(a::GenOrdering)
  return [a]
end
function flat(a::ProdOrdering)
  return vcat(flat(a.a), flat(a.b))
end

@doc Markdown.doc"""
    anti_diagonal(R::Ring, n::Int)

A square matrix with `1` on the anti-diagonal.
"""
function anti_diagonal(R::Ring, n::Int)
  a = zero_matrix(R, n, n)
  for i=1:n
    a[i, n-i+1] = one(R)
  end
  return a
end

#not user facing
function weights(a::GenOrdering)
  if a.ord == :lex || a.ord == Symbol("Singular(lp)")
    return identity_matrix(ZZ, length(a.vars))
  end
  if a.ord == :deglex
    return [matrix(ZZ, 1, length(a.vars), ones(fmpz, length(a.vars)));
            identity_matrix(ZZ, length(a.vars)-1) zero_matrix(ZZ, length(a.vars)-1, 1)]
  end
  if a.ord == :degrevlex || a.ord == Symbol("Singular(dp)")
    return [matrix(ZZ, 1, length(a.vars), ones(fmpz, length(a.vars))) ;
            zero_matrix(ZZ, length(a.vars)-1, 1) anti_diagonal(ZZ, length(a.vars)-1)]
  end              
  if a.ord == Symbol("Singular(ls)")
    return -identity_matrix(ZZ, length(a.vars))
  end
  if a.ord == Symbol("Singular(ds)")
    return [-matrix(ZZ, 1, length(a.vars), ones(fmpz, length(a.vars))) ;
            zero_matrix(ZZ, length(a.vars)-1, 1) anti_diagonal(ZZ, length(a.vars)-1)]
  end              
  if a.ord == Symbol("Singular(a)") || a.ord == Symbol("Singular(M)")
    return a.wgt
  end              
end

#not user facing
function weights(a::AbsOrdering)
  aa = flat(a)
  m = matrix(ZZ, 0, 0, [])
  for o = aa
    w = weights(o)
    if maximum(o.vars) > ncols(m)
      m = hcat(m, zero_matrix(ZZ, nrows(m), maximum(o.vars) - ncols(m)))
    end
    mm = zero_matrix(ZZ, nrows(w), ncols(m))
    for r = 1:nrows(w)
      for c = 1:length(o.vars)
        mm[r, o.vars[c]] = w[r, c]
      end
    end
    m = vcat(m, mm)
  end
  return m
end

"""
Orderings actually applied to polynomial rings (as opposed to variable indices)
"""
mutable struct MonomialOrdering{S}
  R::S
  o::AbsOrdering
end

#not really user facing, not exported
@doc Markdown.doc"""
    ordering(a::Vector{MPolyElem}, s::Symbol)
    ordering(a::Vector{MPolyElem}, m::fmpz_mat)
    ordering(a::Vector{MPolyElem}, s::Symbol, m::fmpz_mat)

Defines an ordering to be applied to the variables in `a`.
In the first form the symbol `s` has to be one of `:lex`, `:deglex` or `:degrevlex`.
In the second form, a weight ordering using the given matrix is used.
In the last version, the symbol if of the form `Singular(..)`.
"""
function ordering(a::AbstractVector{<:MPolyElem}, s...)
  R = parent(first(a))
  g = gens(R)
  aa = [findfirst(x -> x == y, g) for y = a]
  if nothing in aa
    error("only variables allowed")
  end
  return ordering(aa, s...)
end

@doc Markdown.doc"""
    :*(M::MonomialOrdering, N::MonomialOrdering)

For orderings on the same ring, the product ordering obained by concatenation
of the weight matrics.
"""
function Base.:*(M::MonomialOrdering, N::MonomialOrdering)
  M.R == N.R || error("wrong rings")
  return MonomialOrdering(M.R, M.o*N.o)
end

function Base.show(io::IO, M::MonomialOrdering)
  a = flat(M.o)
  if length(a) > 1
    print(io, "Product ordering: ")
    for i=1:length(a)-1
      show(io, M.R, a[i])
      print(io, " \\times ")
    end
  end
  show(io, M.R, a[end])
end

function Base.show(io::IO, R::MPolyRing, o::GenOrdering)
  if o.ord == :weight
    print(io, "weight($(gens(R)[o.vars]) via $(o.wgt))")
  else
    print(io, "$(String(o.ord))($(gens(R)[o.vars]))")
  end
end

@doc Markdown.doc"""
    lex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `lex` (lexicographic) ordering on the variables given.
"""
function lex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :lex))
end
@doc Markdown.doc"""
    deglex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `deglex` ordering on the variables given.
"""
function deglex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :deglex))
end
@doc Markdown.doc"""
    degrevlex(v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines the `degreveex` ordering on the variables given.
"""
function degrevlex(v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, :degrevlex))
end

@doc Markdown.doc"""
    singular(ord::Symbol, v::AbstractVector{<:MPolyElem}) -> MonomialOrdering

Defines an ordering given in terms of Singular primitives on the variables given.
`ord` can be one of `:lp`, `:ls`, `:dp`, `:ds`.
"""
function singular(ord::Symbol, v::AbstractVector{<:MPolyElem})
  return MonomialOrdering(parent(first(v)), ordering(v, Symbol("Singular($(string(ord)))")))
end

@doc Markdown.doc"""
    singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractMatrix{<:IntegerUnion}) -> MonomialOrdering

Defines an ordering given in terms of Singular weight ordering (`M`) with the
matrix given. `ord` has to be `:M` here.
"""
function singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractMatrix{<:IntegerUnion})
  @assert ord == :M
  W = matrix(ZZ, size(w, 1), size(w, 2), w)
  return MonomialOrdering(parent(first(v)), ordering(v, Symbol("Singular($(string(ord)))"), W))
end

@doc Markdown.doc"""
    singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractVector{<:IntegerUnion}) -> MonomialOrdering

Defines an ordering given in terms of Singular weight ordering (`a`) with the
weights given. `ord` has to be `:a` here. The weights will be supplemented by
`0`.
"""
function singular(ord::Symbol, v::AbstractVector{<:MPolyElem}, w::AbstractVector{<:IntegerUnion})
  @assert ord == :a
  W = map(fmpz, w)
  while length(v) > length(W)
    push!(W, 0)
  end

  return MonomialOrdering(parent(first(v)), ordering(v, Symbol("Singular($(string(ord)))"), matrix(ZZ, 1, length(W), W)))

end

@doc Markdown.doc"""
    weights(M::MonomialOrdering)
 
Compute a corresponding weight matrix for the given ordering.
"""
function weights(M::MonomialOrdering)
  return weights(M.o)
end

@doc Markdown.doc"""
    simplify(M::MonomialOrdering) -> MonomialOrdering

Compute a weight ordering with a unique weight matrix.    
"""
function Hecke.simplify(M::MonomialOrdering)
  ww = simplify_weight_matrix(M.o)
  return MonomialOrdering(M.R, ordering(1:ncols(ww), ww))
end

function simplify_weight_matrix(M::AbsOrdering)
  w = weights(M)
  ww = matrix(ZZ, 0, ncols(w), [])
  for i=1:nrows(w)
    if iszero_row(w, i)
      continue
    end
    nw = w[i, :]
    c = content(nw)
    if c != 1
      nw = divexact(nw, c)
    end
    for j=1:nrows(ww)
      h = findfirst(x->ww[j, x] != 0, 1:ncols(w))
      if nw[1, h] != 0
        nw = abs(ww[j, h])*nw - sign(ww[j, h])*nw[1, h]*ww[j, :]
      end
    end
    if !iszero(nw)
      c = content(nw)
      if !isone(c)
        nw = divexact(nw, c)
      end
      ww = vcat(ww, nw)
    end
  end
  return ww
end

function ==(M::MonomialOrdering, N::MonomialOrdering)
  return Hecke.simplify(M).o.wgt == Hecke.simplify(N).o.wgt
end

function Base.hash(M::MonomialOrdering, u::UInt)
  return hash(Hecke.simplify(M).o.wgt, u)
end

end  # module Orderings

using .Orderings
export lex, deglex, degrevlex, weights, MonomialOrdering, singular



##############################################################################
#
# workhorse: BiPolyArray
# ideals are (mostly) generated on the Nemo side, but structural computations
# are in Singular. To avoid permanent conversion, the list of generators = sideal
# is captured in BiPolyArray: for Ocsar this is Vector{MPoly}
#                                 Singular      sideal
#
#TODO/ to think
#  default in Nemo is     :lex
#             Singular is :degrevlex -> better for std
#by default, use different orders???
#make BiPolyArray use different orders for both? make the type depend on it?
#
#for std: abstraction to allow Christian to be used
#
#type for orderings, use this...
#in general: all algos here needs revision: do they benefit from gb or not?

mutable struct BiPolyArray{S}
  O::Vector{S}
  S::Singular.sideal
  Ox #Oscar Poly Ring
  Sx # Singular Poly Ring, poss. with different ordering
  isGB::Bool #if the Singular side (the sideal) will be a GB
  ord :: Orderings.AbsOrdering #for this ordering
  keep_ordering::Bool

  function BiPolyArray(a::Vector{T}; keep_ordering::Bool = true, isGB::Bool = false) where {T <: MPolyElem}
    return BiPolyArray(parent(a[1]), a; keep_ordering = keep_ordering,
                                        isGB = isGB)
  end

  function BiPolyArray(R::MPolyRing, a::Vector{T}; keep_ordering::Bool = true, isGB::Bool = false) where {T <: MPolyElem}
    r = new{T}()
    r.O = a
    r.Ox = R
    r.isGB = isGB
    r.keep_ordering = keep_ordering
    return r
  end

  function BiPolyArray(Ox::T, b::Singular.sideal) where {T <: MPolyRing}
    r = new{elem_type(T)}()
    r.S = b
    r.O = Array{elem_type(T)}(undef, Singular.ngens(b))
    r.Ox = Ox
    r.isGB = b.isGB
    r.Sx = base_ring(b)
    r.keep_ordering = true
    return r
  end
end

function Base.getindex(A::BiPolyArray, ::Val{:S}, i::Int)
  if !isdefined(A, :S)
    A.S = Singular.Ideal(A.Sx, [A.Sx(x) for x = A.O])
  end
  return A.S[i]
end

function Base.getindex(A::BiPolyArray, ::Val{:O}, i::Int)
  if !isassigned(A.O, i)
    A.O[i] = A.Ox(A.S[i])
  end
  return A.O[i]
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
  return A[Val(:O), s], s+1
end

Base.eltype(::BiPolyArray{S}) where S = S

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
# singular_ring(Nemo-Ring) tries to create the appropriate Ring
#

for T in [:MPolyRing, :(AbstractAlgebra.Generic.MPolyRing)]
@eval function (Ox::$T)(f::Singular.spoly)
  O = base_ring(Ox)
  Sx = parent(f)
  @assert ngens(Sx) == ngens(Ox)
  g = MPolyBuildCtx(Ox)
  for (c, e) = Base.Iterators.zip(Singular.coefficients(f), Singular.exponent_vectors(f))
    push_term!(g, O(c), e)
  end
  return finish(g)
end
end


function (S::Singular.Rationals)(a::fmpq)
  b = Base.Rational{BigInt}(a)
  return S(b)
end

(F::Singular.N_ZpField)(a::Nemo.gfp_elem) = F(lift(a))
(F::Singular.N_ZpField)(a::Nemo.nmod) = F(lift(a))
(F::Nemo.GaloisField)(a::Singular.n_Zp) = F(Int(a))
(F::Nemo.NmodRing)(a::Singular.n_Zp) = F(Int(a))

#Note: Singular crashes if it gets Nemo.ZZ instead of Singular.ZZ ((Coeffs(17)) instead of (ZZ))
singular_ring(::Nemo.FlintIntegerRing) = Singular.Integers()
singular_ring(::Nemo.FlintRationalField) = Singular.Rationals()

# if the characteristic overflows an Int, Singular doesn't support it anyways
singular_ring(F::Nemo.GaloisField) = Singular.Fp(Int(characteristic(F)))

function singular_ring(F::Union{Nemo.NmodRing, Nemo.FmpzModRing})
  return Singular.ResidueRing(Singular.Integers(), BigInt(modulus(F)))
end

singular_ring(R::Singular.PolyRing; keep_ordering::Bool = true) = R

# Note: Several Singular functions crash if they get the catch-all
# Singular.CoefficientRing(F) instead of the native Singular equivalent as
# conversions to/from factory are not implemented.
function singular_ring(K::AnticNumberField)
  minpoly = defining_polynomial(K)
  Qa = parent(minpoly)
  a = gen(Qa)
  SQa, (Sa,) = Singular.FunctionField(Singular.QQ, map(String, symbols(Qa)))
  Sminpoly = SQa(coeff(minpoly, 0))
  for i in 1:degree(minpoly)
    Sminpoly += SQa(coeff(minpoly, i))*Sa^i
  end
  SK, _ = Singular.AlgebraicExtensionField(SQa, Sminpoly)
  return SK
end

function singular_ring(F::FqNmodFiniteField)
  # TODO: the Fp(Int(char)) can throw
  minpoly = modulus(F)
  Fa = parent(minpoly)
  SFa, (Sa,) = Singular.FunctionField(Singular.Fp(Int(characteristic(F))),
                                                    map(String, symbols(Fa)))
  Sminpoly = SFa(coeff(minpoly, 0))
  for i in 1:degree(minpoly)
    Sminpoly += SFa(coeff(minpoly, i))*Sa^i
  end
  SF, _ = Singular.AlgebraicExtensionField(SFa, Sminpoly)
  return SF
end

#### TODO stuff to move to singular.jl
function (F::Singular.N_FField)(a::Union{nmod, gfp_elem})
  return F(a.data)
end

function (K::FqNmodFiniteField)(a::Singular.n_algExt)
  SK = parent(a)
  SF = parent(Singular.modulus(SK))
  SFa = SF(a)
  numSa = Singular.n_transExt_to_spoly(numerator(SFa))
  denSa = first(coefficients(Singular.n_transExt_to_spoly(denominator(SFa))))
  @assert isone(denSa)
  res = zero(K)
  Ka = gen(K)
  for (c, e) in zip(coefficients(numSa), exponent_vectors(numSa))
    res += K(Int(c))*Ka^e[1]
  end
  return res
end

function (SF::Singular.N_AlgExtField)(a::fq_nmod)
  F = parent(a)
  SFa = gen(SF)
  res = SF(coeff(a, 0))
  for i in 1:degree(F)-1
    res += SF(coeff(a, i))*SFa^i
  end
  return res
end
#### end stuff to move to singular.jl

function singular_ring(Rx::MPolyRing{T}; keep_ordering::Bool = false) where {T <: RingElem}
  if keep_ordering
    return Singular.PolynomialRing(singular_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ordering(Rx),
              cached = false)[1]
  else
    return Singular.PolynomialRing(singular_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              cached = false)[1]
  end
end

function singular_ring(Rx::MPolyRing{T}, ord::Symbol) where {T <: RingElem}
  return Singular.PolynomialRing(singular_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ord,
              cached = false)[1]
end

function singular(o::Orderings.GenOrdering)
  v = o.vars
  @assert minimum(v)+length(v) == maximum(v)+1
  if o.ord == :lex
    return Singular.ordering_lp(length(v))
  elseif o.ord == :degrevlex
    return Singular.ordering_dp(length(v))
  elseif o.ord == :weights
    return Singular.ordering_M(o.wgt)
  else
    error("not done yet")
  end
end

function singular_ring(Rx::MPolyRing{T}, ord::Orderings.AbsOrdering) where {T <: RingElem}
  #test if it can be mapped directly to singular:
  # - consecutive, non-overlapping variables
  # - covering everything
  # if this fails, create a matrix...
  f = Orderings.flat(ord)
  st = 1
  iseasy = true
  for i = 1:length(f)
    mi = minimum(f[i].vars)
    ma = maximum(f[i].vars)
    if mi == st && length(f[i].vars) + st == ma+1
      st = ma+1
    else
      iseasy = false
      break
    end
  end

  if iseasy
    o = singular(f[1])
    for i=2:length(f)
      o = o*singular(f[i])
    end
  else
    o = Singular.ordering_M(Orderings.simplify_weight_matrix(ord))
  end

  return Singular.PolynomialRing(singular_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = o,
              cached = false)[1]
end


#catch all for generic nemo rings
function Oscar.singular_ring(F::AbstractAlgebra.Ring)
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

@Markdown.doc """
    mutable struct MPolyIdeal{S} <: Ideal{S}

Ideal in a multivariate polynomial ring R with elements of type `S`.

Fields:
  * `gens::BiPolyArray{S}`, a bi-list of generators of the ideal. This is not supposed to be altered ever after assignment of the ideal;
  * `gb::BiPolyArray{S}`, a field used for caching of Groebner basis computations;
  * `dim::Int`, a field used for caching the dimension of the ideal.
"""
mutable struct MPolyIdeal{S} <: Ideal{S}
  gens::BiPolyArray{S}
  gb::BiPolyArray{S}
  dim::Int

  function MPolyIdeal(g::Vector{T}) where {T <: MPolyElem}
    r = new{T}()
    r.dim = -1 #not known
    r.gens = BiPolyArray(g, keep_ordering = false)
    return r
  end
  function MPolyIdeal(Ox::T, s::Singular.sideal) where {T <: MPolyRing}
    r = new{elem_type(T)}()
    r.dim = -1 #not known
    r.gens = BiPolyArray(Ox, s)
    if s.isGB
      r.gb = r.gens
    end
    return r
  end
  function MPolyIdeal(B::BiPolyArray{T}) where T
    r = new{T}()
    r.dim = -1
    r.gens = B
    return r
  end
end

@enable_all_show_via_expressify MPolyIdeal

function AbstractAlgebra.expressify(a::MPolyIdeal; context = nothing)
  return Expr(:call, :ideal, [expressify(g, context = context) for g in collect(a.gens)]...)
end


function ideal(g::Vector{Any})
  return ideal(typeof(g[1])[x for x = g])
end

function ideal(Rx::MPolyRing, s::Singular.sideal)
  return MPolyIdeal(Rx, s)
end

function singular_assure(I::MPolyIdeal)
  singular_assure(I.gens)
end

function singular_assure(I::BiPolyArray)
  if !isdefined(I, :S)
    I.Sx = singular_ring(I.Ox, keep_ordering=I.keep_ordering)
    I.S = Singular.Ideal(I.Sx, elem_type(I.Sx)[I.Sx(x) for x = I.O])
    if I.isGB
      I.S.isGB = true
    end
  end
end

function oscar_assure(I::MPolyIdeal)
  if !isdefined(I.gens, :O)
    I.gens.O = [I.gens.Ox(x) for x = gens(I.gens.S)]
  end
  if isdefined(I, :gb)
    I.gb.O = [I.gb.Ox(x) for x = gens(I.gb.S)]
  end
end

function Base.copy(f::MPolyElem)
    Ox = parent(f)
    g = MPolyBuildCtx(Ox)
    for (c,e) = Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
        push_term!(g, c, e)
    end
    return finish(g)
end

function map_entries(R, M::Singular.smatrix)
  s = nrows(M), ncols(M)
  S = parent(R(zero(base_ring(M))))
  return matrix(S, s[1], s[2], elem_type(S)[R(M[i,j]) for i=1:s[1] for j=1:s[2]])
end

function syzygy_module(a::Vector{MPolyElem})
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
    f = Base.findfirst(x->x==i, pos)
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

@doc Markdown.doc"""
    jacobi_matrix(f::MPolyElem)

Given a polynomial $f$ this function returns the Jacobian matrix ``J_f=(\partial_{x_1}f,...,\partial_{x_n}f)^T`` of $f$.
"""
function jacobi_matrix(f::MPolyElem)
  R = parent(f)
  n = nvars(R)
  return matrix(R, n, 1, [derivative(f, i) for i=1:n])
end

@doc Markdown.doc"""
    jacobi_ideal(f::MPolyElem)

Given a polynomial $f$ this function returns the Jacobian ideal of $f$.
"""
function jacobi_ideal(f::MPolyElem)
  R = parent(f)
  n = nvars(R)
  return ideal(R, [derivative(f, i) for i=1:n])
end

@doc Markdown.doc"""
    jacobi_matrix(g::Vector{<:MPolyElem})

Given an array ``g=[f_1,...,f_m]`` of polynomials over the same base ring,
this function returns the Jacobian matrix ``J=(\partial_{x_i}f_j)_{i,j}`` of ``g``.
"""
function jacobi_matrix(g::Vector{<:MPolyElem})
  R = parent(g[1])
  n = nvars(R)
  @assert all(x->parent(x) == R, g)
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
function im_func(f::MPolyElem, S::MPolyRing, i::Vector{Int})
  O = base_ring(S)
  g = MPolyBuildCtx(S)
  for (c, e) = Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
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
    r = new()
    p = sortperm(i)
    j = Int[]
    for h = 1:length(p)
      if i[p[h]] != 0
        j = p[h:length(p)]
        break
      end
    end
    r.header = MapHeader{T1, T2}(R, S, x -> im_func(x, S, i), y-> im_func(y, R, j))
    r.i = i
    return r
  end

  function MPolyHom_vars{T1, T2}(R::T1, S::T2; type::Symbol = :none) where {T1 <: MPolyRing, T2 <: MPolyRing}

    if type == :names
      i = Int[]
      for h = symbols(R)
        push!(i, findfirst(x -> x == h, symbols(S)))
      end
      return MPolyHom_vars{T1, T2}(R, S, i)
    end
    error("type not supported")
  end
end

(f::MPolyHom_vars)(g::MPolyElem) = image(f, g)

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
@doc Markdown.doc"""
    coordinates(a::Vector{<:MPolyElem}, b::Vector{<:MPolyElem})

Tries to write the entries of `b` as linear combinations of `a`.
"""
function coordinates(a::Vector{<:MPolyElem}, b::Vector{<:MPolyElem})
  ia = ideal(a)
  ib = ideal(b)
  singular_assure(ia)
  singular_assure(ib)
  c = _lift(ia.gens.S, ib.gens.S)
  F = free_module(parent(a[1]), length(a))
  return [F(c[x]) for x = 1:Singular.ngens(c)]
end

function coordinates(a::Vector{<:MPolyElem}, b::MPolyElem)
  return coordinates(a, [b])[1]
end

###################################################

# Some isless functions for orderings:
# _isless_:ord(f, k, l) returns true if the k-th term is lower than the l-th
# term of f in the ordering :ord.

function _isless_lex(f::MPolyElem, k::Int, l::Int)
  n = nvars(parent(f))
  for i = 1:n
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek == el
      continue
    elseif ek > el
      return false
    else
      return true
    end
  end
  return false
end

function _isless_neglex(f::MPolyElem, k::Int, l::Int)
  n = nvars(parent(f))
  for i = 1:n
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek == el
      continue
    elseif ek < el
      return false
    else
      return true
    end
  end
  return false
end

function _isless_revlex(f::MPolyElem, k::Int, l::Int)
  n = nvars(parent(f))
  for i = n:-1:1
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek == el
      continue
    elseif ek > el
      return false
    else
      return true
    end
  end
  return false
end

function _isless_negrevlex(f::MPolyElem, k::Int, l::Int)
  n = nvars(parent(f))
  for i = n:-1:1
    ek = exponent(f, k, i)
    el = exponent(f, l, i)
    if ek == el
      continue
    elseif ek < el
      return false
    else
      return true
    end
  end
  return false
end

function _isless_deglex(f::MPolyElem, k::Int, l::Int)
  tdk = total_degree(term(f, k))
  tdl = total_degree(term(f, l))
  if tdk < tdl
    return true
  elseif tdk > tdl
    return false
  end
  return _isless_lex(f, k, l)
end

function _isless_degrevlex(f::MPolyElem, k::Int, l::Int)
  tdk = total_degree(term(f, k))
  tdl = total_degree(term(f, l))
  if tdk < tdl
    return true
  elseif tdk > tdl
    return false
  end
  return _isless_negrevlex(f, k, l)
end

function _isless_negdeglex(f::MPolyElem, k::Int, l::Int)
  tdk = total_degree(term(f, k))
  tdl = total_degree(term(f, l))
  if tdk > tdl
    return true
  elseif tdk < tdl
    return false
  end
  return _isless_lex(f, k, l)
end

function _isless_negdegrevlex(f::MPolyElem, k::Int, l::Int)
  tdk = total_degree(term(f, k))
  tdl = total_degree(term(f, l))
  if tdk > tdl
    return true
  elseif tdk < tdl
    return false
  end
  return _isless_negrevlex(f, k, l)
end

# Return the degree of the k-th term of f weighted by w,
# that is deg(x^a) = w_1a_1 + ... + w_na_n.
# No sanity checks are performed!
function weighted_degree(f::MPolyElem, k::Int, w::Vector{Int})
  ek = exponent_vector(f, k)
  return dot(ek, w)
end

function _isless_weightlex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
  dk = weighted_degree(f, k, w)
  dl = weighted_degree(f, l, w)
  if dk < dl
    return true
  elseif dk > dl
    return false
  end
  return _isless_lex(f, k, l)
end

function _isless_weightrevlex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
  dk = weighted_degree(f, k, w)
  dl = weighted_degree(f, l, w)
  if dk < dl
    return true
  elseif dk > dl
    return false
  end
  return _isless_negrevlex(f, k, l)
end

function _isless_weightneglex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
  dk = weighted_degree(f, k, w)
  dl = weighted_degree(f, l, w)
  if dk < dl
    return true
  elseif dk > dl
    return false
  end
  return _isless_lex(f, k, l)
end

function _isless_weightnegrevlex(f::MPolyElem, k::Int, l::Int, w::Vector{Int})
  dk = weighted_degree(f, k, w)
  dl = weighted_degree(f, l, w)
  if dk > dl
    return true
  elseif dk < dl
    return false
  end
  return _isless_negrevlex(f, k, l)
end

function _isless_matrix(f::MPolyElem, k::Int, l::Int, M::Union{ Matrix{T}, MatElem{T} }) where T
  ek = exponent_vector(f, k)
  el = exponent_vector(f, l)
  n = nvars(parent(f))
  for i = 1:size(M, 1)
    eki = sum( M[i, j]*ek[j] for j = 1:n )
    eli = sum( M[i, j]*el[j] for j = 1:n )
    if eki == eli
      continue
    elseif eki > eli
      return false
    else
      return true
    end
  end
  return false
end

function _perm_of_terms(f::MPolyElem, ord_lt::Function)
  p = collect(1:length(f))
  sort!(p, lt = (k, l) -> ord_lt(f, k, l), rev = true)
  return p
end

# Requiring R for consistence with the other lt_from_ordering functions
function lt_from_ordering(::MPolyRing, ord::Symbol)
  if ord == :lex || ord == :lp
    return _isless_lex
  elseif ord == :revlex || ord == :rp
    return _isless_revlex
  elseif ord == :deglex || ord == :Dp
    return _isless_deglex
  elseif ord == :degrevlex || ord == :dp
    return _isless_degrevlex
  elseif ord == :neglex || ord == :ls
    return _isless_neglex
  elseif ord == :negrevlex || ord == :rs
    return _isless_negrevlex
  elseif ord == :negdeglex || ord == :Ds
    return _isless_negdeglex
  elseif ord == :negdegrevlex || ord == :ds
    return _isless_negdegrevlex
  else
    error("Ordering $ord not available")
  end
end

function lt_from_ordering(R::MPolyRing, ord::Symbol, w::Vector{Int})
  @assert length(w) == nvars(R) "Number of weights has to match number of variables"

  if ord == :weightlex || ord == :Wp
    @assert all(x -> x > 0, w) "Weights have to be positive"
    return (f, k, l) -> _isless_weightlex(f, k, l, w)
  elseif ord == :weightrevlex || ord == :wp
    @assert all(x -> x > 0, w) "Weights have to be positive"
    return (f, k, l) -> _isless_weightrevlex(f, k, l, w)
  elseif ord == :weightneglex || ord == :Ws
    @assert !iszero(w[1]) "First weight must not be 0"
    return (f, k, l) -> _isless_weightneglex(f, k, l, w)
  elseif ord == :weightnegrevlex || ord == :ws
    @assert !iszero(w[1]) "First weight must not be 0"
    return (f, k, l) -> _isless_weightnegrevlex(f, k, l, w)
  else
    error("Ordering $ord not available")
  end
end

function lt_from_ordering(R::MPolyRing, M::Union{ Matrix{T}, MatElem{T} }) where T
  @assert size(M, 2) == nvars(R) "Matrix dimensions have to match number of variables"

  return (f, k, l) -> _isless_matrix(f, k, l, M)
end

function terms(f::MPolyElem, ord::Function)
  perm = _perm_of_terms(f, ord)
  return ( term(f, perm[i]) for i = 1:length(f) )
end

function coefficients(f::MPolyElem, ord::Function)
  perm = _perm_of_terms(f, ord)
  return ( coeff(f, perm[i]) for i = 1:length(f) )
end

function exponent_vectors(f::MPolyElem, ord::Function)
  perm = _perm_of_terms(f, ord)
  return ( exponent_vector(f, perm[i]) for i = 1:length(f) )
end

function monomials(f::MPolyElem, ord::Function)
  perm = _perm_of_terms(f, ord)
  return ( monomial(f, perm[i]) for i = 1:length(f) )
end

for s in (:terms, :coefficients, :exponent_vectors, :monomials)
  @eval begin
    function ($s)(f::MPolyElem, ord::Symbol)
      R = parent(f)
      if ord == ordering(R)
        return ($s)(f)
      end

      lt = lt_from_ordering(R, ord)
      return ($s)(f, lt)
    end

    function ($s)(f::MPolyElem, M::Union{ Matrix{T}, MatElem{T} }) where T
      R = parent(f)
      lt = lt_from_ordering(R, M)
      return ($s)(f, lt)
    end

    function ($s)(f::MPolyElem, ord::Symbol, weights::Vector{Int})
      R = parent(f)
      lt = lt_from_ordering(R, ord, weights)
      return ($s)(f, lt)
    end

    function ($s)(f::MPolyElem, ord::MonomialOrdering)
      R = parent(f)
      lt = lt_from_ordering(R, weights(ord))
      return ($s)(f, lt)
    end
  end
end

for s in ("term", "coefficient", "monomial")
  @eval begin
    function ($(Symbol("leading_$s")))(args...)
      return first($(Symbol("$(s)s"))(args...))
    end
  end
end

function leading_term(f::MPolyElem)
  return leading_term(f, ordering(parent(f)))
end

function leading_coefficient(f::MPolyElem)
  return leading_coefficient(f, ordering(parent(f)))
end

function leading_monomial(f::MPolyElem)
  return leading_monomial(f, ordering(parent(f)))
end

##############################################################################
#
##############################################################################

#=
function factor(f::MPolyElem)
  I = ideal(parent(f), [f])
  fS = Singular.factor(I.gens[Val(:S), 1])
  R = parent(f)
  return Nemo.Fac(R(fS.unit), Dict(R(k) =>v for (k,v) = fS.fac))
end
=#

# generic fallback since this is not implemented specifically anywhere yet
function isirreducible(a::MPolyElem)
  af = factor(a)
  return !(length(af.fac) > 1 || any(x->x>1, values(af.fac)))
end

################################################################################

@doc Markdown.doc"""
    divrem(a::Vector{T}, b::Vector{T}) where T <: MPolyElem{S} where S <: RingElem
Return an array of tuples (qi, ri) consisting of an array of polynomials qi, one
for each polynomial in b, and a polynomial ri such that
a[i] = sum_i b[i]*qi + ri.
"""
function divrem(a::Vector{T}, b::Vector{T}) where T <: MPolyElem{S} where S <: RingElem
  return [divrem(x, b) for x in a]
end

################################################################################

function Base.:*(f::MPolyElem, I::MPolyIdeal)
  R = base_ring(I)
  R == parent(f) || error("base rings do not match")
  return ideal(R, f.*gens(I))
end

*(I::MPolyIdeal, f::MPolyElem) = f*I


################################################################################

