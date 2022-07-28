#module MPolyModule

export PolynomialRing, total_degree, degree,  MPolyIdeal, MPolyElem, ideal, coordinates,
       jacobi_matrix, jacobi_ideal,  normalize, divrem, is_primary, is_prime

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

function PolynomialRing(R::AbstractAlgebra.Ring, v1::Pair{<:Union{String, Symbol}, <:Any}, v...; cached::Bool = false, ordering::Symbol = :lex)
  w = (v1, v...)
  str = _make_strings(w)
  strings = vcat(str...)
  Rx, c = PolynomialRing(R, strings, cached = cached, ordering = ordering)
  # Now we need to collect the variables
  # We do it recursively to make it type stable
  Rx, _collect_variables(c, w)...
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
    res = Vector{eltype(c)}(undef, length(lv))
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

ngens(F::AbstractAlgebra.Generic.FracField{T}) where {T <: MPolyElem} = ngens(base_ring(F))

function gen(F::AbstractAlgebra.Generic.FracField{T}) where {T <: PolyElem}
  return F(gen(base_ring(F)))
end

function gen(F::AbstractAlgebra.Generic.FracField{T}, i::Int) where {T <: MPolyElem}
  return F(gen(base_ring(F), i))
end

function gens(F::AbstractAlgebra.Generic.FracField{T}) where {T <: Union{PolyElem, MPolyElem}}
  return map(F, gens(base_ring(F)))
end

function Base.getindex(F::AbstractAlgebra.Generic.FracField{T}, i::Int) where {T <: MPolyElem}
  i == 0 && return zero(F)
  return gen(F, i)
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

using .Orderings
export lex, deglex, degrevlex, revlex, neglex, negrevlex, negdeglex,
       negdegrevlex, wdeglex, wdegrevlex, negwdeglex, negwdegrevlex,
       matrix_ordering, weight_matrix, MonomialOrdering, singular,
       is_global, is_local, is_mixed, monomial_ordering,
       permutation_of_terms, weighted_ordering, canonical_weight_matrix


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

default_ordering(R::MPolyRing) = degrevlex(gens(R))

mutable struct BiPolyArray{S}
  Ox::NCRing #Oscar Poly Ring or Algebra
  O::Vector{S}
  Sx # Singular Poly Ring or Algebra, poss. with different ordering
  S::Singular.sideal
  isGB::Bool #if the Singular side (the sideal) will be a GB
  ord::Orderings.AbsOrdering #for this ordering
  keep_ordering::Bool

  function BiPolyArray(O::Vector{T}; keep_ordering::Bool = true, isGB::Bool = false) where {T <: NCRingElem}
    return BiPolyArray(parent(O[1]), O; keep_ordering = keep_ordering,
                                        isGB = isGB)
  end

  function BiPolyArray(Ox::NCRing, O::Vector{T}; keep_ordering::Bool = true, isGB::Bool = false) where {T <: NCRingElem}
    r = new{T}(Ox, O)
    r.isGB = isGB
    r.keep_ordering = keep_ordering
    return r
  end

  function BiPolyArray(Ox::T, S::Singular.sideal) where {T <: NCRing}
    O = Vector{elem_type(T)}(undef, Singular.ngens(S))
    Sx = base_ring(S)
    r = new{elem_type(T)}(Ox, O, Sx, S)
    r.isGB = S.isGB
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
# Singular wrapper of the Oscar type fmpq_poly:
#  singular_coeff_ring(QQ[t]) => Singular.N_Ring{fmpq_poly}
#
# even more wrappings of the immutable Oscar type gfp_fmpz_elem:
#  singular_coeff_ring(GF(fmpz(5))) => Singular.N_Field{Singular.FieldElemWrapper{GaloisFmpzField, gfp_fmpz_elem}}

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
singular_coeff_ring(::Nemo.FlintIntegerRing) = Singular.Integers()
singular_coeff_ring(::Nemo.FlintRationalField) = Singular.Rationals()

# if the characteristic overflows an Int, Singular doesn't support it anyways
singular_coeff_ring(F::Nemo.GaloisField) = Singular.Fp(Int(characteristic(F)))

function singular_coeff_ring(F::Union{Nemo.NmodRing, Nemo.FmpzModRing})
  return Singular.ResidueRing(Singular.Integers(), BigInt(modulus(F)))
end

singular_poly_ring(R::Singular.PolyRing; keep_ordering::Bool = true) = R

# Note: Several Singular functions crash if they get the catch-all
# Singular.CoefficientRing(F) instead of the native Singular equivalent as
# conversions to/from factory are not implemented.
function singular_coeff_ring(K::AnticNumberField)
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

function singular_coeff_ring(F::FqNmodFiniteField)
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

function singular_poly_ring(Rx::MPolyRing{T}; keep_ordering::Bool = false) where {T <: RingElem}
  if keep_ordering
    return Singular.PolynomialRing(singular_coeff_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ordering(Rx),
              cached = false)[1]
  else
    return Singular.PolynomialRing(singular_coeff_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              cached = false)[1]
  end
end

function singular_poly_ring(Rx::MPolyRing{T}, ord::Symbol) where {T <: RingElem}
  return Singular.PolynomialRing(singular_coeff_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ord,
              cached = false)[1]
end

function singular_ring(Rx::MPolyRing{T}, ord::Singular.sordering) where {T <: RingElem}
  return Singular.PolynomialRing(singular_coeff_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ord,
              cached = false)[1]
end

function singular_poly_ring(Rx::MPolyRing{T}, ord::MonomialOrdering) where {T <: RingElem}
  return Singular.PolynomialRing(singular_coeff_ring(base_ring(Rx)),
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = singular(ord),
              cached = false)[1]
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
  gb::Dict{MonomialOrdering, BiPolyArray{S}}
  dim::Int

  function MPolyIdeal(g::Vector{T}) where {T <: MPolyElem}
    return MPolyIdeal(BiPolyArray(g, keep_ordering = false))
  end

  function MPolyIdeal(Ox::T, s::Singular.sideal) where {T <: MPolyRing}
    r = MPolyIdeal(BiPolyArray(Ox, s))
    if s.isGB
      ord = monomial_ordering(Ox, ordering(base_ring(s)))
      r.gb[ord] = r.gens
    end
    return r
  end

  function MPolyIdeal(B::BiPolyArray{T}) where T
    r = new{T}(B, Dict(), -1)
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
    I.Sx = singular_poly_ring(I.Ox, keep_ordering=I.keep_ordering)
    I.S = Singular.Ideal(I.Sx, elem_type(I.Sx)[I.Sx(x) for x = I.O])
  end
  if I.isGB
    I.S.isGB = true
  end
end

function singular_assure(I::MPolyIdeal, ordering::MonomialOrdering)
   singular_assure(I.gens, ordering)
end

function singular_assure(I::BiPolyArray, ordering::MonomialOrdering)
    if !isdefined(I, :S) 
        I.ord = ordering.o
        I.Sx = singular_poly_ring(I.Ox, ordering)
        I.S = Singular.Ideal(I.Sx, elem_type(I.Sx)[I.Sx(x) for x = I.O])
        if I.isGB
            I.S.isGB = true
        end
    else
        #= singular ideal exists, but the singular ring has the wrong ordering
         = attached, thus we have to create a new singular ring and map the ideal. =#
        if !isdefined(I, :ord) || I.ord != ordering.o
            I.ord = ordering.o
            SR    = singular_poly_ring(I.Ox, ordering)
            f     = Singular.AlgebraHomomorphism(I.Sx, SR, gens(SR))
            I.S   = Singular.map_ideal(f, I.S)
            I.Sx  = SR
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

function oscar_assure(B::BiPolyArray)
  if !isdefined(B, :O) || !isassigned(B.O, 1)
    B.O = [B.Ox(x) for x = gens(B.S)]
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

function terms(f::MPolyElem, ord::MonomialOrdering)
  perm = permutation_of_terms(f, ord)
  return ( term(f, perm[i]) for i = 1:length(f) )
end

function coefficients(f::MPolyElem, ord::MonomialOrdering)
  perm = permutation_of_terms(f, ord)
  return ( coeff(f, perm[i]) for i = 1:length(f) )
end

function exponent_vectors(f::MPolyElem, ord::MonomialOrdering)
  perm = permutation_of_terms(f, ord)
  return ( exponent_vector(f, perm[i]) for i = 1:length(f) )
end

function monomials(f::MPolyElem, ord::MonomialOrdering)
  perm = permutation_of_terms(f, ord)
  return ( monomial(f, perm[i]) for i = 1:length(f) )
end

for s in (:terms, :coefficients, :exponent_vectors, :monomials)
  @eval begin
    function ($s)(f::MPolyElem, ord::Symbol)
      R = parent(f)
      if ord == ordering(R)
        return ($s)(f)
      end
      return ($s)(f, monomial_ordering(gens(R), ord))
    end

    function ($s)(f::MPolyElem, M::Union{ Matrix{T}, MatElem{T} }) where T
      R = parent(f)
      return ($s)(f, matrix_ordering(gens(R), M))
    end

    function ($s)(f::MPolyElem, ord::Symbol, weights::Vector{Int})
      R = parent(f)
      return ($s)(f, monomial_ordering(gens(R), ord, weights))
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
function is_irreducible(a::MPolyElem)
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

