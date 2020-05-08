#module MPolyModule

import AbstractAlgebra: PolyRing, PolynomialRing, total_degree, degree, Ideal,
                        MPolyElem, Generic.MPolyExponentVectors, Generic.MPolyCoeffs,
                        Generic.MPolyBuildCtx, Generic.push_term!, Generic.finish, MPolyRing,
                        base_ring, ngens, gens, dim, ordering, SetMap, Map
import Nemo
import Nemo: fmpz, fmpq

import Singular

import Hecke
import Hecke: MapHeader, math_html

export PolynomialRing, total_degree, degree, MPolyElem, ordering, ideal, groebner_basis

#allows
# PolynomialRing(QQ, :a=>1:3, "b"=>1:3, "c=>1:5:10)
# -> QQx, [a1, a2, a3], [b1 ,b2, b3], ....
function PolynomialRing(R::AbstractAlgebra.Ring, v::Pair{<:Union{String, Symbol}, <:Union{StepRange{Int, Int}, UnitRange{Int}}}...; cached::Bool = false)
  s = String[]
  g = []
  j = 1
  for (a, b) = v
    h = []
    for i = b
      if occursin('#', "$a")
        aa = replace("$a", '#' => "$i")
      else
        if Hecke.inNotebook()
          aa = "$(a)_{$i}"
        else
          aa = "$a$i"
        end
      end
      push!(s, aa)
      push!(h, j)
      j += 1
    end
    push!(g, h)
  end
  Rx, c = PolynomialRing(R, s, cached = cached)
  return Rx, [c[x] for x = g]...
end

function Base.getindex(R::MPolyRing, i::Int)
  i == 0 && return zero(R)
  return gen(R, i)
end

function Base.show(io::IO, mime::MIME"text/html", R::MPolyRing)
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

function math_html(io::IO, R::MPolyElem)
  f = "$R"
  f = replace(f, "*" => "")
  f = replace(f, r"\^([0-9]*)" => s"^{\1}")
  print(io, f)
end

function Base.show(io::IO, ::MIME"text/html", R::MPolyElem)
  print(io, "\$")
  math_html(io, R)
  print(io, "\$")
end

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

mutable struct BiPolyArray
  O::Array{MPolyElem, 1} 
  S::Singular.sideal
  Ox #Oscar Poly Ring
  Sx # Singular Poly Ring, poss. with different ordering
  function BiPolyArray(a::Array{T, 1}; keep_ordering::Bool = true) where {T <: MPolyElem}
    r = new()
    r.O = a
    r.Ox = parent(a[1])
    r.Sx = singular_ring(r.Ox, keep_ordering = keep_ordering)
    return r
  end
  function BiPolyArray(Ox, b::Singular.sideal)
    r =new()
    r.S = b
    r.O = Array{MPolyElem}(undef, Singular.ngens(b))
    r.Ox = Ox
    r.Sx = base_ring(r.S)
    return r
  end
end

function Base.getindex(A::BiPolyArray, ::Val{:S}, i::Int)
  if !isdefined(A, :S)
    A.S = Singular.Ideal(A.Sx, [convert(A.Sx, x) for x = A.O])
  end
  return A.S[i]
end

function Base.getindex(A::BiPolyArray, ::Val{:O}, i::Int)
  if !isassigned(A.O, i)
    A.O[i] = convert(A.Ox, A.S[i])
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

function Base.convert(Ox::MPolyRing, f::MPolyElem) 
  O = base_ring(Ox)
  g = MPolyBuildCtx(Ox)
  for (c, e) = Base.Iterators.zip(MPolyCoeffs(f), MPolyExponentVectors(f))
    push_term!(g, O(c), e)
  end
  return finish(g)
end

function Base.convert(::Type{fmpz}, a::Singular.n_Z)
  return fmpz(BigInt(a))
end

function Base.convert(::Type{fmpq}, a::Singular.n_Q)
  return fmpq(Base.Rational{BigInt}(a))
end

function (::Nemo.FlintRationalField)(a::Singular.n_Q)
  return convert(fmpq, a)
end

function (S::Singular.Rationals)(a::fmpq)
  b = Base.Rational{BigInt}(a)
  return S(b)
end
(F::Singular.N_ZpField)(a::Nemo.gfp_elem) = F(lift(a))
(F::Nemo.GaloisField)(a::Singular.n_Zp) = F(Int(a))

singular_ring(::Nemo.FlintRationalField) = Singular.Rationals()
singular_ring(F::Nemo.GaloisField) = Singular.Fp(Int(characteristic(F)))

function singular_ring(Rx::Nemo.FmpqMPolyRing; keep_ordering::Bool = true)
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

function singular_ring(Rx::Generic.MPolyRing{T}; keep_ordering::Bool = true) where {T <: RingElem}
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

function singular_ring(Rx::Nemo.FmpqMPolyRing, ord::Symbol)
  return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ord,
              cached = false)[1]
end

function singular_ring(Rx::Generic.MPolyRing{T}, ord::Symbol) where {T <: RingElem}
  return Singular.PolynomialRing(singular_ring(base_ring(Rx)), 
              [string(x) for x = Nemo.symbols(Rx)],
              ordering = ord,
              cached = false)[1]
end

mutable struct MPolyIdeal <: Ideal{MPolyElem}
  gens::BiPolyArray
  gb::BiPolyArray
  dim::Int

  function MPolyIdeal(g::Array{T, 1}) where {T <: MPolyElem}
    r = new()
    r.dim = -1 #not known
    r.gens = BiPolyArray(g, keep_ordering = false)
    return r
  end
  function MPolyIdeal(Ox::MPolyRing, s::Singular.sideal)
    r = new()
    r.dim = -1 #not known
    r.gens = BiPolyArray(Ox, s)
    if s.isGB
      r.gb = gens
    end
    return r
  end
end

function Base.show(io::IO, I::MPolyIdeal)
  print(io, "ideal generated by: ")
  g = collect(I.gens)
  first = true
  for i = g
    if first
      print(io, i)
      first = false
    else
      print(io, ", ", i)
    end
  end
  print(io, "")
end

function Base.show(io::IO, ::MIME"text/html", I::MPolyIdeal)
  print(io, "\$")
  math_html(io, I)
  print(io, "\$")
end

function math_html(io::IO, I::MPolyIdeal)
  print(io, "\\text{ideal generated by: }")
  g = collect(I.gens)
  first = true
  for i = g
    if first
      math_html(io, i)
      first = false
    else
      print(io, ", ")
      math_html(io, i)
    end
  end
  print(io, "")
end



function ideal(g::Array{T, 1}) where {T <: MPolyElem}
  return MPolyIdeal(g)
end

function ideal(g::Array{Any, 1})
  return ideal(typeof(g[1])[x for x = g])
end

function ideal(Rx::MPolyRing, g::Array{<:Any, 1})
  f = elem_type(Rx)[Rx(f) for f = g]
  return ideal(f)
end

function singular_assure(I::MPolyIdeal)
  if !isdefined(I.gens, :S)
    I.gens.S = Singular.Ideal(I.gens.Sx, [convert(I.gens.Sx, x) for x = I.gens.O])
  end
end

function oscar_assure(I::MPolyIdeal)
  if !isdefined(I.gens, :O)
    I.gens.S = Singular.Ideal(I.gens.Sx, [convert(I.gens.Sx, x) for x = I.gens.O])
  end
end

function Base.:*(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, I.gens.S * J.gens.S)
end

function Base.:+(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, I.gens.S + J.gens.S)
end
Base.:-(I::MPolyIdeal, J::MPolyIdeal) = I+J

function Base.:(==)(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return Singular.equal(I.gens.S, J.gens.S)
end

function Base.:^(I::MPolyIdeal, j::Int)
  singular_assure(I)
  return MPolyIdeal(I.gens.Ox, I.gens.S^j)
end

function Base.intersect(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, Singular.intersection(I.gens.S, J.gens.S))
end

function ngens(I::MPolyIdeal)
  return length(I.gens)
end

function Base.issubset(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return Singular.contains(I.gens.S, J.gens.S)
end

function gens(I::MPolyIdeal)
  return [I.gens[Val(:O), i] for i=1:ngens(I)]
end

function groebner_assure(I::MPolyIdeal)
  if !isdefined(I, :gb)
    singular_assure(I)
    I.gb = BiPolyArray(I.gens.Ox, Singular.std(I.gens.S))
  end
end

function dim(I::MPolyIdeal)
  if I.dim > -1
    return I.dim
  end
  groebner_assure(I)
  I.dim = Singular.dimension(I.gb.S)
  return I.dim
end

function Base.in(f::MPolyElem, I::MPolyIdeal)
  groebner_assure(I)
  Sx = base_ring(I.gb.S)
  return Singular.iszero(reduce(convert(Sx, f), I.gb.S))
end

function base_ring(I::MPolyIdeal)
  return I.gens.Ox
end

function im_func(f::MPolyElem, S::MPolyRing, i::Array{Int, 1})
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
  Hecke.@declare_other
  i::Array{Int, 1}

  function MPolyHom_vars{T1, T2}(R::T1, S::T2, i::Array{Int, 1}) where {T1 <: MPolyRing, T2 <: MPolyRing}
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
      return MPolyHom_vars(R, S, i)
    end
    error("type not supported")
  end
end

(f::MPolyHom_vars)(g::MPolyElem) = image(f, g)

function Hecke.hom(R::MPolyRing, S::MPolyRing, i::Array{Int, 1})
  return MPolyHom_vars{typeof(R), typeof(S)}(R, S, i)
end

function eliminate(I::MPolyIdeal, l::Array{<:MPolyElem, 1})
  #wrong!!! ideal intersection is not eliminate!!! 
  J = ideal(base_ring(I), [x for x in gens(base_ring(I)) if !(x in l)])
  return intersect(I, J)
end

function groebner_basis(I::MPolyIdeal)
  groebner_assure(I)
  return collect(I.gb)
end

function groebner_basis(I::MPolyIdeal, ord::Symbol)
  R = singular_ring(base_ring(I), ord)
  i = Singular.Ideal(R, [convert(R, x) for x = gens(I)])
  return collect(BiPolyArray(base_ring(I), i))
end

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

# Returns the degree of the k-th term of f weighted by w,
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

function _isless_matrix(f::MPolyElem, k::Int, l::Int, M::Union{ Array{T, 2}, MatElem{T} }) where T
  ek = exponent_vector(f, k)
  el = exponent_vector(f, l)
  n = nvars(parent(f))
  for i = 1:n
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
function lt_from_ordering(R::MPolyRing, ord::Symbol)
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

function lt_from_ordering(R::MPolyRing, M::Union{ Array{T, 2}, MatElem{T} }) where T
  @assert size(M, 1) == nvars(R) && size(M, 2) == nvars(R) "Matrix dimensions have to match number of variables"

  return (f, k, l) -> _isless_matrix(f, k, l, M)
end

function terms(f::MPolyElem, ord::Function)
  perm = _perm_of_terms(f, ord)
  return ( term(f, perm[i]) for i = 1:length(f) )
end

function coeffs(f::MPolyElem, ord::Function)
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

for s in (:terms, :coeffs, :exponent_vectors, :monomials)
  @eval begin
    function ($s)(f::MPolyElem, ord::Symbol)
      R = parent(f)
      if ord == ordering(R)
        return ($s)(f)
      end

      lt = lt_from_ordering(R, ord)
      return ($s)(f, lt)
    end

    function ($s)(f::MPolyElem, M::Union{ Array{T, 2}, MatElem{T} }) where T
      R = parent(f)
      lt = lt_from_ordering(R, M)
      return ($s)(f, lt)
    end

    function ($s)(f::MPolyElem, ord::Symbol, weights::Vector{Int})
      R = parent(f)
      lt = lt_from_ordering(R, ord, weights)
      return ($s)(f, lt)
    end
  end
end

for s in ("term", "coeff", "monomial")
  @eval begin
    function ($(Symbol("leading_$s")))(args...)
      return first($(Symbol("$(s)s"))(args...))
    end
  end
end

function leading_term(f::MPolyElem)
  return leading_term(f, ordering(parent(f)))
end

function leading_coeff(f::MPolyElem)
  return leading_coeff(f, ordering(parent(f)))
end

function leading_monomial(f::MPolyElem)
  return leading_monomial(f, ordering(parent(f)))
end

function leading_ideal(g::Array{T, 1}, args...) where { T <: MPolyElem }
  return ideal([ leading_monomial(f, args...) for f in g ])
end

function leading_ideal(g::Array{Any, 1}, args...)
  return leading_ideal(typeof(g[1])[ f for f in g ], args...)
end

function leading_ideal(Rx::MPolyRing, g::Array{Any, 1}, args...)
  h = elem_type(Rx)[ Rx(f) for f in g ]
  return leading_ideal(h, args...)
end

function leading_ideal(I::MPolyIdeal)
  return leading_ideal(groebner_basis(I))
end

function leading_ideal(I::MPolyIdeal, ord::Symbol)
  return leading_ideal(groebner_basis(I, ord), ord)
end


function factor(f::MPolyElem)
  I = ideal(parent(f), [f])
  fS = Singular.factor(I.gens[Val(:S), 1])
  R = parent(f)
  return Nemo.Fac(convert(R, fS.unit), Dict(convert(R, k) =>v for (k,v) = fS.fac))
end
#end #MPolyModule
