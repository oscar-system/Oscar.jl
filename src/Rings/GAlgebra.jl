export g_algebra, build_ctx

mutable struct GAlgRing{T, S} <: NCRing
  sring::Singular.PluralRing{S}
  relations
  ordering::MonomialOrdering
end

mutable struct GAlgElem{T, S} <: NCRingElem
  parent::GAlgRing{T, S}
  sdata::Singular.spluralg{S}
end

mutable struct GAlgIdeal{T, S}
  parent::GAlgRing{T, S}
  sdata::Singular.sideal{Singular.spluralg{S}}
end

####

elem_type(R::GAlgRing{T, S}) where {T, S} = GAlgElem{T, S}

parent_type(a::GAlgElem{T, S}) where {T, S} = GAlgRing{T, S}

parent(a::GAlgElem) = a.parent

symbols(a::GAlgRing) = symbols(a.sring)

coefficient_ring(a::GAlgRing) = coefficient_ring(base_ring(a.ordering))

coefficient_ring(a::GAlgElem) = coefficient_ring(parent(a))

function Base.deepcopy_internal(a::GAlgElem, dict::IdDict)
  return GAlgElem(parent(a), deepcopy_internal(a.sdata, dict))
end

function expressify(a::GAlgElem; context = nothing)
  return expressify(a.sdata; context=context)
end

@enable_all_show_via_expressify GAlgElem

function expressify(a::GAlgRing; context = nothing)
  x = symbols(a)
  n = length(x)
  rel = [Expr(:call, :(==), Expr(:call, :*, x[j], x[i]), expressify(a.relations[i,j]))
         for i in 1:n-1 for j in i+1:n]
  return Expr(:sequence, Expr(:text, "G-algebra over "),
                         expressify(coefficient_ring(a); context=context),
                         Expr(:text, " with relations "),
                         Expr(:series, rel...))
end

@enable_all_show_via_expressify GAlgRing

####

function leading_exponent_vector(a::GAlgElem)
  return leading_exponent_vector(a.sdata)
end

function leading_coefficient(a::GAlgElem{T})::T where T
  return coefficient_ring(a)(leading_coefficient(a.sdata))
end

function trailing_coefficient(a::GAlgElem{T})::T where T
  return coefficient_ring(a)(trailing_coefficient(a.sdata))
end

function constant_coefficient(a::GAlgElem{T})::T where T
  return coefficient_ring(a)(constant_coefficient(a.sdata))
end

function leading_term(a::GAlgElem)
  return GAlgElem(parent(a), leading_term(a.sdata))
end

function leading_monomial(a::GAlgElem)
  return GAlgElem(parent(a), leading_monomial(a.sdata))
end

function tail(a::GAlgElem)
  return GAlgElem(parent(a), tail(a.sdata))
end

struct oitwrap{S, T}
  ring::S
  iter::T
end

function exponent_vectors(a::GAlgElem)
  return exponent_vectors(a.sdata)
end

function terms(a::GAlgElem)
  return oitwrap(parent(a), terms(a.sdata))
end

function Base.iterate(a::oitwrap{<:GAlgRing, <:Singular.SPolyTerms})
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (GAlgElem(a.ring, b[1]), b[2])
end

function Base.iterate(a::oitwrap{<:GAlgRing, <:Singular.SPolyTerms}, state)
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (GAlgElem(a.ring, b[1]), b[2])
end

function monomials(a::GAlgElem)
  return oitwrap(parent(a), monomials(a.sdata))
end

function Base.iterate(a::oitwrap{<:GAlgRing, <:Singular.SPolyMonomials})
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (GAlgElem(a.ring, b[1]), b[2])
end

function Base.iterate(a::oitwrap{<:GAlgRing, <:Singular.SPolyMonomials}, state)
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (GAlgElem(a.ring, b[1]), b[2])
end

function coefficients(a::GAlgElem)
  return oitwrap(parent(a), coefficients(a.sdata))
end

function Base.iterate(a::oitwrap{<:GAlgRing{T}, <:Singular.SPolyCoeffs}) where T
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (coefficient_ring(a.ring)(b[1])::T, b[2])
end

function Base.iterate(a::oitwrap{<:GAlgRing{T}, <:Singular.SPolyCoeffs}, state) where T
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (coefficient_ring(a.ring)(b[1])::T, b[2])
end

function build_ctx(R::GAlgRing)
  return oitwrap(R, MPolyBuildCtx(R.sring))
end

function push_term!(M::oitwrap{<:GAlgRing{T,S}, <:MPolyBuildCtx}, c, e::Vector{Int}) where {T, S}
  c = coefficient_ring(M.ring)(c)::T
  c = base_ring(M.ring.sring)(c)::S
  push_term!(M.iter, c, e)
end

function finish(M::oitwrap{<:GAlgRing{T,S}, <:MPolyBuildCtx}) where {T, S}
  return GAlgElem(M.ring, finish(M.iter))
end

####

function zero(R::GAlgRing)
  return GAlgElem(R, zero(R.sring))
end

function one(R::GAlgRing)
  return GAlgElem(R, one(R.sring))
end

function gens(R::GAlgRing)
  return [GAlgElem(R, x) for x in gens(R.sring)]
end

function Base.:(==)(a::GAlgElem, b::GAlgElem)
  return a.sdata == b.sdata
end

function Base.:+(a::GAlgElem, b::GAlgElem)
  return GAlgElem(parent(a), a.sdata + b.sdata)
end

function Base.:-(a::GAlgElem, b::GAlgElem)
  return GAlgElem(parent(a), a.sdata - b.sdata)
end

function Base.:-(a::GAlgElem)
  return GAlgElem(parent(a), -a.sdata)
end

function Base.:*(a::GAlgElem, b::GAlgElem)
  return GAlgElem(parent(a), a.sdata*b.sdata)
end

function Base.:^(a::GAlgElem, b::Int)
  return GAlgElem(parent(a), a.sdata^b)
end

####

function AbstractAlgebra.promote_rule(::Type{GAlgElem{T, S}}, ::Type{GAlgElem{T}}) where {T, S}
  return GAlgElem{T, S}
end

function AbstractAlgebra.promote_rule(::Type{GAlgElem{T, S}}, ::Type{U}) where {T, S, U}
  a = AbstractAlgebra.promote_rule(T, U)
  return a == T ? GAlgElem{T, S} : Union{}
end

function (R::GAlgRing{T, S})(c::T) where {T, S}
  c = coefficient_ring(R)(c)::T
  c = base_ring(R.sring)(c)::S
  return GAlgElem(R, R.sring(c))
end

function (R::GAlgRing{T, S})(c::Integer) where {T, S}
  c = base_ring(R.sring)(c)::S
  return GAlgElem(R, R.sring(c))
end

####

function g_algebra(r::MPolyRing{T}, rel, ord::MonomialOrdering) where T
  n = nvars(r)
  nrows(rel) == n && ncols(rel) == n || error("oops")
  scr = singular_coeff_ring(coefficient_ring(r))
  S = elem_type(scr)
  sr, _ = Singular.PolynomialRing(scr, [string(x) for x in symbols(r)]; ordering = singular(ord))
  sr::Singular.PolyRing{S}
  C = Singular.zero_matrix(sr, n, n)
  D = Singular.zero_matrix(sr, n, n)
  for i in 1:n-1, j in i+1:n
    t = sr(rel[i,j])
    leading_monomial(t) == gen(sr, i)*gen(sr, j) || error("incorrect leading monomial in relations")
    C[i,j] = sr(leading_coefficient(t))
    D[i,j] = tail(t)
  end
  s, gs = Singular.GAlgebra(sr, C, D)
  R = GAlgRing{T, S}(s, rel, ord)
  return R, [GAlgElem(R, x) for x in gs]
end

