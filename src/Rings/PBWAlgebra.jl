export pbw_algebra, build_ctx

mutable struct PBWAlgRing{T, S} <: NCRing
  sring::Singular.PluralRing{S}
  relations
  ordering::MonomialOrdering
end

mutable struct PBWAlgElem{T, S} <: NCRingElem
  parent::PBWAlgRing{T, S}
  sdata::Singular.spluralg{S}
end

mutable struct PBWAlgIdeal{T, S}
  parent::PBWAlgRing{T, S}
  sdata::Singular.sideal{Singular.spluralg{S}}
end

####

elem_type(R::PBWAlgRing{T, S}) where {T, S} = PBWAlgElem{T, S}

parent_type(a::PBWAlgElem{T, S}) where {T, S} = PBWAlgRing{T, S}

parent(a::PBWAlgElem) = a.parent

symbols(a::PBWAlgRing) = symbols(a.sring)

coefficient_ring(a::PBWAlgRing) = coefficient_ring(base_ring(a.ordering))

coefficient_ring(a::PBWAlgElem) = coefficient_ring(parent(a))

function Base.deepcopy_internal(a::PBWAlgElem, dict::IdDict)
  return PBWAlgElem(parent(a), deepcopy_internal(a.sdata, dict))
end

function expressify(a::PBWAlgElem; context = nothing)
  return expressify(a.sdata; context=context)
end

@enable_all_show_via_expressify PBWAlgElem

function expressify(a::PBWAlgRing; context = nothing)
  x = symbols(a)
  n = length(x)
  rel = [Expr(:call, :(==), Expr(:call, :*, x[j], x[i]), expressify(a.relations[i,j]))
         for i in 1:n-1 for j in i+1:n]
  return Expr(:sequence, Expr(:text, "G-algebra over "),
                         expressify(coefficient_ring(a); context=context),
                         Expr(:text, " with relations "),
                         Expr(:series, rel...))
end

@enable_all_show_via_expressify PBWAlgRing

####

function leading_exponent_vector(a::PBWAlgElem)
  return leading_exponent_vector(a.sdata)
end

function leading_coefficient(a::PBWAlgElem{T})::T where T
  return coefficient_ring(a)(leading_coefficient(a.sdata))
end

function trailing_coefficient(a::PBWAlgElem{T})::T where T
  return coefficient_ring(a)(trailing_coefficient(a.sdata))
end

function constant_coefficient(a::PBWAlgElem{T})::T where T
  return coefficient_ring(a)(constant_coefficient(a.sdata))
end

function leading_term(a::PBWAlgElem)
  return PBWAlgElem(parent(a), leading_term(a.sdata))
end

function leading_monomial(a::PBWAlgElem)
  return PBWAlgElem(parent(a), leading_monomial(a.sdata))
end

function tail(a::PBWAlgElem)
  return PBWAlgElem(parent(a), tail(a.sdata))
end

struct oitwrap{S, T}
  ring::S
  iter::T
end

function exponent_vectors(a::PBWAlgElem)
  return exponent_vectors(a.sdata)
end

function terms(a::PBWAlgElem)
  return oitwrap(parent(a), terms(a.sdata))
end

function Base.iterate(a::oitwrap{<:PBWAlgRing, <:Singular.SPolyTerms})
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function Base.iterate(a::oitwrap{<:PBWAlgRing, <:Singular.SPolyTerms}, state)
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function monomials(a::PBWAlgElem)
  return oitwrap(parent(a), monomials(a.sdata))
end

function Base.iterate(a::oitwrap{<:PBWAlgRing, <:Singular.SPolyMonomials})
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function Base.iterate(a::oitwrap{<:PBWAlgRing, <:Singular.SPolyMonomials}, state)
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function coefficients(a::PBWAlgElem)
  return oitwrap(parent(a), coefficients(a.sdata))
end

function Base.iterate(a::oitwrap{<:PBWAlgRing{T}, <:Singular.SPolyCoeffs}) where T
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (coefficient_ring(a.ring)(b[1])::T, b[2])
end

function Base.iterate(a::oitwrap{<:PBWAlgRing{T}, <:Singular.SPolyCoeffs}, state) where T
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (coefficient_ring(a.ring)(b[1])::T, b[2])
end

function build_ctx(R::PBWAlgRing)
  return oitwrap(R, MPolyBuildCtx(R.sring))
end

function push_term!(M::oitwrap{<:PBWAlgRing{T,S}, <:MPolyBuildCtx}, c, e::Vector{Int}) where {T, S}
  c = coefficient_ring(M.ring)(c)::T
  c = base_ring(M.ring.sring)(c)::S
  push_term!(M.iter, c, e)
end

function finish(M::oitwrap{<:PBWAlgRing{T,S}, <:MPolyBuildCtx}) where {T, S}
  return PBWAlgElem(M.ring, finish(M.iter))
end

####

function zero(R::PBWAlgRing)
  return PBWAlgElem(R, zero(R.sring))
end

function one(R::PBWAlgRing)
  return PBWAlgElem(R, one(R.sring))
end

function gens(R::PBWAlgRing)
  return [PBWAlgElem(R, x) for x in gens(R.sring)]
end

function Base.:(==)(a::PBWAlgElem, b::PBWAlgElem)
  return a.sdata == b.sdata
end

function Base.:+(a::PBWAlgElem, b::PBWAlgElem)
  return PBWAlgElem(parent(a), a.sdata + b.sdata)
end

function Base.:-(a::PBWAlgElem, b::PBWAlgElem)
  return PBWAlgElem(parent(a), a.sdata - b.sdata)
end

function Base.:-(a::PBWAlgElem)
  return PBWAlgElem(parent(a), -a.sdata)
end

function Base.:*(a::PBWAlgElem, b::PBWAlgElem)
  return PBWAlgElem(parent(a), a.sdata*b.sdata)
end

function Base.:^(a::PBWAlgElem, b::Int)
  return PBWAlgElem(parent(a), a.sdata^b)
end

####

function AbstractAlgebra.promote_rule(::Type{PBWAlgElem{T, S}}, ::Type{PBWAlgElem{T}}) where {T, S}
  return PBWAlgElem{T, S}
end

function AbstractAlgebra.promote_rule(::Type{PBWAlgElem{T, S}}, ::Type{U}) where {T, S, U}
  a = AbstractAlgebra.promote_rule(T, U)
  return a == T ? PBWAlgElem{T, S} : Union{}
end

function (R::PBWAlgRing{T, S})(c::T) where {T, S}
  c = coefficient_ring(R)(c)::T
  c = base_ring(R.sring)(c)::S
  return PBWAlgElem(R, R.sring(c))
end

function (R::PBWAlgRing{T, S})(c::Integer) where {T, S}
  c = base_ring(R.sring)(c)::S
  return PBWAlgElem(R, R.sring(c))
end

####

function pbw_algebra(r::MPolyRing{T}, rel, ord::MonomialOrdering) where T
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
  R = PBWAlgRing{T, S}(s, rel, ord)
  return R, [PBWAlgElem(R, x) for x in gs]
end

