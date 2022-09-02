export pbw_algebra, build_ctx, PBWAlgElem, PBWAlgRing, is_two_sided

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
  basering::PBWAlgRing{T, S}
  sdata::Singular.sideal{Singular.spluralg{S}}
  two_sided::Bool
  gb::Singular.sideal{Singular.spluralg{S}}
  # Singular.jl may or may not keep track of two-sidedness correctly
  function PBWAlgIdeal{T, S}(p::PBWAlgRing{T, S},
                             d::Singular.sideal{Singular.spluralg{S}},
                             t::Bool) where {T, S}
    d.isTwoSided = t
    return new(p, d, t)
  end
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
  return Expr(:sequence, Expr(:text, "PBW-algebra over "),
                         expressify(coefficient_ring(a); context=context),
                         Expr(:text, " with relations "),
                         Expr(:series, rel...))
end

@enable_all_show_via_expressify PBWAlgRing

####

function length(a::PBWAlgElem)
  return length(a.sdata)
end

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

struct owrap{S, T}
  ring::S
  iter::T
end

function Base.length(x::owrap)
   return length(x.iter)
end

function exponent_vectors(a::PBWAlgElem)
  return exponent_vectors(a.sdata)
end

function terms(a::PBWAlgElem)
  return owrap(parent(a), terms(a.sdata))
end

function Base.eltype(x::owrap{<:PBWAlgRing{T, S}, <:Singular.SPolyTerms}) where {T, S}
   return PBWAlgElem{T, S}
end

function Base.iterate(a::owrap{<:PBWAlgRing, <:Singular.SPolyTerms})
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function Base.iterate(a::owrap{<:PBWAlgRing, <:Singular.SPolyTerms}, state)
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function monomials(a::PBWAlgElem)
  return owrap(parent(a), monomials(a.sdata))
end

function Base.eltype(x::owrap{<:PBWAlgRing{T, S}, <:Singular.SPolyMonomials}) where {T, S}
   return PBWAlgElem{T, S}
end

function Base.iterate(a::owrap{<:PBWAlgRing, <:Singular.SPolyMonomials})
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function Base.iterate(a::owrap{<:PBWAlgRing, <:Singular.SPolyMonomials}, state)
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (PBWAlgElem(a.ring, b[1]), b[2])
end

function coefficients(a::PBWAlgElem)
  return owrap(parent(a), coefficients(a.sdata))
end

function Base.eltype(x::owrap{<:PBWAlgRing{T, S}, <:Singular.SPolyCoeffs}) where {T, S}
   return T
end

function Base.iterate(a::owrap{<:PBWAlgRing{T}, <:Singular.SPolyCoeffs}) where T
  b = Base.iterate(a.iter)
  b == nothing && return b
  return (coefficient_ring(a.ring)(b[1])::T, b[2])
end

function Base.iterate(a::owrap{<:PBWAlgRing{T}, <:Singular.SPolyCoeffs}, state) where T
  b = Base.iterate(a.iter, state)
  b == nothing && return b
  return (coefficient_ring(a.ring)(b[1])::T, b[2])
end

function build_ctx(R::PBWAlgRing)
  return owrap(R, MPolyBuildCtx(R.sring))
end

function push_term!(M::owrap{<:PBWAlgRing{T,S}, <:MPolyBuildCtx}, c, e::Vector{Int}) where {T, S}
  c = coefficient_ring(M.ring)(c)::T
  c = base_ring(M.ring.sring)(c)::S
  push_term!(M.iter, c, e)
end

function finish(M::owrap{<:PBWAlgRing{T,S}, <:MPolyBuildCtx}) where {T, S}
  return PBWAlgElem(M.ring, finish(M.iter))
end

####

function ngens(R::PBWAlgRing)
  return Singular.nvars(R.sring)
end

function gens(R::PBWAlgRing)
  return elem_type(R)[PBWAlgElem(R, x) for x in gens(R.sring)]
end

function gen(R::PBWAlgRing, i::Int)
  return PBWAlgElem(R, gen(R.sring, i))
end

function Base.getindex(R::PBWAlgRing, i::Int)
  return gen(R, i)
end

function zero(R::PBWAlgRing)
  return PBWAlgElem(R, zero(R.sring))
end

function one(R::PBWAlgRing)
  return PBWAlgElem(R, one(R.sring))
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

function (R::PBWAlgRing)(a::PBWAlgElem)
  parent(a) == R || error("coercion impossible")
  return a
end

function (R::PBWAlgRing)(cs::AbstractVector, es::AbstractVector{Vector{Int}})
  z = build_ctx(R)
  @assert length(cs) == length(es)
  for (c, e) in zip(cs, es)
    push_term!(z, c, e)
  end
  return finish(z)
end

####
@doc Markdown.doc"""
    pbw_algebra(R::MPolyRing{T}, rel, ord::MonomialOrdering) where T

Given a multivariate polynomial ring `R` over a field, say ``R=K[x_1, \dots, x_n]``, given
a strictly upper triangular matrix `rel` with entries in `R` of type ``c_{ij} \cdot x_ix_j+d_{ij}``,
where the ``c_{ij}`` are nonzero scalars and where we think of the ``x_jx_i = c_{ij} \cdot x_ix_j+d_{ij}``
as setting up relations in the free associative algebra ``K\langle x_1, \dots , x_n\rangle``, and given
an ordering `ord` on ``\text{Mon}(x_1, \dots, x_n)``, return the PBW-algebra
```math
A = K\langle x_1, \dots , x_n \mid x_jx_i = c_{ij} \cdot x_ix_j+d_{ij},  \ 1\leq i<j \leq n \rangle.
```

!!! note
    The input data gives indeed rise to  a PBW-algebra if:
    - The ordering `ord` is admissible for `A`.
    - The standard monomials in ``K\langle x_1, \dots , x_n\rangle`` represent a `K`-basis for `A`.
    See the definition of PBW-algebras in the OSCAR documentation for details.

!!! warning
    The `K`-basis condition above is not checked by the function.

# Examples
```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)))
(PBW-algebra over Rational Field with relations y*x = x*y, z*x = x*z, z*y = y*z + 1, PBWAlgElem{fmpq, Singular.n_Q}[x, y, z])
```
"""
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

####

function base_ring(a::PBWAlgIdeal)
  return a.basering
end

function is_two_sided(a::PBWAlgIdeal)
  return a.two_sided
end

function ngens(a::PBWAlgIdeal)
  return ngens(a.sdata)
end

function gens(a::PBWAlgIdeal{T, S}) where {T, S}
  R = base_ring(a)
  return PBWAlgElem{T, S}[PBWAlgElem(R, x) for x in gens(a.sdata)]
end

function AbstractAlgebra.expressify(a::PBWAlgIdeal; context = nothing)
  return Expr(:call, :ideal, [expressify(g, context = context) for g in collect(gens(a))]...)
end

@enable_all_show_via_expressify PBWAlgIdeal

@doc Markdown.doc"""
    ideal(g::Vector{<:PBWAlgElem}; two_sided=false)

Given a vector `g` of elements in a polynomial ring `A`, say, return the left ideal of `A` generated by these elements.

    ideal(A::PBWAlgRing{T, S}, g::AbstractVector; two_sided=false) where {T, S}

Given a vector `g` of elements in `A`, return the left ideal of `A` generated by these elements.

# Examples
```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x,y,z) = pbw_algebra(R, REL, deglex(gens(R)))
(PBW-algebra over Rational Field with relations y*x = x*y, z*x = x*z, z*y = y*z + 1, PBWAlgElem{fmpq, Singular.n_Q}[x, y, z])

julia> I = ideal(A, [x^2*y^2, x*z+y*z])
ideal(x^2*y^2, x*z + y*z)
```
"""
function ideal(g::Vector{<:PBWAlgElem}; two_sided=false)
  @assert length(g) > 0
  R = parent(g[1])
  @assert all(x->parent(x) == R, g)
  return ideal(R, g; two_sided=two_sided)
end

function ideal(R::PBWAlgRing{T, S}, g::AbstractVector; two_sided=false) where {T, S}
  i = Singular.sideal{Singular.spluralg{S}}(R.sring, [R(x).sdata for x in g], two_sided)
  return PBWAlgIdeal{T, S}(R, i, two_sided)
end

function groebner_assure!(a::PBWAlgIdeal)
  if !isdefined(a, :gb)
    a.gb = Singular.std(a.sdata)
  end
  return a
end

function iszero(I::PBWAlgIdeal)
  return iszero(I.sdata)
end

function _one_check(I::Singular.sideal)
  for g in gens(I)
    if is_constant(g) && is_unit(leading_coefficient(g))
      return true
    end
  end
  return false
end

function isone(I::PBWAlgIdeal)
  if iszero(I)
    return false
  end
  if _one_check(I.sdata)
    return true
  end
  groebner_assure!(I)
  return _one_check(I.gb)
end


function Base.:+(a::PBWAlgIdeal{T, S}, b::PBWAlgIdeal{T, S}) where {T, S}
  return PBWAlgIdeal{T, S}(base_ring(a), a.sdata + b.sdata,
                           is_two_sided(a)&&is_two_sided(b))
end

function Base.:*(a::PBWAlgIdeal{T, S}, b::PBWAlgIdeal{T, S}) where {T, S}
  return PBWAlgIdeal{T, S}(base_ring(a), a.sdata + b.sdata, is_two_sided(b))
end

function Base.:^(a::PBWAlgIdeal{T, S}, b::Int) where {T, S}
  return PBWAlgIdeal{T, S}(base_ring(a), a.sdata^b, is_two_sided(a))
end

function Base.intersect(a::PBWAlgIdeal{T, S}, b::PBWAlgIdeal{T, S}) where {T, S}
  return PBWAlgIdeal{T, S}(base_ring(a), Singular.intersection(a.sdata, b.sdata),
                           is_two_sided(a)&&is_two_sided(b))
end

function ideal_membership(f::PBWAlgElem{T, S}, I::PBWAlgIdeal{T, S}) where {T, S}
  parent(f) == base_ring(I) || error("parent mismatch")
  is_two_sided(I) && throw(NotImplementedError(:ideal_membership, f, I))
  groebner_assure!(I)
  return Singular.iszero(Singular.reduce(f.sdata, I.gb))
end

function Base.in(f::PBWAlgElem, I::PBWAlgIdeal)
  return ideal_membership(f, I)
end
