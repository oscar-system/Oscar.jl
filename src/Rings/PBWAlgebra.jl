# Use attribute :is_weyl_algebra to permit better printing (see expressify, below)
@attributes mutable struct PBWAlgRing{T, S} <: NCRing
  sring::Singular.PluralRing{S}
  relations::Singular.smatrix{Singular.spoly{S}}
  coeff_ring
  poly_ring
  opposite::PBWAlgRing{T, S}

  function PBWAlgRing{T, S}(sring, relations, coeff_ring, poly_ring) where {T, S}
    return new{T, S}(sring, relations, coeff_ring, poly_ring)
  end
end

struct PBWAlgOppositeMap{T, S}
  source::PBWAlgRing{T, S}  # target is _opposite(source)
end

mutable struct PBWAlgElem{T, S} <: NCRingElem
  parent::PBWAlgRing{T, S}
  sdata::Singular.spluralg{S}
end

mutable struct PBWAlgIdeal{D, T, S}
  basering::PBWAlgRing{T, S}
  sdata::Singular.sideal{Singular.spluralg{S}}    # the gens of this ideal, always defined
  sopdata::Singular.sideal{Singular.spluralg{S}}  # the gens mapped to the opposite
  gb::Singular.sideal{Singular.spluralg{S}}
  opgb::Singular.sideal{Singular.spluralg{S}}
  # Singular.jl may or may not keep track of two-sidedness correctly
  function PBWAlgIdeal{D, T, S}(p::PBWAlgRing{T, S},
                      d::Singular.sideal{Singular.spluralg{S}}) where {D, T, S}
    d.isTwoSided = (D == 0)
    return new{D, T, S}(p, d)
  end
  function PBWAlgIdeal{D, T, S}(p::PBWAlgRing{T, S},
                      d::Singular.sideal{Singular.spluralg{S}},
                    opd::Singular.sideal{Singular.spluralg{S}}) where {D, T, S}
    d.isTwoSided = (D == 0)
    opd.isTwoSided = (D == 0)
    return new{D, T, S}(p, d, opd)
  end
end
# the meaning of the direction parameter D
is_left(a::PBWAlgIdeal{D}) where D = (D <= 0)
is_right(a::PBWAlgIdeal{D}) where D = (D >= 0)
is_two_sided(a::PBWAlgIdeal{D}) where D = (D == 0)

####

function is_domain_type(a::Type{U}) where {T, U <: PBWAlgElem{T}}
   return is_domain_type(T)
end

function is_exact_type(a::Type{U}) where {T, U <: PBWAlgElem{T}}
   return is_exact_type(T)
end

elem_type(::PBWAlgRing{T, S}) where {T, S} = PBWAlgElem{T, S}

elem_type(::Type{PBWAlgRing{T, S}}) where {T, S} = PBWAlgElem{T, S}

parent_type(::PBWAlgElem{T, S}) where {T, S} = PBWAlgRing{T, S}

parent_type(::Type{PBWAlgElem{T, S}}) where {T, S} = PBWAlgRing{T, S}

parent(a::PBWAlgElem) = a.parent

symbols(a::PBWAlgRing) = symbols(a.sring)

coefficient_ring(a::PBWAlgRing) = a.coeff_ring

coefficient_ring(a::PBWAlgElem) = coefficient_ring(parent(a))

base_ring(a::PBWAlgRing) = a.poly_ring

base_ring(a::PBWAlgElem) = base_ring(parent(a))

function Base.deepcopy_internal(a::PBWAlgElem, dict::IdDict)
  return PBWAlgElem(parent(a), deepcopy_internal(a.sdata, dict))
end

function Base.hash(a::PBWAlgElem, h::UInt)
  return hash(a.sdata, h)
end

function expressify(a::PBWAlgElem; context = nothing)
  return expressify(a.sdata; context=context)
end

@enable_all_show_via_expressify PBWAlgElem

function expressify(a::PBWAlgRing; context = nothing)
  x = symbols(a)
  n = length(x)
  # Next if stmt handles special printing for Weyl algebras
  if get_attribute(a, :is_weyl_algebra) === :true
      return Expr(:sequence, Expr(:text, "Weyl-algebra over "),
                             expressify(coefficient_ring(a); context=context),
                             Expr(:text, " in variables ("),
                             Expr(:series, first(x,div(n,2))...),
                             Expr(:text, ")"))
  end
  rel = [Expr(:call, :(==), Expr(:call, :*, x[j], x[i]), expressify(a.relations[i,j]))
         for i in 1:n-1 for j in i+1:n]
  return Expr(:sequence, Expr(:text, "PBW-algebra over "),
                         expressify(coefficient_ring(a); context=context),
                         Expr(:text, " in "),
                         Expr(:series, x...),
                         Expr(:text, " with relations "),
                         Expr(:series, rel...))
end

@enable_all_show_via_expressify PBWAlgRing

#### AA prefix here because these all use the ordering in the parent

function length(a::PBWAlgElem)
  return length(a.sdata)
end

function AbstractAlgebra.leading_exponent_vector(a::PBWAlgElem)
  return AbstractAlgebra.leading_exponent_vector(a.sdata)
end

function AbstractAlgebra.leading_coefficient(a::PBWAlgElem{T})::T where T
  return coefficient_ring(a)(AbstractAlgebra.leading_coefficient(a.sdata))
end

function AbstractAlgebra.trailing_coefficient(a::PBWAlgElem{T})::T where T
  return coefficient_ring(a)(AbstractAlgebra.trailing_coefficient(a.sdata))
end

function constant_coefficient(a::PBWAlgElem{T})::T where T
  return coefficient_ring(a)(constant_coefficient(a.sdata))
end

function AbstractAlgebra.leading_term(a::PBWAlgElem)
  return PBWAlgElem(parent(a), AbstractAlgebra.leading_term(a.sdata))
end

function AbstractAlgebra.leading_monomial(a::PBWAlgElem)
  return PBWAlgElem(parent(a), AbstractAlgebra.leading_monomial(a.sdata))
end

function AbstractAlgebra.tail(a::PBWAlgElem)
  return PBWAlgElem(parent(a), AbstractAlgebra.tail(a.sdata))
end


function AbstractAlgebra.exponent_vectors(a::PBWAlgElem)
  return AbstractAlgebra.exponent_vectors(a.sdata)
end

function terms(a::PBWAlgElem)
  return OscarPair(parent(a), terms(a.sdata))
end

function Base.length(x::OscarPair{<:PBWAlgRing, <:Singular.SPolyTerms})
   return length(x.second)
end

function Base.eltype(x::OscarPair{<:PBWAlgRing{T, S}, <:Singular.SPolyTerms}) where {T, S}
   return PBWAlgElem{T, S}
end

function Base.iterate(a::OscarPair{<:PBWAlgRing, <:Singular.SPolyTerms})
  b = Base.iterate(a.second)
  b === nothing && return b
  return (PBWAlgElem(a.first, b[1]), b[2])
end

function Base.iterate(a::OscarPair{<:PBWAlgRing, <:Singular.SPolyTerms}, state)
  b = Base.iterate(a.second, state)
  b === nothing && return b
  return (PBWAlgElem(a.first, b[1]), b[2])
end

function AbstractAlgebra.monomials(a::PBWAlgElem)
  return OscarPair(parent(a), AbstractAlgebra.monomials(a.sdata))
end

function Base.length(x::OscarPair{<:PBWAlgRing, <:Singular.SPolyMonomials})
   return length(x.second)
end

function Base.eltype(x::OscarPair{<:PBWAlgRing{T, S}, <:Singular.SPolyMonomials}) where {T, S}
   return PBWAlgElem{T, S}
end

function Base.iterate(a::OscarPair{<:PBWAlgRing, <:Singular.SPolyMonomials})
  b = Base.iterate(a.second)
  b === nothing && return b
  return (PBWAlgElem(a.first, b[1]), b[2])
end

function Base.iterate(a::OscarPair{<:PBWAlgRing, <:Singular.SPolyMonomials}, state)
  b = Base.iterate(a.second, state)
  b === nothing && return b
  return (PBWAlgElem(a.first, b[1]), b[2])
end

function AbstractAlgebra.coefficients(a::PBWAlgElem)
  return OscarPair(parent(a), AbstractAlgebra.coefficients(a.sdata))
end

function Base.length(x::OscarPair{<:PBWAlgRing, <:Singular.SPolyCoeffs})
   return length(x.second)
end

function Base.eltype(x::OscarPair{<:PBWAlgRing{T, S}, <:Singular.SPolyCoeffs}) where {T, S}
   return T
end

function Base.iterate(a::OscarPair{<:PBWAlgRing{T}, <:Singular.SPolyCoeffs}) where T
  b = Base.iterate(a.second)
  b === nothing && return b
  return (coefficient_ring(a.first)(b[1])::T, b[2])
end

function Base.iterate(a::OscarPair{<:PBWAlgRing{T}, <:Singular.SPolyCoeffs}, state) where T
  b = Base.iterate(a.second, state)
  b === nothing && return b
  return (coefficient_ring(a.first)(b[1])::T, b[2])
end

function build_ctx(R::PBWAlgRing)
  return OscarPair(R, MPolyBuildCtx(R.sring))
end

function push_term!(M::OscarPair{<:PBWAlgRing{T,S}, <:MPolyBuildCtx}, c, e::Vector{Int}) where {T, S}
  c = coefficient_ring(M.first)(c)::T
  c = base_ring(M.first.sring)(c)::S
  push_term!(M.second, c, e)
end

function finish(M::OscarPair{<:PBWAlgRing{T,S}, <:MPolyBuildCtx}) where {T, S}
  return PBWAlgElem(M.first, finish(M.second))
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

function var_index(a::PBWAlgElem)
  return Singular.var_index(a.sdata)
end

function is_unit(a::PBWAlgElem)
  return Singular.is_unit(a.sdata)
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

function divexact_left(a::PBWAlgElem, b::PBWAlgElem; check::Bool = true)
  throw(NotImplementedError(:divexact_left, a, b))
end

function divexact_right(a::PBWAlgElem, b::PBWAlgElem; check::Bool = true)
  throw(NotImplementedError(:divexact_right, a, b))
end


####

function AbstractAlgebra.promote_rule(::Type{PBWAlgElem{T, S}}, ::Type{PBWAlgElem{T, S}}) where {T, S}
  return PBWAlgElem{T, S}
end

function AbstractAlgebra.promote_rule(::Type{PBWAlgElem{T, S}}, ::Type{U}) where {T, S, U}
  a = AbstractAlgebra.promote_rule(T, U)
  return a == T ? PBWAlgElem{T, S} : Union{}
end

function (R::PBWAlgRing)()
  return PBWAlgElem(R, R.sring())
end

function (R::PBWAlgRing{T, S})(c::T) where {T, S}
  c = coefficient_ring(R)(c)::T
  c = base_ring(R.sring)(c)::S
  return PBWAlgElem(R, R.sring(c))
end

function (R::PBWAlgRing{T, S})(c::IntegerUnion) where {T, S}
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

function (R::PBWAlgRing)(a::MPolyRingElem)
  @assert parent(a) == R.poly_ring
  z = build_ctx(R)
  for (c, e) in zip(AbstractAlgebra.coefficients(a), AbstractAlgebra.exponent_vectors(a))
    push_term!(z, c, e)
  end
  return finish(z)
end

####

function _unsafe_coerce(R::Union{MPolyRing, Singular.PluralRing}, a::Union{MPolyRingElem, Singular.spluralg}, rev::Bool)
  z = MPolyBuildCtx(R)
  for (c, e) in zip(AbstractAlgebra.coefficients(a), AbstractAlgebra.exponent_vectors(a))
    push_term!(z, base_ring(R)(c), rev ? reverse(e) : e)
  end
  return finish(z)
end

function _unsafe_coerse(R::Singular.PluralRing, I::Singular.sideal, rev::Bool)
  return Singular.Ideal(R, elem_type(R)[_unsafe_coerce(R, a, rev) for a in gens(I)])
end

function is_admissible_ordering(R::PBWAlgRing, o::MonomialOrdering)
  r = base_ring(o)
  n = ngens(R)
  gs = gens(r)
  @assert n == length(gs)
  for i in 1:n-1, j in i+1:n
    t = _unsafe_coerce(r, R.relations[i,j], false)
    if leading_monomial(t; ordering = o) != gs[i]*gs[j]
      return false
    end
  end
  return true
end

function _g_algebra_internal(sr::Singular.PolyRing, rel)
  n = nvars(sr)
  srel = Singular.zero_matrix(sr, n, n)
  C = Singular.zero_matrix(sr, n, n)
  D = Singular.zero_matrix(sr, n, n)
  for i in 1:n-1, j in i+1:n
    t = _unsafe_coerce(sr, rel[i,j], false)
    AbstractAlgebra.leading_monomial(t) == gen(sr, i)*gen(sr, j) ||
                              error("incorrect leading monomial in relations")
    C[i,j] = sr(AbstractAlgebra.leading_coefficient(t))
    D[i,j] = AbstractAlgebra.tail(t)
    srel[i,j] = t
  end
  s, gs = Singular.GAlgebra(sr, C, D)
  return s, gs, srel
end


@doc raw"""
    pbw_algebra(R::MPolyRing{T}, rel, ord::MonomialOrdering; check::Bool = true) where T

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

!!! note
    The `K`-basis condition above is checked by default. This check may be
    skipped by passing `check = false`.

# Examples
```jldoctest
julia> R, (x, y, z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x, y, z) = pbw_algebra(R, REL, deglex(gens(R)))
(PBW-algebra over Rational Field in x, y, z with relations y*x = x*y, z*x = x*z, z*y = y*z + 1, PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, z])
```
"""
function pbw_algebra(r::MPolyRing{T}, rel, ord::MonomialOrdering; check::Bool = true) where T
  n = nvars(r)
  nrows(rel) == n && ncols(rel) == n || error("oops")
  scr = singular_coeff_ring(coefficient_ring(r))
  S = elem_type(scr)
  sr, _ = Singular.polynomial_ring(scr, symbols(r); ordering = singular(ord))
  sr::Singular.PolyRing{S}
  s, gs, srel = _g_algebra_internal(sr, rel)
  if check && !is_zero(Singular.LibNctools.ndcond(s))
    error("PBW-basis condition not satisfied")
  end
  R = PBWAlgRing{T, S}(s, srel, coefficient_ring(r), r)
  return R, [PBWAlgElem(R, x) for x in gs]
end

function pbw_algebra(r::MPolyRing{T}, rel::Vector{Tuple{Int, Int, U}}, ord::MonomialOrdering; check::Bool = true) where {T, U <: MPolyRingElem{T}}
  n = nvars(r)
  gs = gens(r)
  relm = strictly_upper_triangular_matrix([gs[i]*gs[j] for i in 1:n-1 for j in i+1:n])
  for (j, i, p) in rel
    i < j || error("variable indices out of order")
    relm[i, j] = p
  end
  return pbw_algebra(r, relm, ord)
end

function pbw_algebra(r::MPolyRing{T}, rel::Vector{Tuple{U, U, U}}, ord::MonomialOrdering; check::Bool = true) where {T, U <: MPolyRingElem{T}}
  rel2 = Tuple{Int, Int, U}[(var_index(i[1]), var_index(i[2]), i[3]) for i in rel]
  return pbw_algebra(r, rel2, ord)
end

macro pbw_relations(relations...)
  z = Expr(:vect)
  for a in relations
    (a isa Expr) && (a.head == :call) && (length(a.args) == 3) && (a.args[1] == :(==)) ||
        error("bad relation: need ==")
    b = a.args[2]
    (b isa Expr) && (b.head == :call) && (length(b.args) == 3) && (b.args[1] == :*) ||
        error("bad relation: need * on left hand side")
    push!(z.args, :(($(b.args[2]), $(b.args[3]), $(a.args[3]))))
  end
  return esc(z)
end

function weyl_algebra(K::Ring, xs::Vector{Symbol}, dxs::Vector{Symbol})
  n = length(xs)
  n > 0 || error("empty list of variables")
  n == length(dxs) || error("number of differentials should match number of variables")
  r, v = polynomial_ring(K, vcat(xs, dxs))
  rel = elem_type(r)[v[i]*v[j] + (j == i + n) for i in 1:2*n-1 for j in i+1:2*n]
  R,vars = pbw_algebra(r, strictly_upper_triangular_matrix(rel), default_ordering(r); check = false)
  set_attribute!(R, :is_weyl_algebra, :true)  # to activate special printing for Weyl algebras
  return (R,vars)
end

function weyl_algebra(
  K::Ring,
  xs::AbstractVector{<:VarName},
  dxs::AbstractVector{<:VarName}
)
  return weyl_algebra(K, [Symbol(i) for i in xs], [Symbol(i) for i in dxs])
end

@doc raw"""
    weyl_algebra(K::Ring, xs::AbstractVector{<:VarName})

Given a field `K` and a vector `xs` of,  say, $n$ Strings, Symbols, or Characters, return the $n$-th Weyl algebra over `K`.

The generators of the returned algebra print according to the entries of `xs`. See the example below.

# Examples
```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
(Weyl-algebra over Rational Field in variables (x, y), PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, dx, dy])

julia> dx*x
x*dx + 1
```
"""
function weyl_algebra(
  K::Ring,
  xs::AbstractVector{<:VarName}
)
  return weyl_algebra(K, [Symbol(i) for i in xs], [Symbol("d", i) for i in xs])
end

####

function expressify(a::PBWAlgOppositeMap; context = nothing)
  return Expr(:sequence, Expr(:text, "Map to opposite of "),
                         expressify(a.source; context=context))
end

@enable_all_show_via_expressify PBWAlgOppositeMap

function _opposite(a::PBWAlgRing{T, S}) where {T, S}
  if !isdefined(a, :opposite)
    ptr = Singular.libSingular.rOpposite(a.sring.ptr)
    revs = reverse(symbols(a))
    n = length(revs)
    bsring = Singular.PluralRing{S}(ptr, a.sring.base_ring, revs)
    bspolyring, _ = Singular.polynomial_ring(a.sring.base_ring,
                                revs, ordering = ordering(bsring))
    bsrel = Singular.zero_matrix(bspolyring, n, n)
    for i in 1:n-1, j in i+1:n
      bsrel[i,j] = _unsafe_coerce(bspolyring, a.relations[n+1-j,n+1-i], true)
    end
    b = PBWAlgRing{T, S}(bsring, bsrel, a.coeff_ring, polynomial_ring(a.coeff_ring, revs)[1])
    a.opposite = b
    b.opposite = a
  end
  return a.opposite
end

@doc raw"""
    opposite_algebra(A::PBWAlgRing)

Return the opposite algebra of `A`.

# Examples
```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
(Weyl-algebra over Rational Field in variables (x, y), PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, dx, dy])

julia> Dop, opp = opposite_algebra(D);

julia> Dop
PBW-algebra over Rational Field in dy, dx, y, x with relations dx*dy = dy*dx, y*dy = dy*y + 1, x*dy = dy*x, y*dx = dx*y, x*dx = dx*x + 1, x*y = y*x

julia> opp
Map to opposite of Weyl-algebra over Rational Field in variables (x, y)

julia> opp(dx*x)
dx*x + 1
```
"""
function opposite_algebra(a::PBWAlgRing)
  return _opposite(a), PBWAlgOppositeMap(a)
end

function inv(a::PBWAlgOppositeMap)
  return PBWAlgOppositeMap(_opposite(a.source))
end

function _opmap(B::PBWAlgRing{T, S}, a::Singular.spluralg{S}, A::PBWAlgRing{T, S}) where {T, S}
  ptr = GC.@preserve A a B Singular.libSingular.pOppose(A.sring.ptr, a.ptr, B.sring.ptr)
  return B.sring(ptr)
end

function _opmap(B::PBWAlgRing{T, S}, a::Singular.sideal{Singular.spluralg{S}}, A::PBWAlgRing{T, S}) where {T, S}
  ptr = GC.@preserve A a B Singular.libSingular.idOppose(A.sring.ptr, a.ptr, B.sring.ptr)
  return B.sring(ptr)
end

function (M::PBWAlgOppositeMap{T, S})(a::PBWAlgElem{T, S}) where {T, S}
  @assert a.parent === M.source
  opM = _opposite(M.source)
  return PBWAlgElem{T, S}(opM, _opmap(opM, a.sdata, M.source))
end

function Base.broadcasted(M::PBWAlgOppositeMap{T, S}, a::PBWAlgIdeal{D, T, S}) where {D, T, S}
  @assert base_ring(a) === M.source
  opM = _opposite(M.source)
  return PBWAlgIdeal{-D, T, S}(opM, _opmap(opM, a.sdata, M.source))
end

####

function base_ring(a::PBWAlgIdeal)
  return a.basering
end

function ngens(a::PBWAlgIdeal)
  return ngens(a.sdata)
end

function gens(a::PBWAlgIdeal{D, T, S}) where {D, T, S}
  R = base_ring(a)
  return PBWAlgElem{T, S}[PBWAlgElem(R, x) for x in gens(a.sdata)]
end

function gen(a::PBWAlgIdeal, i::Int)
  R = base_ring(a)
  return PBWAlgElem(R, a.sdata[i])
end

getindex(I::PBWAlgIdeal, i::Int) = gen(I, i)

function expressify(a::PBWAlgIdeal{D}; context = nothing) where D
  dir = D < 0 ? :left_ideal : D > 0 ? :right_ideal : :two_sided_ideal
  return Expr(:call, dir, [expressify(g, context = context) for g in gens(a)]...)
end

@enable_all_show_via_expressify PBWAlgIdeal

@doc raw"""
    left_ideal(g::Vector{<:PBWAlgElem})

Given a vector `g` of elements in a PBW-algebra `A`, say, return the left ideal of `A` generated by these elements.

    left_ideal(A::PBWAlgRing, g::AbstractVector)

Given a vector `g` of elements of `A`, return the left ideal of `A` generated by these elements.

# Examples
```jldoctest
julia> R, (x, y, z) = QQ["x", "y", "z"];

julia> L = [x*y, x*z, y*z + 1];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x, y, z) = pbw_algebra(R, REL, deglex(gens(R)))
(PBW-algebra over Rational Field in x, y, z with relations y*x = x*y, z*x = x*z, z*y = y*z + 1, PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, z])

julia> I = left_ideal(A, [x^2*y^2, x*z+y*z])
left_ideal(x^2*y^2, x*z + y*z)
```
"""
function left_ideal(g::Vector{<:PBWAlgElem})
  @assert length(g) > 0
  R = parent(g[1])
  @assert all(x->parent(x) == R, g)
  return left_ideal(R, g)
end

function left_ideal(R::PBWAlgRing{T, S}, g::AbstractVector) where {T, S}
  i = Singular.sideal{Singular.spluralg{S}}(R.sring, [R(x).sdata for x in g], false)
  return PBWAlgIdeal{-1, T, S}(R, i)
end

@doc raw"""
    two_sided_ideal(g::Vector{<:PBWAlgElem})

Given a vector `g` of elements in a PBW-algebra `A`, say, return the two-sided ideal of `A` generated by these elements.

    two_sided_ideal(A::PBWAlgRing, g::AbstractVector)

Given a vector `g` of elements of `A`, return the two-sided ideal of `A` generated by these elements.
"""
function two_sided_ideal(g::Vector{<:PBWAlgElem})
  @assert length(g) > 0
  R = parent(g[1])
  @assert all(x->parent(x) == R, g)
  return two_sided_ideal(R, g)
end

function two_sided_ideal(R::PBWAlgRing{T, S}, g::AbstractVector) where {T, S}
  i = Singular.sideal{Singular.spluralg{S}}(R.sring, [R(x).sdata for x in g], true)
  return PBWAlgIdeal{0, T, S}(R, i)
end

@doc raw"""
    right_ideal(g::Vector{<:PBWAlgElem})

Given a vector `g` of elements in a PBW-algebra `A`, say, return the right ideal of `A` generated by these elements.

    right_ideal(A::PBWAlgRing, g::AbstractVector)

Given a vector `g` of elements of `A`, return the right ideal of `A` generated by these elements.
"""
function right_ideal(g::Vector{<:PBWAlgElem})
  @assert length(g) > 0
  R = parent(g[1])
  @assert all(x->parent(x) == R, g)
  return right_ideal(R, g)
end

function right_ideal(R::PBWAlgRing{T, S}, g::AbstractVector) where {T, S}
  i = Singular.sideal{Singular.spluralg{S}}(R.sring, [R(x).sdata for x in g], true)
  return PBWAlgIdeal{1, T, S}(R, i)
end

# assure a.sopdata is defined
function _sopdata_assure!(a::PBWAlgIdeal)
  if !isdefined(a, :sopdata)
    R = base_ring(a)
    a.sopdata = _opmap(_opposite(R), a.sdata, R)
  end
end


# for D < 0, a.gb is a left gb of left_ideal(a.sdata)
# for D = 0, a.gb is a left gb of two_sided_ideal(a.sdata)
function groebner_assure!(a::PBWAlgIdeal{D}) where D
  @assert D <= 0
  if !isdefined(a, :gb)
    a.gb = Singular.std(a.sdata)
    if D == 0
      a.gb.isTwoSided = false
    end
  end
end

# for D > 0, a.sopdata are gens of the left ideal opposite(right_ideal(a.sdata))
#            a.opgb is a left gb of left_ideal(a.sopdata)
# for D = 0, a.sopdata are gens of the two sided ideal opposite(two_sided_ideal(a.sdata))
#            a.opgb is a left gb of two_sided_ideal(a.sopdata)
function opgroebner_assure!(a::PBWAlgIdeal{D}) where D
  @assert D >= 0
  _sopdata_assure!(a)
  if !isdefined(a, :opgb)
    a.opgb = Singular.std(a.sopdata)
    if D == 0
      a.opgb.isTwoSided = false
    end
  end
end

@doc raw"""
    is_zero(I::PBWAlgIdeal)

Return `true` if `I` is the zero ideal, `false` otherwise.

# Examples
```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
(Weyl-algebra over Rational Field in variables (x, y), PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, dx, dy])

julia> I = left_ideal(D, [x, dx])
left_ideal(x, dx)

julia> is_zero(I)
false
```
"""
function is_zero(a::PBWAlgIdeal)
  return is_zero(a.sdata)
end

function _one_check(I::Singular.sideal)
  for g in gens(I)
    if is_constant(g) && is_unit(AbstractAlgebra.leading_coefficient(g))
      return true
    end
  end
  return false
end

@doc raw"""
    is_one(I::PBWAlgIdeal{D}) where D

Return `true` if `I` is generated by `1`, `false` otherwise.

# Examples
```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
(Weyl-algebra over Rational Field in variables (x, y), PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, dx, dy])

julia> I = left_ideal(D, [x, dx])
left_ideal(x, dx)

julia> is_one(I)
true

julia> J = left_ideal(D, [y*x])
left_ideal(x*y)

julia> is_one(J)
false

julia> K = two_sided_ideal(D, [y*x])
two_sided_ideal(x*y)

julia> is_one(K)
true
```

```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(GF(3), ["x", "y"]);

julia> I = two_sided_ideal(D, [x^3])
two_sided_ideal(x^3)

julia> is_one(I)
false
```
"""
function is_one(a::PBWAlgIdeal{D}) where D
  if is_zero(a.sdata)
    return false
  end
  if _one_check(a.sdata)
    return true
  end
  if D > 0
    opgroebner_assure!(a)
    return _one_check(a.opgb)
  else
    groebner_assure!(a)
    return _one_check(a.gb)
  end
end

@doc raw"""
    +(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}

Return the sum of `I` and `J`.
"""
function Base.:+(a::PBWAlgIdeal{D, T, S}, b::PBWAlgIdeal{D, T, S}) where {D, T, S}
  return PBWAlgIdeal{D, T, S}(base_ring(a), a.sdata + b.sdata)
end


function _as_left_ideal(a::PBWAlgIdeal{D}) where D
  is_left(a) || error("cannot convert to left ideal")
  if D < 0
    return a.sdata
  else
    groebner_assure!(a)
    return a.gb
  end
end

function _as_right_ideal(a::PBWAlgIdeal{D}) where D
  is_right(a) || error("cannot convert to right ideal")
  if D > 0
    return a.sdata
  else
    opgroebner_assure!(a)
    R = base_ring(a)
    return _opmap(R, a.opgb, _opposite(R))
  end
end

@doc raw"""
    *(I::PBWAlgIdeal{DI, T, S}, J::PBWAlgIdeal{DJ, T, S}) where {DI, DJ, T, S}

Given two ideals `I` and `J` such that both `I` and `J` are two-sided ideals
or `I` and `J` are a left and a right ideal, respectively, return the product of `I` and `J`.

# Examples
```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(GF(3), ["x", "y"]);

julia> I = left_ideal(D, [x^3+y^3, x*y^2])
left_ideal(x^3 + y^3, x*y^2)

julia> J = right_ideal(D, [dx^3, dy^5])
right_ideal(dx^3, dy^5)

julia> I*J
two_sided_ideal(x^3*dx^3 + y^3*dx^3, x^3*dy^5 + y^3*dy^5, x*y^2*dx^3, x*y^2*dy^5)
```
"""
function Base.:*(a::PBWAlgIdeal{Da, T, S}, b::PBWAlgIdeal{Db, T, S}) where {Da, Db, T, S}
  @assert base_ring(a) == base_ring(b)
  is_left(a) && is_right(b) || throw(NotImplementedError(:*, a, b))
  # Singular.jl's cartesian product is correct for left*right
  return PBWAlgIdeal{0, T, S}(base_ring(a), _as_left_ideal(a)*_as_right_ideal(b))
end

@doc raw"""
    ^(I::PBWAlgIdeal{D, T, S}, k::Int) where {D, T, S}

Given a two_sided ideal `I`, return the `k`-th power of `I`.

# Examples
```jldoctest
julia> D, (x, dx) = weyl_algebra(GF(3), ["x"]);

julia> I = two_sided_ideal(D, [x^3])
two_sided_ideal(x^3)

julia> I^2
two_sided_ideal(x^6)
```
"""
function Base.:^(a::PBWAlgIdeal{D, T, S}, b::Int) where {D, T, S}
  @assert b >= 0

  if b == 0
    R = base_ring(a)
    return PBWAlgIdeal{D, T, S}(R, Singular.Ideal(R.sring, one(R.sring)))
  elseif b == 1
    return a
  end

  if D == 0
    # Note: repeated mul seems better than nested squaring
    res = a
    while (b -= 1) > 0
      res = res*a
    end
    return res
  else
    throw(NotImplementedError(:^, a, b))
  end
end

@doc raw"""
    intersect(I::PBWAlgIdeal{D, T, S}, Js::PBWAlgIdeal{D, T, S}...) where {D, T, S}
    intersect(V::Vector{PBWAlgIdeal{D, T, S}}) where {D, T, S}

Return the intersection of two or more ideals.

# Examples
```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"]);

julia> I = intersect(left_ideal(D, [x^2, x*dy, dy^2])+left_ideal(D, [dx]), left_ideal(D, [dy^2-x^3+x]))
left_ideal(-x^3 + dy^2 + x)
```
"""
function Base.intersect(a::PBWAlgIdeal{D, T, S}, b::PBWAlgIdeal{D, T, S}...) where {D, T, S}
  R = base_ring(a)
  isempty(b) && return a
  for bi in b
    @assert R === base_ring(bi)
  end
  if D < 0
    res = a.sdata
    res = Singular.intersection(res, [bi.sdata for bi in b]...)
    return PBWAlgIdeal{D, T, S}(R, res)
  elseif D > 0
    _sopdata_assure!(a)
    res = a.sopdata
    for bi in b
      _sopdata_assure!(bi)
    end
    res = Singular.intersection(res, [bi.sopdata for bi in b]...)
    return PBWAlgIdeal{D, T, S}(R, _opmap(R, res, _opposite(R)), res)
  else
    res = _as_left_ideal(a)
    res = Singular.intersection(res, [_as_left_ideal(bi) for bi in b]...)
    return PBWAlgIdeal{D, T, S}(R, res)
  end
end

function Base.intersect(V::Vector{PBWAlgIdeal{D, T, S}}) where {D, T, S}
  @assert length(V) != 0
  length(V) == 1 && return V[1]

  return Base.intersect(V[1], V[2:end]...)
end


@doc raw"""
    ideal_membership(f::PBWAlgElem{T, S}, I::PBWAlgIdeal{D, T, S}) where {D, T, S}

Return `true` if `f` is contained in `I`, `false` otherwise. Alternatively, use `f in I`.

# Examples

```jldoctest
julia> D, (x, dx) = weyl_algebra(QQ, ["x"]);

julia> I = left_ideal(D, [x*dx^4, x^3*dx^2])
left_ideal(x*dx^4, x^3*dx^2)

julia> dx^2 in I
true
```

```jldoctest
julia> D, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"]);

julia> I = two_sided_ideal(D, [x, dx])
two_sided_ideal(x, dx)

julia> one(D) in I
true
```
"""
function ideal_membership(f::PBWAlgElem{T, S}, I::PBWAlgIdeal{D, T, S}) where {D, T, S}
  R = base_ring(I)
  @assert R === parent(f)
  if D <= 0
    # this code works for both D < 0 and D = 0 since:
    #  - groebner_assure! gives a left gb for D = 0 as well (id_TwoStd)
    #  - Singular.reduce is a left normal form
    groebner_assure!(I)
    return Singular.is_zero(Singular.reduce(f.sdata, I.gb))
  else
    opgroebner_assure!(I)
    opf = _opmap(_opposite(R), f.sdata, R)
    return Singular.is_zero(Singular.reduce(opf, I.opgb))
  end
end

function Base.in(f::PBWAlgElem, I::PBWAlgIdeal)
  return ideal_membership(f, I)
end

@doc raw"""
    is_subset(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}

Return `true` if `I` is contained in `J`, `false` otherwise.
# Examples
```jldoctest
julia> D, (x, dx) = weyl_algebra(QQ, ["x"]);

julia> I = left_ideal(D, [dx^2])
left_ideal(dx^2)

julia> J = left_ideal(D, [x*dx^4, x^3*dx^2])
left_ideal(x*dx^4, x^3*dx^2)

julia> is_subset(I, J)
true
```
"""
function is_subset(a::PBWAlgIdeal{D, T, S}, b::PBWAlgIdeal{D, T, S}) where {D, T, S}
    # NB There is an alias is_subset = issubset.
    @assert base_ring(a) === base_ring(b)
    if D <= 0
        # Ditto comment ideal_membership
        groebner_assure!(b)
        return Singular.is_zero(Singular.reduce(a.sdata, b.gb))
    else
        _sopdata_assure!(a)
        opgroebner_assure!(b)
        return Singular.is_zero(Singular.reduce(a.sopdata, b.opgb))
    end
end



@doc raw"""
    ==(I::PBWAlgIdeal{D, T, S}, J::PBWAlgIdeal{D, T, S}) where {D, T, S}

Return `true` if `I` is equal to `J`, `false` otherwise.

# Examples
```jldoctest
julia> D, (x, dx) = weyl_algebra(QQ, ["x"]);

julia> I = left_ideal(D, [dx^2])
left_ideal(dx^2)

julia> J = left_ideal(D, [x*dx^4, x^3*dx^2])
left_ideal(x*dx^4, x^3*dx^2)

julia> I == J
true
```
"""
function Base.:(==)(a::PBWAlgIdeal{D, T, S}, b::PBWAlgIdeal{D, T, S}) where {D, T, S}
  a === b && return true
  gens(a) == gens(b) && return true
  return is_subset(a, b) && is_subset(b, a)
end


#### elimination

function _depends_on_vars(p::Union{Singular.spoly, Singular.spluralg}, sigmaC::Vector{Int})
  for e in AbstractAlgebra.exponent_vectors(p)
    for k in sigmaC
      e[k] == 0 || return true
    end
  end
  return false
end

function is_elimination_subalgebra_admissible(R::PBWAlgRing, sigmaC::Vector{Int})
  n = ngens(R)
  varmap, sigma, sigmaC = Orderings._elimination_data(n, sigmaC)
  for i in sigma, j in sigma
    i < j || continue
    if _depends_on_vars(AbstractAlgebra.tail(R.relations[i,j]), sigmaC)
      return false
    end
  end
  return true
end

function _elimination_ordering_weights(R::PBWAlgRing, sigmaC::Vector{Int})
  n = ngens(R)
  varmap, sigma, sigmaC = Orderings._elimination_data(n, sigmaC)

  # would like
  #  w[i] = 0, for i in sigma
  #  w[i] >= 1, for i in sigma^C
  #  w[i] + w[j] >= e.w, for 1 <= i < j <= n, for e in D[i,j]

  # Build Ax <= b. The set of variables is 1:n
  A = Vector{Vector{Int}}()
  b = Vector{BigInt}()
  for i in sigmaC
    r = zeros(Int, n)
    r[i] = -1
    push!(A, r)
    push!(b, -1)
  end
  for i in 1:n-1, j in i+1:n, e in AbstractAlgebra.exponent_vectors(AbstractAlgebra.tail(R.relations[i,j]))
    e[i] -= 1
    e[j] -= 1
    push!(A, e)
    push!(b, 0)
  end

  # compress variables to sigma^C since w[i] = 0, for i in sigma
  AA = zeros(BigInt, length(A)+1, length(sigmaC))
  for i in 1:length(A), j in sigmaC
    AA[i,varmap[j]] = A[i][j]
  end

  # minimize c.w
  c = ones(Int, length(sigmaC))

  # extra condition c.w <= 2^16 to help polymake
  i = length(A) + 1
  for j in 1:length(sigmaC)
    AA[i,j] = c[j]
  end
  push!(b, 2^16)

  P = polyhedron(AA, b)
  LP = mixed_integer_linear_program(P, c; convention = :min)
  s = optimal_solution(LP)
  if isnothing(s)
    return false, Int[]
  end

  w = Int[varmap[i] > 0 ? ZZ(s[varmap[i]]) : 0 for i in 1:n]
  return true, w
end

# use I's ordering to do the elimination
function _left_eliminate_via_given_ordering(I::Singular.sideal{<:Singular.spluralg}, sigmaC)
  @assert !I.isTwoSided
  J = Singular.std(I)
  @assert J !== I
  for i in 1:ngens(J)
    g = J[i]
    if _depends_on_vars(g, sigmaC)
      J[i] = zero(base_ring(J))
    end
  end
  GC.@preserve J Singular.libSingular.idSkipZeroes(J.ptr)
  return J
end

function _left_eliminate(R::PBWAlgRing, I::Singular.sideal, sigma, sigmaC, ordering)
  r = R.poly_ring
  if !isnothing(ordering)
    @assert is_elimination_ordering(ordering, sigmaC)
    @assert is_admissible_ordering(R, ordering)
    o = singular(ordering)
  else
    # if R's given orderings works, use that
    oo = monomial_ordering(r, Singular.ordering(base_ring(R.relations)))
    if is_elimination_ordering(oo, sigmaC)
      return _left_eliminate_via_given_ordering(I, sigmaC)
    end

    # if degrevlex(sigmaC)*degrevlex(sigma) works, use that
    dpdp = MonomialOrdering(r, Orderings.SymbOrdering(:degrevlex, sigmaC)*
                               Orderings.SymbOrdering(:degrevlex, sigma))
    # This dpdp is fast when sigma and sigmaC are consecutive variables
    # _elimination_data sorts sigma and sigmaC, so consecutive test is easy.
    if length(sigma) == 1 + sigma[end] - sigma[1] &&
       length(sigmaC) == 1 + sigmaC[end] - sigmaC[1] &&
       is_admissible_ordering(R, dpdp)
      o = singular(dpdp)
    else
      # dpdp didn't work, so try to prepend a weight vector to the ordering
      ok, w = _elimination_ordering_weights(R, sigmaC)
      ok || error("could not find elimination ordering")
      o = Singular.ordering_a(w)*Singular.ordering(base_ring(R.relations))
      @assert is_admissible_ordering(R, monomial_ordering(r, o))
    end
  end

  sr, _ = Singular.polynomial_ring(base_ring(R.sring), symbols(r); ordering = o)
  s, gs, _ = _g_algebra_internal(sr, R.relations)
  Io = _unsafe_coerse(s, I, false)
  Io = _left_eliminate_via_given_ordering(Io, sigmaC)
  z = _unsafe_coerse(R.sring, Io, false)
  return z
end

function eliminate(I::PBWAlgIdeal{D, T, S}, sigmaC::Vector{Int}; ordering = nothing) where {D, T, S}
  R = base_ring(I)
  _, sigma, sigmaC = Orderings._elimination_data(ngens(R), sigmaC)

  if is_zero(I)
    return I
  elseif is_empty(sigmaC)
    # eliminating no variables
    return I
  elseif is_empty(sigma)
    # eliminating all variables
    z = is_one(I) ? one(R.sring) : zero(R.sring)
    return PBWAlgIdeal{D, T, S}(R, Singular.Ideal(R.sring, z))
  end

  if !is_elimination_subalgebra_admissible(R, sigmaC)
    error("no elimination is possible: subalgebra is not admissible")
  end

  if D > 0
    Rop = _opposite(R)
    _sopdata_assure!(I)
    n = ngens(R)
    sigmaop = reverse!(n + 1 .- sigma)
    sigmaCop = reverse!(n + 1 .- sigmaC)
    orderingop = isnothing(ordering) ? ordering :
                                     opposite_ordering(Rop.poly_ring, ordering)
    zop = _left_eliminate(Rop, I.sopdata, sigmaop, sigmaCop, orderingop)
    z = _opmap(R, zop, Rop)
    return PBWAlgIdeal{D, T, S}(R, z, zop)
  else
    z = _left_eliminate(R, _as_left_ideal(I), sigma, sigmaC, ordering)
    return PBWAlgIdeal{D, T, S}(R, z)
  end
end

@doc raw"""
    eliminate(I::PBWAlgIdeal, V::Vector{<:PBWAlgElem}; ordering = nothing)

Given a vector `V` of variables, these variables are eliminated from `I`.
That is, return the ideal generated by all polynomials in `I` which only involve the remaining variables.

    eliminate(I::PBWAlgIdeal, V::Vector{Int}; ordering = nothing)

Given a vector `V` of indices which specify variables, these variables are eliminated from `I`.
That is, return the ideal generated by all polynomials in `I` which only involve the remaining variables.


!!! note
    The return value is an ideal of the original algebra.

!!! note
    If provided, the `ordering` must be an admissible elimination ordering (this is checked by the function).   
    If not provided, finding an admissible elimination ordering may involve solving a particular
    linear programming problem. Here, the function is implemented so that
    it searches for solutions in a certain range only. If no solution is found
    in that range, the function will throw an error.

# Examples
```jldoctest
julia> R, (x, y, z, a) = QQ["x", "y", "z", "a"];

julia> L = [x*y-z, x*z+2*x, x*a, y*z-2*y, y*a, z*a];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (x, y, z, a) = pbw_algebra(R, REL, deglex(gens(R)))
(PBW-algebra over Rational Field in x, y, z, a with relations y*x = x*y - z, z*x = x*z + 2*x, a*x = x*a, z*y = y*z - 2*y, a*y = y*a, a*z = z*a, PBWAlgElem{QQFieldElem, Singular.n_Q}[x, y, z, a])

julia> f = 4*x*y+z^2-2*z-a;

julia> I = left_ideal(A, [x^2, y^2, z^2-1, f])
left_ideal(x^2, y^2, z^2 - 1, 4*x*y + z^2 - 2*z - a)

julia> eliminate(I, [x, y, z])
left_ideal(a - 3)

julia> eliminate(I, [1, 2 ,3])
left_ideal(a - 3)

julia> try eliminate(I, [z, a]); catch e; e; end
ErrorException("no elimination is possible: subalgebra is not admissible")
```

```jldoctest
julia> R, (p, q) = QQ["p", "q"];

julia> L = [p*q+q^2];

julia> REL = strictly_upper_triangular_matrix(L);

julia> A, (p, q) = pbw_algebra(R, REL, lex(gens(R)))
(PBW-algebra over Rational Field in p, q with relations q*p = p*q + q^2, PBWAlgElem{QQFieldElem, Singular.n_Q}[p, q])

julia> I = left_ideal(A, [p, q])
left_ideal(p, q)

julia> try eliminate(I, [q]); catch e; e; end   # in fact, no elimination ordering exists
ErrorException("could not find elimination ordering")
```
"""
function eliminate(I::PBWAlgIdeal, sigmaC::Vector{<:PBWAlgElem}; ordering = nothing)
  return eliminate(I, [var_index(i) for i in sigmaC]; ordering = ordering)
end
