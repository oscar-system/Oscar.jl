###############################################################################
# FreeModElem term orderings
###############################################################################

function Orderings.lex(F::ModuleFP)
   return Orderings.ModuleOrdering(F, Orderings.ModOrdering(1:ngens(F), :lex))
end

function Orderings.invlex(F::ModuleFP)
   return Orderings.ModuleOrdering(F, Orderings.ModOrdering(1:ngens(F), :invlex))
end

@doc raw"""
    cmp(ord::ModuleOrdering, a::FreeModElem{T}, b::FreeModElem{T}) where T <: MPolyRingElem

Compare monomials `a` and `b` with regard to the ordering `ord`: Return `-1` for `a < b`
and `1` for `a > b` and `0` for `a == b`. An error is thrown if `ord` is
a partial ordering that does not distinguish `a` from `b`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = FreeMod(R, 3);

julia> cmp(lex(R)*invlex(F), x*F[2], y*F[1])
1

julia> cmp(invlex(F)*lex(R), x*F[2], y*F[1])
-1

julia> cmp(lex(F)*lex(R), x*F[1], x*F[1])
0
```
"""
function Base.cmp(ord::ModuleOrdering, a::FreeModElem{T}, b::FreeModElem{T}) where T <: MPolyRingElem
  A = coordinates(a)
  B = coordinates(b)
  @assert length(A.pos) == 1
  @assert length(B.pos) == 1
  av = A.values[1]
  bv = B.values[1]
  @assert is_monomial(av)
  @assert is_monomial(bv)
  ap = A.pos[1]
  bp = B.pos[1]
  c = Orderings._cmp_vector_monomials(ap, av, 1, bp, bv, 1, ord.o)
  if c == 0 && (ap != bp || av != bv)
    error("$a and $b are incomparable with respect to $ord")
  end
  return c
end

@doc raw"""
    permutation_of_terms(f::FreeModElem{<:MPolyRingElem}, ord::ModuleOrdering)

Return an array of `Tuple{Int, Int}` that puts the terms of `f` in the order
`ord`. The index tuple `(i, j)` corresponds to `term(f[i], j)`.
"""
function Orderings.permutation_of_terms(f::FreeModElem{<:MPolyRingElem}, ord::ModuleOrdering)
  # redispatch to take advantage of the type of ord.o
  return Orderings.__permutation_of_terms(f, ord.o)
end

function Orderings.__permutation_of_terms(f::FreeModElem{<:MPolyRingElem}, ord::Orderings.AbsOrdering)
  ff = coordinates(f)
  p = collect((i, j) for i in ff.pos for j in 1:length(ff[i]))
  sort!(p, lt = (k, l) -> (Orderings._cmp_vector_monomials(k[1], ff[k[1]], k[2],
                                             l[1], ff[l[1]], l[2], ord) > 0))
  return p
end

function Orderings.index_of_leading_term(f::FreeModElem{<:MPolyRingElem}, ord::ModuleOrdering)
  p = Orderings.permutation_of_terms(f, ord)
  @req !isempty(p) "zero element does not have a leading term"
  return p[1]
end


# expressify wrt arbitrary permutation
function expressify(a::OscarPair{<:FreeModElem{<:MPolyRingElem}, Vector{Tuple{Int, Int}}}; context = nothing)
  f = a.first
  x = symbols(base_ring(parent(f)))
  e = generator_symbols(parent(f))
  s = Expr(:call, :+)
  for (i, j) in a.second
    prod = Expr(:call, :*)
    fi = f[i]
    c = coeff(fi, j)
    if !isone(c)
      push!(prod.args, expressify(c, context = context))
    end
    _push_monomial_expr!(prod, x, exponent_vector(fi, j))
    push!(prod.args, e[i])
    push!(s.args, prod)
  end
  return s
end
@enable_all_show_via_expressify OscarPair{<:FreeModElem{<:MPolyRingElem}, <:Vector{Tuple{Int, Int}}}

# expressify wrt ordering
function expressify(a::OscarPair{<:FreeModElem{<:MPolyRingElem}, <:ModuleOrdering}; context = nothing)
  perm = Orderings.permutation_of_terms(a.first, a.second)
  return expressify(OscarPair(a.first, perm); context = context)
end
@enable_all_show_via_expressify OscarPair{<:FreeModElem{<:MPolyRingElem}, <:ModuleOrdering}

@doc raw"""
    coefficients(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return an iterator for the coefficients of `f` with respect to the order `ordering`.
"""
function coefficients(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:coefficients}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:coefficients, <:FreeModElem{<:MPolyRingElem}}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  (i, j) = a.perm[state]
  return coeff(a.elem[i], j), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:coefficients, T}}) where T <: FreeModElem{<:MPolyRingElem{C}} where C
  return C
end

@doc raw"""
    coefficients_and_exponents(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return an iterator whose elements are tuples of coefficients of `f` and exponent
vectors (as `Tuple{Vector{Int}, Int}`) with respect to the order `ordering`.
"""
function coefficients_and_exponents(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:coefficients_and_exponents}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:coefficients_and_exponents, <:FreeModElem{<:MPolyRingElem}}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  (i, j) = a.perm[state]
  return (coeff(a.elem[i], j), (exponent_vector(a.elem[i], j), i)), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:coefficients_and_exponents, T}}) where T <: FreeModElem{<:MPolyRingElem{C}} where C
  return Tuple{C, Tuple{Vector{Int}, Int}}
end

@doc raw"""
    exponents(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return an iterator for the exponents of `f` with respect to the order
`ordering`. Each exponent is a tuple whose first element is a monomial exponent
vectors and whose second element is the index of the module generator.
"""
function exponents(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:exponents}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:exponents, <:FreeModElem{<:MPolyRingElem}}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  (i, j) = a.perm[state]
  return (exponent_vector(a.elem[i], j), i), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:exponents, T}}) where T <: FreeModElem{<:MPolyRingElem}
  return Tuple{Vector{Int}, Int}
end

@doc raw"""
    terms(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return an iterator for the terms of `f` with respect to the order `ordering`.
"""
function terms(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:terms}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:terms, <:FreeModElem{<:MPolyRingElem}}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  (i, j) = a.perm[state]
  return FreeModElem(i, term(a.elem[i], j), parent(a.elem)), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:terms, T}}) where T <: FreeModElem{<:MPolyRingElem}
  return T
end

@doc raw"""
    monomials(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return an iterator for the monomials of `f` with respect to the order `ordering`.
"""
function monomials(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  return GeneralPermutedIterator{:monomials}(f, permutation_of_terms(f, ordering))
end

function Base.iterate(a::GeneralPermutedIterator{:monomials, <:FreeModElem{<:MPolyRingElem}}, state = 0)
  state += 1
  state <= length(a.perm) || return nothing
  (i, j) = a.perm[state]
  return FreeModElem(i, monomial(a.elem[i], j), parent(a.elem)), state
end

function Base.eltype(::Type{<:GeneralPermutedIterator{:monomials, T}}) where T <: FreeModElem{<:MPolyRingElem}
  return T
end

@doc raw"""
    leading_coefficient(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return the leading coefficient and exponent of `f` with respect to the order `ordering`.
"""
function leading_coefficient(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  (i, j) = index_of_leading_term(f, ordering)
  return coeff(f[i], j)
end

@doc raw"""
    leading_coefficient_and_exponent(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return the leading coefficient and exponent of `f` with respect to the order `ordering`.
"""
function leading_coefficient_and_exponent(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  (i, j) = index_of_leading_term(f, ordering)
  return (coeff(f[i], j), (exponent_vector(f[i], j), i))
end

@doc raw"""
    leading_exponent(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return the leading exponent vector (as `Vector{Int}`) of `f` with
respect to the order `ordering`.
"""
function leading_exponent(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  (i, j) = index_of_leading_term(f, ordering)
  return (exponent_vector(f[i], j), i)
end

@doc raw"""
    leading_monomial(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return the leading monomial of `f` with respect to the order `ordering`.
"""
function leading_monomial(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  (i, j) = index_of_leading_term(f, ordering)
  return FreeModElem(i, monomial(f[i], j), parent(f))
end

@doc raw"""
    leading_term(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return the leading term of `f` with respect to the order `ordering`.
"""
function leading_term(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  (i, j) = index_of_leading_term(f, ordering)
  return FreeModElem(i, term(f[i], j), parent(f))
end

@doc raw"""
    tail(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))

Return the tail of `f` with respect to the order `ordering`.
"""
function tail(f::FreeModElem; ordering::ModuleOrdering = default_ordering(parent(f)))
  (i, j) = index_of_leading_term(f, ordering)
  c = collect(k == i ? _delete_index(p, j) : p for (k, p) in coordinates(f))
  return FreeModElem(c, parent(f))
end


