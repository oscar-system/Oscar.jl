###############################################################################
# A place to accumulate code that should eventually be moved to Hecke.jl
###############################################################################

function getindex(r::Hecke.SRow, u::AbstractUnitRange)
  s = sparse_row(base_ring(r))
  shift = 1-first(u)
  for (p,v) = r
    if p in u
      push!(s.pos, p+shift)
      push!(s.values, v)
    end
  end
  return s
end

canonical_unit(x::AbsSimpleNumFieldOrderQuoRingElem) = one(parent(x))

function numerator(f::QQPolyRingElem, parent::ZZPolyRing = Hecke.Globals.Zx)
  g = parent()
  ccall((:fmpq_poly_get_numerator, Nemo.libflint), Cvoid, (Ref{ZZPolyRingElem}, Ref{QQPolyRingElem}), g, f)
  return g
end

Hecke.minpoly(a::QQBarFieldElem) = minpoly(Hecke.Globals.Qx, a)

function extend_domain_to_fraction_field(phi::Map{<:MPolyRing, <:Ring})
  ext_dom = fraction_field(domain(phi))
  return MapFromFunc(ext_dom, codomain(phi), x->phi(numerator(x))*inv(phi(denominator(x))))
end


function iszero(P::EllipticCurvePoint)
  return iszero(P[1]) && isone(P[2]) && iszero(P[3])
end

@doc raw"""
    extended_ade(ADE::Symbol, n::Int)

Return the dual intersection matrix of an extended ade Dynkin diagram
as well as the isotropic vector (with positive coefficients in the roots).
"""
function extended_ade(ADE::Symbol, n::Int)
  R = change_base_ring(ZZ,gram_matrix(root_lattice(ADE,n)))
  G = block_diagonal_matrix([ZZ[2;],R])
  if ADE == :E && n == 8
    G[1,n] = -1
    G[n,1] = -1
  end
  if ADE == :E && n == 7
    G[1,2] = -1
    G[2,1] = -1
  end
  if ADE == :E && n == 6
    G[1,n+1] = -1
    G[n+1,1] = -1
  end
  if ADE == :A && n > 0
    G[1,2] = -1
    G[2,1] = -1
    G[1,n+1] = -1
    G[n+1,1] = -1
  end
  if ADE == :A && n ==1 0
    G[1,2]= -2
    G[2,1] = -2
  end
  if ADE == :D
    G[1,n] = -1
    G[n,1] = -1
  end
  @assert rank(G) == n
  return -G, kernel(G; side = :left)
end

@doc raw"""
    disc_log(b::T, x::T) where {T <: FinFieldElem}

Return an integer `s` such that $b^s = x$.
If no such `x` exists, an exception is thrown.

# Examples
```jldoctest
julia> F = GF(3,4); a = gen(F)^21;

julia> disc_log(gen(F), a)
21
```
"""
function disc_log(b::T, x::T) where {T <: FinFieldElem}
  @assert parent(b) === parent(x)
  return Hecke.disc_log_bs_gs(b, x, order(parent(b)))
end

###############################################################################
# Part of https://github.com/thofma/Hecke.jl/pull/1800, but breaking for OSCAR 1.3.1.
# To allow progress, we include it here until a future breaking Hecke release.
#
function (::Type{T})(G::FinGenAbGroup) where T <: Group
  return codomain(isomorphism(T, G))
end

function (::Type{FinGenAbGroup})(G::Group)
  return codomain(isomorphism(FinGenAbGroup, G))
end

function isomorphism(::Type{FinGenAbGroup}, G::FinGenAbGroup; on_gens::Bool=false)
  # Known isomorphisms are cached in the attribute `:isomorphisms`.
  on_gens = true # we ignore the on_gens flag, the identity will *always* map gens onto gens
  isos = get_attribute!(Dict{Tuple{Type, Bool}, Any}, G, :isomorphisms)::Dict{Tuple{Type, Bool}, Any}
  return get!(isos, (FinGenAbGroup, on_gens)) do
    return id_hom(G)
  end::FinGenAbGroupHom
end

function isomorphism(::Type{T}, G::FinGenAbGroup; on_gens::Bool=false) where T <: Group
  throw(NotImplementedError(:isomorphism, T, G))
end

function isomorphism(::Type{FinGenAbGroup}, G::Group; on_gens::Bool=false)
  throw(NotImplementedError(:isomorphism, FinGenAbGroup, G))
end
#
###############################################################################
