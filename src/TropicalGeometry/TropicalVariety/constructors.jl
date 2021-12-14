######################
# 1: The Julia type for TropicalVarieties
#    To Do (brainstorming):
#    - additional member objects for specialized structs (?)
#      before adding, check whether it is not already included in the polymake object
#      (e.g. a tropical hypersurface in polymake knows its tropical polynomial)
#      and check whether it might not be better as its own struct
#      (e.g. integral structures on abstract tropical varieties as desired by Janko)
#    - tropical linear spaces need not be realizable, but what about the rest?
######################


# We use M to record whether we are in the min or max case
# M is either typeof(min) or typeof(max)
# We use EMB to record whether the variety is embedded or abstract
# EMB is either true or false
abstract type TropicalVarietySupertype{M,EMB} end


# abstract tropical varieties are hypergraphs,
#   i.e., polyhedral complexes whose vertices are merely enumerated
#         and have no coordinates assigned to them
# embedded tropical varieties are polyhedral complexes in euclidean space
# min or max influences:
# - initial ideals for embedded tropical varieties
struct TropicalVariety{M,EMB} <: TropicalVarietySupertype{M,EMB}
    # GapTV::GapObj
    polymakeTV::Polymake.BigObject
    algebraicTV
end
export TropicalVariety
function pm_object(v::TropicalVariety)
  return v.polymakeTV
end


# abstract tropical hypersurfaces do not make any concrete sense to me,
#   so I don't think we should actively support it for now (Yue)
# embedded tropical varieties are one-codimensional polyhedral complexes in euclidean space
# min or max should influence:
# - tropical polynomials and their initial forms
struct TropicalHypersurface{M,EMB} <: TropicalVarietySupertype{M,EMB}
    # GapTV::GapObj
    polymakeTV::Polymake.BigObject
    tropicalPolynomial
    tropicalInitials
    algebraicPolynomial
    algebraicInitials
end
export TropicalHypersurface
function pm_object(v::TropicalHypersurface)
  return v.polymakeTV
end


# abstract tropical curves are abstracts graphs
# embedded tropical curves are one-dimensional polyhedral complexes
# min or max should influence:
# - initial ideals for embedded tropical curves (if realizable)
struct TropicalCurve{M,T} <: TropicalVarietySupertype{M,T}
    # GapTV::GapObj
    polymakeTV::Polymake.BigObject
    algebraicTV
end
export TropicalCurve
function pm_object(v::TropicalCurve)
  return v.polymakeTV
end


# abstract tropical linear spaces are valuated matroids (?)
# embedded tropical linear spaces are degree 1 embedded tropical varieties
# min or max should influence:
# - Pluecker vector
# - initial ideals for embedded tropical linear spaces (if realizable)
struct TropicalLinearSpace{M,T} <: TropicalVarietySupertype{M,T}
    # GapTV::GapObj
    polymakeTV::Polymake.BigObject
    algebraicTV
end
export TropicalLinearSpace


function pm_object(v::TropicalVariety)
  return v.polymakeTV
end


######################
# 2: Generic constructors
######################


######################
# 2.1: Generic constructors for tropical varieties
# - embedded tropical varieties of polynomial ideals
# - abstract tropical varieties from vertices and maximal cells
######################
@doc Markdown.doc"""
    TropicalVariety()

Construct the embedded tropical variety of a polynomial ideal over a (possibly trivially) valued field

# Examples
"""
function TropicalVariety()
    # ...
    return #...
end


@doc Markdown.doc"""
    TropicalVariety()

Construct the abstract tropical variety from a list of vertices and maximal cells

# Examples
"""
function TropicalVariety()
    # ...
    return #...
end



######################
# 2.2: Generic constructors for tropical hypersurfaces
# - embedded tropical hypersurfaces from polynomials
# - embedded tropical hypersurfaces from tropical polynomials
######################
@doc Markdown.doc"""
    TropicalHypersurface()

Construct the tropical hypersurface of a polynomial over a valued field

# Examples
"""
function TropicalHypersurface()
  return #...
end


@doc Markdown.doc"""
    TropicalHypersurface()

Construct the tropical hypersurface of a polynomial over the tropical numbers

# Examples
```jldoctest
julia> T = tropical_ring(min)
Tropical ring (min)

julia> Txy,(x,y) = T["x","y"]
(Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalRingElem{typeof(min)}}[x, y])

julia> f = x+y+1
x + y + (1)

julia> hyp = TropicalHypersurface(f)
TropicalHypersurface{min, true}(Polymake.BigObjectAllocated(Ptr{Nothing} @0x000000001adc7d70), Any[])
```
"""
function TropicalHypersurface(f)
    if total_degree(f) <= 0
        error("Tropical variety of constant polynomials not supported.")
    end
    convention = fun(base_ring(f))
    fstr = Tuple(tropical_polynomial_to_polymake(f))
    pmpoly = Polymake.common.totropicalpolynomial(fstr...)
    pmhyp = Polymake.tropical.Hypersurface{convention}(POLYNOMIAL=pmpoly)
    return TropicalHypersurface{convention, true}(pmhyp, [])
end




######################
# 2.3: Generic constructors for tropical curves
# - embedded tropical curves from one-dimensional polynomial ideals
# - abstract tropical curves from vertices, maximal cells, edge lengths, etc.
######################
@doc Markdown.doc"""
    TropicalCurve()

Construct an abstract tropical curve from a list of edges and a vector of their lengths

# Examples
"""
function TropicalCurve()
  return #...
end


@doc Markdown.doc"""
    TropicalCurve()

Construct an abstract tropical curve from a list of edges and a vector of their lengths

# Examples
"""
function TropicalCurve()
  return #...
end



######################
# 2.4: Generic constructors for tropical linear spaces
# - embedded tropical linear spaces from degree 1 polynomial ideals
# - embedded tropical linear spaces from Pluecker vectors
######################
@doc Markdown.doc"""
    TropicalLinearSpace()

Construct a tropical linear space from a degree 1 polynomial ideal

# Examples
"""
function TropicalLinearSpace()
  return #...
end


@doc Markdown.doc"""
    TropicalLinearSpace()

Construct a tropical linear space from its Pluecker vector

# Examples
"""
function TropicalLinearSpace()
  return #...
end


######################
# 3: Special constructors
######################



######################
# 4: Properties
# TropicalVarietySupertype
# - dimension
# - degree (top-dimensional only?)
# - vertices
# - polyhedra
# - IsPure
# - IsConnected
# - IsIrreducible / IsReducible # I think this is somewhat non-trivial to decide
# TropicalVarieties
# - InitialIdeal (if constructed from ideal)
# TropicalHypersurface
# - NewtonSubdivision
# - TropicalPolynomial
# - Polynomial (if constructed from polynomial over valued field) (?)
# TropicalCurve
# - should we support marked points, piece-wise linear functions, divisors, chip-firing?
# TropicalLinearSpace
# - PlueckerVector
######################


# @doc Markdown.doc"""
#     isnormal( v::AbstractNormalToricVariety )

# Checks if the normal toric variety `v` is normal. (This function is somewhat tautological at this point.)

# # Examples
# ```julia-repl
# julia> isnormal(projective_space( 2 ))
# true
# ```
# """
# function isnormal( v::AbstractNormalToricVariety )
#   return true
# end
# export isnormal




###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
# function Base.show(io::IO, ntv::AbstractNormalToricVariety)
#   # fan = get_polyhedral_fan(ntv)
#   pmntv = pm_ntv(ntv)
#   ambdim = pmntv.FAN_AMBIENT_DIM
#   print(io, "A normal toric variety corresponding to a polyhedral fan in ambient dimension $(ambdim)")
# end
