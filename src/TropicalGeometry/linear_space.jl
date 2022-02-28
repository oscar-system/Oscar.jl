###
# Tropical linear spaces in Oscar
# ===============================
###



###
# 1. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention, affecting initial ideals, Pluecker vector, etc.
# EMB = true or false:
#   embedded or abstract tropical linear space
#   embedded tropical linear space = weighted polyhedral complex in euclidean space
#   abstract tropical linear space = weighted hypergraph with enumerated vertices
###

@attributes mutable struct TropicalLinearSpace{M,T} <: TropicalVarietySupertype{M,T}
    # GapTV::GapObj
    polymakeTV::Polymake.BigObject
end
export TropicalLinearSpace



###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalLinearSpace(ideal::MPolyIdeal{fmpq_poly})

Construct a tropical linear space from a degree 1 polynomial ideal

# Examples
"""
function TropicalLinearSpace(ideal::MPolyIdeal{fmpq_poly})
    error("TODO: Not implemented yet.")
    return #...
end


@doc Markdown.doc"""
    TropicalLinearSpace(plv::Vector)

Construct a tropical linear space from its Pluecker vector

# Examples
"""
function TropicalLinearSpace(plv::Vector)
    error("TODO: Not implemented yet.")
    return #...
end

#needs Oscar type as entry 
function TropicalLinearSpace(tropicalmatrix::MatElem)
  plv = Nemo.minors(tropicalmatrix, min(size(tropicalmatrix)[1], size(tropicalmatrix)[2])) 
  rk = rank(tropicalmatrix)
  nelement = size(tropicalmatrix)[2]
  return (plv,nelement, rk)
  ### TropicalLinearSpace(plv, nelement, rk)
end

function TropicalLinearSpace(tropicalmatrix::Matrix{Int})
  return TropicalLinearSpace(matrix(ZZ, tropicalmatrix))
end


function TropicalLinearSpace(tropicalmatrix::Matrix{fmpq})
  return TropicalLinearSpace(matrix(base_ring(tropicalmatrix), tropicalmatrix))  
end


function TropicalLinearSpace(tropicalmatrix::Matrix{fmpz})
  return TropicalLinearSpace(matrix(base_ring(tropicalmatrix), tropicalmatrix))  
end

function TropicalLinearSpace(tropicalmatrix::Matrix{Rational})
  return TropicalLinearSpace(matrix(QQ, tropicalmatrix))  
end


###
# 3. Basic properties
# -------------------
###
