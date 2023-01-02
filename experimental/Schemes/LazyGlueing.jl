@attributes mutable struct LazyGlueing{
                                       LeftSpecType<:AbsSpec, 
                                       RightSpecType<:AbsSpec,
                                       GlueingDataType
                                      }<:AbsGlueing{
                                                    LeftSpecType,
                                                    RightSpecType, 
                                                    Scheme, Scheme,
                                                    Hecke.Map, Hecke.Map
                                                   }

  X::LeftSpecType
  Y::RightSpecType
  GD::GlueingDataType
  compute_function::Function
  G::AbsGlueing

  function LazyGlueing(X::AbsSpec, Y::AbsSpec, 
      compute_function::Function, GD::GlueingDataType
    ) where {GlueingDataType}
    return new{typeof(X), typeof(Y), GlueingDataType}(X, Y, GD, compute_function)
  end
end

patches(G::LazyGlueing) = (G.X, G.Y)

# The underlying_glueing triggers the computation of the glueing on request.
function underlying_glueing(G::LazyGlueing)
  if !isdefined(G, :G) 
    G.G = G.compute_function(G.GD)
  end
  return G.G
end


### Preparations for some sample use cases
struct RestrictionDataIsomorphism
  G::AbsGlueing
  i::AbsSpecMor
  j::AbsSpecMor
end

struct RestrictionDataClosedEmbedding
  G::AbsGlueing
  X::AbsSpec
  Y::AbsSpec
end

function _compute_restriction(GD::RestrictionDataClosedEmbedding)
  return restrict(GD.G, GD.X, GD.Y, check=false)
end

function _compute_restriction(GD::RestrictionDataIsomorphism)
  return restrict(GD.G, GD.i, GD.j, check=false)
end

# Method for compatibility with internal methods of the glueings
  
