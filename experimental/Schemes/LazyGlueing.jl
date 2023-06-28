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
  compute_glueing_domains::Function
  G::AbsGlueing

  function LazyGlueing(X::AbsSpec, Y::AbsSpec, 
      compute_function::Function, GD::GlueingDataType
    ) where {GlueingDataType}
    return new{typeof(X), typeof(Y), GlueingDataType}(X, Y, GD, compute_function)
  end
  function LazyGlueing(X::AbsSpec, Y::AbsSpec, 
      compute_function::Function, compute_glueing_domains::Function, 
      GD::GlueingDataType
    ) where {GlueingDataType}
    return new{typeof(X), typeof(Y), GlueingDataType}(X, Y, GD, compute_function, 
                                                      compute_glueing_domains)
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

@attr function glueing_domains(G::LazyGlueing) # TODO: Type annotation
  if isdefined(G, :compute_glueing_domains)
    return compute_glueing_domains(G)
  end
  # If no extra function was provided, fall back to computing the full glueing.
  return glueing_domains(underlying_glueing(G))
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
  
