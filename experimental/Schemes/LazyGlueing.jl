@attributes mutable struct LazyGlueing{
                                       LeftSpecType<:AbsSpec, 
                                       RightSpecType<:AbsSpec,
                                       GlueingDataType
                                      }<:AbsGlueing{
                                                    LeftSpecType,
                                                    RightSpecType, 
                                                    Scheme, Scheme,
                                                    Map, Map
                                                   }

  X::LeftSpecType
  Y::RightSpecType
  GD::GlueingDataType
  compute_function::Function
  compute_glueing_domains::Function
  glueing_domains::Union{Tuple{PrincipalOpenSubset, PrincipalOpenSubset},
                         Tuple{SpecOpen, SpecOpen}}
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

function is_computed(G::LazyGlueing)
  return isdefined(G, :G)
end

function glueing_domains(G::LazyGlueing) # TODO: Type annotation
  if !isdefined(G, :glueing_domains)
    # If there was a function provided to do the extra computation, use it.
    if isdefined(G, :compute_glueing_domains)
      G.glueing_domains = G.compute_glueing_domains(G.GD)
    else
      # Otherwise default to computation of the whole glueing.
      G.glueing_domains = glueing_domains(underlying_glueing(G))
    end
  end
  return G.glueing_domains
end


### Preparations for some sample use cases
mutable struct RestrictionDataIsomorphism
  G::AbsGlueing
  i::AbsSpecMor
  j::AbsSpecMor
  i_res::SchemeMor
  j_res::SchemeMor
  function RestrictionDataIsomorphism(G::AbsGlueing, i::AbsSpecMor, j::AbsSpecMor)
    return new(G, i, j)
  end
end

mutable struct RestrictionDataClosedEmbedding
  G::AbsGlueing
  X::AbsSpec
  Y::AbsSpec
  UX::Scheme
  VY::Scheme
  function RestrictionDataClosedEmbedding(G::AbsGlueing, X::AbsSpec, Y::AbsSpec)
    return new(G, X, Y)
  end
end

function _compute_restriction(GD::RestrictionDataClosedEmbedding)
  _compute_domains(GD) # Fills in the empty fields for the glueing domains
  return restrict(GD.G, GD.X, GD.Y, GD.UX, GD.VY, check=false)
end

function _compute_domains(GD::RestrictionDataClosedEmbedding)
  if isdefined(GD, :UX) && isdefined(GD, :VY)
    return GD.UX, GD.UY
  end
  (U, V) = glueing_domains(GD.G)
  if U isa PrincipalOpenSubset
    GD.UX = PrincipalOpenSubset(GD.X, OO(GD.X)(lifted_numerator(complement_equation(U))))
    GD.VY = PrincipalOpenSubset(GD.Y, OO(GD.Y)(lifted_numerator(complement_equation(V))))
    return GD.UX, GD.VY
  elseif U isa SpecOpen
    GD.UX = intersect(GD.X, U, check=false)
    GD.VY = intersect(GD.Y, V, check=false)
    return GD.UX, GD.VY
  end
end

function _compute_restriction(GD::RestrictionDataIsomorphism)
  _compute_domains(GD) # Fills in the empty fields for the glueing domains
  return restrict(GD.G, GD.i, GD.j, GD.i_res, GD.j_res, check=false)
end

function _compute_domains(GD::RestrictionDataIsomorphism)
  if !(isdefined(GD, :i_res) && isdefined(GD, :j_res))
    U, V = glueing_domains(GD.G)

    GD.i_res = restrict(GD.i, U, preimage(inverse(GD.i), U, check=false), check=false)
    i_res_inv = restrict(inverse(GD.i), codomain(GD.i_res), U, check=false)
    set_attribute!(GD.i_res, :inverse, i_res_inv)
    set_attribute!(i_res_inv, :inverse, GD.i_res)

    GD.j_res = restrict(GD.j, V, preimage(inverse(GD.j), V, check=false), check=false)
    j_res_inv = restrict(inverse(GD.j), codomain(GD.j_res), V, check=false)
    set_attribute!(GD.j_res, :inverse, j_res_inv)
    set_attribute!(j_res_inv, :inverse, GD.j_res)
  end
  return codomain(GD.i_res), codomain(GD.j_res)
end

# Method for compatibility with internal methods of the glueings
  
Glueing(g::LazyGlueing) = Glueing(underlying_glueing(g))

