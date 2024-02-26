@attributes mutable struct LazyGluing{
                                       LeftAffineSchemeType<:AbsAffineScheme, 
                                       RightAffineSchemeType<:AbsAffineScheme,
                                       GluingDataType
                                      }<:AbsGluing{
                                                    LeftAffineSchemeType,
                                                    RightAffineSchemeType, 
                                                    Scheme, Scheme,
                                                    Map, Map
                                                   }

  X::LeftAffineSchemeType
  Y::RightAffineSchemeType
  GD::GluingDataType
  compute_function::Function
  compute_gluing_domains::Function
  gluing_domains::Union{Tuple{PrincipalOpenSubset, PrincipalOpenSubset},
                         Tuple{AffineSchemeOpenSubscheme, AffineSchemeOpenSubscheme}}
  G::AbsGluing

  function LazyGluing(X::AbsAffineScheme, Y::AbsAffineScheme, 
      compute_function::Function, GD::GluingDataType
    ) where {GluingDataType}
    return new{typeof(X), typeof(Y), GluingDataType}(X, Y, GD, compute_function)
  end
  function LazyGluing(X::AbsAffineScheme, Y::AbsAffineScheme, 
      compute_function::Function, compute_gluing_domains::Function, 
      GD::GluingDataType
    ) where {GluingDataType}
    return new{typeof(X), typeof(Y), GluingDataType}(X, Y, GD, compute_function, 
                                                      compute_gluing_domains)
  end
end

patches(G::LazyGluing) = (G.X, G.Y)

# The underlying_gluing triggers the computation of the gluing on request.
function underlying_gluing(G::LazyGluing)
  if !isdefined(G, :G) 
    G.G = G.compute_function(G.GD)
  end
  return G.G
end

function is_computed(G::LazyGluing)
  return isdefined(G, :G)
end

function gluing_domains(G::LazyGluing) # TODO: Type annotation
  if !isdefined(G, :gluing_domains)
    # If there was a function provided to do the extra computation, use it.
    if isdefined(G, :compute_gluing_domains)
      G.gluing_domains = G.compute_gluing_domains(G.GD)
    else
      # Otherwise default to computation of the whole gluing.
      G.gluing_domains = gluing_domains(underlying_gluing(G))
    end
  end
  return G.gluing_domains
end


### Preparations for some sample use cases
mutable struct RestrictionDataIsomorphism
  G::AbsGluing
  i::AbsAffineSchemeMor
  j::AbsAffineSchemeMor
  i_res::SchemeMor
  j_res::SchemeMor
  function RestrictionDataIsomorphism(G::AbsGluing, i::AbsAffineSchemeMor, j::AbsAffineSchemeMor)
    return new(G, i, j)
  end
end

mutable struct RestrictionDataClosedEmbedding
  G::AbsGluing
  X::AbsAffineScheme
  Y::AbsAffineScheme
  UX::Scheme
  VY::Scheme
  function RestrictionDataClosedEmbedding(G::AbsGluing, X::AbsAffineScheme, Y::AbsAffineScheme)
    return new(G, X, Y)
  end
end

function _compute_restriction(GD::RestrictionDataClosedEmbedding)
  _compute_domains(GD) # Fills in the empty fields for the gluing domains
  return restrict(GD.G, GD.X, GD.Y, GD.UX, GD.VY, check=false)
end

function _compute_domains(GD::RestrictionDataClosedEmbedding)
  if isdefined(GD, :UX) && isdefined(GD, :VY)
    return GD.UX, GD.UY
  end
  (U, V) = gluing_domains(GD.G)
  if U isa PrincipalOpenSubset
    GD.UX = PrincipalOpenSubset(GD.X, OO(GD.X)(lifted_numerator(complement_equation(U))))
    GD.VY = PrincipalOpenSubset(GD.Y, OO(GD.Y)(lifted_numerator(complement_equation(V))))
    return GD.UX, GD.VY
  elseif U isa AffineSchemeOpenSubscheme
    GD.UX = intersect(GD.X, U, check=false)
    GD.VY = intersect(GD.Y, V, check=false)
    return GD.UX, GD.VY
  end
end

function _compute_restriction(GD::RestrictionDataIsomorphism)
  _compute_domains(GD) # Fills in the empty fields for the gluing domains
  return restrict(GD.G, GD.i, GD.j, GD.i_res, GD.j_res, check=false)
end

function _compute_domains(GD::RestrictionDataIsomorphism)
  if !(isdefined(GD, :i_res) && isdefined(GD, :j_res))
    U, V = gluing_domains(GD.G)

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

# Method for compatibility with internal methods of the gluings
  
Gluing(g::LazyGluing) = Gluing(underlying_gluing(g))

