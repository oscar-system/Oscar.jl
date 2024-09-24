patches(G::LazyGluing) = (G.X, G.Y)

# The underlying_gluing triggers the computation of the gluing on request.
#
# For the latter, the user has to implement a method of `_compute_gluing` 
# for their specific type of gluing data. 
function underlying_gluing(G::LazyGluing)
  if !isdefined(G, :G)
    @vprintln :Gluing 5 "computing gluing domains"
    G.G = _compute_gluing(G.GD)
  end
  return G.G
end

function is_computed(G::LazyGluing)
  U, V = patches(G)
  U === V && return true
  return isdefined(G, :G)
end

# We have offered the possibility to compute the gluing domains only, possibly 
# avoiding computations of the gluing morphisms. To this end, the user 
# must implement _gluing_domains for their specific type of gluing data. 
#
# However, we decided to abandon this feature now, since it did not really 
# pay off in performance and it made the things unnecessarily complicated 
# for the serialization of gluings.
gluing_domains(G::LazyGluing) = gluing_domains(underlying_gluing(G))

# Fallback method for the deprecated feature described above. 
# This should not be called, unless the user really has intended this to be used 
# for their particular type of lazy gluing. Therefore we didn't put a fallback 
# which would just be `gluing_domains(_compute_gluing)`. That wouldn't cache the 
# gluing and hence encourages the user to use it in very bad ways. 
_gluing_domains(gd::Any) = error("no shortcut for computation of gluing domains for gluing data of type $(typeof(gd))")

function _compute_gluing(GD::RestrictionDataClosedEmbedding)
  _compute_domains(GD) # Fills in the empty fields for the gluing domains
  return restrict(GD.G, GD.X, GD.Y, GD.UX, GD.VY, check=false)
end

function _compute_domains(GD::RestrictionDataClosedEmbedding)
  if isdefined(GD, :UX) && isdefined(GD, :VY)
    return GD.UX, GD.VY
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

function _compute_gluing(GD::RestrictionDataIsomorphism)
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

