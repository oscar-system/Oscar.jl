
########################################################################
# Properties of AbsCoveredSchemes                                      #
#    by referral to the patches                                        #
########################################################################

########################################################################
# (1) Check and store, whether a covered scheme is empty               #
########################################################################
@doc raw"""
    is_empty(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is empty.

"""
is_empty(X::AbsCoveredScheme) = is_empty(underlying_scheme(X))
@attr Bool function is_empty(X::CoveredScheme)
  if !isdefined(X, :coverings)
    return true
  end
  return all(isempty, all_patches(default_covering(X)))
end

########################################################################
# (2) Check and store, whether a covered scheme is smooth              #
########################################################################
@doc raw"""
    is_smooth(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is smooth.
"""
is_smooth(X::AbsCoveredScheme) = is_smooth(underlying_scheme(X))

@attr Bool function is_smooth(X::CoveredScheme)
  if !isdefined(X, :coverings)
    return true
  end
  # TODO: Can this be optimized? We only need to check the rank of
  # the Jacobian matrices where we haven't checked in another chart before.
  # This is a tradeoff between cost of Jacobian criterion and cost of
  # optimizing the use of the covering
  return all(is_smooth, affine_charts(X))
end

function _jacobian_criterion(X::CoveredScheme{<:Field})
  if !isdefined(X, :coverings)
    return true
  end

  if !has_decomposition_info(default_covering(X))
    throw(NotImplementedError(:_jacobian_criterion, "only implemented when decomposition info is available"))
  end

  dec_info = decomposition_info(default_covering(X))
  for (V, fs) in dec_info
    R_quotient = OO(V)
    R = R_quotient isa MPolyRing ? R_quotient : base_ring(R_quotient)
    I = saturated_ideal(defining_ideal(V))
    mat = jacobian_matrix(R, gens(I))
    sing_locus = ideal(R, fs) + ideal(R, minors(mat, codim(V)))
    sing_subscheme = subscheme(V, sing_locus)
    if !isempty(sing_subscheme)
      return false
    end
  end
  return true
end

########################################################################
# (3) Check and store, whether a covered scheme is integral            #
#     i.e. irreducible and reduced                                     #
########################################################################
@doc raw"""
    is_integral(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is integral.

"""
@attr Bool function is_integral(X::AbsCoveredScheme)
  !is_empty(X) || return false
  return is_reduced(X) && is_irreducible(X)
end

# auxiliary function for connectedness of gluing graph
#      do not confuse with connectedness of the scheme
# Note: This does not work with gluing_graph, because empty patches
#      need to be ignored without changing the covering
@doc raw"""
    is_connected_gluing(X::AbsCoveredScheme)

Return the boolean value whether the gluing graph of the default
covering of the scheme X is connected.

!!! note
    This function is designed to ignore empty patches, which may arise e.g. upon creation of subschemes of covered schemes,

"""
@attr Bool function is_connected_gluing(X::AbsCoveredScheme)
  return is_connected(pruned_gluing_graph(default_covering(X)))
end

########################################################################
# (4) Check and store, whether a covered scheme is connected           #
########################################################################
@doc raw"""
    is_connected(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is connected.

"""
@attr Bool function is_connected(X::AbsCoveredScheme)
  is_connected_gluing(X) || return false
  # note for future implementation: expensive property
  # 1) do primary decomposition
  # 2) check connectedness of lowest two layers of the intersection lattice
  error("not implemented yet")
end

##############################################################################
# (5) Check and store, whether a scheme is reduced
##############################################################################
@doc raw"""
    is_reduced(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is reduced.

"""
@attr Bool function is_reduced(X::AbsCoveredScheme)
  return all(is_reduced, affine_charts(X))
end

##############################################################################
# (6) Check and store, whether a covered scheme is irreducible
##############################################################################
@doc raw"""
    is_irreducible(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is irreducible.

"""
@attr Bool function is_irreducible(X::AbsCoveredScheme)
  is_connected_gluing(X) || return false
  !is_empty(X) || return false
  v=findall(!is_empty, affine_charts(X))  ## only check non-empty patches
  return all(is_irreducible(affine_charts(X)[v[i]]) for i in 1:length(v) )
end
