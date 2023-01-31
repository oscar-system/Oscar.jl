export is_empty

########################################################################
# Properties of AbsCoveredSchemes                                      #
########################################################################
@attr function is_empty(X::AbsCoveredScheme)
  if !isdefined(X, :coverings) 
    return true
  end
  return all(x->isempty(x), all_patches(default_covering(X)))
end

@attr function is_smooth(X::AbsCoveredScheme)
  if !isdefined(X, :coverings) 
    return true
  end
  # TODO: This can be optimized! We only need to check the rank of 
  # the jacobian matrices where we haven't checked in another chart before. 
  return all(x->is_smooth(x), affine_charts(X))
end

@attr function is_integral(X::AbsCoveredScheme)
  return is_connected(X) && all(U->(is_integral(U)), affine_charts(X))
end

@attr function is_connected(X::AbsCoveredScheme)
  return is_connected(glueing_graph(default_covering(X)))
end

