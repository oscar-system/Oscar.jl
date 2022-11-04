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


