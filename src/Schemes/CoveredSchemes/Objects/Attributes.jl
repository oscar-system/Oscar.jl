export base_ring_type, coverings, default_covering, patches, glueings, name_of

########################################################################
# Attributes of AbsCoveredScheme                                       #
########################################################################

########################################################################
# Type getters                                                         #
########################################################################
base_ring_type(::Type{T}) where {BRT, T<:AbsCoveredScheme{BRT}} = BRT
base_ring_type(X::AbsCoveredScheme) = base_ring_type(typeof(X))

########################################################################
# Basic getters                                                        #
########################################################################
base_ring(X::AbsCoveredScheme) = base_ring(underlying_scheme(X))

@Markdown.doc """
    coverings(X::AbsCoveredScheme)

Returns the available coverings for ``X``.
"""
function coverings(X::AbsCoveredScheme) ::Vector{<:Covering}
  return coverings(underlying_scheme(X))
end

@Markdown.doc """
    default_covering(X::AbsCoveredScheme)::Covering

Returns the default covering for ``X``.
"""
function default_covering(X::AbsCoveredScheme)
  return default_covering(underlying_scheme(X))::Covering
end

@Markdown.doc """
    patches(X::AbsCoveredScheme) = patches(default_covering(X))

Returns the affine patches in the `default_covering` of ``X``.
"""
patches(X::AbsCoveredScheme) = patches(default_covering(X))


########################################################################
# Attributes of CoveredScheme                                          #
########################################################################

########################################################################
# Basic getters                                                        #
########################################################################
coverings(X::CoveredScheme) = X.coverings
default_covering(X::CoveredScheme) = X.default_covering
patches(X::CoveredScheme) = patches(default_covering(X))
glueings(X::CoveredScheme) = glueings(default_covering(X))
base_ring(X::CoveredScheme) = X.kk

########################################################################
# Names of CoveredSchemes                                              #
########################################################################
set_name!(X::AbsCoveredScheme, name::String) = set_attribute!(X, :name, name)
name_of(X::AbsCoveredScheme) = get_attribute(X, :name)::String
has_name(X::AbsCoveredScheme) = has_attribute(X, :name)

########################################################################
# Auxiliary attributes                                                 #
########################################################################
function dim(X::AbsCoveredScheme) 
  if !has_attribute(X, :dim)
    d = -1
    is_equidimensional=true
    for U in patches(default_covering(X))
      e = dim(U)
      if e > d
        d == -1 || (is_equidimensional=false)
        d = e
      end
    end
    set_attribute!(X, :dim, d)
    if !is_equidimensional
      # the above is not an honest check for equidimensionality,
      # because in each chart the output of `dim` is only the 
      # supremum of all components. Thus we can only infer 
      # non-equidimensionality in case this is already visible
      # from comparing the diffent charts
      set_attribute(X, :is_equidimensional, false)
    end
  end
  return get_attribute(X, :dim)::Int
end

