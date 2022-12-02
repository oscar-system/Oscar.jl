export base_ring_type, coverings, default_covering, patches, glueings,  affine_charts, name

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

Return the list of internally stored `Covering`s of ``X``.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> Pcov = covered_scheme(P);

julia> coverings(Pcov)
1-element Vector{Covering{FlintRationalField}}:
 Covering with 3 patches
```
"""
function coverings(X::AbsCoveredScheme) ::Vector{<:Covering}
  return coverings(underlying_scheme(X))
end

@Markdown.doc """
    default_covering(X::AbsCoveredScheme)

Return the default covering for ``X``.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> S = ambient_coordinate_ring(P);

julia> I = ideal(S, [S[1]*S[2]-S[3]^2]);

julia> X = subscheme(P, I);

julia> Xcov = covered_scheme(X)
covered scheme with 3 affine patches in its default covering

julia> default_covering(Xcov)
Covering with 3 patches

```
"""
function default_covering(X::AbsCoveredScheme)
  return default_covering(underlying_scheme(X))::Covering
end

@Markdown.doc """
    patches(X::AbsCoveredScheme) = patches(default_covering(X))

Return the affine patches in the `default_covering` of ``X``.
"""
patches(X::AbsCoveredScheme) = patches(default_covering(X))

@Markdown.doc """
    affine_charts(X::AbsCoveredScheme)

Return the affine charts in the `default_covering` of ``X``.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> S = ambient_coordinate_ring(P);

julia> I = ideal(S, [S[1]*S[2]-S[3]^2]);

julia> X = subscheme(P, I);

julia> Xcov = covered_scheme(X)
covered scheme with 3 affine patches in its default covering

julia> affine_charts(Xcov)
3-element Vector{AbsSpec}:
 Spec of Quotient of Multivariate Polynomial Ring in (s1//s0), (s2//s0) over Rational Field by ideal((s1//s0) - (s2//s0)^2, 0)
 Spec of Quotient of Multivariate Polynomial Ring in (s0//s1), (s2//s1) over Rational Field by ideal((s0//s1) - (s2//s1)^2, 0)
 Spec of Quotient of Multivariate Polynomial Ring in (s0//s2), (s1//s2) over Rational Field by ideal((s0//s2)*(s1//s2) - 1, 0)

```
"""
affine_charts(X::AbsCoveredScheme) = basic_patches(default_covering(X))


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
name(X::AbsCoveredScheme) = get_attribute(X, :name)::String
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

