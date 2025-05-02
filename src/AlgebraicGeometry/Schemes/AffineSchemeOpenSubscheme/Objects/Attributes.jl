
########################################################################
# (1) Type getters for AffineSchemeOpenSubscheme                                        #
########################################################################
open_subset_type(::Type{AffineSchemeType}) where {BRT, RT, AffineSchemeType<:AbsAffineScheme{BRT, RT}} = AffineSchemeOpenSubscheme{AffineSchemeType, BRT}
open_subset_type(X::AffineScheme) = open_subset_type(typeof(X))

ambient_type(U::AffineSchemeOpenSubscheme{AffineSchemeType, BRT}) where {AffineSchemeType<:AffineScheme, BRT} = AffineSchemeType
ambient_type(::Type{AffineSchemeOpenSubscheme{AffineSchemeType, BRT}}) where {AffineSchemeType<:AffineScheme, BRT} = AffineSchemeType

poly_type(::Type{AffineSchemeOpenSubschemeType}) where {AffineSchemeOpenSubschemeType<:AffineSchemeOpenSubscheme} = poly_type(ambient_type(AffineSchemeOpenSubschemeType))
poly_type(U::AffineSchemeOpenSubscheme) = poly_type(typeof(U))

########################################################################
# (2) Getter methods for the internally stored data                    #
########################################################################
@doc raw"""
    affine_patches(U::AffineSchemeOpenSubscheme)

Return a list of principal affine open subschemes covering ``U``.
TODO: Add example!
"""
function affine_patches(U::AffineSchemeOpenSubscheme)
  if !isdefined(U, :patches)
    X = ambient_scheme(U)
    U.patches = [PrincipalOpenSubset(X, OO(X)(f)) for f in complement_equations(U)]
  end
  return U.patches
end

@doc raw"""
    intersections(U::AffineSchemeOpenSubscheme)

Return a list of pairwise intersections of the
principal open subschemes covering ``U``.
TODO: Add example!
"""
function intersections(U::AffineSchemeOpenSubscheme)
  if !isdefined(U, :intersections)
    X = ambient_scheme(U)
    V = affine_patches(U)
    for i in 1:length(V)
      for j in 1:i-1
        U.intersections[(i,j)] = U.intersections[(j,i)] = intersect(V[i], V[j])
      end
    end
  end
  return U.intersections
end

@doc raw"""
    ambient_scheme(U::AffineSchemeOpenSubscheme)

Return the ambient scheme ``X`` of a Zariski open subset ``U ⊂ X``.
TODO: Add example!
"""
ambient_scheme(U::AffineSchemeOpenSubscheme) = U.X

@doc raw"""
    ambient_coordinate_ring(U::AffineSchemeOpenSubscheme)

For the open set `U = X \ V ` return the ambient coordinate ring of `X`.
TODO: Add example!
"""
ambient_coordinate_ring(U::AffineSchemeOpenSubscheme) = ambient_coordinate_ring(ambient_scheme(U))

@doc raw"""
    ambient_space(U::AffineSchemeOpenSubscheme) -> AffineScheme

For ``U ⊆ X \subseteq 𝔸 ⁿ`` return the affine space``𝔸 ⁿ``.
"""
ambient_space(U::AffineSchemeOpenSubscheme) = ambient_space(ambient_scheme(U))

@doc raw"""
    ambient_coordinates(U::AffineSchemeOpenSubscheme)

Return the coordinates of the ambient affine space of ``U``.
"""
ambient_coordinates(U::AffineSchemeOpenSubscheme) = coordinates(ambient_space(U))

@doc raw"""
    number_of_patches(U::AffineSchemeOpenSubscheme)

Return the number of generators stored for describing the complement of ``U``.
"""
number_of_patches(U::AffineSchemeOpenSubscheme) = length(U.gens)

@doc raw"""
    complement_equations(U::AffineSchemeOpenSubscheme)

Return the generators ``[f₁,…,fᵣ]`` stored for the description
of the complement of ``U``.
"""
complement_equations(U::AffineSchemeOpenSubscheme) = U.gens::Vector{elem_type(ambient_coordinate_ring(ambient_scheme(U)))}
number_of_complement_equations(U::AffineSchemeOpenSubscheme) = length(U.gens)

@doc raw"""
    affine_patch(U::AffineSchemeOpenSubscheme, i::Int)

Return the hypersurface complement of ``fᵢ`` in the
ambient scheme ``X`` of ``U`` where ``f₁,…,fᵣ`` are
the generators stored for the description of the complement
of ``U``. This function can also be called using the
`getindex` method or simply via `U[i]`.
"""
affine_patch(U::AffineSchemeOpenSubscheme, i::Int) = affine_patches(U)[i]
gens(U::AffineSchemeOpenSubscheme) = affine_patches(U)
gen(U::AffineSchemeOpenSubscheme, i::Int) = affine_patches(U)[i]
getindex(U::AffineSchemeOpenSubscheme, i::Int) = affine_patches(U)[i]
number_of_generators(U::AffineSchemeOpenSubscheme) = number_of_patches(U)

function getindex(U::AffineSchemeOpenSubscheme, X::AbsAffineScheme)
  for i in 1:n_patches(U)
    X === U[i] && return i
  end
  error("scheme $X not found among the open patches in $U")
end

function getindex(U::AffineSchemeOpenSubscheme, i::Int, j::Int)
  if !haskey(intersections(U), (i, j))
    intersections(U)[(i, j)] = hypersurface_complement(U[i], complement_equations(U)[j])
    intersections(U)[(j, i)] = intersections(U)[(i, j)]
  end
  return intersections(U)[(i,j)]
end

#TODO: Add docstring.
function complement_ideal(U::AffineSchemeOpenSubscheme)
  if !isdefined(U, :complement_ideal)
    I = ideal(OO(ambient_scheme(U)), complement_equations(U))
    U.complement_ideal = I
  end
  return U.complement_ideal::Ideal
end

# TODO: Add docstring.
function complement(U::AffineSchemeOpenSubscheme)
  if !isdefined(U, :complement)
    #I = radical(saturated_ideal(ideal(localized_ring(OO(ambient_scheme(U))), complement_equations(U))))
    #U.complement = subscheme(ambient_scheme(U), I)
    U.complement = subscheme(ambient_scheme(U), complement_equations(U))
  end
  return U.complement
end
function set_name!(U::AffineSchemeOpenSubscheme, name::String)
  U.name = name
end

function name(U::AffineSchemeOpenSubscheme)
  if isdefined(U, :name)
    return U.name
  end
  return "open subset of $(ambient_scheme(U))"
end

