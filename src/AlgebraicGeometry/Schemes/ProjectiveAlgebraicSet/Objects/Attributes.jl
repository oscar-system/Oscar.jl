########################################################################
#
# (1) AbsProjectiveScheme interface
#
########################################################################

fat_scheme(X::ProjectiveAlgebraicSet) = X.X

function underlying_scheme(X::ProjectiveAlgebraicSet)
  if isdefined(X, :Xred)
    return X.Xred
  end
  X.Xred = reduced_scheme(fat_scheme(X))
  return X.Xred
end

########################################################################
#
# (2) Further attributes
#
########################################################################

@doc raw"""
    vanishing_ideal(X::AbsProjectiveAlgebraicSet) -> Ideal

Return the ideal of all homogeneous polynomials vanishing in ``X``.
"""
vanishing_ideal(X::AbsProjectiveAlgebraicSet) = defining_ideal(X)


@doc raw"""
    fat_ideal(X::AbsProjectiveAlgebraicSet) -> Ideal

Return a homogeneous ideal whose radical is the vanishing ideal of `X`.
"""
fat_ideal(X::AbsProjectiveAlgebraicSet) = defining_ideal(fat_scheme(X))


base_ring(X::AbsProjectiveAlgebraicSet) = base_ring(fat_scheme(X))

dim(X::AbsProjectiveAlgebraicSet) = dim(fat_scheme(X))

codim(X::AbsProjectiveAlgebraicSet) = codim(fat_scheme(X))

homogeneous_coordinate_ring(X::AbsProjectiveAlgebraicSet) = homogeneous_coordinate_ring(fat_scheme(X))

relative_ambient_dimension(X::AbsProjectiveAlgebraicSet) = relative_ambient_dimension(fat_scheme(X))

