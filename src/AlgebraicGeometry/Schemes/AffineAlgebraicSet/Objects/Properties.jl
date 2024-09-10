################################################################################
# Avoid computing the reduced structure by using fat_scheme
################################################################################

@attr Any dim(X::AbsAffineAlgebraicSet) = dim(fat_scheme(X))
