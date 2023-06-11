################################################################################
# Avoid computing the reduced structure by using fat_scheme
################################################################################

dim(X::AbsAffineAlgebraicSet) = dim(fat_scheme(X))
