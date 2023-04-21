################################################################################
# Avoid computing the reduced structure by using overlying_scheme
################################################################################

dim(X::AbsAffineAlgebraicSet) = dim(overlying_scheme(X))
