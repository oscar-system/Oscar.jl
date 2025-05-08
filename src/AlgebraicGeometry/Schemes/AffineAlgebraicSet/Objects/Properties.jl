################################################################################
# Avoid computing the reduced structure by using fat_scheme
################################################################################

@attr Union{Int, NegInf} dim(X::AbsAffineAlgebraicSet) = dim(fat_scheme(X))
