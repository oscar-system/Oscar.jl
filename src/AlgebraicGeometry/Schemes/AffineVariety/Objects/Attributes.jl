########################################################################
#
# (1) AbsSpec interface
#
########################################################################

underlying_scheme(X::AffineVariety) = X.X

########################################################################
#
# (1) AbsAffineAlgebaicSet interface
#
########################################################################

overlying_scheme(X::AffineVariety) = overlying_scheme(underlying_scheme(X))
