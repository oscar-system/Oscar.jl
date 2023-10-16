########################################################################
#
# (1) AbsProjectiveScheme interface
#
########################################################################

underlying_scheme(X::ProjectiveVariety) = underlying_scheme(X.X)

fat_scheme(X::ProjectiveVariety) = fat_scheme(X.X)
