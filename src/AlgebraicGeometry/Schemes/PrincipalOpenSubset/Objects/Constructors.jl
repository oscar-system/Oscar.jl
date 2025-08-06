########################################################################
# Constructors for PrincipalOpenSubsets                                #
########################################################################

### Conversion of an arbitrary AbsAffineScheme
PrincipalOpenSubset(X::AbsAffineScheme) = PrincipalOpenSubset(X, one(OO(X)))
