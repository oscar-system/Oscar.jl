########################################################################
# Constructors for PrincipalOpenSubsets                                #
########################################################################

### Conversion of an arbitrary AbsSpec
PrincipalOpenSubset(X::AbsSpec) = PrincipalOpenSubset(X, one(OO(X)))
