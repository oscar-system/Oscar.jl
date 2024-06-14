########################################################################
# Constructors for PrincipalOpenSubsets                                #
########################################################################

### Conversion of an arbitrary AbsAffineScheme
PrincipalOpenSubset(X::AbsAffineScheme) = PrincipalOpenSubset(X, one(OO(X)))

function EmptyPrincipalOpenSubset(X::AbsAffineScheme)
  Ue = subscheme(X, OO(X)(1))
  U = PrincipalOpenSubset(X, Ue, OO(X)(0))
  return U
end 

