export vanishing_ideal

########################################################################
#
# (1) AbsSpec interface
#
########################################################################

underlying_scheme(X::AffineAlgebraicSet) = X.X

@doc Markdown.doc"""
    vanishing_ideal(X::AbsAffineAlgebraicSet) -> Ideal

Return the ideal of all polynomials vanishing in ``X``.
"""
vanishing_ideal(X::AbsAffineAlgebraicSet) = ambient_closure_ideal(X)

@doc Markdown.doc"""
    ideal(X::AbsAffineAlgebraicSet{}) -> Ideal

Return the ideal of all polynomials vanishing in ``X``.
"""
ideal(X::AbsAffineAlgebraicSet) = vanishing_ideal(X)
