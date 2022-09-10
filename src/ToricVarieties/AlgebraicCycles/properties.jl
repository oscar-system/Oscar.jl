@attr Bool is_trivial(ac::RationalEquivalenceClass) = iszero(polynomial(ac))
export is_trivial
