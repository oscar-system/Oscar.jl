@attr Bool is_trivial(cc::CohomologyClass) = iszero(polynomial(cc))
export is_trivial
