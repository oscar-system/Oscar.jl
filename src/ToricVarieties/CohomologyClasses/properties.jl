@attr Bool istrivial(cc::CohomologyClass) = iszero(polynomial(cc))
export istrivial
