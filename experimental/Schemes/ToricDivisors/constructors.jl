# Method ambiguity requires the following
Base.:*(c::Int, td::ToricDivisor) = toric_divisor(toric_variety(td), [ZZRingElem(c)*x for x in coefficients(td)])
