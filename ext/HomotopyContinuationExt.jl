module HomotopyContinuationExt

import HomotopyContinuation: HomotopyContinuation, Variable

"""
    poly_to_expr(f::MPolyRingElem)

Takes a polynomial `f` an `MPolyRingElem` and converts it into a
`HomotopyContinuation.Expression`. The `HomotopyContinuation.Variable`
objects used in the expression have the same names as the generators
in the `parent(f)` (including brackets as in `x[1]`).
"""
function poly_to_expr(f::MPolyRingElem)
  # Get a list of HomotopyContinuation variables whose names are the
  # same as the ones in the Oscar polynomial f.
  v = map(x -> Variable(symbol(x)), gens(parent(f)))
  # Make the HomotopyContinuation expression
  sum([
    prod([Rational(c), [v[i]^e for (i,e) in enumerate(a)]...]...)
    for (c,a) in coefficients_and_exponents(f) ]...)
end

end
