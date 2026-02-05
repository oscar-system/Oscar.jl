module HomotopyContinuationExt

using Oscar, HomotopyContinuation

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
  v = Variable.(symbols(parent(f)))
  # Make the HomotopyContinuation expression
  +([
    *([Float64(c), [v[i]^e for (i,e) in enumerate(a)]...]...)
    for (c,a) in coefficients_and_exponents(f) ]...)
end

"""
    System(I::MPolyIdeal; args...)

Takes an `MPolyIdeal` and turns it into a `HomotopyContinuation.System`
containing the ideal generators. It forwards all `args` to `HomotopyContinuation.System`.
"""
System(I::MPolyIdeal; args...) = HomotopyContinuation.System(poly_to_expr.(gens(I)); args...)

"""
    nsolve(I::MPolyIdeal; args...)

Call `HomotopyContinuation.solve` on the `HomotopyContinuation.System` derived
from `I` forwarding all `args`.
"""
nsolve(I::MPolyIdeal; show_progress=false, args...) = HomotopyContinuation.solve(System(I); show_progress, args...)

witness_set(I::MPolyIdeal; show_progress=false, args...) = HomotopyContinuation.witness_set(System(I); show_progress, args...)

"""
    ndim(I::MPolyIdeal)

Compute the dimension of `I` numerically.
"""
function ndim(I::MPolyIdeal)
  # This is provided by HomotopyContinuation.jl and computes the
  # dimension based on the Jacobian rank at a random point.
  F = System(I)
  HomotopyContinuation.nvariables(F) - LinearAlgebra.rank(HomotopyContinuation.fixed(F))
end

end
