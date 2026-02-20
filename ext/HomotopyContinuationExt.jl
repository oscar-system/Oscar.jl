module HomotopyContinuationExt

using Oscar, HomotopyContinuation

import HomotopyContinuation: HomotopyContinuation, Variable

"""
    Expression(f::MPolyRingElem)

Takes a polynomial `f` an `MPolyRingElem` and converts it into a
`HomotopyContinuation.Expression`. The `HomotopyContinuation.Variable`
objects used in the expression have the same names as the generators
in the `parent(f)` (including brackets as in `x[1]`).
"""
function HomotopyContinuation.Expression(f::MPolyRingElem)
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
HomotopyContinuation.System(I::MPolyIdeal; args...) = HomotopyContinuation.System(Expression.(gens(I)); args...)

"""
    solve(I::MPolyIdeal; args...)
    solve(I::Vector{MPolyRingElem})

Call `HomotopyContinuation.solve` on the `HomotopyContinuation.System` derived
from `I` forwarding all `args`.
"""
function HomotopyContinuation.solve(I::Vector{MPolyRingElem}; show_progress=false, args...)
  HomotopyContinuation.solve(Expression.(I); show_progress, args...)
end

function HomotopyContinuation.solve(I::MPolyIdeal; args...)
  return HomotopyContinuation.solve(HomotopyContinuation.System(gens(I)); args...)
end

HomotopyContinuation.witness_set(I::MPolyIdeal; show_progress=false, args...) = HomotopyContinuation.witness_set(HomotopyContinuation.System(I); show_progress, args...)

"""
    dim_numerical(I::MPolyIdeal)

Compute the dimension of `I` numerically.
"""
function Oscar.dim_numerical(I::MPolyIdeal)
  # This is provided by HomotopyContinuation.jl and computes the
  # dimension based on the Jacobian rank at a random point.
  F = HomotopyContinuation.System(I)
  HomotopyContinuation.nvariables(F) - rank(HomotopyContinuation.fixed(F))
end

end
