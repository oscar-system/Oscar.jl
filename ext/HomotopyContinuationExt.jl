module HomotopyContinuationExt

using Oscar, HomotopyContinuation

import HomotopyContinuation: HomotopyContinuation, Variable, Expression, System, witness_set

"""
    Expression(f::QQMPolyRingElem)

Takes a polynomial `f` an `MPolyRingElem` and converts it into a
`HomotopyContinuation.Expression`. The `HomotopyContinuation.Variable`
objects used in the expression have the same names as the generators
in the `parent(f)` (including brackets as in `x[1]`).
"""
function Expression(f::QQMPolyRingElem)
  # Get a list of HomotopyContinuation variables whose names are the
  # same as the ones in the Oscar polynomial f.
  v = Variable.(symbols(parent(f)))
  # Make the HomotopyContinuation expression
  return sum(Float64(c) * prod(v[i]^e for (i,e) in enumerate(a)) for (c,a) in coefficients_and_exponents(f))
end

"""
    System(I::MPolyIdeal{QQMPolyRingElem}; args...)
    System(I::Vector{QQMPolyRingElem}; args...)

Takes an `MPolyIdeal` and turns it into a `HomotopyContinuation.System`
containing the ideal generators. It forwards all `args` to `HomotopyContinuation.System`.
"""
System(I::Vector{QQMPolyRingElem}; args...) = System(Expression.(I); args...)
System(I::MPolyIdeal{QQMPolyRingElem}; args...) = System(gens(I); args...)

"""
    solve(I::MPolyIdeal{QQMPolyRingElem}; show_progress=false, args...)
    solve(I::Vector{MPolyRingElem}; show_progress=false, args...)

Call `HomotopyContinuation.solve` on the `HomotopyContinuation.System` derived
from `I` forwarding all `args`.
"""
function HomotopyContinuation.solve(I::Vector{QQMPolyRingElem}; show_progress=false, args...)
  return HomotopyContinuation.solve(System(Expression.(I), args...), show_progress=show_progress)
end
                                    
function HomotopyContinuation.solve(I::MPolyIdeal{QQMPolyRingElem}; args...)
  return HomotopyContinuation.solve(gens(I); args...)
end

witness_set(I::Vector{QQMPolyRingElem}; show_progress=false, args...) = witness_set(System(I); show_progress, args...)
witness_set(I::MPolyIdeal{QQMPolyRingElem}; show_progress=false, args...) = witness_set(gens(I); show_progress, args...)

"""
    dim_numerical(I::MPolyIdeal)

Compute the dimension of `I` numerically.
"""
function Oscar.dim_numerical(I::MPolyIdeal{QQMPolyRingElem})
  # This is provided by HomotopyContinuation.jl and computes the
  # dimension based on the Jacobian rank at a random point.
  F = System(I)
  return HomotopyContinuation.nvariables(F) - rank(HomotopyContinuation.fixed(F))
end

end
