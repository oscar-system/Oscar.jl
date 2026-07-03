## conversions of Oscar objects to GAP objects
## (extends the conversions from GAP.jl's `src/julia_to_gap.jl`,
## where low level Julia objects are treated)


## `QQAbFieldElem` to GAP cyclotomic
GAP.@install function GapObj(elm::QQAbFieldElem)
    coeffs = [Nemo.coeff(elm.data, i) for i in 0:(elm.c-1)]  # QQFieldElem
    return GAPWrap.CycList(GapObj(coeffs; recursive = true))
end


function has_GapObj_with_GapObj(input; recursive::Bool = false)
  try
    res = GapObj(input, recursive = recursive)
    return true, res
  catch e
    return false, GAP.Globals.fail
  end
end
