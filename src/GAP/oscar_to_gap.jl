## conversions of Oscar objects to GAP objects
## (extends the conversions from GAP.jl's `src/julia_to_gap.jl`,
## where low level Julia objects are treated)

import GAP: julia_to_gap

## `fmpz` to GAP integer
function GAP.julia_to_gap(obj::fmpz)
  Nemo._fmpz_is_small(obj) && return GAP.julia_to_gap(Int(obj))
  GC.@preserve obj begin
    x = Nemo._as_bigint(obj)
    return ccall((:MakeObjInt, GAP.libgap), GapObj, (Ptr{UInt64}, Cint), x.d, x.size)
  end
end

## `fmpq` to GAP rational
GAP.julia_to_gap(obj::fmpq) = GAP.Globals.QUO(GAP.julia_to_gap(numerator(obj)), GAP.julia_to_gap(denominator(obj)))

## `fmpz_mat` to matrix of GAP integers
GAP.julia_to_gap(obj::fmpz_mat) = GAP.julia_to_gap(Matrix(obj), recursive = true)

## `fmpq_mat` to matrix of GAP rationals or integers
GAP.julia_to_gap(obj::fmpq_mat) = GAP.julia_to_gap(Matrix(obj), recursive = true)

## `GapGroup` to GAP group
GAP.GapObj(obj::GAPGroup) = return obj.X

convert(::Type{GAP.GapObj}, obj::GAPGroup) = GAP.GapObj(obj)

## `GapGroupElem` to GAP group element
GAP.GapObj(obj::GAPGroupElem) = return obj.X

convert(::Type{GAP.GapObj}, obj::GAPGroupElem) = GAP.GapObj(obj)
