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

## element of cyclotomic field to GAP cyclotomic
function GAP.julia_to_gap(obj::nf_elem)
    F = parent(obj)
    Nemo.iscyclo_type(F) || throw(ArgumentError("the element does not lie in a cyclotomic field"))
    N = Oscar.get_special(F, :cyclo)
    v = zeros(fmpq, N)
    coeffs = coefficients(obj)
    v[1:length(coeffs)] = coeffs
    return GAP.Globals.CycList(GAP.julia_to_gap(v, recursive = true))
end

## matrix of elements of cyclotomic field to GAP matrix of cyclotomics
function GAP.julia_to_gap(obj::AbstractAlgebra.Generic.MatSpaceElem{nf_elem})
    F = base_ring(obj)
    Nemo.iscyclo_type(F) || throw(ArgumentError("the matrix entries do not lie in a cyclotomic field"))
    mat = [GAP.julia_to_gap(obj[i,j]) for i in 1:nrows(obj), j in 1:ncols(obj)]
    return GAP.julia_to_gap(mat)
end
