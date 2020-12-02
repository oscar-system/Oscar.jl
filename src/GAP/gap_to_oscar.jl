## conversions of GAP objects to Oscar objects
## (extends the conversions from GAP.jl's `src/gap_to_julia.jl` and
## `src/constructors.jl`, where low level Julia objects are treated)

import GAP: gap_to_julia
import Nemo: FlintIntegerRing, FlintRationalField, MatrixSpace, fmpz, fmpq, fmpz_mat, fmpq_mat

## GAP integer to `fmpz`
GAP.gap_to_julia(::Type{fmpz}, obj::Union{GAP.GapObj, Int64}) = fmpz(BigInt(obj))

fmpz(obj::GAP.GapObj) = gap_to_julia(fmpz, obj)

(::FlintIntegerRing)(obj::GAP.GapObj) = gap_to_julia(fmpz, obj)

## large GAP rational or integer to `fmpq`
function GAP.gap_to_julia(::Type{fmpq}, obj::Union{GAP.GapObj, Int64})
    GAP.Globals.IsRat(obj) || throw(GAP.ConversionError(obj, fmpq))
    return fmpq(BigInt(GAP.Globals.NumeratorRat(obj)), BigInt(GAP.Globals.DenominatorRat(obj)))
end

fmpq(obj::GAP.GapObj) = gap_to_julia(fmpq, obj)

(::FlintRationalField)(obj::GAP.GapObj) = gap_to_julia(fmpq, obj)

## matrix of GAP integers to `fmpz_mat`
function GAP.gap_to_julia(::Type{fmpz_mat}, obj::GAP.GapObj)
    m = GAP.gap_to_julia(Matrix{fmpz}, obj)
    return fmpz_mat(size(m, 1), size(m, 2), m)
end

fmpz_mat(obj::GAP.GAP.GapObj) = gap_to_julia(fmpz_mat, obj)

## matrix of GAP rationals or integers to `fmpq_mat`
function GAP.gap_to_julia(::Type{fmpq_mat}, obj::GAP.GapObj)
    m = GAP.gap_to_julia(Matrix{fmpq}, obj)
#   return fmpq_mat(size(m, 1), size(m, 2), m) # corrupted object!
    S = MatrixSpace(FlintRationalField(), size(m, 1), size(m, 2))
    return S(m)
end

fmpq_mat(obj::GAP.GAP.GapObj) = gap_to_julia(fmpq_mat, obj)
