## conversions of GAP objects to Oscar objects
## (extends the conversions from GAP.jl's `src/gap_to_julia.jl` and
## `src/constructors.jl`, where low level Julia objects are treated)

import GAP: gap_to_julia
import Nemo: FlintIntegerRing, FlintRationalField, MatrixSpace, fmpz, fmpq, fmpz_mat, fmpq_mat

## GAP integer to `fmpz`
GAP.gap_to_julia(::Type{fmpz}, obj::Int64) = fmpz(obj)
function GAP.gap_to_julia(::Type{fmpz}, obj::GAP.GapObj)
    GAP.Globals.IsInt(obj) || throw(GAP.ConversionError(obj, fmpz))
    GC.@preserve obj fmpz(GAP.ADDR_OBJ(obj), div(GAP.SIZE_OBJ(obj), sizeof(Int)))
end

fmpz(obj::GAP.GapObj) = gap_to_julia(fmpz, obj)

(::FlintIntegerRing)(obj::GAP.GapObj) = gap_to_julia(fmpz, obj)

## large GAP rational or integer to `fmpq`
GAP.gap_to_julia(::Type{fmpq}, obj::Int64) = fmpq(obj)
function GAP.gap_to_julia(::Type{fmpq}, obj::GAP.GapObj)
    GAP.Globals.IsRat(obj) || throw(GAP.ConversionError(obj, fmpq))
    return fmpq(fmpz(GAP.Globals.NumeratorRat(obj)), fmpz(GAP.Globals.DenominatorRat(obj)))
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
    return matrix(FlintQQ, m)
end

fmpq_mat(obj::GAP.GAP.GapObj) = gap_to_julia(fmpq_mat, obj)

## nonempty list of GAP matrices over the same cyclotomic field
function matrices_over_cyclotomic_field(gapmats::GAP.GapObj)
    GAP.Globals.IsList(gapmats) || throw(ArgumentError("gapmats is not a GAP list"))
    GAP.Globals.IsEmpty(gapmats) && throw(ArgumentError("gapmats is empty"))
    GAP.Globals.ForAll(gapmats, GAP.Globals.IsCyclotomicCollColl) ||
      throw(ArgumentError("gapmats is not a GAP list of matrices of cyclotomics"))

    gapF = GAP.Globals.FieldOfMatrixList(gapmats)
    N = GAP.Globals.Conductor(gapF)
    root = GAP.Globals.E(N)
    elms = [root^i for i in 0:(euler_phi(N)-1)]
    B = GAP.Globals.Basis(gapF, GAP.GapObj(elms))
    F, z = CyclotomicField(N)

    result = []
    for mat in gapmats
      m = GAP.Globals.NumberRows(mat)
      n = GAP.Globals.NumberColumns(mat)
      entries = [F(Vector{fmpz}(GAP.Globals.Coefficients(B, mat[i, j])))
                 for i in 1:m for j in 1:n]
      push!(result, matrix(F, m, n, entries))
    end

    return convert(Vector{typeof(result[1])}, result), F, z
end
