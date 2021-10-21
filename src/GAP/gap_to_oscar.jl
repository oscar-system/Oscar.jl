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
    return matrix(FlintQQ, m)
end

fmpq_mat(obj::GAP.GAP.GapObj) = gap_to_julia(fmpq_mat, obj)

## cache default bases of GAP's cyclotomic fields
const default_bases_GAP_cyclotomic_fields = AbstractAlgebra.WeakValueDict{Int64, GAP.GapObj}()

function default_basis_GAP_cyclotomic_field(N::Int64, F::GAP.GapObj = GAP.Globals.CyclotomicField(N))
    get!(default_bases_GAP_cyclotomic_fields, N) do
      root = GAP.Globals.E(N)
      elms = [root^i for i in 0:(euler_phi(N)-1)]
      return GAP.Globals.Basis(F, GAP.GapObj(elms))
    end
end

## single GAP cyclotomic to element of cyclotomic field
function (F::AnticNumberField)(obj::GAP.GapObj)
    Nemo.iscyclo_type(F) || throw(ArgumentError("F is not a cyclotomic field"))
    GAP.Globals.IsCyclotomic(obj) || throw(ArgumentError("input is not a GAP cyclotomic"))
    N = Oscar.get_special(F, :cyclo)
    mod(N, GAP.Globals.Conductor(obj)) == 0 || throw(ArgumentError("obj does not embed into F"))

    B = default_basis_GAP_cyclotomic_field(N)
    return F(Vector{fmpz}(GAP.Globals.Coefficients(B, obj)))
end

## nonempty list of GAP matrices over the same cyclotomic field
function matrices_over_cyclotomic_field(gapmats::GAP.GapObj)
    GAP.Globals.IsList(gapmats) || throw(ArgumentError("gapmats is not a GAP list"))
    GAP.Globals.IsEmpty(gapmats) && throw(ArgumentError("gapmats is empty"))
    GAP.Globals.ForAll(gapmats, GAP.Globals.IsCyclotomicCollColl) ||
      throw(ArgumentError("gapmats is not a GAP list of matrices of cyclotomics"))

    gapF = GAP.Globals.FieldOfMatrixList(gapmats)
    N = GAP.Globals.Conductor(gapF)
    B = default_basis_GAP_cyclotomic_field(N, gapF)
    F, z = CyclotomicField(N)

    result = dense_matrix_type(F)[]
    for mat in gapmats
      m = GAP.Globals.NumberRows(mat)
      n = GAP.Globals.NumberColumns(mat)
      entries = [F(Vector{fmpz}(GAP.Globals.Coefficients(B, mat[i, j])))
                 for i in 1:m, j in 1:n]
      push!(result, matrix(F, entries))
    end

    return result, F, z
end

## single GAP matrix of cyclotomics
function matrix(F::AnticNumberField, mat::GAP.GapObj)
    Nemo.iscyclo_type(F) || throw(ArgumentError("F is not a cyclotomic field"))
    (GAP.Globals.IsCyclotomicCollColl(mat) && GAP.Globals.IsMatrixOrMatrixObj(mat)) || throw(ArgumentError("mat is not a GAP matrix of cyclotomics"))
    N = Oscar.get_special(F, :cyclo)
    B = default_basis_GAP_cyclotomic_field(N)
    m = GAP.Globals.NumberRows(mat)
    n = GAP.Globals.NumberColumns(mat)
    entries = [F(Vector{fmpz}(GAP.Globals.Coefficients(B, mat[i, j])))
               for i in 1:m, j in 1:n]
    return matrix(F, entries)
end
