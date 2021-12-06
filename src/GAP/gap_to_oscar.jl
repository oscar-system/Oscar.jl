## conversions of GAP objects to Oscar objects
## (extends the conversions from GAP.jl's `src/gap_to_julia.jl` and
## `src/constructors.jl`, where low level Julia objects are treated)

import GAP: gap_to_julia
import GAP: GapObj, GapInt
import Nemo: FlintIntegerRing, FlintRationalField, MatrixSpace, fmpz, fmpq, fmpz_mat, fmpq_mat
import Oscar.AbelianClosure: QabElem

##
## GAP integer to `fmpz`
##
function fmpz(obj::GapObj)
    GAP.GAP_IS_INT(obj) || throw(GAP.ConversionError(obj, fmpz))
    GC.@preserve obj fmpz(GAP.ADDR_OBJ(obj), div(GAP.SIZE_OBJ(obj), sizeof(Int)))
end

GAP.gap_to_julia(::Type{fmpz}, obj::GapInt) = fmpz(obj)
(::FlintIntegerRing)(obj::GapObj) = gap_to_julia(fmpz, obj)

##
## large GAP rational or integer to `fmpq`
##
function fmpq(obj::GapObj)
    GAP.GAP_IS_INT(obj) && return fmpq(fmpz(obj))
    GAP.GAP_IS_RAT(obj) || throw(GAP.ConversionError(obj, fmpq))
    return fmpq(fmpz(GAPWrap.NumeratorRat(obj)), fmpz(GAPWrap.DenominatorRat(obj)))
end

GAP.gap_to_julia(::Type{fmpq}, obj::GapInt) = fmpq(obj)
(::FlintRationalField)(obj::GapObj) = gap_to_julia(fmpq, obj)

##
## matrix of GAP integers to `fmpz_mat`
##
function fmpz_mat(obj::GapObj)
    m = Matrix{fmpz}(obj)
    return matrix(ZZ, m)
end

GAP.gap_to_julia(::Type{fmpz_mat}, obj::GapObj) = fmpz_mat(obj)
matrix(R::FlintIntegerRing, obj::GapObj) = fmpz_mat(obj)

##
## matrix of GAP rationals or integers to `fmpq_mat`
##
function fmpq_mat(obj::GapObj)
    m = Matrix{fmpq}(obj)
    return matrix(QQ, m)
end

GAP.gap_to_julia(::Type{fmpq_mat}, obj::GapObj) = fmpq_mat(obj)
matrix(R::FlintRationalField, obj::GapObj) = fmpq_mat(obj)

##
## cache default bases of GAP's cyclotomic fields
##
const default_bases_GAP_cyclotomic_fields = AbstractAlgebra.WeakValueDict{Int64, GapObj}()

function default_basis_GAP_cyclotomic_field(N::Int64, F::GapObj = GAP.Globals.CyclotomicField(N))
    get!(default_bases_GAP_cyclotomic_fields, N) do
      root = GAPWrap.E(N)
      elms = [root^i for i in 0:(euler_phi(N)-1)]
      return GAP.Globals.Basis(F, GapObj(elms))
    end
end

## single GAP cyclotomic to element of cyclotomic field
function (F::AnticNumberField)(obj::GapObj)
    Nemo.iscyclo_type(F) || throw(ArgumentError("F is not a cyclotomic field"))
    GAPWrap.IsCyclotomic(obj) || throw(ArgumentError("input is not a GAP cyclotomic"))
    N = Oscar.get_special(F, :cyclo)
    mod(N, GAPWrap.Conductor(obj)) == 0 || throw(ArgumentError("obj does not embed into F"))

    B = default_basis_GAP_cyclotomic_field(N)
    return F(Vector{fmpz}(GAPWrap.Coefficients(B, obj)))
end

## single GAP cyclotomic to `QabElem`
function QabElem(cyc::GapInt)
    GAPWrap.IsCyc(cyc) || error("cyc must be a GAP cyclotomic")
    denom = GAPWrap.DenominatorCyc(cyc)
    n = GAPWrap.Conductor(cyc)
    coeffs = GAP.Globals.ExtRepOfObj(cyc * denom)
    cycpol = GAP.Globals.CyclotomicPol(n)
    dim = length(cycpol)-1
    GAP.Globals.ReduceCoeffs(coeffs, cycpol)
    coeffs = Vector{fmpz}(coeffs)
    coeffs = coeffs[1:dim]
    denom = fmpz(denom)
    FF = abelian_closure(QQ)[1]
    F, z = Oscar.AbelianClosure.cyclotomic_field(FF, n)
    val = Nemo.elem_from_mat_row(F, Nemo.matrix(Nemo.ZZ, 1, dim, coeffs), 1, denom)
    return QabElem(val, n)
end

## nonempty list of GAP matrices over the same cyclotomic field
function matrices_over_cyclotomic_field(gapmats::GapObj)
    GAPWrap.IsList(gapmats) || throw(ArgumentError("gapmats is not a GAP list"))
    GAPWrap.IsEmpty(gapmats) && throw(ArgumentError("gapmats is empty"))
    GAP.Globals.ForAll(gapmats, GAPWrap.IsCyclotomicCollColl) ||
      throw(ArgumentError("gapmats is not a GAP list of matrices of cyclotomics"))

    gapF = GAP.Globals.FieldOfMatrixList(gapmats)
    N = GAPWrap.Conductor(gapF)
    B = default_basis_GAP_cyclotomic_field(N, gapF)
    F, z = CyclotomicField(N)

    result = dense_matrix_type(F)[]
    for mat in gapmats
      m = GAPWrap.NrRows(mat)
      n = GAPWrap.NrCols(mat)
      entries = [F(Vector{fmpz}(GAPWrap.Coefficients(B, mat[i, j])))
                 for i in 1:m, j in 1:n]
      push!(result, matrix(F, entries))
    end

    return result, F, z
end

## single GAP matrix of cyclotomics
function matrix(F::AnticNumberField, mat::GapObj)
    Nemo.iscyclo_type(F) || throw(ArgumentError("F is not a cyclotomic field"))
    (GAPWrap.IsCyclotomicCollColl(mat) && GAPWrap.IsMatrixOrMatrixObj(mat)) || throw(ArgumentError("mat is not a GAP matrix of cyclotomics"))
    N = Oscar.get_special(F, :cyclo)
    B = default_basis_GAP_cyclotomic_field(N)
    m = GAPWrap.NrRows(mat)
    n = GAPWrap.NrCols(mat)
    entries = [F(Vector{fmpz}(GAPWrap.Coefficients(B, mat[i, j])))
               for i in 1:m, j in 1:n]
    return matrix(F, entries)
end
