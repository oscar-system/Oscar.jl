## conversions of GAP objects to Oscar objects
## (extends the conversions from GAP.jl's `src/gap_to_julia.jl` and
## `src/constructors.jl`, where low level Julia objects are treated)

##
## GAP integer to `fmpz`
##
function fmpz(obj::GapObj)
    GAP.GAP_IS_INT(obj) || throw(GAP.ConversionError(obj, fmpz))
    result = GC.@preserve obj fmpz(GAP.ADDR_OBJ(obj), div(GAP.SIZE_OBJ(obj), sizeof(Int)))
    obj < 0 && Nemo.neg!(result, result)
    return result
end

GAP.gap_to_julia(::Type{fmpz}, obj::GapInt) = fmpz(obj)
(::FlintIntegerRing)(obj::GapObj) = fmpz(obj)

##
## large GAP rational or integer to `fmpq`
##
function fmpq(obj::GapObj)
    GAP.GAP_IS_INT(obj) && return fmpq(fmpz(obj))
    GAP.GAP_IS_RAT(obj) || throw(GAP.ConversionError(obj, fmpq))
    return fmpq(fmpz(GAPWrap.NumeratorRat(obj)), fmpz(GAPWrap.DenominatorRat(obj)))
end

GAP.gap_to_julia(::Type{fmpq}, obj::GapInt) = fmpq(obj)
(::FlintRationalField)(obj::GapObj) = fmpq(obj)

###
### GAP finite field elements to Oscar (generically)
###

# TODO: the following could be made faster for GAP.FFE by extracting the
# characteristic directly from the GAP FFE
characteristic(x::GAP.FFE) = fmpz(GAPWrap.CHAR_FFE_DEFAULT(x))
characteristic(x::GapObj) = fmpz(GAPWrap.Characteristic(x))

# test code for producing an FFE:  `GAP.Globals.Z(5)`
function (F::FinField)(x::GAP.FFE)
    characteristic(x) == characteristic(F) || error("characteristic does not match")

    if GAPWrap.DegreeFFE(x) == 1
        # FFE in GAP only exist for "small" characteristic, so we know the int
        # value fits into an Int; telling Julia about this via a type assertion
        # results in slightly better code
        val = GAPWrap.INT_FFE_DEFAULT(x)::Int
        return F(val)
    end

    # HACK: use `iso_oscar_gap` for now, until `iso_gap_oscar` becomes available
    iso = iso_oscar_gap(F)
    return preimage(iso, x)
end

# test code for producing a gap finite field element not stored as FFE:  `GAP.Globals.Z(65537)`
function (F::FinField)(x::GapObj)
    GAP.GAP_IS_INT(x) && return F(fmpz(x))
    GAPWrap.IsFFE(x) || error("<x> must be a GAP large integer or a GAP finite field element")
    characteristic(x) == characteristic(F) || error("characteristic does not match")

    if GAPWrap.DegreeFFE(x) == 1
        val = GAPWrap.IntFFE(x)
        return F(val)
    end

    # HACK: use `iso_oscar_gap` for now, until `iso_gap_oscar` becomes available
    iso = iso_oscar_gap(F)
    return preimage(iso, x)
end

##
## matrix conversion
##

function __ensure_gap_matrix(obj::GapObj)
    GAPWrap.IsMatrixOrMatrixObj(obj) || throw(ArgumentError("<obj> is not a GAP matrix"))
end

##
## matrix of GAP integers to `fmpz_mat`
##
function fmpz_mat(obj::GapObj)
    __ensure_gap_matrix(obj)
    nrows = GAPWrap.NrRows(obj)
    ncols = GAPWrap.NrCols(obj)
    m = zero_matrix(ZZ, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        x = obj[i,j]
        @inbounds m[i,j] = x isa Int ? x : fmpz(x::GapObj)
    end
    return m
end

GAP.gap_to_julia(::Type{fmpz_mat}, obj::GapObj) = fmpz_mat(obj) # TODO: deprecate/remove this

##
## matrix of GAP rationals or integers to `fmpq_mat`
##
function fmpq_mat(obj::GapObj)
    __ensure_gap_matrix(obj)
    nrows = GAPWrap.NrRows(obj)
    ncols = GAPWrap.NrCols(obj)
    m = zero_matrix(QQ, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        x = obj[i,j]
        @inbounds m[i,j] = x isa Int ? x : fmpq(x::GapObj)
    end
    return m
end

GAP.gap_to_julia(::Type{fmpq_mat}, obj::GapObj) = fmpq_mat(obj) # TODO: deprecate/remove this

##
## generic matrix() method for GAP matrices which converts each element on its
## own: this is inefficient but almost always works, so we use it as our base
## case
##
function matrix(R::Ring, obj::GapObj)
    # TODO: add special code for compressed matrices, resp. MatrixObj, so that
    # we can perform the characteristic check once, instead of nrows*ncols
    # times
    __ensure_gap_matrix(obj)
    nrows = GAPWrap.NrRows(obj)
    ncols = GAPWrap.NrCols(obj)
    m = zero_matrix(R, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        x = obj[i,j]::Union{Int,GapObj,GAP.FFE} # type annotation so Julia generates better code
        @inbounds m[i,j] = x isa Int ? x : R(x)
    end
    return m
end

# also allow map_entries to make Claus happy ;-)
map_entries(R::Ring, obj::GapObj) = matrix(R, obj)

# TODO: cache conversion tables


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
    N = get_attribute(F, :cyclo)
    mod(N, GAPWrap.Conductor(obj)) == 0 || throw(ArgumentError("obj does not embed into F"))

    B = default_basis_GAP_cyclotomic_field(N)
    return F(Vector{fmpz}(GAPWrap.Coefficients(B, obj)))
end

## single GAP cyclotomic to `QabElem`
function QabElem(a::GapInt)
  c = GAPWrap.Conductor(a)
  E = abelian_closure(QQ)[2](c)
  z = parent(E)(0)
  co = GAP.Globals.CoeffsCyc(a, c)
  for i=1:c
    if !iszero(co[i])
      z += fmpq(co[i])*E^(i-1)
    end
  end
  return z
end

GAP.gap_to_julia(::Type{QabElem}, a::GapInt) = QabElem(a)

(::QabField)(a::GAP.GapObj) = GAP.gap_to_julia(QabElem, a)

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
    N = get_attribute(F, :cyclo)
    B = default_basis_GAP_cyclotomic_field(N)
    m = GAPWrap.NrRows(mat)
    n = GAPWrap.NrCols(mat)
    entries = [F(Vector{fmpz}(GAPWrap.Coefficients(B, mat[i, j])))
               for i in 1:m, j in 1:n]
    return matrix(F, entries)
end
