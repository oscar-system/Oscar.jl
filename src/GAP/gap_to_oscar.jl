## conversions of GAP objects to Oscar objects
## (extends the conversions from GAP.jl's `src/gap_to_julia.jl` and
## `src/constructors.jl`, where low level Julia objects are treated)

##
## GAP integer to `ZZRingElem`
##
function ZZRingElem(obj::GapObj)
    GAP.GAP_IS_INT(obj) || throw(GAP.ConversionError(obj, ZZRingElem))
    result = GC.@preserve obj ZZRingElem(GAP.ADDR_OBJ(obj), div(GAP.SIZE_OBJ(obj), sizeof(Int)))
    obj < 0 && Nemo.neg!(result, result)
    return result
end

GAP.gap_to_julia(::Type{ZZRingElem}, obj::GapInt) = ZZRingElem(obj)
(::ZZRing)(obj::GapObj) = ZZRingElem(obj)

##
## large GAP rational or integer to `QQFieldElem`
##
function QQFieldElem(obj::GapObj)
    GAP.GAP_IS_INT(obj) && return QQFieldElem(ZZRingElem(obj))
    GAP.GAP_IS_RAT(obj) || throw(GAP.ConversionError(obj, QQFieldElem))
    return QQFieldElem(ZZRingElem(GAPWrap.NumeratorRat(obj)), ZZRingElem(GAPWrap.DenominatorRat(obj)))
end

GAP.gap_to_julia(::Type{QQFieldElem}, obj::GapInt) = QQFieldElem(obj)
(::QQField)(obj::GapObj) = QQFieldElem(obj)

###
### GAP finite field elements to Oscar (generically)
###

# TODO: the following could be made faster for GAP.FFE by extracting the
# characteristic directly from the GAP FFE
characteristic(x::GAP.FFE) = ZZRingElem(GAPWrap.CHAR_FFE_DEFAULT(x))
characteristic(x::GapObj) = ZZRingElem(GAPWrap.Characteristic(x))

# test code for producing an FFE:  `GAP.Globals.Z(5)`
function (F::FinField)(x::GAP.FFE)
    characteristic(x) == characteristic(F) || error("characteristic does not match")

    if GAPWrap.DegreeFFE(x) == 1
        # FFE in GAP only exist for "small" characteristic, so we know the int
        # value fits into an Int; telling Julia about this via a type assertion
        # results in slightly better code
        val = GAPWrap.INT_FFE_DEFAULT(x)
        return F(val)
    end

    # Use `iso_oscar_gap` not `iso_gap_oscar` in order to make sure
    # that the result is in `F`, and in order to cache the isomorphism in `F`.
    iso = iso_oscar_gap(F)
    return preimage(iso, x)
end

# test code for producing a gap finite field element not stored as FFE:  `GAP.Globals.Z(65537)`
function (F::FinField)(x::GapObj)
    GAP.GAP_IS_INT(x) && return F(ZZRingElem(x))
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
    @req GAPWrap.IsMatrixOrMatrixObj(obj) "<obj> is not a GAP matrix"
end

##
## matrix of GAP integers to `ZZMatrix`
##
function ZZMatrix(obj::GapObj)
    __ensure_gap_matrix(obj)
    nrows = GAPWrap.NrRows(obj)
    ncols = GAPWrap.NrCols(obj)
    m = zero_matrix(ZZ, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        x = obj[i,j]
        @inbounds m[i,j] = x isa Int ? x : ZZRingElem(x::GapObj)
    end
    return m
end

GAP.gap_to_julia(::Type{ZZMatrix}, obj::GapObj) = ZZMatrix(obj) # TODO: deprecate/remove this

##
## matrix of GAP rationals or integers to `QQMatrix`
##
function QQMatrix(obj::GapObj)
    __ensure_gap_matrix(obj)
    nrows = GAPWrap.NrRows(obj)
    ncols = GAPWrap.NrCols(obj)
    m = zero_matrix(QQ, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        x = obj[i,j]
        @inbounds m[i,j] = x isa Int ? x : QQFieldElem(x::GapObj)
    end
    return m
end

GAP.gap_to_julia(::Type{QQMatrix}, obj::GapObj) = QQMatrix(obj) # TODO: deprecate/remove this

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

## single GAP cyclotomic to element of cyclotomic field
function (F::AnticNumberField)(obj::GapInt)
    @req Nemo.is_cyclo_type(F) "F is not a cyclotomic field"
    @req GAPWrap.IsCyc(obj) "input is not a GAP cyclotomic"
    N = get_attribute(F, :cyclo)
    @req mod(N, GAPWrap.Conductor(obj)) == 0 "obj does not embed into F"
    return preimage(iso_oscar_gap(F), obj)
end

## single GAP cyclotomic to `QQAbElem`
function QQAbElem(a::GapInt)
  c = GAPWrap.Conductor(a)
  E = abelian_closure(QQ)[2](c)
  z = parent(E)(0)
  co = GAP.Globals.CoeffsCyc(a, c)
  for i=1:c
    if !iszero(co[i])
      z += QQFieldElem(co[i])*E^(i-1)
    end
  end
  return z
end

GAP.gap_to_julia(::Type{QQAbElem}, a::GapInt) = QQAbElem(a)

(::QQAbField)(a::GAP.GapObj) = GAP.gap_to_julia(QQAbElem, a)

## nonempty list of GAP matrices over a given cyclotomic field
function matrices_over_cyclotomic_field(F::AnticNumberField, gapmats::GapObj)
    @req Nemo.is_cyclo_type(F) "F is not a cyclotomic field"
    @req GAPWrap.IsList(gapmats) "gapmats is not a GAP list"
    @req !GAPWrap.IsEmpty(gapmats) "gapmats is empty"
    @req GAPWrap.IsCyclotomicCollCollColl(gapmats) "gapmats is not a GAP list of matrices of cyclotomics"

    iso = iso_oscar_gap(F)
    result = dense_matrix_type(F)[]
    for mat in gapmats
      m = GAPWrap.NrRows(mat)
      n = GAPWrap.NrCols(mat)
      entries = [preimage(iso, mat[i, j]) for i in 1:m, j in 1:n]
      push!(result, matrix(F, entries))
    end

    return result, F, gen(F)
end

## nonempty list of GAP matrices over the same cyclotomic field
## (backwards compatibility)
matrices_over_cyclotomic_field(gapmats::GapObj) = matrices_over_field(gapmats)

## single GAP matrix of cyclotomics
function matrix(F::AnticNumberField, mat::GapObj)
    __ensure_gap_matrix(mat)
    @req Nemo.is_cyclo_type(F) "F is not a cyclotomic field"
    @req GAPWrap.IsCyclotomicCollColl(mat) "mat is not a GAP matrix of cyclotomics"
    m = GAPWrap.NrRows(mat)
    n = GAPWrap.NrCols(mat)
    iso = iso_oscar_gap(F)
    entries = [preimage(iso, mat[i, j]) for i in 1:m, j in 1:n]
    return matrix(F, entries)
end

## nonempty list of GAP matrices over the same field
function matrices_over_field(gapmats::GapObj)
    @req GAPWrap.IsList(gapmats) "gapmats is not a GAP list"
    @req !GAPWrap.IsEmpty(gapmats)  "gapmats is empty"

    if GAPWrap.IsFFECollCollColl(gapmats) ||
       GAPWrap.IsCyclotomicCollCollColl(gapmats)
      # If the entries of the matrices are FFEs or cyclotomics then
      # we construct the field generated by all entries.
      gapF = GAP.Globals.FieldOfMatrixList(gapmats)
    elseif GAPWrap.IsAlgebraicElementCollCollColl(gapmats) &&
           GAPWrap.Characteristic(gapmats) > 0
#T omit this restriction as soon as we have an iso for general alg. ext.
      # If the entries of the matrices are 'IsAlgebraicElement' then
      # we take the "parent field".
      # (This is more Oscar-like, but the reason why we use this field
      # is that subfields are not supported on the GAP side;
      # for example, already 'FieldOfMatrixList' would not work.)
      gapF = GAP.getbangproperty(GAP.Globals.FamilyObj(gapmats[1][1,1]), :wholeField)
    else
      throw(ArgumentError("gapmats is not a GAP list of matrices over a finite field"))
    end

    f = iso_gap_oscar(gapF)
    F = codomain(f)

    result = dense_matrix_type(F)[]
    for mat in gapmats
      m = GAPWrap.NrRows(mat)
      n = GAPWrap.NrCols(mat)
      entries = [f(mat[i, j]) for i in 1:m, j in 1:n]
      push!(result, matrix(F, entries))
    end

    z = gen(F)
    return result, F, z
end
