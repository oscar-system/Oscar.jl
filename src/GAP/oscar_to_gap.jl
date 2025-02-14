## conversions of Oscar objects to GAP objects
## (extends the conversions from GAP.jl's `src/julia_to_gap.jl`,
## where low level Julia objects are treated)

## `ZZRingElem` to GAP integer
GAP.@install function GapObj(obj::ZZRingElem)
  Nemo._fmpz_is_small(obj) && return GapObj(Int(obj))
  GC.@preserve obj begin
    x = Nemo._as_bigint(obj)
    return ccall((:MakeObjInt, GAP.libgap), GapObj, (Ptr{UInt64}, Cint), x.d, x.size)
  end
end

## `QQFieldElem` to GAP rational
GAP.@install GapObj(obj::QQFieldElem) = GAPWrap.QUO(GapObj(numerator(obj)), GapObj(denominator(obj)))

## `PosInf` to GAP infinity
GAP.@install GapObj(obj::PosInf) = GAP.Globals.infinity

## `ZZMatrix` to matrix of GAP integers
GAP.@install GapObj(obj::ZZMatrix) = GAP.GapObj(Matrix(obj); recursive = true)

## `QQMatrix` to matrix of GAP rationals or integers
GAP.@install GapObj(obj::QQMatrix) = GAP.GapObj(Matrix(obj); recursive = true)

## element of cyclotomic field to GAP cyclotomic
GAP.@install function GapObj(obj::AbsSimpleNumFieldElem)
    F = parent(obj)
    @req Nemo.is_cyclo_type(F) "the element does not lie in a cyclotomic field"
    N = get_attribute(F, :cyclo)
    v = zeros(QQFieldElem, N)
    coeffs = coefficients(obj)
    v[1:length(coeffs)] = coeffs
    return GAPWrap.CycList(GapObj(v; recursive = true))
end

## `QQAbFieldElem` to GAP cyclotomic
GAP.@install function GapObj(elm::QQAbFieldElem)
    coeffs = [Nemo.coeff(elm.data, i) for i in 0:(elm.c-1)]  # QQFieldElem
    return GAPWrap.CycList(GapObj(coeffs; recursive = true))
end

## matrix of elements of cyclotomic field to GAP matrix of cyclotomics
GAP.@install function GapObj(obj::AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem})
    F = base_ring(obj)
    @req Nemo.is_cyclo_type(F) "the matrix entries do not lie in a cyclotomic field"
    mat = [GapObj(obj[i,j]) for i in 1:nrows(obj), j in 1:ncols(obj)]
    return GapObj(mat)
end
