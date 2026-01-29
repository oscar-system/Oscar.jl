## conversions of Oscar objects to GAP objects
## (extends the conversions from GAP.jl's `src/julia_to_gap.jl`,
## where low level Julia objects are treated)

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
