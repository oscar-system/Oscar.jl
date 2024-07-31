## conversions of Oscar objects to GAP objects
## (extends the conversions from GAP.jl's `src/julia_to_gap.jl`,
## where low level Julia objects are treated)

## `ZZRingElem` to GAP integer
function GAP.julia_to_gap(obj::ZZRingElem)
  Nemo._fmpz_is_small(obj) && return GAP.julia_to_gap(Int(obj))
  GC.@preserve obj begin
    x = Nemo._as_bigint(obj)
    return ccall((:MakeObjInt, GAP.libgap), GapObj, (Ptr{UInt64}, Cint), x.d, x.size)
  end
end

## `QQFieldElem` to GAP rational
GAP.julia_to_gap(obj::QQFieldElem) = GAP.Globals.QUO(GAP.julia_to_gap(numerator(obj)), GAP.julia_to_gap(denominator(obj)))

## `PosInf` to GAP infinity
GAP.julia_to_gap(obj::PosInf) = GAP.Globals.infinity

## `ZZMatrix` to matrix of GAP integers
GAP.julia_to_gap(obj::ZZMatrix) = GAP.julia_to_gap(Matrix(obj); recursive = true)

## `QQMatrix` to matrix of GAP rationals or integers
GAP.julia_to_gap(obj::QQMatrix) = GAP.julia_to_gap(Matrix(obj); recursive = true)

## element of cyclotomic field to GAP cyclotomic
function GAP.julia_to_gap(obj::AbsSimpleNumFieldElem)
    F = parent(obj)
    @req Nemo.is_cyclo_type(F) "the element does not lie in a cyclotomic field"
    N = get_attribute(F, :cyclo)
    v = zeros(QQFieldElem, N)
    coeffs = coefficients(obj)
    v[1:length(coeffs)] = coeffs
    return GAPWrap.CycList(GAP.julia_to_gap(v; recursive = true))
end

## `QQAbFieldElem` to GAP cyclotomic
function GAP.julia_to_gap(elm::QQAbFieldElem)
    coeffs = [Nemo.coeff(elm.data, i) for i in 0:(elm.c-1)]  # QQFieldElem
    return GAPWrap.CycList(GapObj(coeffs; recursive = true))
end

## matrix of elements of cyclotomic field to GAP matrix of cyclotomics
function GAP.julia_to_gap(obj::AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem})
    F = base_ring(obj)
    @req Nemo.is_cyclo_type(F) "the matrix entries do not lie in a cyclotomic field"
    mat = [GAP.julia_to_gap(obj[i,j]) for i in 1:nrows(obj), j in 1:ncols(obj)]
    return GAP.julia_to_gap(mat)
end

## TODO: remove the following once GAP.jl has it
function GAP.julia_to_gap(
    obj::Set{T},
    recursion_dict::IdDict{Any,Any} = IdDict();
    recursive::Bool = false,
) where {T}

    gapset = GAP.NewPlist(length(obj))
    if recursive
        recursion_dict[obj] = gapset
    end
    for x in obj
        if recursive
            x = get!(recursion_dict, x) do
                GAP.julia_to_gap(x, recursion_dict; recursive)
            end
        end
        GAP.Globals.Add(gapset, x)
    end
    GAP.Globals.Sort(gapset)
    @assert GAPWrap.IsSet(gapset)

    return gapset
end

## TODO: remove the following once GAP.jl has it
## (This will be the case when the change from
## https://github.com/oscar-system/GAP.jl/pull/989
## will be available.)
using JSON3

function GAP.julia_to_gap(
    obj::JSON3.Array,
    recursion_dict::IdDict{Any,Any} = IdDict();
    recursive::Bool = false)

    return GAP.julia_to_gap(copy(obj), recursion_dict; recursive)
end
