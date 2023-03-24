# manages the basisvectors and its linear independence

using Oscar
using SparseArrays

ZZ = Int
TVec = SparseVector{ZZ, Int} # values ZZ, indices Int (TVec is datatype of basisvectors basisLieHighestWeight)
Short = UInt8 # for exponents of monomials; max. 255

struct VSBasis
    A::Vector{TVec} # vector of basisvectors
    pivot::Vector{Int} # vector of pivotelements, i.e. pivot[i] is first nonzero element of A[i]
    dim::Vector{Int} # dimension
end


nullSpace() = VSBasis([], [], []) # empty Vektorraum


function normalize(v::TVec)
    """
    divides vector by gcd of nonzero entries, returns vector and first nonzero index
    used: addAndReduce!
    """
    dropzeros!(v)
    if isempty(v.nzind)
        return (v, 0)
    end

    pivot = v.nzind[1]  # first nonzero element of vector

    return v .÷ gcd(v.nzval), pivot
end


reduceCol(a, b, i::Int) = a*b[i] - b*a[i] # create zero entry in a


function addAndReduce!(sp::VSBasis, v::TVec)
    """
    for each pivot of sp.A we make entry of v zero and return the result and insert it into sp
    0 => linear dependent
    * => linear independent, new column element of sp.A since it increases basis
    invariants: the row of a pivotelement in any column in A is 0 (except the pivotelement)
                elements of A are integers, gcd of each column is 1
    """
    A = sp.A
    pivot = sp.pivot
    dim = sp.dim

    if dim == [] # update dimension
        insert!(dim, 1, length(v))
    end

    if length(A) == sp.dim[1] # space already full => linear dependence guaranteed
        v = spzeros(ZZ, sp.dim[1])
        return v
    end

    v, newPivot = normalize(v) 
    if newPivot == 0 # v zero vector
        return v
    end
 
    for j = 1:length(A)
        i = pivot[j]
        if i != newPivot
            continue
        end
        v = reduceCol(v, A[j], i)
        v, newPivot = normalize(v)
        if newPivot == 0
            return v
        end
    end

    pos = findfirst(pivot .> newPivot)
    if (pos === nothing)
        pos = length(pivot) + 1
    end

    insert!(A, pos, v)
    insert!(pivot, pos, newPivot)
    return v
end



