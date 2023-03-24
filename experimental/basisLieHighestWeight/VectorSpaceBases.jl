# manages the basisvectors and its linear independence

using Oscar
# SparseArrays

TVec = SRow{ZZRingElem} # values ZZ, indices Int (TVec is datatype of basisvectors basisLieHighestWeight)
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
    if is_empty(v)
        return (v, 0)
    end

    pivot = first(v)[1]  # first nonzero element of vector

    return divexact(v, gcd(map(y->y[2], union(v)))), pivot
end


reduceCol(a, b, i::Int) = b[i]*a - a[i]*b # create zero entry in a


function addAndReduce!(sp::VSBasis, v::TVec)
    """
    for each pivot of sp.A we make entry of v zero and return the result and insert it into sp
    0 => linear dependent
    * => linear independent, new column element of sp.A since it increases basis
    invariants: the row of a pivotelement in any column in A is 0 (except the pivotelement)
                elements of A are integers, gcd of each column is 1
    """
    #println("sp: ", sp)
    #println("v: ", v)
    A = sp.A
    pivot = sp.pivot

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
            #return 0
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



