# manages the linear (in-)dependence of integer vectors
# this file is only of use to basis_lie_highest_weight for the basis vectors

using Oscar

TVec = SRow{ZZRingElem} # TVec is datatype of basisvectors in basisLieHighestWeight
Short = UInt8 # for exponents of monomials; max. 255

struct VSBasis
    basis_vectors::Vector{TVec} # vector of basisvectors
    pivot::Vector{Int} # vector of pivotelements, i.e. pivot[i] is first nonzero element of basis_vectors[i]
end

nullSpace() = VSBasis([], []) # empty Vektorraum

reduce_col(a::TVec, b::TVec, i::Int) = (b[i]*a - a[i]*b)::TVec # create zero entry in a

function normalize(v::TVec)::Tuple{TVec, Int64}
    """
    divides vector by gcd of nonzero entries, returns vector and first nonzero index
    used: add_and_reduce!
    """
    if is_empty(v)
        return (v, 0)
    end
    pivot = first(v)[1]  # first nonzero element of vector
    return divexact(v, gcd(map(y->y[2], union(v)))), pivot
end

function add_and_reduce!(sp::VSBasis, v::TVec)::TVec
    """
    for each pivot of sp.basis_vectors we make entry of v zero and return the result and insert it into sp
    0 => linear dependent
    * => linear independent, new column element of sp.basis_vectors since it increases basis
    invariants: the row of a pivotelement in any column in basis_vectors is 0 (except the pivotelement)
                elements of basis_vectors are integers, gcd of each column is 1
    """
    # initialize objects
    basis_vectors = sp.basis_vectors
    pivot = sp.pivot
    v, newPivot = normalize(v)

    # case v zero vector
    if newPivot == 0
        return v
    end

    # use pivots of basis basis_vectors to create zeros in v
    for j in 1:length(basis_vectors)
        i = pivot[j]
        if i != newPivot
            continue
        end
        v = reduce_col(v, basis_vectors[j], i)
        v, newPivot = normalize(v)
        if newPivot == 0
            #return 0
            return v
        end
    end

    # new pivot element of v
    pos = findfirst(pivot .> newPivot)
    if (pos === nothing)
        pos = length(pivot) + 1
    end

    # save result
    insert!(basis_vectors, pos, v)
    insert!(pivot, pos, newPivot)
    return v
end



