using Oscar

mutable struct SAGBICandidate
    ring::MPolyRing
    elements::Vector{<:MPolyRingElem}
    leading_monomials::Dict{MPolyRingElem, Vector{<:MPolyRingElem}}
    sagbi_degree::Integer
end

# function leading_monomial(
#     f::MPolyRingElem
# )
#     R = f.parent
#     lm = R(1)
#     for m in monomials(f)
#         if m > lm
#             lm = m
#         end
#     end
#     return lm
# end

function _initialize_sagbi_candidate(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    ring = B[1].parent
    leading_monomials = Dict{MPolyRingElem, Vector{<:MPolyRingElem}}()

    for b in B
        lm = leading_monomial(b; ordering)
        push!(get!(leading_monomials, lm, MPolyRingElem[]), b)
    end

    return SAGBICandidate(
        ring, B, leading_monomials, 0
    )
end

@doc raw"""
    monoid_representation(
        m::MPolyRingElem, 
        generators::Vector{<:MPolyRingElem}
    )-> Vector{Int}

    For a monomial $m$ and a generating set of monomials $\{g_i\}$, find the
    exponent vector $v$ such that $m = g_i^{v_i}$.
    Return `nothing` if no such representation exists.
"""
function monoid_representation(
    m::MPolyRingElem,
    generators::Vector{<:MPolyRingElem},
)
    isempty(generators) && return isone(m) ? Int[] : nothing

    R = parent(m)
    # Check that everything lives in the same polynomial ring.
    for g in generators
        parent(g) === R || throw(ArgumentError(
            "all monomials must belong to the same polynomial ring"
        ))
    end

    target = exponent_vector(m, 1)

    A = hcat(
        [exponent_vector(g, 1) for g in generators]...
    )

    # solve_non_negative solves Ax = target with
    # x ∈ Z_{\geq 0}^n.
    Azz = ZZMatrix(A)
    bzz = ZZMatrix(reshape(target, :, 1))

    solutions = solve_non_negative(Azz, bzz)

    # No solution.
    isempty(solutions) && return nothing

    # `solve_non_negative` returns the solutions as rows.
    # We only need one representation.
    return Vector{Int}(solutions[1, :])
end

@doc raw"""
    monoid_representation(
        m::MPolyRingElem, 
        B::SAGBICandidate,
    )-> Vector{Int}

    For a monomial $m$ and generators `B.elements={g_i}` find the exponent 
    vector $v$ such that $m = LM(g_i)^{v_i}$.
    Return `nothing` if no such representation exists.
"""
function monoid_representation(
    m::MPolyRingElem,
    B::SAGBICandidate,
)
    generators = collect(keys(B.leading_monomials))
    return monoid_representation(m, generators)
end

"""
Columns of the returned `ZZMatrix` are the exponent vectors of the leading
monomials of `B`.
"""
function _exponent_matrix(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    isempty(B) && return zero_matrix(ZZ, 0, 0)
    n = ngens(parent(B[1]))
    A = zero_matrix(ZZ, n, length(B))
    for (j, b) in enumerate(B)
        v = exponent_vector(leading_monomial(b; ordering), 1)
        for i in 1:n
            A[i, j] = v[i]
        end
    end
    return A
end

_exponent_matrix(B::SAGBICandidate) = _exponent_matrix(B.elements)

"""
    toric_ideal_of_leading_monomials(B::Vector{<:MPolyRingElem}) -> ideal

The toric ideal `ker(ZZ[y₁,…,yₙ] → ZZ[x₁,…,xₘ])` of the monomial map
`yᵢ ↦ leading_monomial(bᵢ)`.
"""
function toric_ideal_of_leading_monomials(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    isempty(B) && throw(ArgumentError("B must be non-empty"))
    R = parent(B[1])
    lms = [leading_monomial(b; ordering) for b in B]
    S, _ = polynomial_ring(base_ring(R), length(B))
    return kernel(hom(S, R, lms))
end

toric_ideal_of_leading_monomials(B::SAGBICandidate) =
    toric_ideal_of_leading_monomials(B.elements)

"""
    tete_a_tetes(B)

Translate each Gröbner basis binomial `y^α - y^β` of the toric ideal of
`LM(B)` into the polynomial tête-à-tête `B^α - B^β` in the original ring.
These are the relations that have to be tested by subduction.

Returns a vector of named tuples `(polynomial, a, b)` where `polynomial`
is the tête-à-tête `B^α - B^β`, and `a`/`b` are the exponent vectors
`α`/`β` of length `length(B)`.
"""
function tete_a_tetes(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    isempty(B) && return NamedTuple{
        (:polynomial, :a, :b),
        Tuple{MPolyRingElem,
        Vector{Int},
        Vector{Int}}
    }[]
    I = toric_ideal_of_leading_monomials(B; ordering)
    R = parent(B[1])
    tetes = NamedTuple{
        (:polynomial, :a, :b),
        Tuple{MPolyRingElem, Vector{Int},
        Vector{Int}}
    }[]
    for g in elements(groebner_basis(I))
        T = zero(R)
        exponents = Vector{Int}[]
        for (c, term) in zip(coefficients(g), terms(g))
            ev = Vector{Int}(exponent_vector(term, 1))
            push!(exponents, ev)
            u = one(R)
            for (j, e) in enumerate(ev)
                e > 0 && (u *= B[j]^e)
            end
            T += c * u
        end
        a = length(exponents) >= 1 ? exponents[1] : Int[]
        b = length(exponents) >= 2 ? exponents[2] : Int[]
        push!(tetes, (polynomial=T, a=a, b=b))
    end
    return tetes
end

tete_a_tetes(
    B::SAGBICandidate;
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
) = tete_a_tetes(B.elements; ordering)

@doc raw"""
    subduct(f::MPolyRingElem, B::SAGBICandidate) -> MPolyRingElem

    Compute the remainder of `f` after subduction by `B.elements`.
"""
function subduct(
    f::MPolyRingElem,
    B::SAGBICandidate;
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    lms = collect(keys(B.leading_monomials))
    R = parent(f)
    while !iszero(f)
        rep = monoid_representation(leading_monomial(f; ordering), lms)
        # Leading monomial not in the monoid generated by LM(B): stuck.
        rep === nothing && break
        h = one(R)
        for (i, lm) in enumerate(lms)
            rep[i] > 0 && (h *= B.leading_monomials[lm][1]^rep[i])
        end
        f = f - divexact(leading_coefficient(f), leading_coefficient(h)) * h
    end
    return f
end


@doc raw"""
    subduct(f::MPolyRingElem, B::Vector{<:MPolyRingElem}) -> MPolyRingElem

    Compute the remainder of `f` after subduction by `B`.

    # Examples
    ```jldoctest
    R, (x, y) = polynomial_ring(QQ, ['x','y'])
    B = [x^2 - x, y+1]
    f = x^2*y + x*y - 1
    
    subduct(f, B)
    
    #output
    2*x*y - 1
    '''
"""
subduct(
    f::MPolyRingElem, B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
) =
    isempty(B) ? f : subduct(f, _initialize_sagbi_candidate(B; ordering); ordering)

"""
    is_sagbi(B::SAGBICandidate) -> Bool

    Check if `B.elements` satisfies the SAGBI criterion.
"""
function is_sagbi(
    B::SAGBICandidate;
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    for t in tete_a_tetes(B; ordering)
        iszero(subduct(t.polynomial, B; ordering)) || return false
    end
    return true
end

"""
    is_sagbi(B::Vector{<:MPolyRingElem}) -> Bool

    Check if `B` satisfies the SAGBI criterion.

    # Examples
    '''jldoctest    
    Qx, x = QQ["x"];
    K, a = number_field(x^2-2, "a")
    R, (x,y,z) = polynomial_ring(K, ['x','y','z'])

    is_sagbi([x^2x, y+1])
    
    #output
    true

    is_sagbi([x+y, x^2+y^2, a*z])
    
    #output
    false
"""
is_sagbi(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
) =
    isempty(B) ? true : is_sagbi(_initialize_sagbi_candidate(B; ordering); ordering)

function _monic_if_nonzero(f::MPolyRingElem)
    if !iszero(f)
        return f/leading_coefficient(f)
    else
        return f
    end
end

function _compute_all_subductions(F::Vector{<:MPolyRingElem})
    R = parent(F[1])
    
    binom_kernel = toric_ideal_of_leading_monomials(F)
    S = base_ring(binom_kernel)
    phi = hom(S, R, F)
        
    r = [
        subduct(phi(beta), F) for beta in gens(binom_kernel)
    ]
    r = _monic_if_nonzero.(r)

    additional_terms = filter(!iszero, r)  
    return additional_terms
end


# Following Bruns & Conca 
# https://www.sciencedirect.com/science/article/pii/S0747717123000512
@doc raw"""
    sagbi(generating_set::Vector{<:MPolyRingElem}; degree_bound::Integer)
    -> SAGBICandidate

    Given `generating_set` generating a subalgebra, return a SAGBICandidate that
    completes the generating set to a SAGBI basis of the subalgebra.
    As SAGBI bases may be countably infinite, the algorithm will stop when
    all leading monomials of total degree < `degree_bound` are included in the
    SAGBI basis. 
"""
function sagbi(
        generating_set::Vector{<:MPolyRingElem};
        degree_bound::Integer,
        ordering::MonomialOrdering = default_ordering(parent(B[1]))
)

    F = generating_set
    
    min_total_deg = 0
    sagbi_degree = 0
    while true
        additional_terms = _compute_all_subductions(F)
        
        if isempty(additional_terms)
            # F is a SAGBI basis
            sagbi_degree = -1
            break
        end

        min_total_deg = minimum(
            total_degree.(additional_terms)
        )
        if min_total_deg > degree_bound
            sagbi_degree = min_total_deg - 1
            break
        end
        
        F = union(F, additional_terms)
    end
    
    leading_monomials = Dict{MPolyRingElem, Vector{<:MPolyRingElem}}()
    for f in F
        lm = leading_monomial(f; ordering)
        push!(get!(leading_monomials, lm, MPolyRingElem[]), f)
    end
    
    return SAGBICandidate(
        R, F, leading_monomials, sagbi_degree 
    )
end

@doc raw"""
    compute_sagbi_degree!(B:SAGBICandidate)

    Computes the minimum $d$ for which all leading monomials of total degree
    < $d$ are in the subalgebra generated by the leading terms of $B$, and 
    updates `B.sagbi_degree` to equal $d$. 
    
    If `B` is a SAGBI basis, then `B.sagbi_degree` will be set to -1.
"""
function compute_sagbi_degree!(B::SAGBICandidate)
    F = B.elements
    
    additional_terms = _compute_all_subductions(F)
    
    if isempty(additional_terms)
        # F is a SAGBI basis
        sagbi_degree = -1
    else
        min_total_deg = minimum(
            total_degree.(additional_terms)
        )
        sagbi_degree = min_total_deg - 1
    end

    B.sagbi_degree = sagbi_degree
    return B
end