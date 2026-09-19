using Oscar

mutable struct SAGBICandidate
    ring::ZZMPolyRing
    elements::Vector{ZZMPolyRingElem}
    leading_monomials::Dict{ZZMPolyRingElem, Vector{ZZMPolyRingElem}}
    sagbi_degree::Integer
end

function leading_monomial(
    f::ZZMPolyRingElem
)
    R = f.parent
    lm = R(1)
    for m in monomials(f)
        if m > lm
            lm = m
        end
    end
    return lm
end

function _initialize_sagbi_candidate(
    B::Vector{ZZMPolyRingElem}
)
    ring = B[1].parent
    leading_monomials = Dict{ZZMPolyRingElem, Vector{ZZMPolyRingElem}}()

    for b in B
        lm = leading_monomial(b)
        push!(get!(leading_monomials, lm, ZZMPolyRingElem[]), b)
    end

    return SAGBICandidate(
        ring, B, leading_monomials, 0
    )
end

@doc raw"""
    monoid_reprsentation(
        m::ZZMPolyRingElem, 
        generators::Vector{ZZMPolyRingElem}
    )-> Vector{Int}

    For a monomial $m$ and generators $\{g_i\}$, find the exponent vector
    $v$ such that $m = LM(g_i)^{v_i}$.
    Return `nothing` if no such representation exists.
"""
function monoid_representation(
    m::ZZMPolyRingElem,
    generators::Vector{ZZMPolyRingElem},
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


"""
    is_monomial_in_leading_monoid(m, generators)

Return `true` iff `m` is a product of nonnegative powers
of the monomials in `generators`.
"""
function is_monomial_in_leading_monoid(
    m::ZZMPolyRingElem,
    generators::Vector{ZZMPolyRingElem},
)
    return monoid_representation(m, generators) !== nothing
end

function monoid_representation(
    m::ZZMPolyRingElem,
    B::SAGBICandidate,
)
    generators = collect(keys(B.leading_monomials))
    return monoid_representation(m, generators)
end

function is_monomial_in_leading_monoid(
    m::ZZMPolyRingElem,
    B::SAGBICandidate,
)
    return monoid_representation(m, B) !== nothing
end

# --- Gröbner-based SAGBI testing (PLAN.md) ---

"""
Columns of the returned `ZZMatrix` are the exponent vectors of the leading
monomials of `B`.
"""
function _exponent_matrix(B::Vector{ZZMPolyRingElem})
    isempty(B) && return zero_matrix(ZZ, 0, 0)
    n = ngens(parent(B[1]))
    A = zero_matrix(ZZ, n, length(B))
    for (j, b) in enumerate(B)
        v = exponent_vector(leading_monomial(b), 1)
        for i in 1:n
            A[i, j] = v[i]
        end
    end
    return A
end

_exponent_matrix(B::SAGBICandidate) = _exponent_matrix(B.elements)

"""
    toric_ideal_of_leading_monomials(B)

The toric ideal `ker(ZZ[y₁,…,yₙ] → ZZ[x₁,…,xₘ])` of the monomial map
`yᵢ ↦ leading_monomial(bᵢ)`.
"""
function toric_ideal_of_leading_monomials(B::Vector{ZZMPolyRingElem})
    isempty(B) && throw(ArgumentError("B must be non-empty"))
    R = parent(B[1])
    lms = [leading_monomial(b) for b in B]
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
function tete_a_tetes(B::Vector{ZZMPolyRingElem})
    isempty(B) && return NamedTuple{(:polynomial, :a, :b), Tuple{ZZMPolyRingElem, Vector{Int}, Vector{Int}}}[]
    I = toric_ideal_of_leading_monomials(B)
    R = parent(B[1])
    ts = NamedTuple{(:polynomial, :a, :b), Tuple{ZZMPolyRingElem, Vector{Int}, Vector{Int}}}[]
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
        push!(ts, (polynomial=T, a=a, b=b))
    end
    return ts
end

tete_a_tetes(B::SAGBICandidate) = tete_a_tetes(B.elements)

@doc raw"""
    subduct(f::ZZMPolyRingElem, B::SAGBICandidate) -> ZZMPolyRingElem

    Compute the remainder of `f` after subduction by `B.elements`.
"""
function subduct(
    f::ZZMPolyRingElem,
    B::SAGBICandidate,
)
    lms = collect(keys(B.leading_monomials))
    R = parent(f)
    while !iszero(f)
        rep = monoid_representation(leading_monomial(f), lms)
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
    subduct(f::ZZMPolyRingElem, B::Vector{ZZMPolyRingElem}) -> ZZMPolyRingElem

    Compute the remainder of `f` after subduction by `B`.

    # Examples
    ```jldoctest
    R, (x, y) = polynomial_ring(ZZ, ['x','y'])
    B = [x^2 - x, y+1]
    f = x^2*y + x*y - 1
    
    subduct(f, B)
    
    #output

    2*x*y - 1
    '''
"""
subduct(f::ZZMPolyRingElem, B::Vector{ZZMPolyRingElem}) =
    isempty(B) ? f : subduct(f, _initialize_sagbi_candidate(B))

"""
    is_sagbi(B)

Gröbner-based SAGBI criterion: `B` is a SAGBI basis iff every tête-à-tête
coming from a Gröbner basis of the toric ideal of `LM(B)` subducts to zero.
"""
function is_sagbi(B::SAGBICandidate)
    for t in tete_a_tetes(B)
        iszero(subduct(t.polynomial, B)) || return false
    end
    return true
end

is_sagbi(B::Vector{ZZMPolyRingElem}) =
    isempty(B) ? true : is_sagbi(_initialize_sagbi_candidate(B))

function _monic_if_nonzero(f::ZZMPolyRingElem)
    if !iszero(f)
        return f/leading_coefficient(f)
    else
        return f
    end
end

function _compute_all_subductions(F::Vector{ZZMPolyRingElem})
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
function sagbi(
        generating_set::Vector{ZZMPolyRingElem};
        degree_bound::Integer
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
    
    leading_monomials = Dict{ZZMPolyRingElem, Vector{ZZMPolyRingElem}}()
    for f in F
        lm = leading_monomial(f)
        push!(get!(leading_monomials, lm, ZZMPolyRingElem[]), f)
    end
    
    return SAGBICandidate(
        R, F, leading_monomials, sagbi_degree 
    )
end

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