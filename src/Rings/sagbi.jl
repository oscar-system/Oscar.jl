using Oscar

mutable struct SAGBICandidate
    ring::MPolyRing
    elements::Vector{<:MPolyRingElem}
    leading_monomials::Dict{MPolyRingElem, Vector{<:MPolyRingElem}}
    least_terms_cache::Vector{<:MPolyRingElem}
    sagbi_degree::Integer
end

function _initialize_sagbi_candidate(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    ring = B[1].parent
    @req all(
        x -> parent(x) === ring, B
    ) "All polynomials must belong to the same ring"
    @req(
        coefficient_ring(ring) isa AbstractAlgebra.Field,
        "The coefficients of the polynomial ring must be a field."
    )

    leading_monomials = Dict{MPolyRingElem, Vector{<:MPolyRingElem}}()

    for b in B
        lm = leading_monomial(b; ordering)
        push!(get!(leading_monomials, lm, MPolyRingElem[]), b)
    end

    for v in values(leading_monomials)
        # this is so that when subducting we default to using the polynomial
        # with the least terms for a given leading monomial
        sort!(v, by = p -> length(p)) 
    end
    least_terms_cache = [
        leading_monomials[lm][1] for lm in keys(leading_monomials)
    ]

    return SAGBICandidate(
        ring, B, leading_monomials, least_terms_cache, 0
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
    @req all(x -> parent(x) === R, generators) "All polynomials must belong to the same ring"
    @req is_monomial(m) "m must be a monomial."
    @req all(is_monomial, generators) "The monoid generators must be monomials."

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
    _toric_ideal_of_leading_monomials(B::Vector{<:MPolyRingElem}) -> ideal

The toric ideal `ker(ZZ[y₁,…,yₙ] → ZZ[x₁,…,xₘ])` of the monomial map
`yᵢ ↦ leading_monomial(bᵢ)`.
"""
function _toric_ideal_of_leading_monomials(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    isempty(B) && throw(ArgumentError("B must be non-empty"))
    R = parent(B[1])
    @req all(x -> parent(x) === R, B) "All polynomials must belong to the same ring"
    lms = [leading_monomial(b; ordering) for b in B]
    S, _ = polynomial_ring(base_ring(R), length(B))
    return kernel(hom(S, R, lms))
end

_toric_ideal_of_leading_monomials(B::SAGBICandidate) =
    _toric_ideal_of_leading_monomials(B.elements)

@doc raw"""
    tete_a_tetes(
        B::Vector{<:MPolyRingElem};
        ordering::MonomialOrdering = default_ordering(parent(B[1]))
    ) -> Vector{NamedTuple}

Returns a vector of named tuples `(polynomial, a, b)` where `polynomial`
is the tête-à-tête 
$\prod_{i=1}^{|B|}(b_i^\alpha_i) - prod_{i=1}^{|B|}(b_i^\beta_i)$,
and `a`/`b` are the exponent vectors $\alpha$/$\beta$ of length `length(B)`.
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
    I = _toric_ideal_of_leading_monomials(B; ordering)
    R = parent(B[1])
    @req all(x -> parent(x) === R, B) "All polynomials must belong to the same ring"

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

@doc raw"""
    tete_a_tetes(
        B::SAGBICandidate;
        ordering::MonomialOrdering = default_ordering(parent(B[1]))
    ) -> Vector{NamedTuple}

Returns a vector of named tuples `(polynomial, a, b)` where `polynomial`
is the tête-à-tête 
$\prod_{i=1}^{|B|}(b_i^\alpha_i) - prod_{i=1}^{|B|}(b_i^\beta_i)$,
and `a`/`b` are the exponent vectors $\alpha$/$\beta$ of length `length(B)`.
"""
tete_a_tetes(
    B::SAGBICandidate;
    ordering::MonomialOrdering = default_ordering(parent(B.elements[1]))
) = tete_a_tetes(B.elements; ordering)

@doc raw"""
    subduct(f::MPolyRingElem, B::SAGBICandidate) -> MPolyRingElem

Compute the remainder of `f` after subduction by `B.elements`.
"""
function subduct(
    f::MPolyRingElem,
    B::SAGBICandidate;
    ordering::MonomialOrdering = default_ordering(parent(B.elements[1]))
)
    lms = collect(keys(B.leading_monomials))
    @req parent(f) === B.ring "The polynomial to subduct must be in the same ring as the basis."

    least_terms = B.least_terms_cache
    while !iszero(f)
        rep = monoid_representation(leading_monomial(f; ordering), lms)
        # Leading monomial not in the monoid generated by LM(B): stuck.
        rep === nothing && break
        h = one(B.ring)
        for (i, _) in enumerate(lms)
            rep[i] > 0 && (h *= least_terms[i]^rep[i])
        end
        c = divexact(leading_coefficient(f), leading_coefficient(h))
        f = sub!(f, f, c*h)
    end
    return f
end


@doc raw"""
    subduct(f::MPolyRingElem, B::Vector{<:MPolyRingElem}) -> MPolyRingElem

Compute the remainder of `f` after subduction by `B`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ['x','y']);

julia> B = [x^2 - x, y+1];

julia> f = x^2*y + x*y - 1;

julia> subduct(f, B);
2*x*y - 1
```
"""
function subduct(
    f::MPolyRingElem, B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    if isempty(B)
        return f
    end
    @req all(x -> parent(x) === parent(B[1]), B) "All polynomials must belong to the same ring"
    return subduct(f, _initialize_sagbi_candidate(B; ordering); ordering)
end

@doc raw"""
    is_sagbi(
        B::SAGBICandidate;
        ordering::MonomialOrdering = default_ordering(parent(B.elements[1]))
    ) -> Bool

Check if `B.elements` satisfies the SAGBI criterion with respect to
the monomial ordering `ordering`.
"""
function is_sagbi(
    B::SAGBICandidate;
    ordering::MonomialOrdering = default_ordering(parent(B.elements[1]))
)
    for t in tete_a_tetes(B; ordering)
        iszero(subduct(t.polynomial, B; ordering)) || return false
    end
    return true
end

"""
    is_sagbi(
        B::Vector{<:MPolyRingElem};
        ordering::MonomialOrdering = default_ordering(parent(B[1]))
    ) -> Bool

Check if `B` satisfies the SAGBI criterion.

# Examples
```jldoctest    
julia> Qx, x = QQ["x"];

julia> K, a = number_field(x^2-2, "a");

julia> R, (x,y,z) = polynomial_ring(K, ['x','y','z']);

julia> is_sagbi([x^2x, y+1])
true

julia> is_sagbi([x+y, x^2+y^2, a*z])
false

julia> R, (x,y) = polynomial_ring(QQ, ['x','y']);

julia> B = [x+y^2, x*y+y^3];

julia> is_sagbi(B, ordering=degrevlex(R))
false

julia> is_sagbi(B, ordering=lex(R))
true
```
"""
function is_sagbi(
    B::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
)
    if isempty(B)
        return true
    end
    @req all(x -> parent(x) === parent(B[1]), B) "All polynomials must belong to the same ring"
    return is_sagbi(_initialize_sagbi_candidate(B; ordering); ordering)
end

function _monic_if_nonzero(f::MPolyRingElem)
    if !iszero(f)
        return f/leading_coefficient(f)
    else
        return f
    end
end

function _compute_all_subductions(
    F::Vector{<:MPolyRingElem};
    ordering=default_ordering(parent(F[1]))
)
    R = parent(F[1])
    
    binom_kernel = _toric_ideal_of_leading_monomials(F; ordering)
    S = base_ring(binom_kernel)
    phi = hom(S, R, F)
        
    r = [
        subduct(phi(beta), F; ordering) for beta in gens(binom_kernel)
    ]
    r = _monic_if_nonzero.(r)

    additional_terms = filter(!iszero, r)  
    return additional_terms
end


# This is the simple algorithm called "SABGI" in Bruns & Conca (Bottom of pg 4) 
# https://doi.org/10.1016/j.jsc.2023.102237
@doc raw"""
    sagbi(
        generating_set::Vector{<:MPolyRingElem};
        degree_bound::Integer,
        ordering::MonomialOrdering = default_ordering(parent(generating_set[1]))
    )
    -> SAGBICandidate

Given `generating_set` generating a subalgebra, return a SAGBICandidate that
completes the generating set to a SAGBI basis of the subalgebra.
As SAGBI bases may be countably infinite, the algorithm will stop when
all leading monomials of total degree < `degree_bound` are included in the
SAGBI basis. 

# Examples
```jldoctest
julia> R, (x,y) = polynomial_ring(QQ, ['x','y']);

julia> B = sagbi([x+y^2, x*y+y^3]; degree_bound=10);

julia> B.elements
[x+y^2, x*y+y^3, x^3 + 2*x^2*y^2 + x*y^4]
```
"""
function sagbi(
        generating_set::Vector{<:MPolyRingElem};
        degree_bound::Integer,
        ordering::MonomialOrdering = default_ordering(parent(generating_set[1]))
)
    @req !isempty(generating_set) "The empty subalgebra has no nonempty SAGBI basis."
    @req all(
        x -> parent(x) === parent(generating_set[1]), generating_set
    ) "All polynomials must belong to the same ring"

    F = generating_set
    
    min_total_deg = 0
    sagbi_degree = 0
    while true
        additional_terms = _compute_all_subductions(F; ordering)
        
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

    B = _initialize_sagbi_candidate(F)
    B.sagbi_degree = sagbi_degree
    return B
end

@doc raw"""
    compute_sagbi_degree!(
        B::SAGBICandidate;
        ordering=default_ordering(parent(B.elements[1]))
    ) -> SAGBICandidate

Computes the minimum $d$ for which all leading monomials of total degree
< $d$ are in the subalgebra generated by the leading terms of $B$, and 
updates `B.sagbi_degree` to equal $d$. 

If `B` is a SAGBI basis, then `B.sagbi_degree` will be set to -1.
"""
function compute_sagbi_degree!(
    B::SAGBICandidate;
    ordering=default_ordering(parent(B.elements[1]))
)
    F = B.elements
    
    additional_terms = _compute_all_subductions(F; ordering)
    
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

@doc raw"""
    compute_sagbi_degree(
        B::Vector{<:MPolyRingElem};
        ordering=default_ordering(parent(B[1]))
    ) -> Integer

Computes the minimum $d$ for which all leading monomials of total degree
< $d$ are in the subalgebra generated by the leading terms of $B$, and 
updates `B.sagbi_degree` to equal $d$. 

If `B` is a SAGBI basis, then `B.sagbi_degree` will be set to -1.

# Examples
```jldoctest
julia> R, (x,y) = polynomial_ring(QQ, ['x','y'])
julia> B = [x, x*y - y^2, x*y^2, x*y^3 - 1//2*y^4, x*y^5 - 1//3*y^6,  x*y^4]
julia> compute_sagbi_degree(B)
6
    ```
"""
function compute_sagbi_degree(
    B::Vector{<:MPolyRingElem};
    ordering=default_ordering(parent(B[1]))
)
    sagbi_candidate = _initialize_sagbi_candidate(B; ordering)
    return compute_sagbi_degree!(sagbi_candidate).sagbi_degree
end