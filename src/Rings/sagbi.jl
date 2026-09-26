using Oscar

mutable struct SAGBICandidate
    ring::MPolyRing
    elements::Vector{<:MPolyRingElem}
    leading_monomials::Dict{MPolyRingElem, Vector{<:MPolyRingElem}}
    ordering::MonomialOrdering
    least_terms_cache::Vector{<:MPolyRingElem}
    sagbi_degree::Integer
end

function _initialize_sagbi_candidate(
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
    ring = generators[1].parent
    @req all(
        x -> parent(x) === ring, generators
    ) "All polynomials must belong to the same ring"
    @req(
        coefficient_ring(ring) isa AbstractAlgebra.Field,
        "The coefficients of the polynomial ring must be a field."
    )

    leading_monomials = Dict{MPolyRingElem, Vector{<:MPolyRingElem}}()

    for b in generators
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
        ring, generators, leading_monomials, ordering, least_terms_cache, 0
    )
end

const zsolve_cmd = Oscar.lib4ti2_jll.zsolve()
function _4ti2_solve(A::ZZMatrix, b::ZZMatrix)
    n, m = nrows(A), ncols(A)

    # mktempdir with a do-block guarantees cleanup, even on errors
    return mktempdir() do dir
        proj = joinpath(dir, "zsolve")

        open("$proj.mat", "w") do f
            write(f, "$n $m\n")
            for i in 1:n
                write(f, join((string(A[i, j]) for j in 1:m), ' '), "\n")
            end
        end
        
        open("$proj.rhs", "w") do f
            write(f, "1 $n\n", join((string(b[i, 1]) for i in 1:n), ' '), "\n")
        end
        
        open("$proj.rel", "w") do f
            write(f, "1 $n\n", join(("=" for _ in 1:n), ' '), "\n")
        end
        
        open("$proj.sign", "w") do f
            write(f, "1 $m\n", join(("1" for _ in 1:m), ' '), "\n")
        end

        success(`$zsolve_cmd -p gmp $proj`) ||
            error("Error running 4ti2 zsolve")
            
        lines = readlines("$proj.zinhom")
        nr = parse(Int, split(lines[1])[1])        
        return nr == 0 ? nothing : parse.(Int, split(lines[2]))
    end
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
    if isempty(generators)
        return isone(m) ? Int[] : nothing
    end 

    R = parent(m)
    @req(all(x -> parent(x) === R, generators),
    "All polynomials must belong to the same ring")
    @req is_monomial(m) "m must be a monomial."
    @req all(is_monomial, generators) "The monoid generators must be monomials."

    target = exponent_vector(m, 1)
    if iszero(target)
        return zeros(Int, length(generators))
    end
    if total_degree(m) < minimum(total_degree, generators)
        return nothing
    end

    gen_exps = [exponent_vector(g, 1) for g in generators]

    gen_vars = Set{Int}()
    for e in gen_exps
        union!(gen_vars, findall(!iszero, e))
    end
    if !issubset(findall(!iszero, target), gen_vars)
        return nothing # m has a variable that isn't in the generators
    end

    A = hcat(gen_exps...)

    # We must solve A*x = target with x >= 0. 
    # We will try just solving A*x = target first, but if the solution
    # is negative then we need to do integer programming.
    Azz = ZZMatrix(A)
    bzz = ZZMatrix(reshape(target, :, 1))
    H, U = hnf_with_transform(Azz)
    solvable, x = can_solve_with_solution(H, U*bzz; side=:right)
    if !solvable
        # Not solvable over ZZ, hence not over Z_{>=0} either.
        return nothing
    end
    if all(i -> x[i, 1] >= 0, 1:nrows(x))
        # The HNF solution is already non-negative.
        return [Int(x[i, 1]) for i in 1:nrows(x)]
    end

    # Fall back to integer programming.
    solution = _4ti2_solve(Azz, bzz)
    return solution
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
function monoid_representation(m::MPolyRingElem, B::SAGBICandidate)
    generators = collect(keys(B.leading_monomials))
    return monoid_representation(m, generators)
end

"""
Columns of the returned `ZZMatrix` are the exponent vectors of the leading
monomials of `B`.
"""
function _exponent_matrix(
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
    isempty(generators) && return zero_matrix(ZZ, 0, 0)
    n = ngens(parent(generators[1]))
    A = zero_matrix(ZZ, n, length(generators))
    for (j, b) in enumerate(generators)
        v = exponent_vector(leading_monomial(b; ordering), 1)
        for i in 1:n
            A[i, j] = v[i]
        end
    end
    return A
end

_exponent_matrix(B::SAGBICandidate) = _exponent_matrix(B.elements)

"""
    _toric_ideal_of_leading_monomials(
        generators::Vector{<:MPolyRingElem}
    ) -> ideal

The toric ideal `ker(ZZ[y₁,…,yₙ] → ZZ[x₁,…,xₘ])` of the monomial map
`yᵢ ↦ leading_monomial(bᵢ)`.
"""
function _toric_ideal_of_leading_monomials(
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
    isempty(generators) && throw(ArgumentError("B must be non-empty"))

    R = parent(generators[1])
    @req(all(x -> parent(x) === R, generators),
    "All polynomials must belong to the same ring")

    lms = [leading_monomial(b; ordering) for b in generators]
    S, _ = polynomial_ring(base_ring(R), length(generators))
    return kernel(hom(S, R, lms))
end

_toric_ideal_of_leading_monomials(B::SAGBICandidate) =
    _toric_ideal_of_leading_monomials(B.elements)

@doc raw"""
    _toric_ideal_lattice(B::SAGBICandidate)

Compute a set of lattice generators for `ker(ZZ[y₁,…,yₙ] → ZZ[x₁,…,xₘ])` of
the monomial map `yᵢ ↦ leading_monomial(bᵢ)`.
"""
function _toric_ideal_lattice(B::SAGBICandidate)
    if isempty(B.elements)
        return zero_matrix(ZZ, 0, 0)
    end

    A = ZZMatrix(hcat([
            exponent_vector(leading_monomial(b; ordering=B.ordering), 1) 
            for b in B.elements]...
        ))

    H, U = hnf_with_transform(transpose(A))

    kernel_basis = Vector{Vector{Int}}()
    for i in 1:nrows(H)
        if is_zero_row(H, i)
            push!(kernel_basis, Vector{Int}(U[i, :]))
        end
    end

    if isempty(kernel_basis)
        return zero_matrix(ZZ, 0, length(B.elements))
    end
    
    return matrix(ZZ, kernel_basis)
end
    

@doc raw"""
    tete_a_tetes(
        generators::Vector{<:MPolyRingElem};
        ordering::MonomialOrdering = default_ordering(parent(B[1]))
    ) -> Vector{NamedTuple}

Returns a vector of named tuples `(polynomial, a, b)` where `polynomial`
is the tête-à-tête 
$\prod_{i=1}^{|B|}(b_i^\alpha_i) - prod_{i=1}^{|B|}(b_i^\beta_i)$,
and `a`/`b` are the exponent vectors $\alpha$/$\beta$ of length `length(B)`.
"""
function tete_a_tetes(
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
    isempty(generators) && return NamedTuple{
        (:polynomial, :a, :b), Tuple{MPolyRingElem, Vector{Int}, Vector{Int}}
    }[]

    R = parent(generators[1])
    @req(all(x -> parent(x) === R, generators),
    "All polynomials must belong to the same ring")
    
    B = _initialize_sagbi_candidate(generators; ordering)
    L = _toric_ideal_lattice(B)
    
    if nrows(L) == 0
        return NamedTuple{
            (:polynomial, :a, :b),
            Tuple{MPolyRingElem, Vector{Int}, Vector{Int}}
        }[]
    end
    
    # Each row r of the Markov basis matrix indexes a binomial generator
    # (Π_i f_i^{a_i} - Π_i f_i^{b_i})
    # where a_i = max(0, r_i) and b_i = max(0, -r_i)
    M = Oscar.markov4ti2(L)

    tetes = NamedTuple{
        (:polynomial, :a, :b), Tuple{MPolyRingElem, Vector{Int}, Vector{Int}}
    }[]
    for i in 1:nrows(M)
        v = Vector{Int}(M[i, :])
        
        # Split into positive and negative parts
        a = [max(0, x) for x in v]
        b = [max(0, -x) for x in v]
        
        # Construct the tête-à-tête polynomial
        # ∏ B.elements[j]^a[j] - ∏ B.elements[j]^b[j]
        term1 = one(R)
        term2 = one(R)
        for (j, (aj, bj)) in enumerate(zip(a, b))
            if aj > 0
                term1 *= B.elements[j]^aj
            end
            if bj > 0
                term2 *= B.elements[j]^bj
            end
        end
        
        push!(tetes, (polynomial=term1 - term2, a=a, b=b))
    end
    
    return tetes
end

@doc raw"""
    tete_a_tetes(B::SAGBICandidate) -> Vector{NamedTuple}

Returns a vector of named tuples `(polynomial, a, b)` where `polynomial`
is the tête-à-tête 
$\prod_{i=1}^{|B|}(b_i^\alpha_i) - prod_{i=1}^{|B|}(b_i^\beta_i)$,
and `a`/`b` are the exponent vectors $\alpha$/$\beta$ of length `length(B)`.
"""
tete_a_tetes(B::SAGBICandidate) = tete_a_tetes(B.elements; ordering=B.ordering)

@doc raw"""
    subduct(f::MPolyRingElem, B::SAGBICandidate) -> MPolyRingElem

Compute the remainder of `f` after subduction by `B.elements`.
"""
function subduct(f::MPolyRingElem, B::SAGBICandidate)
    lms = collect(keys(B.leading_monomials))
    @req(parent(f) === B.ring,
    "The polynomial to subduct must be in the same ring as the basis.")

    least_terms = B.least_terms_cache
    while !iszero(f)
        rep = monoid_representation(leading_monomial(f; ordering=B.ordering), lms)
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
    subduct(
        f::MPolyRingElem, generators::Vector{<:MPolyRingElem}
    )-> MPolyRingElem

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
    f::MPolyRingElem, generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
    if isempty(generators)
        return f
    end
    @req(all(x -> parent(x) === parent(generators[1]), generators),
    "All polynomials must belong to the same ring")

    return subduct(f, _initialize_sagbi_candidate(generators; ordering))
end

@doc raw"""
    is_sagbi(B::SAGBICandidate) -> Bool

Check if `B.elements` satisfies the SAGBI criterion with respect to
the monomial ordering `B.ordering`.
"""
function is_sagbi(
    B::SAGBICandidate
)
    for t in tete_a_tetes(B)
        iszero(subduct(t.polynomial, B)) || return false
    end
    return true
end

"""
    is_sagbi(
        generators::Vector{<:MPolyRingElem};
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
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
    if isempty(generators)
        return true
    end
    @req(all(x -> parent(x) === parent(generators[1]), generators),
    "All polynomials must belong to the same ring")
    return is_sagbi(_initialize_sagbi_candidate(generators; ordering))
end

function _monic_if_nonzero(f::MPolyRingElem)
    if !iszero(f)
        return f/leading_coefficient(f)
    else
        return f
    end
end

function _compute_all_subductions(
    B::SAGBICandidate
)
    tetes = tete_a_tetes(B)

    r = [
        subduct(t.polynomial, B) for t in tetes
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
        additional_terms = _compute_all_subductions(
            _initialize_sagbi_candidate(F; ordering)
        )
        
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
    compute_sagbi_degree!(B::SAGBICandidate) -> SAGBICandidate

Computes the minimum $d$ for which all leading monomials of total degree
< $d$ are in the subalgebra generated by the leading terms of $B$, and 
updates `B.sagbi_degree` to equal $d$. 

If `B` is a SAGBI basis, then `B.sagbi_degree` will be set to -1.
"""
function compute_sagbi_degree!(B::SAGBICandidate)
    additional_terms = _compute_all_subductions(B)
    
    if isempty(additional_terms)
        # B is a SAGBI basis
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
        generators::Vector{<:MPolyRingElem};
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
    generators::Vector{<:MPolyRingElem};
    ordering=default_ordering(parent(generators[1]))
)
    sagbi_candidate = _initialize_sagbi_candidate(generators; ordering)
    return compute_sagbi_degree!(sagbi_candidate).sagbi_degree
end