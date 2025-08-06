@doc raw"""
    chow_ring(M::Matroid; ring::MPolyRing=nothing, extended::Bool=false)

Return the Chow ring of a matroid, optionally also with the simplicial generators and the polynomial ring.

See [AHK18](@cite) and [BES23](@cite).

# Examples
The following computes the Chow ring of the Fano matroid.
```jldoctest
julia> M = fano_matroid();

julia> R = chow_ring(M);

julia> R[1]*R[8]
-x_{3,4,7}^2
```

The following computes the Chow ring of the Fano matroid including variables for the simplicial generators.
```jldoctest
julia> M = fano_matroid();

julia> R = chow_ring(M, extended=true);

julia> f = R[22] + R[8] - R[29]
x_{1,2,3} + h_{1,2,3} - h_{1,2,3,4,5,6,7}

julia> f==0
true
```

The following computes the Chow ring of the free matroid on three elements in a given graded polynomial ring.
```jldoctest
julia> M = uniform_matroid(3,3);

julia> GR, _ = graded_polynomial_ring(QQ,[:a,:b,:c,:d,:e,:f]);

julia> R = chow_ring(M, ring=GR);

julia> hilbert_series_reduced(R)
(t^2 + 4*t + 1, 1) 

```
"""
function chow_ring(M::Matroid; ring::Union{MPolyRing,Nothing}=nothing, extended::Bool=false, graded::Bool=false)
    is_loopless(M) || error("Matroid has loops")
    Flats = flats(M)
    number_flats = length(Flats)
    number_flats >= 2 || error("matroid has to few flats")
    proper_flats = Flats[2:number_flats-1]

    #construct polynomial ring and extract variables
    if ring===nothing
        # create variable names, indexed by the proper flats of M
        var_names = ["x_{" * join(S, ",") * "}" for S in proper_flats]
        if extended
            #add the variables for the simplicial generators
            var_names = [var_names; ["h_{" * join(S, ",") * "}" for S in [proper_flats;[Flats[number_flats]]]]]
        end

        @req length(var_names) > 0 "Chow ring is empty"
        if graded
            ring, vars = graded_polynomial_ring(QQ, var_names; cached=false)
        else
            ring, vars = polynomial_ring(QQ, var_names; cached=false)
        end
    else
        if extended
            nvars(ring) == 2*length(proper_flats)+1 || error("the ring has the wrong number of variables")
        else
            nvars(ring) == length(proper_flats) || error("the ring has the wrong number of variables")
        end
        vars = gens(ring)
    end

    #construct ideal and quotient
    I = linear_relations(ring, proper_flats, vars, M)
    J = quadratic_relations(ring, proper_flats, vars)
    Ex = elem_type(ring)[]
    if extended
        Ex = relations_extended_ring(ring, proper_flats, vars)
    end
    chow_modulus = ideal(ring, vcat(I, J, Ex))
    chow_ring, projection = quo(ring, chow_modulus)
    return chow_ring
end

function linear_relations(ring::MPolyRing, proper_flats::Vector{Vector{T}}, vars::Vector, M::Matroid) where T <: ElementType
    alpha = zero(ring)
    relations = elem_type(ring)[]
    for i in M.groundset
        poly = ring(0)
        for index in findall(F->issubset([i], F), proper_flats)
            poly+= vars[index]
        end
        if i==M.groundset[1]
            alpha = poly
        else
            push!(relations, alpha-poly)
        end
    end
    return relations
end

function quadratic_relations(ring::MPolyRing, proper_flats::Vector{Vector{T}}, vars::Vector) where T <: ElementType
    relations = elem_type(ring)[]
    for i in 1:length(proper_flats)
        F = proper_flats[i]
        for j in 1:i-1
            G = proper_flats[j]
            if !issubset(F,G) && !issubset(G,F)
                push!(relations, vars[i]*vars[j])
            end
        end
    end
    return relations
end

function relations_extended_ring(ring::MPolyRing, proper_flats::Vector{Vector{T}}, vars::Vector)  where T <: ElementType
    relations = elem_type(ring)[]
    s = length(proper_flats)
    # h_E = alpha = -x_E
    poly = ring(0)
    for index in findall(F->issubset([1], F), proper_flats)
        poly+= vars[index]
    end
    push!(relations, poly-vars[2*s+1])

    # h_F = h_E - sum x_G where G is a proper flat containing F
    for i in 1:s
        F = proper_flats[i]
        poly = ring(0)
        for index in findall(G->issubset(F,G), proper_flats)
            poly+= vars[index]
        end
        push!(relations, poly+vars[s+i]-vars[2*s+1]) #add h_F and subtract h_E
    end
    return relations
end

@doc raw"""
    augmented_chow_ring(M::Matroid)

Return an augmented Chow ring of a matroid. As described in [BHMPW22](@cite).

# Examples
```jldoctest
julia> M = fano_matroid();

julia> R = augmented_chow_ring(M);
```
"""
function augmented_chow_ring(M::Matroid)
    #This function was implemented by Fedor Glazov
    Flats = flats(M)
    sizeFlats = length(Flats)
    n = length(M)

    is_loopless(M) || error("Matroid has loops")
    sizeFlats>3 || error("Matroid has too few flats")

    proper_flats = Flats[1:sizeFlats-1]
    element_var_names = [string("y_", S) for S in M.groundset]

    flat_var_names = ["x_{" * join(S, ",") * "}" for S in proper_flats]
    var_names = vcat(element_var_names, flat_var_names)
    s = length(var_names)

    ring, vars = polynomial_ring(QQ, var_names; cached=false)
    element_vars = vars[1:n]
    flat_vars = vars[n+1:s]

    I = augmented_linear_relations(ring, proper_flats, element_vars, flat_vars, M)
    J = augmented_quadratic_relations(ring, proper_flats, element_vars, flat_vars, M)

    chow_modulus = ideal(ring, vcat(I, J))
    chow_ring, projection = quo(ring, chow_modulus)

    return chow_ring
end

function augmented_linear_relations(ring::MPolyRing, proper_flats::Vector{Vector{T}}, element_vars::Vector, flat_vars::Vector, M::Matroid)  where T <: ElementType
    n = length(M)
    relations = Vector{elem_type(ring)}(undef,n)
    i = 1
    for element in M.groundset
        relations[i] = element_vars[i]
        j = 1
        for proper_flat in proper_flats
            if !(element in proper_flat)
                relations[i] -= flat_vars[j]
            end
            j += 1
        end
        i += 1
    end
    return relations
end

function augmented_quadratic_relations(ring::MPolyRing, proper_flats::Vector{Vector{T}}, element_vars::Vector, flat_vars::Vector, M::Matroid)  where T <: ElementType
    incomparable_polynomials = quadratic_relations(ring, proper_flats, flat_vars)
    xy_polynomials = elem_type(ring)[]

    i = 1
    for element in M.groundset
        j = 1
        for proper_flat in proper_flats
            if !(element in proper_flat)
                push!(xy_polynomials, element_vars[i] * flat_vars[j])
            end
            j += 1
        end
        i += 1
    end

    return vcat(incomparable_polynomials, xy_polynomials)
end

@doc raw"""
    volume_map(M::Matroid; A::MPolyQuoRing)

Return A (normalized) function that maps the top degree component of the Chow ring to the base ring.

# Examples
The following computes the volume map of the Chow ring of the Fano matroid.
```jldoctest
julia> M = fano_matroid();

julia> R = chow_ring(M);

julia> f = volume_map(M,R);

julia> f(R[1]*R[8])
1

julia> f(R[1]^2)
-2
```
"""
function volume_map(M::Matroid, A::MPolyQuoRing)
    mflats = flats(M)
    flat = mflats[1]
    prod = one(A)
    for i in 2:length(mflats)-1
        if is_subset(flat, mflats[i])
            flat = mflats[i]
            prod*= A[i-1]
        end
    end
    return f -> is_zero(f) ? 0//coeff(lift(prod),1) : coeff(lift(f),1)//coeff(lift(prod),1)
end
