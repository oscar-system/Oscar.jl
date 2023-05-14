###############################################################
# Not tested version of the Generic Walk with different degree
# of pertubation. This implementation is not verified.
###############################################################

#=
###############################################################
# Generic Walk with choosable dergree of pertubation.
# This version is not tested yet.
###############################################################

function pgeneric_walk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    p::Int,
    infoLevel::Int,
)
    R = base_ring(G)
    Rn = change_order(G.base_ring, T)
    v = next_gamma(G, [0], S, T, p)
    Lm = [change_ring(Singular.leading_term(g), Rn) for g in gens(G)]
    G = Singular.Ideal(Rn, [change_ring(x, Rn) for x in gens(G)])

    println("generic_walk results")
    println("Crossed Cones with facetNormal: ")
    while !isempty(v)
        global counter = getCounter() + 1
        println(v)
        G, Lm = generic_step(G, Lm, v, T, R)
        v = next_gamma(G, Lm, v, S, T, p)
    end
    return Singular.interreduce(G)
end

function generic_step(
    G::Oscar.IdealGens,
    Lm::Vector{T},
    v::Vector{Int},
    T::Matrix{Int},
    R::Singular.PolyRing,
) where {T<:MPolyRingElem}

    Rn = Singular.base_ring(G)
    facet_Generators = facet_initials(G, Lm, v)
    H = Singular.std(
        Singular.Ideal(Rn, facet_Generators),
        complete_reduction = true,
    )
    H, Lm = lift_generic(G, Lm, H)
    G = interreduce(H, Lm)
    G = Singular.Ideal(Rn, G)
    G.isGB = true
    return G, Lm
end
=#

###############################################################
# Procedures
###############################################################

#=
function facet_initials(
    G::Oscar.IdealGens,
    lm::Vector{T},
    v::Vector{Int},
) where {T<:MPolyRingElem}
    Rn = base_ring(G)
    initials = Array{Singular.elem_type(Rn),1}(undef, 0)
    count = 1
    for g in Singular.gens(G)
        inw = Singular.MPolyBuildCtx(Rn)
        el = first(Singular.exponent_vectors(lm[count]))
        for (e, c) in
            zip(Singular.exponent_vectors(g), Singular.coefficients(g))
            if el == e || isParallel(el - e, v)
                Singular.push_term!(inw, c, e)
            end
        end
        h = finish(inw)
        push!(initials, h)
        count += 1
    end
    return initials
end

function difference_lead_tail(
    I::Oscar.IdealGens,
    Lm::Vector{T},
) where {T<:MPolyRingElem}
    v = Vector{Int}[]
    for i = 1:ngens(I)
        ltu = Singular.leading_exponent_vector(Lm[i])
        for e in Singular.exponent_vectors(gens(I)[i])
            if ltu != e
                push!(v, ltu .- e)
            end
        end
    end
    return unique!(v)
end

function isParallel(u::Vector{Int}, v::Vector{Int})
    count = 1
    x = 0
    for i = 1:length(u)
        if u[i] == 0
            if v[count] == 0
                count += +1
            else
                return false
            end
        else
            x = v[count] // u[i]
            count += 1
            break
        end
    end
    if count > length(v)
        return true
    end
    for i = count:length(v)
        @inbounds if v[i] != x * u[i]
            return false
        end
    end
    return true
end

function lift_generic(
    G::Oscar.IdealGens,
    Lm::Vector{T},
    H::Oscar.IdealGens,
) where {T<:MPolyRingElem}
    Rn = base_ring(G)
    Newlm = Array{Singular.elem_type(Rn),1}(undef, 0)
    liftPolys = Array{Singular.elem_type(Rn),1}(undef, 0)
    for g in Singular.gens(H)
        r = modulo(g, gens(G), Lm)
        diff = g - r
        if diff != 0
            push!(Newlm, Singular.leading_term(g))
            push!(liftPolys, diff)
        end
    end
    return liftPolys, Newlm
end

function filter_btz(S::Matrix{Int}, V::Vector{Vector{Int}})
    btz = Set{Vector{Int}}()
    for v in V
        if bigger_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end

function filter_ltz(S::Matrix{Int}, V::Set{Vector{Int}})
    btz = Set{Vector{Int}}()
    for v in V
        if less_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end
function filter_lf(
    w::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
    V::Set{Vector{Int}},
)
    btz = Set{Vector{Int}}()
    for v in V
        if less_facet(w, v, S, T)
            push!(btz, v)
        end
    end
    return btz
end

function next_gamma(
    G::Oscar.IdealGens,
    Lm::Vector{T},
    w::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
    p::Int,
) where {T<:MPolyRingElem}
    V = filter_btz(S, difference_lead_tail(G, Lm))
    V = filter_ltz(T, V)
    if (w != [0])
        V = filter_lf(w, S, T, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V
        if less_facet(v, minV, S, T, p)
            minV = v
        end
    end
    return minV
end

function next_gamma(
    G::Oscar.IdealGens,
    w::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
    p::Int,
)
    V = filter_btz(S, difference_lead_tail(G))
    V = filter_ltz(T, V)
    if (w != [0])
        V = filter_lf(w, S, T, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V
        if less_facet(v, minV, S, T, p)
            minV = v
        end
    end
    return minV
end

function bigger_than_zero(M::Matrix{Int}, v::Vector{Int})
    nrows, ncols = size(M)
    for i = 1:nrows
        d = 0
        for j = 1:ncols
            @inbounds d += M[i, j] * v[j]
        end
        if d != 0
            return d > 0
        end
    end
    return false
end

function less_than_zero(M::Matrix{Int}, v::Vector{Int})
    nrows, ncols = size(M)
    for i = 1:nrows
        d = 0
        for j = 1:ncols
            @inbounds d += M[i, j] * v[j]
        end
        if d != 0
            return d < 0
        end
    end
    return false
end

function less_facet(
    u::Vector{Int},
    v::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
    d::Int,
)
    for i = 1:d
        for j = 1:d
            @inbounds Tuv = dot(T[i, :], u) * dot(S[j, :], v)
            @inbounds Tvu = dot(T[i, :], v) * dot(S[j, :], u)
            if Tuv != Tvu
                return Tuv < Tvu
            end
        end
    end
    return false
end

function modulo(
    p::Singular.spoly,
    G::Vector{T},
    Lm::Vector{T},
) where {T<:MPolyRingElem}
    for i = 1:length(G)
        (q, b) = dividesGW(p, Lm[i], parent(p))

        if b
            return modulo(p - (q * G[i]), G, Lm)
        end
    end
    return p
end

function interreduce(
    G::Vector{T},
    Lm::Vector{T},
) where {T<:MPolyRingElem}
    Rn = parent(first(G))
    for i = 1:Singular.length(G)
        gensrest = Array{Singular.elem_type(Rn),1}(undef, 0)
        Lmrest = Array{Singular.elem_type(Rn),1}(undef, 0)
        for j = 1:length(G)
            if i != j
                push!(gensrest, G[j])
                push!(Lmrest, Lm[j])
            end
        end
        G[i] = modulo(G[i], gensrest, Lmrest)
    end
    return G
end
=#
