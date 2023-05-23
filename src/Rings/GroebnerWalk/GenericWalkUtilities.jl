#################################################################
# Procedures of the generic walk.
# The generic walk is proposed by Fukuda, Lauritzen & Thomas (2005).
#################################################################

# returns the initials of the polynomials w.r.t. the vector v.
function facet_initials(
    G::Oscar.IdealGens,
    lm::Vector{T},
    v::Vector{Int},
) where {T<:MPolyRingElem}
    Rn = parent(first(G))
    initials = Array{Singular.elem_type(Rn),1}(undef, 0)
    count = 1
    for g in G
        inw = Singular.MPolyBuildCtx(Rn)
        el = first(Singular.exponent_vectors(lm[count]))
        for (e, c) in
            zip(Singular.exponent_vectors(g), Singular.coefficients(g))
            if el == e || isparallel(el - e, v)
                Singular.push_term!(inw, c, e)
            end
        end
        h = finish(inw)
        push!(initials, h)
        count += 1
    end
    return initials
end

# returns the differences of the exponent vectors of the leading terms and the polynomials of the generators of I.
function difference_lead_tail(
    G::Oscar.IdealGens,
    Lm::Vector{L},
    T::Matrix{Int}
) where {L<:MPolyRingElem}
    v = Vector{Int}[]
    for i = 1:length(G)
        ltu = collect(exponent_vectors(leading_term(Lm[i],ordering= matrix_ordering(base_ring(G), T))))[1]
        for e in exponent_vectors(G[i])
            if ltu != e
                push!(v, ltu .- e)
            end
        end
    end
    return unique!(v)
end

# returns true if the vector u is parallel to the vector v.
function isparallel(u::Vector{Int}, v::Vector{Int})
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

# performs the lifting in the generic Walk like it´s proposed by Fukuda et al. (2005).
function lift_generic(
    G::Vector{T},
    Lm::Vector{T},
    H::Vector{T},
    ord::MonomialOrdering
) where {T<:MPolyRingElem}
    R = parent(first(G))
    Newlm = Array{elem_type(R),1}(undef, 0)
    liftPolys = Array{elem_type(R),1}(undef, 0)
    for g in H
        push!(Newlm, leading_term(g, ordering=ord))
        push!(liftPolys, g - reduce_walk(g, G, Lm,ord))
    end
    return liftPolys, Newlm
end

# returns all v \in V if v<0 w.r.t. the ordering represented by T and v>0 w.r.t the ordering represented by S.
function filter_by_ordering(S::Matrix{Int},T::Matrix{Int}, V::Vector{Vector{Int}})
    btz = Set{Vector{Int}}()
    for v in V
        if less_than_zero(T, v) && bigger_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end

# returns all v \in V if w<v w.r.t. the facet-preorder.
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

# computes the next vector in the generic walk.
function next_gamma(
    G::Oscar.IdealGens,
    Lm::Vector{L},
    w::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
) where {L<:MPolyRingElem}
    V = filter_by_ordering(S, T, difference_lead_tail(G, Lm, T))
    if (w != [0])
        V = filter_lf(w, S, T, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V
        if less_facet(v, minV, S, T)
            minV = v
        end
    end
    return minV
end

# tests if v>0 w.r.t. the ordering M.
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

# tests if v<0 w.r.t. the ordering M.
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

# tests if u<v w.r.t. the facet-preorder represented by the matrices S and T.
function less_facet(
    u::Vector{Int},
    v::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
)
    for i = 1:size(T, 1)
        for j = 1:size(S, 1)
            @inbounds Tuv = dot(T[i, :], u) * dot(S[j, :], v)
            @inbounds Tvu = dot(T[i, :], v) * dot(S[j, :], u)
            if Tuv != Tvu
                return Tuv < Tvu
            end
        end
    end
    return false
end

# returns the multiple m for all terms q in p with lm * m = q.
function divides_walk(
    p::MPolyRingElem,
    lm::MPolyRingElem,
    S::MPolyRing,
)
    div = false
    newpoly = MPolyBuildCtx(S)
    for term in terms(p)
        (b, c) = divides(term, lm)
        if b
            push_term!(
                newpoly,
                first(coefficients(c)),
                first(exponent_vectors(c)),
            )
            div = true
        end
    end
    return finish(newpoly), div
end

# returns p reduced by the Groebner basis G w.r.t. the leading terms Lm.
function reduce_walk(
    p::MPolyRingElem,
    G::Vector{T},
    Lm::Vector{T},
    ord::MonomialOrdering
) where {T<:MPolyRingElem}
    for i = 1:length(G)
        (q, b) = divides_walk(p, Lm[i], parent(p))
        if b
            return reduce_walk(p - (q * G[i]), G, Lm, ord)
            #return reduce_walk(submult(p,q,G[i],nvars(parent(p))),G,Lm)
        end
    end
    return p
end

# this function interreduces the Groebner basis G w.r.t. the leading terms Lm with tail-reduction.
function interreduce(
    G::Vector{T},
    Lm::Vector{T},
    ord::MonomialOrdering
) where {T<:MPolyRingElem}
    for i = 1:Singular.length(G)
        G[i] = reduce_walk(G[i], G[1:end.!=i], Lm[1:end.!=i], ord)
    end
    return G
end

#=
helping function for the reduction
function submult(
    p::Singular.T,
    q::Singular.T,
    c::Singular.T,
    dim::Int64,
) where {T<:MPolyRingElem}
    sol = MPolyBuildCtx(parent(p))
    ep = collect(Singular.exponent_vectors(p))
    eq = collect(Singular.exponent_vectors(q))
    cp = collect(Singular.coefficients(p))
    cq = collect(Singular.coefficients(q))
    ec = collect(Singular.exponent_vectors(c))
    cc = collect(Singular.coefficients(c))
    multc = Array{RingElem,1}(undef, 0)
    multe = Array{Vector{Int},1}(undef, 0)

    for m = 1:length(cq)
        for n = 1:length(cc)
            skip = false
            so = eq[m] + ec[n]
            for i = 1:length(multe)
                if multe[i] == so
                    multc[i] = multc[i] + (cq[m] * cc[n])
                    skip = true
                end
            end
            if !skip
                push!(multe, so)
                push!(multc, cq[m] * cc[n])
            end
        end
    end
    for j = 1:length(ep)
        fin = true
        for k = 1:length(multe)
            equals = true
            if ep[j] != multe[k]
                equals = false
            end
            if equals
                diff = cp[j] - multc[k]
                if diff != 0
                    Singular.push_term!(sol, diff, ep[j])
                    diff = nothing
                end
                fin = false
                multe[k][1] = -1 #delete Vector if equal
                break
            end
        end
        if fin
            Singular.push_term!(sol, cp[j], ep[j])
        end
    end
    for i = 1:length(multe)
        if multe[i][1] != -1
            Singular.push_term!(sol, -multc[i], multe[i])
        end
    end
    ep = nothing
    eq = nothing
    cp = nothing
    cq = nothing
    ec = nothing
    cc = nothing
    multc = nothing
    multe = nothing
    return Singular.finish(sol)
end
=#
