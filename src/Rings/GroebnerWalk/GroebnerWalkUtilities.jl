###############################################################
# Several Procedures for the Groebner Walk
###############################################################

# returns the next intermediate weight vector.
function next_weight(
    G::Oscar.IdealGens,
    cw::Vector{Int},
    tw::Vector{Int},
) where {K<:Number}
    tmin = 1
    for v in difference_lead_tail(G)
        tdotw = dot(tw, v)
        if tdotw < 0
            cdotw = dot(cw, v)
            t = cdotw // (cdotw - tdotw)
            if t < tmin
                tmin = t
            end
        end
    end
    return convert_bounding_vector(cw + BigInt(numerator(tmin))//BigInt(denominator(tmin)) * (tw - cw))
end

# multiplys every entry of the given weight w with 0.1 as long as it stays on the same halfspace as w.
function truncw(
    G::Oscar.IdealGens,
    w::Vector{Int},
    inw::Vector{T},
) where {T<:MPolyRingElem}
    while !checkInt32(w)
        for i = 1:length(w)
            w[i] = round(w[i] * 0.10)
        end
        w = convert_bounding_vector(w)
        if inw != initials(G, w)
            # initials are different - return unrounded weight
            return w, false
        end
    end
    # converted to Vector w of the same face
    return w, true
end

# returns the initialform of G w.r.t. the given weight vector.
function initials(
    G::Oscar.IdealGens,
    w::Vector{Int},
)
R= base_ring(G)
    inits = elem_type(R)[]
    for g in gens(G)
        inw = MPolyBuildCtx(R)
        maxw = 0
        eczip = zip(exponent_vectors(g), coefficients(g))
        for (e, c) in eczip
            tmpw = dot(w, e)
            if maxw == tmpw
                push_term!(inw, c, e)
            elseif maxw < tmpw
                inw = MPolyBuildCtx(R)
                push_term!(inw, c,e)
                maxw = tmpw
            end
        end
        h = finish(inw)
        push!(inits, h)
    end
    return inits
end

function difference_lead_tail(I::Oscar.IdealGens)
    v = Vector{Int}[]
    for g in gens(I)
        leadexpv = first(exponent_vectors(leading_term(g, ordering=I.ord)))
        tailexpvs = exponent_vectors(tail(g, ordering=I.ord))
        for e in tailexpvs
            push!(v, leadexpv .- e)
        end
    end
    return unique!(v)
end

# computes a p-pertubed vector from the matrix M.
function pertubed_vector(G::Oscar.IdealGens, M::Matrix{Int}, p::Integer)
    m = Int[]
    n = size(M, 1)
    for i = 1:p
        max = M[i, 1]
        for j = 1:n
            temp = abs(M[i, j])
            if temp > max
                max = temp
            end
        end
        push!(m, max)
    end
    msum = 0
    for i = 2:p
        msum += m[i]
    end
    maxdeg = 0
    for g in gens(G)
        td = deg(g, n)
        if (td > maxdeg)
            maxdeg = td
        end
    end
    e = maxdeg * msum + 1
    w = M[1, :] * e^(p - 1)
    for i = 2:p
        w += e^(p - i) * M[i, :]
    end
    return convert_bounding_vector(w)
end

# returns 'true' if the leading terms of G w.r.t the matrixorder T are the same as the leading terms of G w.r.t the weighted monomial order with weight vector t and matrix T.
#function inCone(G::Oscar.IdealGens, T::Matrix{Int}, t::Vector{Int})
#    R = change_order(G.base_ring, T)
#    I = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
#    cvzip = zip(Singular.gens(I), initials(R, Singular.gens(I), t))
#    for (g, ing) in cvzip
#        if !isequal(Singular.leading_exponent_vector(g), Singular.leading_exponent_vector(ing))
#            return false
#        end
#    end
#    return true
#end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G w.r.t the weighted monomial order with weight vector t and the matrix order T.
#function inCone(G::Oscar.IdealGens, t::Vector{Int})
#    cvzip = zip(Singular.gens(G), initials(base_ring(G), gens(G), t))
#    for (g, ing) in cvzip
#        if !isequal(leading_exponent_vector(g), leading_exponent_vector(ing))
#            return false
#        end
#    end
#    return true
#end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G with the current ordering.
function same_cone(G::Oscar.IdealGens, T::Matrix{Int})
    R = base_ring(G)
    ord = matrix_ordering(R, T)
    for g in gens(G)
        if leading_term(g, ordering=ord) != leading_term(g, ordering=G.ord)
            return false
        end
    end
    return true
end

# lifts the Groebner basis G to the Groebner basis w.r.t. the Ring Rn like it´s presented in Fukuda et al. (2005).
function lift(
    G::Oscar.IdealGens,
    orderingAlt::MonomialOrdering,
    H::Oscar.IdealGens,
    ordering::MonomialOrdering
)
    G = Oscar.IdealGens(
        [
            gen - Oscar.IdealGens([reduce(gen, gens(G), ordering=orderingAlt)], ordering)[1]
            for gen in gens(H)
                ],
        ordering,
        isGB = true
    )
     return G
end

# lifts the Groebner basis G to the Groebner basis w.r.t. the Ring Rn like it´s done in Collart et al. (1997).
function liftGW2(
    G::Oscar.IdealGens,
    orderingAlt::MonomialOrdering,
    inG::MPolyIdeal{T},
    H::Oscar.IdealGens,
    ordering::MonomialOrdering
) where {T<:MPolyRingElem}
    for i = 1:length(gens(H))
        q = divrem(gens(Oscar.IdealGens([H[i]], orderingAlt)), gens(inG))
        gens(H)[i] = base_ring(G)(0)
        for j = 1:length(gens(inG))
            gens(H)[i] = Oscar.IdealGens([gens(H)[i] +q[j] * gens(G)[j]], ordering)
        end
    end
    G = Oscar.IdealGens(H, ordering)
    G.isGB = true
    return G
end

# divisionalgorithm that returns q with q_1*f_1 + ... + q_s *f_s=p.
function division_algorithm(
    p::T,
    f::Vector{T},
    R::MPolyRing,
) where {T<:MPolyRingElem}
    q = Array{Singular.elem_type(R),1}(undef, length(f))
    for i = 1:length(f)
        q[i] = R(0)
    end
    while !isequal(p, R(0))
        i = 1
        div = false
        while (div == false && i <= length(f))
            b, m = divides(leading_term(p), leading_term(f[i]))
            if b
                q[i] = q[i] + m
                p = p - (m * f[i])
                div = true
            else
                i = i + 1
            end
        end
        if div == false
            p = p - leading_term(p)
        end
    end
    return q
end

# converts a vector wtemp by dividing the entries with gcd(wtemp).
function convert_bounding_vector(wtemp::Vector{T}) where {T<:Rational{BigInt}}
    w = Vector{Int}()
    g = gcd(wtemp)
    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], g)))
    end
    return w
end

# converts a vector wtemp by dividing the entries with gcd(wtemp).
function convert_bounding_vector(wtemp::Vector{T}) where {T<:Number}
    w = Vector{Int}()
    g = gcd(wtemp)
    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], g)))
    end
    return w
end

# returns a copy of the PolynomialRing I, equipped with the ordering a(cw)*ordering_M(T).
function change_order(
    R::MPolyRing,
    cw::Array{L,1},
    T::Matrix{Int},
) where {L<:Number,K<:Number}
    s = size(T)
    if s[1] == s[2]
        S = singular_poly_ring(R, weight_ordering(cw, matrix_ordering(R,T)))
    elseif s[1] - s[2] == 1
       S= singular_poly_ring(R,weight_ordering(cw,
        weight_ordering(T[1, :],
                       matrix_ordering(R,T[2:end, :]))))
    else
        S = singular_poly_ring(R, weight_ordering(cw,
        weight_ordering(R,T[1, :] ,
                       weight_ordering(T[2, :] ,
                       matrix_ordering(R,T[3:end, :])))),
        )
    end
    return S
end

# returns a copy of the PolynomialRing I, equipped with the ordering a(cw)*ordering_M(T).
function create_order(
    R::MPolyRing,
    cw::Array{L,1},
    T::Matrix{Int},
) where {L<:Number,K<:Number}
    s = size(T)
    if s[1] == s[2]
        ord = weight_ordering(cw, matrix_ordering(R,T))
    elseif s[1] - s[2] == 1
       ord= weight_ordering(cw,
        weight_ordering(T[1, :],
                       matrix_ordering(R,T[2:end, :])))
    else
        ord = weight_ordering(cw,
        weight_ordering(R,T[1, :] ,
                       weight_ordering(T[2, :] ,
                       matrix_ordering(R,T[3:end, :]))),
        )
    end
    return ord
end

# returns a copy of the PolynomialRing I, equipped with the ordering a(w)*a(t)*ordering_M(T).
function change_order(
    R::MPolyRing,
    w::Vector{Int},
    t::Vector{Int},
    T::Matrix{Int},
) where {}
    S = singular_poly_ring(R,  weight_ordering(w, weight_ordering(t,
                   matrix_ordering(R, T))))
    return S
end

# returns a copy of the PolynomialRing I, equipped with the ordering ordering_M(T).
function change_order(
    R::MPolyRing{T},
    M::Matrix{Int},
) where {T<:RingElem}
singular_poly_ring(R, matrix_ordering(R,M))
end

# recreates the polynomials p equipped with ring R.
function change_ring(p::RingElement, R::MPolyRing)
    cvzip = zip(Singular.coefficients(p), Singular.exponent_vectors(p))
    M = MPolyBuildCtx(R)
    for (c, v) in cvzip
        push_term!(M, c, v)
    end
    return finish(M)
end

# interreduces the Groebner basis G.
function interreduce_walk(G::Oscar.IdealGens) where {T<:MPolyRingElem}
    Rn = base_ring(G)
    Generator = collect(gens(G))
    I = 0
    for i = 1:length(gens(G))
        Generator[i] = reduce(Generator[i], Generator[1:end.!=i], ordering=G.ord)
    end
    return Oscar.IdealGens(Generator, G.ord)
end
