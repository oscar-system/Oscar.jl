#################################################################
# Procedures of the fractal walk.
# The fractal walk is proposed by Amrhein & Gloor (1998).
#################################################################

# lifts the Groebner basis G to the Groebner basis w.r.t. in the Fractal Walk like it´s done in Fukuda et. al (2005).
function lift_fractal_walk(
    G::Oscar.IdealGens,
    H::Oscar.IdealGens,
    Rn::MPolyRing,
)
    R = base_ring(G)
    G.isGB = true
    G = Singular.Ideal(
        Rn,
        [
            change_ring(gen, Rn) -
            change_ring(Singular.reduce(change_ring(gen, R), G), Rn) for
            gen in Singular.gens(H)
        ],
    )
    G.isGB = true
    return G
end

# returns ´true´ if all polynomials of the given array are monomials.
function ismonomial(Gw::Vector{T}) where {T<:MPolyRingElem}
    for g in Gw
        if length(Singular.coefficients(g)) > 1
            return false
        end
    end
    return true
end

# returns ´true´ if all polynomials of the given array are binomials or less.
function isbinomial(Gw::Vector{T}) where {T<:MPolyRingElem}
    for g in Gw
        if length(Singular.coefficients(g)) > 2
            return false
        end
    end
    return true
end

# returns the next t to compute the next weight vector w(t) = w + t * (tw - w) like it´s done in Amrhein & Gloor (1998). This Method is NOT tested sufficiently.
function nextT(
    G::Oscar.IdealGens,
    w::Array{T,1},
    tw::Array{K,1},
) where {T<:Number,K<:Number}
    if (w == tw)
        return [0]
    end
    tmin = 2
    t = 0
    for g in gens(G)
        a = Singular.leading_exponent_vector(g)
        d = Singular.exponent_vectors(tail(g))
        for v in d
            frac = (dot(w, a) - dot(w, v) + dot(tw, v) - dot(tw, a))
            if frac != 0
                t = (dot(w, a) - dot(w, v)) // frac
            end
            if t > 0 && t < tmin
                tmin = t
            end
        end
    end
    if tmin <= 1
        return tmin
    else
        return [0]
    end
end

# returns the next t to compute the next weight vector w(t) = w + t * (tw - w) like it´s done in Cox, Little & O'Sheao (2005).
function next_weightfr(
    G::Oscar.IdealGens,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
    if (cweight == tweight)
        return [0]
    end
    tmin = 1
    for v in difference_lead_tail(G)
        tw = dot(tweight, v)
        if tw < 0
            cw = dot(cweight, v)
            t = cw // (cw - tw)
            if t < tmin
                tmin = t
            end
        end
    end

    # BigInt is needed to prevent overflows in the conversion of the weight vectors.
    return BigInt(numerator(tmin)) // BigInt(denominator(tmin))
end

# returns 'true' if the leading terms of G w.r.t the matrixordering T are the same as the leading terms of G w.r.t the weighted monomial ordering with weight vector of pvecs[p] (pvecs[p-1]) and the matrixordering T.
function inCone(
    G::Oscar.IdealGens,
    T::Matrix{Int},
    pvecs::Vector{Vector{Int}},
    p::Int,
)
    if p == 1
        return true
    end
    R = change_order(G.base_ring, T)
    cvzip = zip(
        Singular.gens(G),
        initials(R, Singular.gens(G), pvecs[p-1]),
        initials(R, Singular.gens(G), pvecs[p]),
    )
    for (g, in, in2) in cvzip
        if !isequal(
            Singular.leading_exponent_vector(change_ring(g, R)),
            Singular.leading_exponent_vector(in),
        ) ||
           !isequal(
            Singular.leading_exponent_vector(change_ring(g, R)),
            Singular.leading_exponent_vector(in2),
        )
            return false
        end
    end
    return true
end
