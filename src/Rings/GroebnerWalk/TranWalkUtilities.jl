# computes the representation of the matrixorder defined by T.
function representation_vector(G::Oscar.IdealGens, T::Matrix{Int})
    n = size(T)[1]
    M = 0
    for i = 1:n
        for j = 1:n
            temp = T[i, j]
            if M < temp
                M = temp
            end
        end
    end
    d0 = 0
    for g in Singular.gens(G)
        temp = deg(g, n)
        if d0 < temp
            d0 = temp
        end
    end
    d = M * (2 * d0^2 + (n + 1) * d0)
    w = d^(n - 1) * T[1, :]
    for i = 2:n
        w = w + d^(n - i) * T[i, :]
    end
    return w
end

# checks if Gw is an initialform of an ideal corresponding to a face of the Groebner fan with dimension less than n-1.
function inSeveralCones(Gw::Vector{T}) where {T<:MPolyRingElem}
    counter = 0
    for g in Gw
        if size(collect(Singular.coefficients(g)))[1] > 2
            return true
        end
        if size(collect(Singular.coefficients(g)))[1] == 2
            counter = counter + 1
        end
    end
    if counter > 1
        return true
    end
    return false
end
