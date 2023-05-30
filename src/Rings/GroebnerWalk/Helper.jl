#############################################
# unspecific help functions
#############################################

function ident_matrix(n::Int)
    M = zeros(Int, n, n)
    for i = 1:n
        M[i, i] = 1
    end
    return M
end

function anti_diagonal_matrix(n::Int)
    M = zeros(Int, n, n)
    for i = 1:n
        M[i, n+1-i] = -1
    end
    return M
end

# Singular.isequal depends on order of generators
function equalitytest(G::Oscar.IdealGens, K::Oscar.IdealGens)
    if length(gens(G)) != length(gens(K))
        return false
    end
    generators = Singular.gens(G)
    count = 0
    for gen in generators
        for r in Singular.gens(K)
            if gen * first(coefficients(leading_term(r, ordering=G.ord))) - r * first(coefficients(leading_term(gen, ordering=G.ord))) == 0
                count += 1
                break
            end
        end
    end
    if count == length(gens(G))
        return true
    end
    return false
end

function dot(v::Vector{Int}, w::Vector{Int})
    sum = 0
    for i = 1:length(v)
        @inbounds sum += v[i] * w[i]
    end
    return sum
end

function ordering_as_matrix(w::Vector{Int}, ord::Symbol)
    if length(w) > 2
        if ord == :lex
            return [
                w'
                ident_matrix(length(w))[1:length(w)-1, :]
            ]
        end
        if ord == :deglex
            return [
                w'
                ones(Int, length(w))'
                ident_matrix(length(w))[1:length(w)-2, :]
            ]
        end
        if ord == :degrevlex
            return [
                w'
                ones(Int, length(w))'
                anti_diagonal_matrix(length(w))[1:length(w)-2, :]
            ]
        end
        if ord == :revlex
            return [
                w'
                anti_diagonal_matrix(length(w))[1:length(w)-1, :]
            ]
        end
    else
        error("not implemented")
    end
end

function change_weight_vector(w::Vector{Int}, M::Matrix{Int})
    return [
        w'
        M[2:length(w), :]
    ]
end

function insert_weight_vector(w::Vector{Int}, M::Matrix{Int})
    return [
        w'
        M[1:length(w)-1, :]
    ]
end

function add_weight_vector(w::Vector{Int}, M::Matrix{Int})
    return [
        w'
        M
    ]
end

function ordering_as_matrix(ord::Symbol, nvars::Int)
    if ord == :lex
        return ident_matrix(nvars)
    end
    if ord == :deglex
        return [
            ones(Int, nvars)'
            ident_matrix(nvars)[1:nvars-1, :]
        ]
    end
    if ord == :degrevlex
        return [
            ones(Int, nvars)'
            anti_diagonal_matrix(nvars)[1:nvars-1, :]
        ]
    end
    if ord == :revlex
        return [
            w'
            anti_diagonal_matrix(length(w))[1:length(w)-1, :]
        ]
    else
        error("not implemented")
    end
end

function deg(p::MPolyRingElem, n::Int)
    max = 0
    for mon in Singular.monomials(p)
        ev = Singular.exponent_vectors(mon)
        sum = 0
        for e in ev
            for i = 1:n
                sum += e[i]
            end
            if (max < sum)
                max = sum
            end
        end
    end
    return max
end

function check_order_M(S::Matrix{Int}, T::Matrix{Int}, G::Oscar.IdealGens)
    (nrows, ncols) = size(T)
    return (
        (nrows, ncols) == size(S) &&
        nrows == ncols &&
        nrows == nvars(base_ring(G)) &&
        rank(T) == nvars(base_ring(G)) &&
        rank(S) == nvars(base_ring(G))
    )
end

function checkInt32(w::Vector{Int})
    for i = 1:length(w)
        if tryparse(Int32, string(w[i])) == nothing
            return false
        end
    end
    return true
end
