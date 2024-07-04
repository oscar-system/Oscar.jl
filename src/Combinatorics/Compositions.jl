
function weak_compositions(n, k)
    if k == 0
        return n == 0 ? [Int[]] : Vector{Int}[]
    end
    compos = Vector{Int}[]
    for i in 0:n
        for p in weak_compositions(n - i, k - 1)
            push!(compos, [i; p])
        end
    end
    return compos
end

function compositions(n, k)
    if k == 0
        return n == 0 ? [Int[]] : Vector{Int}[]
    end
    if n < k
        return []
    end
    compos = Vector{Int}[]
    for i in 1:n
        for p in compositions(n - i, k - 1)
            push!(compos, [i; p])
        end
    end
    return compos
end
