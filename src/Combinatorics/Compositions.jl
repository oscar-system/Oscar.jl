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
