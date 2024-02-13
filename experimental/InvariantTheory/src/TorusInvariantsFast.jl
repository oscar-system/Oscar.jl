mutable struct ReductiveGroupFastTorus
    field::Field
    rank::Int
    weights::Vector{Vector{ZZRingElem}}

    function ReductiveGroupFastTorus(F::Field, n::Int, W::Vector{Vector{ZZRingElem}})
        z = new()
        z.field = F
        z.rank = n
        z.weights = W
        return z
    end
end

#unsure of this here.
function reductive_group(sym::Symbol, n::Int, F::Field, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})
    if sym == :torus

        return ReuctiveGroupFastTorus(F,n,V)
    end
end

function weights_from_matrix(n::Int, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})
    V = Vector{Vector{ZZRingElem}}()
    if W isa Vector
        G.group[2] == 1 || error("Incompatible weights")
        for i in 1:length(W)
            push!(V, [ZZRingElem(W[i])])
        end
    else
        G.group[2] == ncols(W) || error("Incompatible weights")
        #assume columns = G.group[2]
        for i in 1:nrows(W)
            push!(V, [ZZRingElem(W[i,j]) for j in 1:ncols(W)])
        end
    end
    return V
end