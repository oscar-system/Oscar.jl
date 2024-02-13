#####################
#Setting up reductive groups for fast torus algorithm
#####################

mutable struct ReductiveGroupFastTorus
    field::Field
    rank::Int
    #weights::Vector{Vector{ZZRingElem}}

    function ReductiveGroupFastTorus(F::Field, n::Int)
        z = new()
        z.field = F
        z.rank = n
        #z.weights = W
        return z
    end
end

rank(G::ReductiveGroupFastTorus) = G.rank
field(G::ReductiveGroupFastTorus) = G.field

#unsure of this here.
function reductive_group(sym::Symbol, n::Int, F::Field, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})
    if sym == :torus
        return ReuctiveGroupFastTorus(F,n)
    end
end

function Base.show(io::IO, G::ReductiveGroupFastTorus)
    io = pretty(io)
    println(io, "Torus of rank ", rank(G))
    print(IOContext(io, :supercompact => true), Indent(), "over ", Lowercase(), field(G))
    print(io, Dedent())
end

#####################
#Setting up weights for fast torus algorithm
#####################

mutable struct RepresentationReductiveGroupFastTorus
    group::ReductiveGroupFastTorus
    weights::Vector{Vector{ZZRingElem}}

    function RepresentationReductiveGroupFastTorus(G::ReductiveGroupFastTorus, W::Vector{Vector{ZZRingElem}})
        z = new()
        z.group = G
        z.weights = W
        return z
    end
end

function representation_from_weights(G::ReductiveGroupFastTorus, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})
    V = weights_from_matrix(W)
    return RepresentationReductiveGroupFastTorus(G,V)
end

function weights_from_matrix(n::Int, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})
    V = Vector{Vector{ZZRingElem}}()
    if W isa Vector
        n == 1 || error("Incompatible weights")
        for i in 1:length(W)
            push!(V, [ZZRingElem(W[i])])
        end
    else
        n == ncols(W) || error("Incompatible weights")
        #assume columns = G.group[2]
        for i in 1:nrows(W)
            push!(V, [ZZRingElem(W[i,j]) for j in 1:ncols(W)])
        end
    end
    return V
end

weights(R::RepresentationReductiveGroupFastTorus) = R.weights
group(R::RepresentationReductiveGroupFastTorus) = R.group

function Base.show(io::IO, R::RepresentationReductiveGroupFastTorus)
    io = pretty(io)
        println(io, "Representation of torus of rank ", rank(group(R)))
        println(IOContext(io, :supercompact => true), Indent(), "over ", Lowercase(), field(group(R)), " and weights ")
        print(io, R.weights)
        print(io, Dedent())
end

#####################
#Setting up invariant ring for fast torus algorithm. 
#####################

mutable struct InvariantRingFastTorus
    field::Field
    poly_ring::MPolyDecRing #graded
    
    group::ReductiveGroupFastTorus
    representation::RepresentationReductiveGroupFastTorus
    
    fundamental::Vector{MPolyDecRingElem}
    
    #Invariant ring of reductive group G (in representation R), no other input.
    function InvariantRingFastTorus(R::RepresentationReductiveGroupFastTorus) #here G already contains information n and rep_mat
        n = length(weights(R))
        super_ring, __ = graded_polynomial_ring(field(group(R)), "X"=>1:n)
        return InvariantRingFastTorus(R, super_ring)
    end
    
    #to compute invariant ring ring^G where G is the reductive group of R. 
    function InvariantRingFastTorus(R::RepresentationReductiveGroupFastTorus, ring_::MPolyDecRing)
        z = new()
        n = length(weights(R))
        z.field = field(group(R))
        z.poly_ring = ring_
        z.representation = R
        z.group = group(R)
    end
end

poly_ring(R::InvariantRingFastTorus) = R.poly_ring
group(R::InvariantRingFastTorus) = R.group
representation(R::InvariantRingFastTorus) = R.representation

function fundamental_invariants(z::InvariantRingFastTorus)
    if isdefined(z, :fundamental)
        return z.fundamental
    else
        R = z.representation
        z.fundamental = torus_invariants_fast(weights(R), poly_ring(z))
        return z.fundamental
    end
end

function Base.show(io::IO, R::InvariantRingFastTorus) 
    io = pretty(io)
    println(io, "Invariant Ring of")
    print(io, Lowercase(), R.poly_ring)
    print(io, Indent(),  " under group action of torus of rank", rank(group(R)))
    print(io, Dedent())
end

##########################
#fast algorithm for invariants of tori
##########################
#Algorithm 4.3.1 from Derksen and Kemper. Computes Torus invariants without Reynolds operator.
function torus_invariants_fast(W::Vector{Vector{ZZRingElem}}, R::MPolyRing)
    #no check that length(W[i]) for all i is the same
    length(W) == ngens(R) || error("number of weights must be equal to the number of generators of the polynomial ring")
    n = length(W)
    r = length(W[1])
    #step 2
    if length(W[1]) == 1
        M = zero_matrix(ZZ, n, 1)
        for i in 1:n
            M[i,1] = W[i][1]
        end
        C1 = lattice_points(convex_hull(M))
    else
        M = zero_matrix(ZZ, 2*n, r)
        for i in 1:n
            M[i, 1:r] = 2*r*W[i]
            M[n + i, 1:r] = -2*r*W[i]   
        end
        C1 = lattice_points(convex_hull(M))
    end
    
    #get a Vector{Vector{ZZRingElem}} from Vector{PontVector{ZZRingElem}}
    C = map(Vector{ZZRingElem}, C1)
    #step 3
    S = Vector{Vector{elem_type(R)}}()
    U = Vector{Vector{elem_type(R)}}()
    index_0 = 0
    for point in C
        if is_zero(point)
            index_0 = findfirst(item -> item == point, C)
        end
        c = true
        for i in 1:n
            if point == W[i]
                push!(S, [gen(R,i)])
                push!(U, [gen(R,i)])
                c = false
                break
            end
        end
        if c == true
            push!(S, elem_type(R)[])
            push!(U, elem_type(R)[])
        end
    end
    #step 4
    count = 0
    while true
        for j in 1:length(U)
            if length(U[j]) != 0
                m = U[j][1]
                w = C[j] #weight_of_monomial(m, W)
                #step 5 - 7
                for i in 1:n
                    u = m*gen(R,i)
                    v = w + W[i]
                    if v in C
                        index = findfirst(item -> item == v, C)
                        c = true
                        for elem in S[index]
                            if is_divisible_by(u, elem)
                                c = false
                                break
                            end
                        end
                        if c == true
                            push!(S[index], u)
                            push!(U[index], u)
                        end
                    end
                end
                deleteat!(U[j], findall(item -> item == m, U[j]))
            else
                count += 1
            end
        end
        if count == length(U)
            return S[index_0]
        else
            count = 0
        end
    end
end
