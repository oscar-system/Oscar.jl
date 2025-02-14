export field
export representation_from_weights
export torus_group

#####################
#Setting up tori for fast torus algorithm
#####################

struct TorusGroup
    field::Field
    rank::Int
    #weights::Vector{Vector{ZZRingElem}}
end

@doc raw"""
    torus_group(K::Field, m::Int)

Return the torus $(K^{\ast})^m$.

!!! note
    In the context of computing invariant rings, there is no need to deal with the group structure of a torus: The torus $(K^{\ast})^m$ is specified by just giving $K$ and $m$.
# Examples
```jldoctest
julia> T = torus_group(QQ,2)
Torus of rank 2
  over QQ
```
"""
torus_group(F::Field, n::Int) = TorusGroup(F,n)

@doc raw"""
    rank(T::TorusGroup)

Return the rank of `T`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> rank(T)
2
```
"""
rank(G::TorusGroup) = G.rank

@doc raw"""
    field(T::TorusGroup)

Return the field over which `T` is defined.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> field(T)
Rational field
```
"""
field(G::TorusGroup) = G.field

function Base.show(io::IO, G::TorusGroup)
    io = pretty(io)
    println(io, "Torus of rank ", rank(G))
    print(terse(io), Indent(), "over ", Lowercase(), field(G))
    print(io, Dedent())
end

#####################
#Setting up weights for fast torus algorithm
#####################

struct RepresentationTorusGroup
    group::TorusGroup
    weights::Vector{Vector{ZZRingElem}}
end

@doc raw"""
    representation_from_weights(T::TorusGroup, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})

Return the diagonal action of `T` with weights given by `W`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1])
Representation of torus of rank 2
  over QQ and weights 
  Vector{ZZRingElem}[[-1, 1], [-1, 1], [2, -2], [0, -1]]
```
"""
function representation_from_weights(G::TorusGroup, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})
    n = rank(G)
    V = weights_from_matrix(n,W)
    return RepresentationTorusGroup(G,V)
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

@doc raw"""
    weights(r::RepresentationTorusGroup)

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> weights(r)
4-element Vector{Vector{ZZRingElem}}:
 [-1, 1]
 [-1, 1]
 [2, -2]
 [0, -1]
```
"""
weights(R::RepresentationTorusGroup) = R.weights

@doc raw"""
    group(r::RepresentationTorusGroup)

Return the torus group represented by `r`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> group(r)
Torus of rank 2
  over QQ
```
"""
group(R::RepresentationTorusGroup) = R.group

function Base.show(io::IO, R::RepresentationTorusGroup)
    io = pretty(io)
        println(io, "Representation of torus of rank ", rank(group(R)))
        println(terse(io), Indent(), "over ", Lowercase(), field(group(R)), " and weights ")
        print(io, R.weights)
        print(io, Dedent())
end

#####################
#Setting up invariant ring for fast torus algorithm. 
#####################

@attributes mutable struct TorGroupInvarRing{FldT, PolyRingElemT, PolyRingT}
    field::FldT
    poly_ring::PolyRingT #graded

    group::TorusGroup
    representation::RepresentationTorusGroup

    fundamental::Vector{PolyRingElemT}
    presentation::MPolyAnyMap{MPolyQuoRing{PolyRingElemT}, PolyRingT, Nothing, PolyRingElemT}


    #Invariant ring of reductive group G (in representation R), no other input.
    function TorGroupInvarRing(R::RepresentationTorusGroup) #here G already contains information n and rep_mat
        n = length(weights(R))
        super_ring, _ = graded_polynomial_ring(field(group(R)), :X=>1:n)
        return TorGroupInvarRing(R, super_ring)
    end

    #to compute invariant ring ring^G where G is the reductive group of R. 
    function TorGroupInvarRing(R::RepresentationTorusGroup, ring_::MPolyDecRing)
        K = field(group(R))
        z = new{typeof(K), elem_type(ring_), typeof(ring_)}()
        n = length(weights(R))
        z.field = K
        z.poly_ring = ring_
        z.representation = R
        z.group = group(R)
        return z
    end
end

@doc raw"""
    invariant_ring(r::RepresentationTorusGroup)

Return the invariant ring of the torus group represented by `r`.

!!! note
    The creation of invariant rings is lazy in the sense that no explicit computations are done until specifically invoked (for example, by the `fundamental_invariants` function).

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r)
Invariant Ring of
graded multivariate polynomial ring in 4 variables over QQ under group action of torus of rank2
```
"""
invariant_ring(R::RepresentationTorusGroup) = TorGroupInvarRing(R)

@doc raw"""
    polynomial_ring(RT::TorGroupInvarRing)

# Examples
```jldoctest
julia> T = torus_group(QQ,2)
Torus of rank 2
  over QQ
```
"""
polynomial_ring(R::TorGroupInvarRing) = R.poly_ring

@doc raw"""
    group(RT::TorGroupInvarRing)

# Examples
```jldoctest
julia> T = torus_group(QQ,2)
Torus of rank 2
  over QQ
```
"""
group(R::TorGroupInvarRing) = R.group

@doc raw"""
    representation(RT::TorGroupInvarRing)

# Examples
```jldoctest
julia> T = torus_group(QQ,2)
Torus of rank 2
  over QQ
```
"""
representation(R::TorGroupInvarRing) = R.representation

@doc raw"""
    fundamental_invariants(RT::TorGroupInvarRing)

Return a system of fundamental invariants for `RT`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r);

julia> fundamental_invariants(RT)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 X[1]^2*X[3]
 X[1]*X[2]*X[3]
 X[2]^2*X[3]
```
"""
function fundamental_invariants(z::TorGroupInvarRing)
    if !isdefined(z, :fundamental)
        R = z.representation
        z.fundamental = torus_invariants_fast(weights(R), polynomial_ring(z))
    end
    return copy(z.fundamental)
end

function Base.show(io::IO, R::TorGroupInvarRing) 
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
            index_0 = findfirst(==(point), C)
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
                        index = findfirst(==(v), C)
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
                deleteat!(U[j], findall(==(m), U[j]))
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

#####################Invariant rings as affine algebras

@doc raw"""
    affine_algebra(RT::TorGroupInvarRing)

Return the invariant ring `RT` as an affine algebra (this amounts to compute the algebra syzygies among the fundamental invariants of `RT`).

In addition, if `A` is this algebra, and `R` is the polynomial ring of which `RT` is a subalgebra,
return the inclusion homomorphism  `A` $\hookrightarrow$ `R` whose image is `RT`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r);

julia> fundamental_invariants(RT)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 X[1]^2*X[3]
 X[1]*X[2]*X[3]
 X[2]^2*X[3]

julia> affine_algebra(RT)
(Quotient of multivariate polynomial ring by ideal (-t[1]*t[3] + t[2]^2), Hom: quotient of multivariate polynomial ring -> graded multivariate polynomial ring)
```
"""
function affine_algebra(R::TorGroupInvarRing)
  if !isdefined(R, :presentation)
    V = fundamental_invariants(R)
    s = length(V)
    weights_ = zeros(Int, s)
    for i in 1:s
        weights_[i] = total_degree(V[i])
    end
    S,_ = graded_polynomial_ring(field(group(representation(R))), :t=>1:s; weights = weights_)
    R_ = polynomial_ring(R)
    StoR = hom(S,R_,V)
    I = kernel(StoR)
    Q, StoQ = quo(S,I)
    QtoR = hom(Q,R_,V)
    R.presentation = QtoR
  end
  return domain(R.presentation), R.presentation
end
