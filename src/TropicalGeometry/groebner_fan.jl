###
# Computing Groebner fans in Oscar
# ================================
#
# For a detailed writeup of the algorithm please see Anders Jensen's PhD dissertation
#   "Algorithmic Aspects of Groebner Fans and Tropical Varieties"
#
# In tropical terminology: Groebner fans contain tropicalizations in the constant
#   coefficient case. The convention is always max, as in OSCAR leading monomials
#   under weighted orderings are the ones with maximal weighted degree.
#
###

# Computes the space of weight vectors w.r.t. which G is weighted homogeneous
@doc raw"""
    homogeneity_space(G::Vector{<:MPolyRingElem})

Return a `Cone` that is the space of all weight vectors under which all elements of `G` are weighted homogeneous.  The cone will be a linear subspace without any rays.

!!! note
    If `G` is a reduced Groebner basis of an ideal `I` with respect a global
    ordering, then the returned cone is the space of all weight vectors with
    respect to which `I` is weighted homogeneous.

# Examples
```jldoctest
julia> Qx,(x1,x2,x3) = polynomial_ring(QQ,[:x1,:x2,:x3])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> G1 = [x1+x2,x2+x3,x3]
3-element Vector{QQMPolyRingElem}:
 x1 + x2
 x2 + x3
 x3

julia> rays_modulo_lineality(homogeneity_space(G1))
(rays_modulo_lineality = RayVector{QQFieldElem}[], lineality_basis = RayVector{QQFieldElem}[[1, 1, 1]])

julia> G2 = [x1^2+x2,x2^2+x3,x3]
3-element Vector{QQMPolyRingElem}:
 x1^2 + x2
 x2^2 + x3
 x3

julia> rays_modulo_lineality(homogeneity_space(G2))
(rays_modulo_lineality = RayVector{QQFieldElem}[], lineality_basis = RayVector{QQFieldElem}[[1//4, 1//2, 1]])

julia> G3 = [x1^2+x2,x2^2+x3,x3+1]
3-element Vector{QQMPolyRingElem}:
 x1^2 + x2
 x2^2 + x3
 x3 + 1

julia> rays_modulo_lineality(homogeneity_space(G3))
(rays_modulo_lineality = RayVector{QQFieldElem}[], lineality_basis = RayVector{QQFieldElem}[])

```
"""
function homogeneity_space(G::Vector{<:MPolyRingElem})
    inequalities = QQMatrix(0,0)
    equations = Vector{Vector{Int}}()
    for g in G
        alpha,tail_exponents = Iterators.peel(exponents(g))
        for beta in tail_exponents
            push!(equations,beta-alpha)
        end
    end
    return cone_from_inequalities(inequalities,matrix(QQ,equations))
end



@doc raw"""
    homogeneity_vector(G::Vector{<:MPolyRingElem})

Return a positive `Vector{ZZRingElem}` under which every element of `G` is weighted
homogeneous. If no such vector exists, return `nothing`.

# Note
Suppose `G` is the reduced Groebner basis of an ideal `I` with respect to a global ordering. If a `Vector{ZZRingElem}` is returned, then `I` is weighted homogeneous with respect to it.  If `nothing` is returned, then `I` is not weighted homogeneous with respect to any positive weight vector.

# Examples
```jldoctest
julia> Qx,(x1,x2,x3) = polynomial_ring(QQ,[:x1,:x2,:x3])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> G1 = [x1+x2,x2+x3,x3]
3-element Vector{QQMPolyRingElem}:
 x1 + x2
 x2 + x3
 x3

julia> homogeneity_vector(G1)
3-element Vector{ZZRingElem}:
 1
 1
 1

julia> G2 = [x1^2+x2,x2^2+x3,x3]
3-element Vector{QQMPolyRingElem}:
 x1^2 + x2
 x2^2 + x3
 x3

julia> homogeneity_vector(G2)
3-element Vector{ZZRingElem}:
 1
 2
 4

julia> G3 = [x1^2+x2,x2^2+x3,x3+1]
3-element Vector{QQMPolyRingElem}:
 x1^2 + x2
 x2^2 + x3
 x3 + 1

julia> homogeneity_vector(G3) === nothing
true

```
"""
function homogeneity_vector(G::Vector{<:MPolyRingElem})
    homogeneitySpace = homogeneity_space(G)

    # return nothing if homogeneitySpace is the origin
    # (which will often be the case)
    if dim(homogeneitySpace)<1
        return nothing
    end

    # test if homogeneitySpace contains the ones vector
    n = ambient_dim(homogeneitySpace)
    homogeneityVector = ones(ZZRingElem,n)
    if homogeneityVector in homogeneitySpace
        return homogeneityVector
    end

    # try finding any positive vector inside of it
    positiveOrthant = positive_hull(identity_matrix(ZZ,n))
    posHomogeneitySpace = intersect(homogeneitySpace,positiveOrthant)
    homogeneityVector = unique_identifying_point(posHomogeneitySpace)
    if isnothing(findfirst(<=(0), homogeneityVector))
        return homogeneityVector
    end

    # return nothing if everything fails
    return nothing
end



@doc raw"""
    maximal_groebner_cone(G::Oscar.IdealGens{<:MPolyRingElem}; homogeneityWeight::Union{Nothing,Vector{ZZRingElem}}=nothing)

Return the maximal Groebner cone of a Groebner basis `G`, i.e., the closure of all weight vectors with respect to whose weighted ordering `G` is a Groebner basis (independent of tie-breaker).

If `homogeneityWeight==nothing`, assumes that `G` is not quasi-homogeneous, i.e. not homogeneous with respect to any positive weight vector, and returns a cone inside the positive orthant.

If `homogeneityWeight!=nothing`, assumes that `G` is quasi-homogeneous with respect to it and returns a cone whose lineality space contains `homogeneityWeight`.

If `homogeneityWeight` is unspecified, will test if `G` is quasi-homogeneous, and behave accordingly.

# Examples
```jldoctest
julia> Qx,(x1,x2,x3) = polynomial_ring(QQ,[:x1,:x2,:x3])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = ideal([x1,x2+x3])
Ideal generated by
  x1
  x2 + x3

julia> G = groebner_basis(I,ordering=lex(Qx))
Gröbner basis with elements
  1: x2 + x3
  2: x1
with respect to the ordering
  lex([x1, x2, x3])

julia> rays_modulo_lineality(maximal_groebner_cone(G))
(rays_modulo_lineality = RayVector{QQFieldElem}[[0, 0, -1]], lineality_basis = RayVector{QQFieldElem}[[1, 0, 0], [0, 1, 1]])

```
"""
function maximal_groebner_cone(G::Oscar.IdealGens{<:MPolyRingElem})
    @req is_groebner_basis(G) "input not a Groebner basis"
    ord = ordering(G)
    G = collect(G)
    homogeneityWeight = homogeneity_vector(G)
    return maximal_groebner_cone(G,ord,homogeneityWeight)
end

function maximal_groebner_cone(G::Oscar.IdealGens{<:MPolyRingElem}, homogeneityWeight::Union{Vector{ZZRingElem},Nothing})
    @req is_groebner_basis(G) "input not a Groebner basis"
    ord = ordering(G)
    G = collect(G)
    return maximal_groebner_cone(G,ord,homogeneityWeight)
end

function maximal_groebner_cone(G::Vector{<:MPolyRingElem}, ord::MonomialOrdering, homogeneityWeight::Nothing)
    # calling with `nothing` as homogeneityWeight
    # means that C needs to be restricted to the positive orthant
    C = maximal_groebner_cone_extended(G,ord)
    positiveOrthant = positive_hull(identity_matrix(ZZ,ambient_dim(C)))
    return intersect(C,positiveOrthant)
end

function maximal_groebner_cone(G::Vector{<:MPolyRingElem}, ord::MonomialOrdering, homogeneityWeight::Vector{ZZRingElem})
    # calling with a positive homogeneityWeight
    # means that C does not be restricted to the positive orthant
    return maximal_groebner_cone_extended(G,ord)
end

function maximal_groebner_cone_extended(G::Vector{<:MPolyRingElem}, ord::MonomialOrdering)
    # iterate over all elements of G and construct the inequalities
    inequalities = Vector{Vector{Int}}()
    for g in G
        alpha,tail_exponents = Iterators.peel(exponents(g,ordering=ord))
        for beta in tail_exponents
            push!(inequalities,beta-alpha)
        end
    end
    # if there are none, which means G is monomial, return the entire ambient space
    if isempty(inequalities)
        return cone_from_inequalities(QQMatrix(0,ngens(base_ring(ord))))
    end
    return cone_from_inequalities(matrix(QQ,inequalities))
end



# Computes the unique primitive vector that lies on the sum of all rays and is orthogonal to the lineality space.
# If the cone has no rays, returns the zero vector.
# Serves as a unique identifier for cones in a polyhedral fan that is compatible with respect to symmetries which permute the coordinates.
function unique_identifying_point(C::Cone)
    # compute rays modulo lineality space
    R,_ = rays_modulo_lineality(C)
    # if there are none, return zero vector
    if isempty(R)
        return zeros(ZZRingElem,ambient_dim(C))
    end
    # otherwise compute the sum of rays, cast it to a rational vector
    # and return a primitive integer vector
    pt = Vector{QQFieldElem}(sum(R))
    return numerator.(pt .* lcm(denominator.(pt)))
end



# Returns an interior point on each facet, respectively.
# We use unique_identifying_points because we also use them to track which facets have been traversed
function interior_facet_points(C::Cone)
    return unique_identifying_point.(faces(C,dim(C)-1))
end



# Returns a primitive outer normal vector for each facet.
function outer_normal_vectors(C::Cone)
    outerNormalMatrix = linear_inequality_matrix(facets(C))
    normals = Vector{Vector{ZZRingElem}}()
    for i in 1:nrows(outerNormalMatrix)
        v = [outerNormalMatrix[i,:]...]
        push!(normals,numerator.(v .* lcm(denominator.(v))))
    end
    return normals
end



# Returns the initial form of polynomial g with respect to weight vector u.
# Requires a monomial ordering that is compatible with respect to u, i.e., the leading monomial is of highest weighted degree.
function initial(g::MPolyRingElem, ordering::MonomialOrdering, u::Vector{ZZRingElem})
    lt,tail = Iterators.peel(terms(g,ordering=ordering));
    d = dot(u,ZZRingElem.(leading_exponent_vector(lt)))
    initial_terms = [s for s in tail if dot(u,ZZRingElem.(leading_exponent_vector(s)))==d]
    push!(initial_terms,lt)
    return sum(initial_terms)
end

# Applies the function above to each polynomial in a list of polynomials.
function initial(G::Vector{<:MPolyRingElem}, ord::MonomialOrdering, u::Vector{ZZRingElem})
    return initial.(G,Ref(ord),Ref(u))
end



# Assumes:
# - G is a Groebner basis with respect to ord
# Returns:
# - a reduced Groebner basis with respect to ord
function interreduce(G::Vector{<:MPolyRingElem}, ord::MonomialOrdering)
    # sort G by its leading monomials (result is smallest first)
    # then reduce G[i] by G[1:i-1] for i=length(G),...,2
    sort!(G,
          by=g->leading_monomial(g,ordering=ord),
          order=ord)
    for i in 2:length(G)
        G[i] = reduce(G[i],G[1:i-1],ordering=ord,complete_reduction=true)
    end
    return G
end



# Assumes:
# - Hnew is a Groebner basis of the initial ideal w.r.t. ordH
# - G is a Groebner basis of the original ideal w.r.t. ordG
# - ordH and ordG are adjacent orderings
# Returns:
# - Gnew, a Groebner basis of the original ideal w.r.t. ordH
function groebner_lift(Hnew::Vector{<:MPolyRingElem},
                       ordH::MonomialOrdering,
                       G::Vector{<:MPolyRingElem},
                       ordG::MonomialOrdering)
    # lift to groebner basis
    Gnew = Hnew - reduce(Hnew,G,ordering=ordG,complete_reduction=true)
    # interreduce to make reduced
    Gnew = interreduce(Gnew,ordH)
    return Gnew
end



# Assumes:
# - G is a Groebner basis with respect to ordG
# - u is an interior facet point of the maximal Groebner cone around ordG
# - v is an outer normal vector of the facet containing u
# Returns:
# - the reduced Groebner basis with respect to the adjacent ordering
#     in direction v
function groebner_flip(G::Vector{<:MPolyRingElem},
                       ordG::MonomialOrdering,
                       homogeneityVector::Union{Vector{ZZRingElem},Nothing},
                       interior_facet_point::Vector{ZZRingElem},
                       outer_normal_vector::Vector{ZZRingElem})
    R = parent(first(G))
    n = ngens(R)
    adjacentOrdering = groebner_flip_adjacent_ordering(R,homogeneityVector,interior_facet_point,outer_normal_vector)
    H = initial(G,ordG,interior_facet_point)
    Hnew = collect(groebner_basis(ideal(H),ordering=adjacentOrdering,complete_reduction=true))
    Gnew = groebner_lift(Hnew,adjacentOrdering,G,ordG)
    return Gnew,adjacentOrdering
end

# helper functions for groebner_flip
function groebner_flip_adjacent_ordering(R::MPolyRing,
                                         homogeneityVector::Vector{ZZRingElem},
                                         interior_facet_point::Vector{ZZRingElem},
                                         outer_normal_vector::Vector{ZZRingElem})
    return weight_ordering(Int.(homogeneityVector),
                           weight_ordering(Int.(interior_facet_point),
                                           weight_ordering(Int.(outer_normal_vector),
                                                           invlex(R))))
end

function groebner_flip_adjacent_ordering(R::MPolyRing,
                                         homogeneityVector::Nothing,
                                         interior_facet_point::Vector{ZZRingElem},
                                         outer_normal_vector::Vector{ZZRingElem})
    return weight_ordering(Int.(interior_facet_point),
                           weight_ordering(Int.(outer_normal_vector),
                                           invlex(R)))
end

@doc raw"""
    groebner_fan(I::MPolyIdeal; return_groebner_bases::Bool=false, return_orderings::Bool=false, return_initial_ideals::Bool=false, marked_groebner_bases::Bool=false, verbose_level::Int=0)

Return a `PolyhedralFan` representing the Groebner fan of `I`, where is a multivariate polynomial ideal over a field.

If `verbose_level` is positive, also print how many cones have been computed during the traversal.

If `I` is not weighted homogeneous with respect to a positive weight vector, the Groebner fan will be restricted to the positive orthant.  Otherwise, the Groebner fan will span the entire space and have a non-trivial lineality space.

If any of `return_interior_points`, `return_groebner_bases`, `return_orderings` or `return_initial_ideals`
is set to true, also return a list with the corresponding information for each cone.

If `return_interior_points==true`, this includes the interior points of the maximal cones.

If `return_groebner_bases==true`, above list includes the Groebner bases for those cones.  Their union will be a universal Groebner basis.
If additionally `marked_groebner_bases==true`, the values the Groebner bases for those cones together with the leading term for each generator.  

If `return_orderings==true`, above list includes the monomial orderings for those cones.  These orderings are suboptimal and hence it is generally recommended to create new orderings with the interior points.  However they do contain information on how the fan was traversed.

If `return_initial_ideals==true`, above list includes the initial ideals for those cones.

# Examples
```jldoctest
julia> Qx,(x1,x2,x3) = polynomial_ring(QQ,[:x1,:x2,:x3])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x1, x2, x3])

julia> I = ideal([x1,x2+x3])
Ideal generated by
  x1
  x2 + x3

julia> SigmaI = groebner_fan(I)
Polyhedral fan in ambient dimension 3

julia> SigmaI,output = groebner_fan(I,return_groebner_bases=true,return_orderings=true)
2-element Vector{Any}:
 Polyhedral fan in ambient dimension 3
 Any[Any[QQMPolyRingElem[x1, x2 + x3], matrix_ordering([x1, x2, x3], [1 1 1])*matrix_ordering([x1, x2, x3], [0 0 0])*matrix_ordering([x1, x2, x3], [0 -1 1])*invlex([x1, x2, x3])], Any[QQMPolyRingElem[x2 + x3, x1], degrevlex([x1, x2, x3])]]

```
"""
function groebner_fan(I::MPolyIdeal; 
                      return_interior_points::Bool=false,
                      return_groebner_bases::Bool=false, 
                      return_orderings::Bool=false, 
                      return_initial_ideals::Bool=false, 
                      marked_groebner_bases::Bool=false, 
                      verbose_level::Int=0)
    ###
    # Preparation:
    #   Test whether the ideal is weighted homogeneous with respect to a positive weight vector
    #   Preference is given to the ones vector for the sake of optimisation
    ###
    G = groebner_basis(I,complete_reduction=true)
    ord = ordering(G)
    G = collect(G)
    homogeneityWeight = homogeneity_vector(G)


    ###
    # Starting cone
    ###
    C = maximal_groebner_cone(G,ord,homogeneityWeight)
    workingList = [(G,ord,C,unique_identifying_point(C))] # list of Groebner cones whose neighbors may be unknown
    finishedList = typeof(workingList)()     # list of Groebner cones whose neighbors are known
    finishedFacets = Vector{Vector{ZZRingElem}}()  # list of interior facet points whose facet has been traversed


    ###
    # Fan traversal
    ###
    while !isempty(workingList)
        if verbose_level>0
            println("#workingList / #finishedList: ",length(workingList)," / ",length(finishedList))
        end
        # take Groebner cone from workingList
        G,ord,C,pt = pop!(workingList)

        # enumerate its facets
        for (u,v) in zip(interior_facet_points(C),outer_normal_vectors(C))
            # check whether facet is supposed to be traversed
            # (i.e., if ideal is not quasi-homogeneous,
            #  check whether facet lies inside the positive orthant)
            if homogeneityWeight === nothing && !reduce(&,isless.(0,u))
                continue
            end

            # check whether facet has already been traversed.
            i = searchsortedfirst(finishedFacets,u)
            if i<=length(finishedFacets) && finishedFacets[i]==u
                continue
            end
            insert!(finishedFacets,i,u)

            # compute adjacent Groebner cone
            Gprime,ordPrime = groebner_flip(G,ord,homogeneityWeight,u,v)
            Cprime = maximal_groebner_cone(Gprime,ordPrime,homogeneityWeight)
            ptPrime = unique_identifying_point(Cprime)

            # check whether adjacent cone already in finishedList
            if insorted((0,0,0,ptPrime),finishedList,by=entry->entry[4])
                continue
            end

            # check whether adjacent cone already in workingList
            i = searchsortedfirst(workingList,(0,0,0,ptPrime),by=entry->entry[4])
            if i <= length(workingList) && workingList[i][4] == ptPrime
                continue
            end

            # insert into workingList in previously computed position
            insert!(workingList,i,(Gprime,ordPrime,Cprime,ptPrime))
        end

        # put Groebner cone into finishedList
        i = searchsortedfirst(finishedList,(0,0,0,pt),by=entry->entry[4])
        insert!(finishedList,i,(G,ord,C,pt))
    end


    ###
    # Post processing
    ###
    # construct polyhedral fan and return it if nothing else was required
    Sigma = polyhedral_fan(getindex.(finishedList,3);non_redundant=true)
    if !return_interior_points && !return_groebner_bases && !return_orderings && !return_initial_ideals
        return Sigma
    end

    # otherwise return what was required
    output = []

    for (GB,ord,_,pt) in finishedList
        output_C = []

        if return_interior_points == true
            push!(output_C,pt)
        end
        if return_groebner_bases == true
            if marked_groebner_bases == true
                push!(output_C,(f,leading_term(f;ordering=ord)))
            else
                push!(output_C,GB)
            end
        end
        if return_orderings == true
            push!(output_C,ord)
        end
        if return_initial_ideals == true
            initial_ideal = leading_ideal(I;ordering=ord)
            push!(output_C,initial_ideal)
        end

        push!(output, output_C)
    end
    return [Sigma, output]
end
