################################################################################
#
#  Stable intersections
#
################################################################################

@doc raw"""
    stable_intersection(TropV1::TropicalVariety, TropV2::TropicalVariety)

Return the stable intersection of `TropV1` and `TropV2`.
"""
function stable_intersection(TropV1::TropicalVarietySupertype{minOrMax,true}, TropV2::TropicalVarietySupertype{minOrMax,true},
                             perturbation::Union{Vector{Int},Nothing}=nothing) where minOrMax

    @req ambient_dim(TropV1)==ambient_dim(TropV2) "different ambient dimenions of tropical varieties"
    if isnothing(perturbation)
        perturbation = rand(Int,ambient_dim(TropV1))
    end

    n = ambient_dim(TropV1)
    Sigma12 = Polyhedron[]
    mults12 = ZZRingElem[]

    # todo: check that perturbation is actually generic
    for (sigma1,m1) in maximal_polyhedra_and_multiplicities(TropV1)
        for (sigma2,m2) in maximal_polyhedra_and_multiplicities(TropV2)
            if dim(sigma1+sigma2)==n && intersect_after_perturbation(sigma1,sigma2,perturbation)
                sigma12 = intersect(sigma1, sigma2)
                i = findfirst(isequal(sigma12), Sigma12)
                if isnothing(i)
                    push!(Sigma12, intersect(sigma1, sigma2))
                    push!(mults12, m1*m2*tropical_intersection_multiplicity(sigma1,sigma2))
                else
                    mults12[i] += m1*m2*tropical_intersection_multiplicity(sigma1,sigma2)
                end
            end
        end
    end

    return tropical_variety(Sigma12,mults12,minOrMax)
end

function intersect_after_perturbation(sigma1::Polyhedron{QQFieldElem}, sigma2::Polyhedron{QQFieldElem}, perturbation::Vector{Int})

    if dim(intersect(sigma1,sigma2))<0
        return false
    end

    # construct sigma1Prime = sigma1 x RR in RR^(n+1)
    M = affine_inequality_matrix(facets(sigma1))
    M = hcat(M, zero_matrix(QQ, nrows(M), 1))
    N = affine_equation_matrix(affine_hull(sigma1))
    N = hcat(N, zero_matrix(QQ, nrows(N), 1))
    sigma1Prime = polyhedron((M[:,2:end],[-M[:,1]...]), (N[:,2:end],[-N[:,1]...]))

    # construct sigma2Prime = (sigma2 x {0}) + RR_{>=0} * (perturbation,1)
    M = affine_inequality_matrix(facets(sigma2))
    M = hcat(M, zero_matrix(QQ, nrows(M), 1))
    N = affine_equation_matrix(affine_hull(sigma2))
    N = hcat(N, zero_matrix(QQ, nrows(N), 1))
    v = zero_matrix(QQ, 1, ncols(N))
    v[1,end] = 1 # v = last unit vector
    N = vcat(N, v)
    sigma2Prime = polyhedron((M[:,2:end],[-M[:,1]...]), (N[:,2:end],[-N[:,1]...])) + convex_hull(zero_matrix(QQ, 1, ncols(N) - 1), matrix(QQ, [vcat(perturbation, [1])]))

    sigma12Prime = intersect(sigma1Prime, sigma2Prime)
    # @req codim(sigma1Prime)+codim(sigma2Prime)==codim(sigma12Prime) "perturbation "*string(perturbation)*" not generic"

    return dim(sigma12Prime)>0 && last(relative_interior_point(sigma12Prime))>0
end


# Input: B1, B2 matrices whose rows are generating sets of two euclidean linear spaces,
#               whose sum is the entire space
# Output: the tropical intersection number as defined in [Maclagan-Sturmfels, Definition 3.6.5]
# todo: rewrite function below so that it takes polyhedra as input
function tropical_intersection_multiplicity(sigma1::Polyhedron,sigma2::Polyhedron)
    B1 = kernel_basis(affine_equation_matrix(affine_hull(sigma1))[:,2:end])
    B1 = matrix(ZZ,[ numerator.(b .* lcm(denominator.(b))) for b in B1 ])
    B1 = saturate(B1)

    B2 = kernel_basis(affine_equation_matrix(affine_hull(sigma2))[:,2:end])
    B2 = matrix(ZZ,[ numerator.(b .* lcm(denominator.(b))) for b in B2 ])
    B2 = saturate(B2)

    @req ncols(B1) == ncols(B2) && nrows(B1)+nrows(B2) >= ncols(B1) "polyhedra do not span ambient space"

    snfB12 = snf(vcat(B1,B2))
    return abs(prod([snfB12[i,i] for i in 1:ncols(snfB12)]))
end



################################################################################
#
#  Stable intersections of tropical linear spaces
#
#  Todo: make this use algebraic_pluecker_vector if available
#  Problem: the +/- in the formula of speyers paper
#
################################################################################

function stable_intersection(TropL1::TropicalLinearSpace{minOrMax,true}, TropL2::TropicalLinearSpace{minOrMax,true}) where minOrMax

    plueckerIndices12 = Vector{Int}[]
    plueckerVector12 = TropicalSemiringElem{typeof(minOrMax)}[]
    d = codim(TropL1)+codim(TropL2)

    for (I1,p1) in zip(pluecker_indices(TropL1),tropical_pluecker_vector(TropL1))
        for (I2,p2) in zip(pluecker_indices(TropL2),tropical_pluecker_vector(TropL2))
            I12 = intersect(I1,I2)
            if length(I12)>d
                continue
            end

            i = searchsortedfirst(plueckerIndices12,I12)
            if i<=length(plueckerIndices12) && plueckerIndices12[i]==I12
                plueckerVector12[i] += p1*p2
            else
                insert!(plueckerIndices12,i,I12)
                insert!(plueckerVector12,i,p1*p2)
            end
        end
    end

    return tropical_linear_space(plueckerIndices12,plueckerVector12)
end
