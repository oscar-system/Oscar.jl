################################################################################
#
#  Stable intersections
#
################################################################################

@doc raw"""
    stable_intersection(TropV1::TropicalVariety, TropV2::TropicalVariety)

Return the stable intersection of `TropV1` and `TropV2`.

# Examples
```jldoctest
julia> R,(x,y) = tropical_semiring()[:x,:y];

julia> f1 = x+y+0;

julia> TropH1 = tropical_hypersurface(f1) # tropical line with apex (0,0)
Min tropical hypersurface

julia> f2 = x^2+2*y^2+2;

julia> TropH2 = tropical_hypersurface(f2) # tropical double line with apex (1,0)
Min tropical hypersurface

julia> TropV = stable_intersection(TropH1, TropH2)
Min tropical variety

julia> dim(TropV)
0

julia> vertices(TropV)
1-element SubObjectIterator{PointVector{QQFieldElem}}:
 [1, 0]

julia> multiplicities(TropV)
1-element Vector{ZZRingElem}:
 2
```
"""
function stable_intersection(TropV1::TropicalVarietySupertype{minOrMax,true}, TropV2::TropicalVarietySupertype{minOrMax,true}) where {minOrMax<:Union{typeof(min),typeof(max)}}

    @req ambient_dim(TropV1)==ambient_dim(TropV2) "tropical varieties must have the same ambient dimension"

    ###
    # Compute intersection and intersection multiplicities of all maximal
    # polyhedra which intersect after a generic perturbation.  To avoid the
    # need for genericity, we will use a lexicographic perturbation (see below)
    ###
    Sigma12 = Polyhedron{QQFieldElem}[]
    mults12 = ZZRingElem[]
    for (sigma1,m1) in maximal_polyhedra_and_multiplicities(TropV1)
        for (sigma2,m2) in maximal_polyhedra_and_multiplicities(TropV2)
            intersectAfterLexPerturbation, sigma12 = intersect_after_lex_perturbation_with_intersection(sigma1,sigma2)
            if intersectAfterLexPerturbation
                mult12 = m1*m2*tropical_intersection_multiplicity(sigma1,sigma2)
                i = findfirst(isequal(sigma12), Sigma12)
                if isnothing(i)
                    push!(Sigma12, sigma12)
                    push!(mults12, mult12)
                else
                    mults12[i] += mult12
                end
            end
        end
    end

    # construct the correct polyhedral complex
    if isempty(Sigma12)
        # empty polyhedral complex in the correct ambient dimension
        Sigma12 = polyhedral_complex(IncidenceMatrix(),zero_matrix(QQ,0,ambient_dim(TropV1)))
    else
        # polyhedral complex from maximal polyhedra
        Sigma12 = polyhedral_complex(Sigma12; non_redundant=true)
    end

    return tropical_variety(Sigma12,mults12,convention(TropV1))
end


# Input: sigma1 and sigma2, two lowerdimensional polyhedra in the same ambient
#   space RR^n (not checked!!!)
# Output: (true, sigma1 \cap sigma2) if sigma1 intersects sigma2 +
#   ε¹·e₁+...+εⁿ·eₙ for ε>0 sufficiently small and where e₁,...,eₙ is the basis
#   of RR^n, (false, sigma1 \cap sigma2) otherwise.
function intersect_after_lex_perturbation_with_intersection(sigma1::Polyhedron{QQFieldElem},sigma2::Polyhedron{QQFieldElem})

    # Quick test 1: return false if sigma1 and sigma2 do not intersect in the first place
    sigma12 = intersect(sigma1,sigma2)
    if dim(sigma12) < 0
        return false, sigma12
    end

    # Quick test 2: return false if the affine hulls of sigma1 and sigma2 do not
    #   span the entire space
    n = ambient_dim(sigma1)
    if n > dim(sigma1)+dim(sigma2)-dim(sigma12)
        return false, sigma12
    end

    # Quick test 3: compute relative interior point of sigma12 (also
    #   necessary for later) and return true, if it lies in the interior of
    #   sigma1 and sigma2 (the correctness of this test requires sigma1 and
    #   sigma2 spannig the entire space)
    p = relative_interior_point(sigma12)
    # if contains_in_interior(sigma1, p) && contains_in_interior(sigma2, p)
    #     return true, sigma12
    # end

    ###
    # Actual test: check that u := ε¹·e₁+...+εⁿ·eₙ is in the tangent cone
    #   C_0(sigma1-sigma2) = C_p(sigma1)-C_p(sigma2) for p some relative
    #   interior point of sigma12.
    ###
    V1, L1 = minimal_faces(sigma1)
    V2, L2 = minimal_faces(sigma2)
    R1, _ = rays_modulo_lineality(sigma1)
    R2, _ = rays_modulo_lineality(sigma2)
    V1p = V1 .- Ref(p)
    V2p = V2 .- Ref(p)
    Cp12 = positive_hull(vcat(V1p,-V2p,R1,-R2), vcat(L1,L2))

    # and check whether u is contained in C_0(sigma1-sigma2) by checking that it
    # is full-dimensional and that the sign of the first nonzero entry of its
    # inequalities is positive.
    if dim(Cp12) < n
        return false, sigma12
    end

    F12 = linear_inequality_matrix(facets(Cp12))
    for i in 1:nrows(F12)
        j = findfirst(!iszero, F12[i, :])
        if is_positive(F12[i, j])
            return false, sigma12
        end
    end

    return true, sigma12
end


# Input: B1, B2 matrices whose rows are generating sets of two euclidean linear spaces,
#               whose sum is the entire space
# Output: the tropical intersection number as defined in [Maclagan-Sturmfels, Definition 3.6.5]
# todo: rewrite function below so that it takes polyhedra as input
function tropical_intersection_multiplicity(sigma1::Polyhedron,sigma2::Polyhedron)
    B1 = kernel(affine_equation_matrix(affine_hull(sigma1))[:,2:end], side = :right)
    B1 = matrix(ZZ,[ numerator.(B1[:, i] .* lcm(denominator.(B1[:, i]))) for i in 1:ncols(B1) ])
    B1 = saturate(B1)

    B2 = kernel(affine_equation_matrix(affine_hull(sigma2))[:,2:end], side = :right)
    B2 = matrix(ZZ,[ numerator.(B2[:, i] .* lcm(denominator.(B2[:, i]))) for i in 1:ncols(B2) ])
    B2 = saturate(B2)

    @req ncols(B1) == ncols(B2) && nrows(B1)+nrows(B2) >= ncols(B1) "polyhedra do not span ambient space"

    return abs(prod(elementary_divisors(vcat(B1,B2))))
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

    d = dim(TropL1)-codim(TropL2)
    if d<=0
        # return empty polyhedral complex
        Sigma = polyhedral_complex(IncidenceMatrix(),zero_matrix(QQ,0,ambient_dim(TropL1)))
        mults = ZZRingElem[]
        return tropical_linear_space(Sigma,mults,convention(TropL1))
    end

    plueckerIndices12 = Vector{Int}[]
    plueckerVector12 = TropicalSemiringElem{minOrMax}[]
    for (I1,p1) in zip(pluecker_indices(TropL1),tropical_pluecker_vector(TropL1))
        for (I2,p2) in zip(pluecker_indices(TropL2),tropical_pluecker_vector(TropL2))
            I12 = intersect(I1,I2)
            if length(I12)!=d
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
