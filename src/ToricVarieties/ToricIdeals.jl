@doc Markdown.doc"""
    toric_ideal_binomial_generators(pts::Matrix{Int})

Get the exponent vectors corresponding to the generators of the toric ideal
coming from the affine relations between the point `pts`.
"""
function toric_ideal_binomial_generators(pts::AbstractMatrix)
    eq = homogenize(transpose(pts), 0)
    ineq = Polymake.common.zero_matrix(0,ncols(eq))
    pmPolytope = Polymake.polytope.Polytope(INEQUALITIES=ineq, EQUATIONS=eq)
    result = (pmPolytope.LATTICE_POINTS_GENERATORS)[3]
    return result[:, 2:ncols(result)]
end



@doc Markdown.doc"""
    toric_ideal_binomial_generators(antv::AffineNormalToricVariety)

Get the exponent vectors corresponding to the generators of the toric ideal
associated to the affine normal toric variety `antv`.

# Examples
Take the cyclic quotient singularity corresponding to the pair of integers
`(2,5)`.
```jldoctest
julia> C = Oscar.positive_hull([-2 5; 1 0])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> toric_ideal_binomial_generators(antv)
pm::Matrix<long>
-1 -1 2 1
-1 0 3 -1
0 -1 -1 2
```
"""
function toric_ideal_binomial_generators(antv::AffineNormalToricVariety)
    cone = Cone(pm_ntv(antv).WEIGHT_CONE)
    return toric_ideal_binomial_generators(hilbert_basis(cone).m)
end
export toric_ideal_binomial_generators
toric_ideal_binomial_generators(ntv::NormalToricVariety) = toric_ideal_binomial_generators(AffineNormalToricVariety(ntv))


function binomial_exponents_to_ideal(binoms::AbstractMatrix)
    nvars = ncols(binoms)
    R, x = PolynomialRing(QQ, "x" => 1:nvars)
    terms = Vector{fmpq_mpoly}(undef, nrows(binoms))
    for i in 1:nrows(binoms)
        binom = binoms[i, :]
        xpos = one(R)
        xneg = one(R)
        for j in 1:nvars
            if binom[j] < 0
                xneg = xneg * x[j]^(-binom[j])
            elseif binom[j] > 0
                xpos = xpos * x[j]^(binom[j])
            end
        end
        terms[i] = xpos-xneg
    end
    return ideal(terms)
end


@doc Markdown.doc"""
    toric_ideal(pts::AbstractMatrix)

Return the toric ideal generated from the affine relations between the points
`pts`.
"""
function toric_ideal(pts::AbstractMatrix)
    binoms = toric_ideal_binomial_generators(pts)
    return binomial_exponents_to_ideal(binoms)
end


@doc Markdown.doc"""
    toric_ideal(antv::AffineNormalToricVariety)

Return the toric ideal defining the affine normal toric variety.
"""
function toric_ideal(antv::AffineNormalToricVariety)
    binoms = toric_ideal_binomial_generators(antv)
    return binomial_exponents_to_ideal(binoms)
end
export toric_ideal
toric_ideal(ntv::NormalToricVariety) = toric_ideal(AffineNormalToricVariety(ntv))




