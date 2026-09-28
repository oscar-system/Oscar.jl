################################################################################
#
#  Homogenization routines
#  =======================
#  (used wherever mathematically necessary)
#
################################################################################

function homogenize_pre_tropicalization(I::MPolyIdeal)
    # Compute reduced Groebner basis w.r.t. degrevlex
    Kx = base_ring(I)
    G = groebner_basis(I, ordering=degrevlex(Kx), complete_reduction=true)

    # Construct homogenization
    Kxhx,_ = polynomial_ring(coefficient_ring(Kx),vcat([:xh],symbols(Kx)); cached=false)
    Gh = Vector{elem_type(Kxhx)}(undef,length(G))
    for (i,g) in enumerate(G)
        gh = MPolyBuildCtx(Kxhx)
        d = total_degree(g)
        for (c,alpha) in coefficients_and_exponents(g)
            pushfirst!(alpha,d-sum(alpha)) # homogenize exponent vector
            push_term!(gh,c,alpha)
        end
        Gh[i] = finish(gh)
    end
    Ih = ideal(Kxhx,Gh)

    # Construct dehomogenization map
    dehomogenizationMap = hom(Kxhx,Kx,vcat(1,gens(Kx)))

    return Ih, Kxhx, dehomogenizationMap
end


function dehomogenize_post_tropicalization(TropV::TropicalVarietySupertype)

    @req lineality_dim(TropV)>0 "dehomogenizing polyhedral complex without lineality"

    ###
    # Construct hyperplane {first coord = 0}
    ###
    n = ambient_dim(TropV)
    zerothUnitRowVector = zeros(Int,1,n)
    zerothUnitRowVector[1,1] = 1
    dehomogenisingHyperplane = polyhedron((zeros(Int,0,n),zeros(Int,0)), (zerothUnitRowVector,[0]))

    ###
    # Construct matrix and incidence matrix of vertices and rays
    ###
    incidenceMatrixVertices = Vector{Int}[]
    dehomogenizedVertices = Vector{QQFieldElem}[]
    incidenceMatrixRays = Vector{Int}[]
    dehomogenizedRays = Vector{QQFieldElem}[]
    for sigma in maximal_polyhedra(TropV)
        sigmaDehomogenized = intersect(sigma,dehomogenisingHyperplane)
        incidenceVectorVertices = Int[]
        V,_ = minimal_faces(sigmaDehomogenized)
        for vertex in V
            vertex = vertex[2:end]
            i = findfirst(isequal(vertex),dehomogenizedVertices)
            if i === nothing
                push!(dehomogenizedVertices,vertex)
                push!(incidenceVectorVertices,length(dehomogenizedVertices))
            else
                push!(incidenceVectorVertices,i)
            end
        end
        push!(incidenceMatrixVertices,incidenceVectorVertices)

        incidenceVectorRays = Int[]
        R,_ = rays_modulo_lineality(sigmaDehomogenized)
        for ray in R
            ray = ray[2:end]
            i = findfirst(isequal(ray),dehomogenizedRays)
            if i === nothing
                push!(dehomogenizedRays,ray)
                push!(incidenceVectorRays,length(dehomogenizedRays))
            else
                push!(incidenceVectorRays,i)
            end
        end
        push!(incidenceMatrixRays,incidenceVectorRays)
    end

    ###
    # Concatenate vertically matrixes of vertices and rays,
    # shift incidence matrix of rays and concatenate it horizontally to incicende matrix of vertices,
    # dehomogenize generators of lineality space
    ###
    dehomogenizedVerticesAndRays = matrix(QQ,vcat(dehomogenizedVertices,dehomogenizedRays))
    incidenceMatrixRaysShifted = (x -> x .+length(dehomogenizedVertices)).(incidenceMatrixRays)
    incidenceMatrixVerticesAndRays = IncidenceMatrix([vcat(iv,ir) for (iv,ir) in zip(incidenceMatrixVertices,incidenceMatrixRaysShifted)])

    ###
    # Dehomogenize lineality space
    ###
    sigma = first(maximal_polyhedra(TropV))
    sigmaDehomogenized = intersect(sigma,dehomogenisingHyperplane)
    dehomogenizedLineality = [linealityVector[2:end] for linealityVector in lineality_space(sigmaDehomogenized)]

    SigmaDehom =  polyhedral_complex(incidenceMatrixVerticesAndRays,
                                     dehomogenizedVerticesAndRays,
                                     collect(length(dehomogenizedVertices)+1:length(dehomogenizedVertices)+length(dehomogenizedRays)),
                                     dehomogenizedLineality)

    TropVDehom = tropical_variety(SigmaDehom,multiplicities(TropV),convention(TropV))
    copy_tropical_attributes!(TropVDehom,TropV)
    return TropVDehom

end


@doc raw"""
    homogenize_pre_tropicalization(ord::MonomialOrdering, Rh::MPolyRing)

Given monomial ordering `ord` on polynomial ring `R` with variables x1, ..., xn,
and polynomial ring `Rh` with variables xh, x1, ..., xn, return an extended
monomial ordering on `Rh` such that:
- named monomial orderings are preserved (e.g. `lex(R))` becomes `lex(Rh)`),
- matrix and weight orderings are extended by prepending a zero column or
  entry for the new variable `Rh`.
"""
homogenize_pre_tropicalization(ord::MonomialOrdering, Rh::MPolyRing) = homogenize_pre_tropicalization(ord.o, Rh)

function homogenize_pre_tropicalization(o::Orderings.ProdOrdering, Rh::MPolyRing)
    return homogenize_pre_tropicalization(o.a, Rh) * homogenize_pre_tropicalization(o.b, Rh)
end

function homogenize_pre_tropicalization(o::Orderings.SymbOrdering, Rh::MPolyRing)
    oldSymbol = typeof(o).parameters[1]
    return monomial_ordering(Rh, oldSymbol)
end

function homogenize_pre_tropicalization(o::Orderings.MatrixOrdering, Rh::MPolyRing)
    newWeights = zero_matrix(ZZ, nrows(o.matrix), nvars(Rh))
    newWeights[:, 2:end] = o.matrix
    return matrix_ordering(Rh, newWeights; check=false)
end

function homogenize_pre_tropicalization(o::Orderings.WSymbOrdering, Rh::MPolyRing)
    oldSymbol = typeof(o).parameters[1]
    newWeights = vcat([0], o.weights)
    return MonomialOrdering(Rh, Orderings.WSymbOrdering(oldSymbol, collect(1:nvars(Rh)), newWeights))
end
