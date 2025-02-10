################################################################################
#
#  Homogenization routines
#  =======================
#  (used wherever mathematically necessary)
#
################################################################################

function homogenize_pre_tropicalization(I::MPolyIdeal)
    # Compute reduced Groebner basis
    G = groebner_basis(I,complete_reduction=true)

    # Construct polynomial ring for homogenization
    Kx = base_ring(I)
    Kxhx,_ = polynomial_ring(coefficient_ring(Kx),vcat([:xh],symbols(Kx)); cached=false)

    # Construct homogenization
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


function dehomogenize_post_tropicalization(Sigma::PolyhedralComplex)
    @req lineality_dim(Sigma)>0 "dehomogenizing polyhedral complex without lineality"

    ###
    # Construct hyperplane {first coord = 0}
    ###
    n = ambient_dim(Sigma)
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
    for sigma in maximal_polyhedra(Sigma)
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
    sigma = first(maximal_polyhedra(Sigma))
    sigmaDehomogenized = intersect(sigma,dehomogenisingHyperplane)
    dehomogenizedLineality = [linealityVector[2:end] for linealityVector in lineality_space(sigmaDehomogenized)]

    return polyhedral_complex(incidenceMatrixVerticesAndRays,
                              dehomogenizedVerticesAndRays,
                              collect(length(dehomogenizedVertices)+1:length(dehomogenizedVertices)+length(dehomogenizedRays)),
                              dehomogenizedLineality)
end
