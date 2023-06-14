NF, sr2 = quadratic_field(2)
Qx, x = QQ["x"]
K, (a1, a2) = embedded_number_field([x^2 - 2, x^3 - 5], [(0, 2), (0, 2)])

# add K to the loop when hash is defined for OscarNumber
for f in (QQ, NF)

    T = elem_type(f)
    @testset "PolyhedralComplex{$T}" begin
    
        I = IncidenceMatrix([[1, 2, 3], [2, 4]])
        P = [0 0; 1 0; 0 1; 1 1]
        P2 = [0 0 0; 1 0 0; 0 1 0; 1 1 0]
        F = [4]
        L = [0 0 1]
    
        @testset "constructors" begin
        
            @test polyhedral_complex(f, I, P) isa PolyhedralComplex{T}
            @test polyhedral_complex(f, I, P, F) isa PolyhedralComplex{T}
            @test polyhedral_complex(f, I, P2, F, L) isa PolyhedralComplex{T}
            @test polyhedral_complex(f, I, P2, F, L; non_redundant = true) isa PolyhedralComplex{T}
            @test polyhedral_complex(f, I, P2, nothing, L) isa PolyhedralComplex{T}
        
        end
    
        PC = polyhedral_complex(f, I, P)
        PCF = polyhedral_complex(f, I, -P, F)
        PCFL = polyhedral_complex(f, I, P2, F, L)
        PCFLN = polyhedral_complex(f, I, P2, F, L; non_redundant = true)
        PCL = polyhedral_complex(f, I, P2, nothing, L)
     
        @test common_refinement(PC, PCF) isa PolyhedralComplex{T}
        PCR = common_refinement(PC, PCF)
     
        @test k_skeleton(PC, 1) isa PolyhedralComplex{T}
        PCK = k_skeleton(PC, 1)
     
        @testset "core functionality" begin
         
            @test ambient_dim(PC) == 2
            @test vertices(PC) isa SubObjectIterator{PointVector{T}}
            @test length(vertices(PC)) == 4
            @test point_matrix(vertices(PC)) == matrix(f, P)
            @test vertices(PC) == [[0, 0], [1, 0], [0, 1], [1, 1]]
            @test rays(PCF) isa SubObjectIterator{RayVector{T}}
            @test length(rays(PCF)) == 1
            @test rays(PCF) == [[-1, -1]]
            @test vector_matrix(rays(PCF)) == matrix(f, [-1 -1])
            @test vertices_and_rays(PCFL) isa SubObjectIterator{Union{RayVector{T}, PointVector{T}}}
            @test length(vertices_and_rays(PCFL)) == 4
            @test vertices_and_rays(PCFL) == [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
            @test typeof.(vertices_and_rays(PCFL)) == [PointVector{T}, PointVector{T}, PointVector{T}, RayVector{T}]
            @test vector_matrix(vertices_and_rays(PCFL)) == matrix(f, P2)
            @test maximal_polyhedra(PC) isa SubObjectIterator{Polyhedron{T}}
            @test length(maximal_polyhedra(PC)) == 2
            @test maximal_polyhedra(PC) == convex_hull.([f], [P[1:3, :], P[[2, 4], :]])
            @test n_maximal_polyhedra(PC) == 2
            @test is_simplicial(PC)
            @test !is_pure(PCL)
            @test dim(PCL) == 3
            @test polyhedra_of_dim(PC, 1) isa SubObjectIterator{Polyhedron{T}}
            @test length(polyhedra_of_dim(PC, 1)) == 4
            if T == QQFieldElem
                @test polyhedra_of_dim(PC, 1) == convex_hull.(T, [P[[2, 4], :], P[[1, 3], :], P[[1, 2], :], P[[2, 3], :]])
            else
                @test polyhedra_of_dim(PC, 1) == convex_hull.([f], [P[[2, 4], :], P[[1, 3], :], P[[2, 3], :], P[[1, 2], :]])
            end
            @test lineality_space(PCL) isa SubObjectIterator{RayVector{T}}
            @test length(lineality_space(PCL)) == 1
            @test lineality_space(PCL) == [L[:]]
            @test generator_matrix(lineality_space(PCL)) == matrix(QQ, L)

            @test lineality_dim(PCFL) == 1
            @test f_vector(PCL) == [0, 4, 4, 1]
            # Since there is lineality, there are no rays or vertices
            @test nrays(PCFL) == 0
            @test nvertices(PCFL) == 0
            @test npolyhedra(PCL) == 9
            @test codim(PCF) == 0
            @test is_embedded(PC)
            
            mfPCFLN = minimal_faces(PCFLN)
            @test mfPCFLN.base_points == [P2[i, :] for i in 1:3]
            rmlPCFLN = rays_modulo_lineality(PCFLN)
            @test rmlPCFLN.rays_modulo_lineality == [P2[4, :]]
            @test lineality_space(PCFLN) == [L[1, :]]
            @test vertex_indices(maximal_polyhedra(PCFLN)) == I[:, 1:3]
            @test ray_indices(maximal_polyhedra(PCFLN)) == I[:, 4:4]
            @test vertex_and_ray_indices(maximal_polyhedra(PCFLN)) == I
            @test IncidenceMatrix(maximal_polyhedra(PCFLN)) == I
            @test maximal_polyhedra(IncidenceMatrix, PCFLN) == I
         
        end

    end
end
