@testset "PolyhedralFan{$T}" for T in [QQFieldElem, nf_elem]
    C0 = cube(T, 2)
    @test normal_fan(C0) isa PolyhedralFan{T}
    NFsquare = normal_fan(C0)
    R = [1 0 0; 0 0 1]
    L = [0 1 0]
    Cone4 = positive_hull(T, R)
    Cone5 = positive_hull(T, [1 0 0; 0 1 0])

    @test polyhedral_fan([Cone4, Cone5]) isa PolyhedralFan{T}
    F0 = polyhedral_fan([Cone4, Cone5])
    I3 = [1 0 0; 0 1 0; 0 0 1]
    incidence1 = IncidenceMatrix([[1,2],[2,3]])
    incidence2 = IncidenceMatrix([[1,2]])
    @test polyhedral_fan(T, I3, incidence1) isa PolyhedralFan{T}
    F1 = polyhedral_fan(T, I3, incidence1)
    F1NR = polyhedral_fan(T, I3, incidence1; non_redundant = true)
    @test polyhedral_fan(T, I3, incidence1) isa PolyhedralFan{T}
    F2 = polyhedral_fan(T, R, L, incidence2)
    F2NR = polyhedral_fan(T, R, L, incidence2; non_redundant = true)

    @testset "core functionality" begin
        if T == QQFieldElem
            @test is_smooth(NFsquare)
            @test vector_matrix(rays(NFsquare)) == matrix(QQ, [1 0; -1 0; 0 1; 0 -1])
        else
            @test vector_matrix(rays(NFsquare)) == [1 0; -1 0; 0 1; 0 -1]
        end
        @test rays(NFsquare) isa SubObjectIterator{RayVector{T}}
        @test length(rays(NFsquare)) == 4
        @test rays(NFsquare) == [[1, 0], [-1, 0], [0, 1], [0, -1]]
        @test is_regular(NFsquare)
        @test is_complete(NFsquare)
        @test !is_complete(F0)
        @test length(rays(F0)) == 3
        @test nrays(F1) == 3
        @test dim(F1) == 2
        @test ambient_dim(F1) == 3
        @test nrays(F2) == 0
        @test lineality_dim(F2) == 1
        RMLF2 = rays_modulo_lineality(F2)
        @test length(RMLF2[:rays_modulo_lineality]) == 2
        @test maximal_cones(F1) isa SubObjectIterator{Cone{T}}
        @test dim.(maximal_cones(F1)) == [2,2]
        @test ray_indices(maximal_cones(F1)) == incidence1
        @test IncidenceMatrix(maximal_cones(F1)) == incidence1
        @test maximal_cones(IncidenceMatrix, F1) == incidence1
        @test n_maximal_cones(F1) == 2
        @test lineality_space(F2) isa SubObjectIterator{RayVector{T}}
        if T == QQFieldElem
            @test generator_matrix(lineality_space(F2)) == matrix(QQ, L)
            @test matrix(QQ, lineality_space(F2)) == matrix(QQ, L)
        else
            @test generator_matrix(lineality_space(F2)) == L
        end
        @test length(lineality_space(F2)) == 1
        @test lineality_space(F2) == [L[:]]
        @test cones(F2, 2) isa SubObjectIterator{Cone{T}}
        @test size(cones(F2, 2)) == (2,)
        @test lineality_space(cones(F2, 2)[1]) == [[0, 1, 0]]
        @test rays.(cones(F2, 2)) == [[], []]
        @test isnothing(cones(F2, 1))
        @test ray_indices(cones(F1, 2)) == incidence1
        @test IncidenceMatrix(cones(F1, 2)) == incidence1
        @test cones(IncidenceMatrix, F1, 2) == incidence1

        II = ray_indices(maximal_cones(NFsquare))
        NF0 = polyhedral_fan(rays(NFsquare), II)
        @test nrays(NF0) == 4
        FF0 = face_fan(C0)
        @test nrays(FF0) == 4
        if T == QQFieldElem
            @test !is_smooth(FF0)
        end
        @test f_vector(NFsquare) == [4, 4]
        @test rays(F1NR) == collect(eachrow(I3))
        @test ray_indices(maximal_cones(F1NR)) == incidence1
        @test IncidenceMatrix(maximal_cones(F1NR)) == incidence1
        @test nrays(F2NR) == 0
        @test lineality_dim(F2NR) == 1
        RMLF2NR = rays_modulo_lineality(F2NR)
        @test length(RMLF2NR[:rays_modulo_lineality]) == 2
        @test RMLF2NR[:rays_modulo_lineality] == collect(eachrow(R))
        @test lineality_space(F2NR) == collect(eachrow(L))
        @test ray_indices(maximal_cones(F2NR)) == incidence2
        @test IncidenceMatrix(maximal_cones(F2NR)) == incidence2
    end

end
