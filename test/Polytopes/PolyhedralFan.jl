@testset "PolyhedralFan" begin
    C0 = cube(2)
    NFsquare = normal_fan(C0)
    R = [1 0 0; 0 0 1]
    L = [0 1 0]
    Cone4 = positive_hull(R)
    Cone5 = positive_hull([1 0 0; 0 1 0])

    F0 = PolyhedralFan([Cone4, Cone5])
    I3 = [1 0 0; 0 1 0; 0 0 1]
    incidence1 = IncidenceMatrix([[1,2],[2,3]])
    incidence2 = IncidenceMatrix([[1,2]])
    F1 = PolyhedralFan(I3, incidence1)
    F2 = PolyhedralFan(R, L, incidence2)

    @testset "core functionality" begin
        @test issmooth(NFsquare)
        @test rays(NFsquare).m == [1 0; -1 0; 0 1; 0 -1]
        @test isregular(NFsquare)
        @test iscomplete(NFsquare)
        @test !iscomplete(F0)
        @test length(rays(F0)) == 3
        @test nrays(F1) == 3
        @test dim(F1) == 2
        @test ambient_dim(F1) == 3
        @test nrays(F2) == 2
        @test dim.(maximal_cones(F1)) == [2,2]
        @test nmaximal_cones(F1) == 2
        @test lineality_space(F2).m == L
        @test length(collect(rays(F0))) == 3

        II = maximal_cones_as_incidence_matrix(NFsquare)
        NF0 = PolyhedralFan(rays(NFsquare).m, II)
        @test nrays(NF0) == 4
        FF0 = face_fan(C0)
        @test nrays(FF0) == 4
        @test !issmooth(FF0)
        @test f_vector(NFsquare) == [4, 4]
    end

end
