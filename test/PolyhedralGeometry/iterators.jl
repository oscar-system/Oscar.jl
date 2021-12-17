@testset "iterators" begin
    
    @testset "SubObjectIterator" begin
        # SubObjectIterator is flexible due to the `Function` object it stores, which is used for `getindex`.
        # For the tests we want this iterator to return a `Pair` consisting of a vertex and a halfspace.
        # This context is very specific and probably not meaningful from a mathematical point of view,
        # but it serves well for abstractly testing this type.
        c = cube(2)
        function _testSOI(::Type{Pair{PointVector{Polymake.Rational}, Halfspace}}, obj::Polymake.BigObject, i::Base.Integer)
            x = PointVector{Polymake.Rational}(obj.VERTICES[i, 2:end])
            a, b = Oscar.decompose_hdata(obj.FACETS)
            y = Halfspace(a[i, :], b[i])
            return Pair{PointVector{Polymake.Rational}, Halfspace}(x, y)
        end
        soi = SubObjectIterator{Pair{PointVector{Polymake.Rational}, Halfspace}}(Oscar.pm_object(c), _testSOI, 4)
        @test soi isa AbstractVector{Pair{PointVector{Polymake.Rational}, Halfspace}}
        @test length(soi) == 4
        @test firstindex(soi) == 1
        @test lastindex(soi) == 4
        for i in 1:4
            p = soi[i]
            @test p[1] == PointVector{Polymake.Rational}(Oscar.pm_object(c).VERTICES[i, 2:end])
            a, b = Oscar.decompose_hdata(Oscar.pm_object(c).FACETS)
            @test p[2] == Halfspace(a[i, :], b[i])
        end
        @test_throws ArgumentError ray_incidences(soi)
        @test_throws ArgumentError linear_inequality_matrix(soi)
        @test_throws ArgumentError halfspace_matrix_pair(soi)
        @test_throws ArgumentError Oscar.linear_matrix_for_polymake(soi)
        @test_throws ArgumentError Oscar.affine_matrix_for_polymake(soi)
        @test_throws ArgumentError Oscar.matrix_for_polymake(soi)
        soi2 = SubObjectIterator{PointVector{Polymake.Rational}}(Oscar.pm_object(c), _testSOI, 4)
        @test_throws ArgumentError point_matrix(soi2)
    end

end
