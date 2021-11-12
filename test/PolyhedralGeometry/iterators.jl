@testset "iterators" begin

    @testset "VectorIterator" begin
        vecmat_old = [1 0 0; 0 1 0; 0 0 1; 1 1 1]
        vecmat = [1 0 0; 1 2 3; 0 0 1; 1 1 1]
        @testset "VectorIterator $T{$U}" for T in (PointVector, RayVector), U in (Polymake.Rational, Polymake.Integer)
            vi = VectorIterator{T{U}}(vecmat_old)
            @test vi isa VectorIterator
            @test vi isa VectorIterator{T{U}}
            @test eltype(vi) == T{U}
            @test size(vi) == (4,)
            @test length(vi) == 4
            @test firstindex(vi) == 1
            @test lastindex(vi) == 4
            vi[2] = [1, 2, 3]
            for i in 1:4
                @test vi[i] isa T{U}
                @test vi[i] == vecmat[i, :]
            end
            @test vi.m isa Polymake.Matrix{U}
            @test matrix(vi) isa (U == Polymake.Rational ? fmpq_mat : fmpz_mat)
            @test matrix(vi) == matrix(QQ, vecmat)
            @test Oscar.matrix_for_polymake(vi) == vecmat
            @test VectorIterator{T{U}}(matrix(vi)).m == vecmat
        end
        @test VectorIterator(vecmat) isa VectorIterator{PointVector{Polymake.Rational}}
    end
    
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
        for i in 1:4
            p = soi[i]
            @test p[1] == PointVector{Polymake.Rational}(Oscar.pm_object(c).VERTICES[i, 2:end])
            a, b = Oscar.decompose_hdata(Oscar.pm_object(c).FACETS)
            @test p[2] == Halfspace(a[i, :], b[i])
        end
    end

end
