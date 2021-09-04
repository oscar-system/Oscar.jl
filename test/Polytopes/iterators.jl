@testset "iterators" begin

    @testset "VectorIterator" begin
        vecmat = [1 0 0; 0 1 0; 0 0 1; 1 1 1]
        for (T, U) in Base.product([PointVector, RayVector], [Polymake.Rational, Polymake.Integer])
            vi = VectorIterator{T{U}}(vecmat)
            @test vi isa VectorIterator
            @test vi isa VectorIterator{T{U}}
            @test eltype(vi) == T{U}
            @test size(vi) == (4,)
            @test length(vi) == 4
            @test firstindex(vi) == 1
            @test lastindex(vi) == 4
            for i in 1:4
                @test vi[i] isa T{U}
                @test vi[i] == vecmat[i, :]
            end
            @test point_matrix(vi) isa Polymake.Matrix
            @test point_matrix(vi) == vecmat
            @test Oscar.matrix_for_polymake(vi) == vecmat
        end
        @test VectorIterator(vecmat) isa VectorIterator{PointVector{Polymake.Rational}}
    end

    @testset "HalfspaceIterator" begin
        A = [-1 0 0; 0 -1 0; 0 0 -1; 1 1 1]
        b = [0, 1, 2, 4]
        for T in [Halfspace, Polyhedron]
            hi = HalfspaceIterator{T}(A, b)
            @test hi isa HalfspaceIterator
            @test hi isa HalfspaceIterator{T}
            @test eltype(hi) == T
            @test size(hi) == (4,)
            @test length(hi) == 4
            @test firstindex(hi) == 1
            @test lastindex(hi) == 4
            for i in 1:4
                @test hi[i] isa T
                @test hi[i] == T(reshape(A[i, :], 1, :), b[i])
            end
            @test halfspace_matrix_pair(hi) isa NamedTuple
            @test halfspace_matrix_pair(hi).A == A
            @test halfspace_matrix_pair(hi).b == b
            @test_throws ArgumentError point_matrix(hi)
        end
        @test HalfspaceIterator(A, b) isa HalfspaceIterator{Halfspace}
    end

    @testset "PolyhedronOrConeIterator" begin
        V = [1 -1 0 0 0 0; 1 0 -1 0 0 0; 1 0 0 -1 0 0; 1 1 1 1 0 0]
        F = Polymake.Set{Polymake.to_cxx_type(Int64)}.([[1 2 3], [1 2 4], [1 3 4], [2 3 4]])
        L = [0 0 0 0 1 0; 0 0 0 0 0 1]
        for T in [Polyhedron, Cone]
            pci = PolyhedronOrConeIterator{T}(V, F, L)
            @test pci isa PolyhedronOrConeIterator
            @test pci isa PolyhedronOrConeIterator{T}
            @test eltype(pci) == T
            @test size(pci) == (4,)
            @test length(pci) == 4
            @test firstindex(pci) == 1
            @test lastindex(pci) == 4
            for i in 1:4
                @test pci[i] isa T
                if T == Polyhedron
                    @test pci[i] == convex_hull(V[[f for f in F[i]], 2:end], nothing, L[:,2:end])
                else
                    @test pci[i] == T(V[[f for f in F[i]], :], L)
                end
            end
            # test incidence_matrix(iter)
        end
    end

end
