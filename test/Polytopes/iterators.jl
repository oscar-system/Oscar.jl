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

    @testset "HalfspaceIterator" begin
        A_old = [-1 0 0; 0 -2 0; 0 0 -1; 1 1 1]
        b_old = [0, 3, 2, 4]
        @testset "HalfspaceIterator{$T}" for T in (Halfspace, Polyhedron, Pair{Polymake.Matrix{Polymake.Rational}, Polymake.Rational}, Hyperplane)
            A = A_old
            b = b_old
            hi = HalfspaceIterator{T}(A, b)
            @test hi isa HalfspaceIterator
            @test hi isa HalfspaceIterator{T}
            @test eltype(hi) == T
            @test size(hi) == (4,)
            @test length(hi) == 4
            @test firstindex(hi) == 1
            @test lastindex(hi) == 4
            A[2, 2] = -1
            b[2] = 1
            hi[2] = Halfspace([0 -1 0], 1)
            for i in 1:4
                @test hi[i] isa T
                @test hi[i] == T(reshape(A[i, :], 1, :), b[i])
            end
            @test halfspace_matrix_pair(hi) isa NamedTuple{(:A, :b), Tuple{fmpq_mat, Vector{fmpq}}}
            @test halfspace_matrix_pair(hi).A == matrix(QQ, A)
            @test halfspace_matrix_pair(hi).b == b
            @test_throws ArgumentError matrix(hi)
            hai = HalfspaceIterator{T}(A)
            @test hai.A == A
            @test hai.b == zeros(4)
            @test matrix(hai) == matrix(QQ, A)
        end
        @test HalfspaceIterator(A_old, b_old) isa HalfspaceIterator{Halfspace}
        hai = HalfspaceIterator{Cone}(A_old)
        @test hai isa HalfspaceIterator{Cone}
        for i in 1:4
            @test hai[i] isa Cone
            @test hai[i] == cone_from_inequalities(reshape(A_old[i, :], 1, :))
        end

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
