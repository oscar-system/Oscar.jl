@testset "types" begin
    
    @testset "IncidenceMatrix" begin
        
        im = IncidenceMatrix([[1,2,3],[4,5,6]])
        @test nrows(im) == 2
        @test ncols(im) == 6
        @test row(im, 1) isa Set{Int}
        @test row(im, 1) == Set{Int}([1, 2, 3])
        @test column(im, 2) isa Set{Int}
        @test column(im, 2) == Set{Int}([1])
        
    end

    a = [1, 2, 3]
    b = [8, 6, 4]

    NF, sr2 = quadratic_field(2)

    #TODO: RayVector
    @testset "$T" for T in (PointVector, RayVector)

        @test T(a) isa T{QQFieldElem}

        for f in (ZZ, QQ, NF)

            U = elem_type(f)
            @testset "$T{$U}" begin

                @test T{U} <: AbstractVector
                @test T{U} <: AbstractVector{U}

                @test T{U}(f.(a)) isa T
                @test T{U}(f.(a)) isa T{U}

                @test T{U}(f, 7) == f.(zeros(Int, 7))

                A = T{U}(f.(a))

                @test A[1] isa U
                @test A[2] == 2
                @test A[begin] == 1
                @test A[end] == 3
                A[3] = f(7)
                @test A[end] == 7
                A[3] = f(3)

                @test size(A) == (3,)

                @test_throws BoundsError A[0]
                @test_throws BoundsError A[4]

                for g in (ZZ, QQ, NF)
                V = elem_type(g)

                    B = T{V}(g.(b))

                    for op in [+, -]
                        @test op(A, B) isa T

                        @test op(A, B) == op(a, b)
                    end

                    @test *(g(3), A) isa T

                    @test *(g(3), A) == 3 * a

                    @test [A; B] == [a; b]

                end

            end
        end

    end





    @testset "$T" for (T, f) in ((AffineHalfspace, affine_halfspace), (AffineHyperplane, affine_hyperplane))
        
        for p in [QQ, NF]
            U = elem_type(p)
            @test f(p, a, 0) isa T{U}
            @test f(p, permutedims(a), 0) isa T{U}

            @test f(p, a, 0) == f(p, permutedims(a), 0) == f(p, a) == f(p, permutedims(a))

            A = f(p, a, 0)
            B = f(p, b, 2)

            @test A != B

            @test normal_vector(A) isa Vector{U}
            @test normal_vector(A) == a
            @test negbias(A) isa U
            @test negbias(A) == 0
            
            
            @test normal_vector(B) == b
            @test negbias(B) == 2
        end
        
        @test f(a, 0) isa T{QQFieldElem}
        
    end
    
    @testset "$T" for (T, f) in ((LinearHalfspace, linear_halfspace), (LinearHyperplane, linear_hyperplane))
        
        for p in [QQ, NF]
            U = elem_type(p)
            @test f(p, a) isa T{U}
            @test f(p, permutedims(a)) isa T{U}

            @test f(p, a) == f(p, permutedims(a))

            A = f(p, a)
            B = f(p, b)

            @test A != B

            @test normal_vector(A) isa Vector{U}
            @test normal_vector(A) == a
            @test negbias(A) isa U
            @test negbias(A) == 0
            
            
            @test normal_vector(B) == b
            @test negbias(B) == 0
        end
        
        @test f(a) isa T{QQFieldElem}
        
    end
    
    for p in [QQ, NF]
        U = elem_type(p)
        let A = linear_halfspace(p, a)
            @test invert(A) isa LinearHalfspace{U}
            Ai = invert(A)
            @test normal_vector(Ai) == -a
        end
        let A = affine_halfspace(p, a, 8)
            @test invert(A) isa AffineHalfspace{U}
            Ai = invert(A)
            @test normal_vector(Ai) == -a
            @test negbias(Ai) == -8
        end
    end
    
    @test halfspace(a) isa LinearHalfspace{QQFieldElem}
    @test hyperplane(a) isa LinearHyperplane{QQFieldElem}
    @test halfspace(a, 0) isa AffineHalfspace{QQFieldElem}
    @test hyperplane(a, 0) isa AffineHyperplane{QQFieldElem}

end
