@testset "types" begin

    a = [1, 2, 3]
    b = [8, 6, 4]

    #TODO: RayVector
    @testset "$T" for T in (PointVector, RayVector)

        @test T(a) isa T{fmpq}

        @testset "$T{$U}" for U in (fmpz, fmpq)

            @test T{U} <: AbstractVector
            @test T{U} <: AbstractVector{U}

            @test T{U}(a) isa T
            @test T{U}(a) isa T{U}

            @test T{U}(7) == zeros(7)

            A = T{U}(a)

            @test A[1] isa U
            @test A[2] == 2
            @test A[begin] == 1
            @test A[end] == 3
            A[3] = 7
            @test A[end] == 7
            A[3] = 3

            @test size(A) == (3,)

            @test_throws BoundsError A[0]
            @test_throws BoundsError A[4]

            for V in [fmpz, fmpq]
            
                B = T{V}(b)
            
                for op in [+, -]
                    @test op(A, B) isa T
                    @test op(A, B) isa T{promote_type(U, V)}
            
                    @test op(A, B) == op(a, b)
                end
            
                @test *(V(3), A) isa T
                @test *(V(3), A) isa T{promote_type(U, V)}
                
                @test *(V(3), A) == 3 * a
                
                @test [A; B] isa T
                @test [A; B] isa T{promote_type(U, V)}
                @test [A; B] == [a; b]
            
            end

        end

    end





    @testset "$T" for T in (AffineHalfspace, AffineHyperplane)
        
        for U in (fmpz, fmpq)
            @test T{U}(a, 0) isa T{U}
            @test T{U}(a', 0) isa T{U}

            @test T{U}(a, 0) == T{U}(a', 0) == T{U}(a) == T{U}(a')

            A = T{U}(a, 0)
            B = T{U}(b, 2)

            @test A != B

            @test A.a == a
            @test A.b == Oscar.negbias(A) == 0

            @test B.a == b
            @test B.b == Oscar.negbias(B) == 2
        end
        
        @test T(a, 0) isa T{fmpq}
        
    end
    
    @test Halfspace(a) isa LinearHalfspace{fmpq}
    @test Hyperplane(a) isa LinearHyperplane{fmpq}
    @test Halfspace(a, 0) isa AffineHalfspace{fmpq}
    @test Hyperplane(a, 0) isa AffineHyperplane{fmpq}

end
