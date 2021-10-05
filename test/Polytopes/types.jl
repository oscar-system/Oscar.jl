@testset "types" begin

    a = [1, 2, 3]
    b = [8, 6, 4]

    @testset "$T" for T in (PointVector, RayVector)

        @test T(a) isa T{Polymake.Rational}

        @testset "$T{$U}" for U in (Int64, Polymake.Integer, Polymake.Rational, Float64)

            @test T{U} <: AbstractVector
            @test T{U} <: AbstractVector{U}

            @test T{U}(a) isa T
            @test T{U}(a) isa T{U}

            @test T{U}(7) == zeros(7)

            A = T{U}(a)

            @test A.p isa Polymake.Vector{Polymake.to_cxx_type(U)}

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

            for V in [Int64, Polymake.Integer, Polymake.Rational, Float64]

                B = T{V}(b)

                for op in [+, -]
                    @test op(A, B) isa T
                    @test op(A, B) isa T{promote_type(U, V)}

                    @test op(A, B) == op(a, b)
                end

                @test *(V(3), A) isa T
                @test *(V(3), A) isa T{promote_type(U, V)}

                @test *(V(3), A) == *(V(3), A)

                @test [A; B] isa T
                @test [A; B] isa T{promote_type(U, V)}
                @test [A; B] == [a; b]

            end

        end

    end





    @testset "$T" for T in (Halfspace, Hyperplane)

        @test T(a, 0) isa T
        @test T(a', 0) isa T

        @test T(a, 0) == T(a', 0) == T(a) == T(a')

        A = T(a, 0)
        B = T(b, 2)

        @test A != B

        @test A.a == a
        @test A.b == 0

        @test B.a == b
        @test B.b == 2

    end

end
