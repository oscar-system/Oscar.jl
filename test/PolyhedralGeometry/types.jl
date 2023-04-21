@testset "types" begin
    
    @testset "IncidenceMatrix" begin
        
        im = IncidenceMatrix([[1,2,3],[4,5,6]])
        @test nrows(im) == 2
        @test ncols(im) == 6
        
    end
    
    @testset "nf_scalar" begin
        
        qe = Polymake.QuadraticExtension{Polymake.Rational}(123, 456, 789)
        @test convert(Oscar.nf_scalar, qe) isa nf_elem
        nfe = convert(Oscar.nf_scalar, qe)
        @test coordinates(nfe) == [123, 456]
        p = defining_polynomial(parent(nfe))
        @test p == parent(p)([-789, 0, 1])
        @test convert(Polymake.QuadraticExtension{Polymake.Rational}, nfe) == qe
        
        qet = Polymake.QuadraticExtension{Polymake.Rational}(9)
        @test convert(Oscar.nf_scalar, qet) isa QQFieldElem
        @test convert(Oscar.nf_scalar, qet) == 9
        @test convert(Polymake.QuadraticExtension{Polymake.Rational}, convert(Oscar.nf_scalar, qet)) == qet
        
        @test convert(Oscar.nf_scalar, 5) isa QQFieldElem
        @test convert(Oscar.nf_scalar, 5) == 5
        
    end

    a = [1, 2, 3]
    b = [8, 6, 4]

    #TODO: RayVector
    @testset "$T" for T in (PointVector, RayVector)

        @test T(a) isa T{QQFieldElem}

        @testset "$T{$U}" for U in (ZZRingElem, QQFieldElem, Oscar.nf_scalar)

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

            for V in [ZZRingElem, QQFieldElem, Oscar.nf_scalar]
            
                B = T{V}(b)
            
                for op in [+, -]
                    @test op(A, B) isa T
            
                    @test op(A, B) == op(a, b)
                end
            
                @test *(V(3), A) isa T
                
                @test *(V(3), A) == 3 * a
                
                @test [A; B] == [a; b]
            
            end

        end

    end





    @testset "$T" for T in (AffineHalfspace, AffineHyperplane)
        
        for U in [QQFieldElem]
            @test T{U}(a, 0) isa T{U}
            @test T{U}(permutedims(a), 0) isa T{U}

            @test T{U}(a, 0) == T{U}(permutedims(a), 0) == T{U}(a) == T{U}(permutedims(a))

            A = T{U}(a, 0)
            B = T{U}(b, 2)

            @test A != B

            @test normal_vector(A) isa Vector{U}
            @test normal_vector(A) == a
            @test negbias(A) isa U
            @test negbias(A) == 0
            
            
            @test normal_vector(B) == b
            @test negbias(B) == 2
        end
        
        @test T(a, 0) isa T{QQFieldElem}
        
    end
    
    @testset "$T" for T in (LinearHalfspace, LinearHyperplane)
        
        for U in [QQFieldElem]
            @test T{U}(a) isa T{U}
            @test T{U}(permutedims(a)) isa T{U}

            @test T{U}(a) == T{U}(permutedims(a))

            A = T{U}(a)
            B = T{U}(b)

            @test A != B

            @test normal_vector(A) isa Vector{U}
            @test normal_vector(A) == a
            @test negbias(A) isa U
            @test negbias(A) == 0
            
            
            @test normal_vector(B) == b
            @test negbias(B) == 0
        end
        
        @test T(a) isa T{QQFieldElem}
        
    end
    
    for U in [QQFieldElem]
        let A = LinearHalfspace{U}(a)
            @test invert(A) isa LinearHalfspace{U}
            Ai = invert(A)
            @test normal_vector(Ai) == -a
        end
        let A = AffineHalfspace{U}(a, 8)
            @test invert(A) isa AffineHalfspace{U}
            Ai = invert(A)
            @test normal_vector(Ai) == -a
            @test negbias(Ai) == -8
        end
    end
    
    @test Halfspace(a) isa LinearHalfspace{QQFieldElem}
    @test Hyperplane(a) isa LinearHyperplane{QQFieldElem}
    @test Halfspace(a, 0) isa AffineHalfspace{QQFieldElem}
    @test Hyperplane(a, 0) isa AffineHyperplane{QQFieldElem}

end
