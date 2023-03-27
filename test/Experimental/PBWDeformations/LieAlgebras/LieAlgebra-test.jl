
function liealgebra_conformance_test(L::LieAlgebra{C}, parentT::DataType, elemT::DataType) where {C <: RingElement}
    @testset "basic manipulation" begin
        x = L(rand(-10:10, dim(L)))

        @test parentT <: LieAlgebra{C}
        @test elemT <: LieAlgebraElem{C}
        @test L isa parentT
        @test x isa elemT

        @test parent_type(elemT) == parentT
        @test elem_type(parentT) == elemT

        @test parent(x) == L

        @test base_ring(x) == base_ring(L)
        @test elem_type(base_ring(L)) == C

        @test dim(L) == ngens(L)
        @test length(gens(L)) == ngens(L)
        @test all(gen(L, i) == gens(L)[i] for i in 1:ngens(L))

        @test isempty(rels(L))
    end

    @testset "parent object call overload" begin
        @test L() == zero(L) == L(zeros(base_ring(L), dim(L)))

        for _ in 1:num_random_tests
            coeffs = rand(-10:10, dim(L))
            x1 = L(coeffs)
            x2 = L(base_ring(L).(coeffs))
            x3 = L(matrix(base_ring(L), 1, dim(L), coeffs))
            x4 = L(sparse_row(matrix(base_ring(L), 1, dim(L), coeffs)))
            x5 = L(x1)
            @test x1 == x2
            @test x1 == x3
            @test x1 == x4
            @test x1 == x5
        end
    end

    @testset "lie algebra axioms" begin
        for _ in 1:num_random_tests
            x = L(rand(-10:10, dim(L)))
            y = L(rand(-10:10, dim(L)))
            z = L(rand(-10:10, dim(L)))

            @test bracket(x + y, z) == bracket(x, z) + bracket(y, z)
            @test bracket(x, y + z) == bracket(x, y) + bracket(x, z)

            @test bracket(x, x) == zero(L)
            @test bracket(x, y) == -bracket(y, x)

            @test bracket(x, bracket(y, z)) + bracket(y, bracket(z, x)) + bracket(z, bracket(x, y)) == zero(L)
        end
    end
end

include("AbstractLieAlgebra-test.jl")
include("LinearLieAlgebra-test.jl")

@testset ExtendedTestSet "All LieAlgebra.jl tests" begin
    # nothing here yet
end
