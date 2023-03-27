
function liealgebra_module_conformance_test(
    L::LieAlgebra{C},
    V::LieAlgebraModule{C},
    parentT::DataType,
    elemT::DataType,
) where {C <: RingElement}
    @testset "basic manipulation" begin
        v = V(rand(-10:10, dim(V)))

        @test parentT <: LieAlgebraModule{C}
        @test elemT <: LieAlgebraModuleElem{C}
        @test V isa parentT
        @test v isa elemT

        @test parent_type(elemT) == parentT
        @test elem_type(parentT) == elemT

        @test parent(v) == V

        @test base_ring(v) == base_ring(V)
        @test elem_type(base_ring(V)) == C

        @test dim(V) == ngens(V)
        @test length(gens(V)) == ngens(V)
        @test all(gen(V, i) == gens(V)[i] for i in 1:ngens(V))

        @test isempty(rels(V))
    end

    @testset "parent object call overload" begin
        @test V() == zero(V) == V(zeros(base_ring(V), dim(V)))

        for _ in 1:num_random_tests
            coeffs = rand(-10:10, dim(V))
            v1 = V(coeffs)
            v2 = V(base_ring(V).(coeffs))
            v3 = V(matrix(base_ring(V), 1, dim(V), coeffs))
            v4 = V(sparse_row(matrix(base_ring(V), 1, dim(V), coeffs)))
            v5 = V(v1)
            @test v1 == v2
            @test v1 == v3
            @test v1 == v4
            @test v1 == v5
        end
    end

    @testset "lie algebra action axioms" begin
        for _ in 1:num_random_tests
            x = L(rand(-10:10, dim(L)))
            y = L(rand(-10:10, dim(L)))
            v = V(rand(-10:10, dim(V)))
            w = V(rand(-10:10, dim(V)))

            @test (x * v) isa elemT
            @test parent(x * v) == parent(v)

            @test (x + y) * v == x * v + y * v
            @test x * (v + w) == x * v + x * w

            @test bracket(x, y) * v == x * (y * v) - y * (x * v)
        end
    end
end

include("LieAlgebraAbstractModule-test.jl")
include("LieAlgebraStdModule-test.jl")
include("LieAlgebraExteriorPowerModule-test.jl")
include("LieAlgebraSymmetricPowerModule-test.jl")
include("LieAlgebraTensorPowerModule-test.jl")

@testset ExtendedTestSet "All LieAlgebraModule.jl tests" begin
    # nothing here yet
end
