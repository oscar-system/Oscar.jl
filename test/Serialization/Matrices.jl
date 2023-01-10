R, x = PolynomialRing(QQ, "x")
q = x^2 + 3//4
K, a = NumberField(q)
Z7 = ResidueRing(ZZ, 7)
Z7t, t = PolynomialRing(Z7, "t")
Fin, d = FiniteField(t^2 + t + 1)
Frac = FractionField(R)

cases = [
    (QQ, [1 2; 3 4//5]),
    (Z7, [1 3; 4 2]),
    (K, [a a + 1; a - 1 0]),
    (Fin, [d d^3; d^2 0]),
    (Frac, [1 // x x^2; 3 0])
]

function test_equality(m::MatrixElem{T}, l::MatrixElem{T}) where T <: Union{
    fmpz, fmpq, nmod, fq_nmod}
    return change_base_ring(base_ring(l), m) == l
end

function test_equality(m::MatrixElem{T}, l:: MatrixElem{T}) where T <: AbstractAlgebra.Generic.Frac{fmpq_poly}
    x = gen(base_ring(m))
    return map(i -> evaluate(i, x), l) == m
end

function test_equality(m::MatrixElem{nf_elem}, l:: MatrixElem{nf_elem}) 
    R = base_ring(m)
    S = base_ring(l)
    h = hom(R, S, gen(S))
    return map(h, m) == l
end

@testset "Matrices" begin
    mktempdir() do path
        for case in cases
            @testset "Matrices over $(case[1])" begin
                m = matrix(case[1], case[2])
                test_save_load_roundtrip(path, m) do loaded
                    @test test_equality(m, loaded)
                end

                test_save_load_roundtrip(path, m; parent=parent(m)) do loaded
                    @test loaded == m
                end
            end
        end
    end
end
