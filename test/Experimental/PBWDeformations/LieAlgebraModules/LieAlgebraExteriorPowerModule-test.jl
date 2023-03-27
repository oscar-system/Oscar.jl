@testset ExtendedTestSet "All LieAlgebraExteriorPowerModule.jl tests" begin

    @testset "constructors for R=$R, n=$n" for R in [QQ, cyclotomic_field(4)[1]],
        n in 1:5,
        L in [general_linear_liealgebra(R, n), special_linear_liealgebra(R, n), special_orthogonal_liealgebra(R, n)]

        for m in 1:min(3, n)
            V1 = standard_module(L)
            V2 = exterior_power(V1, m)
            @test dim(V2) == binomial(dim(V1), m)

            V1 = symmetric_power(standard_module(L), 3)
            V2 = exterior_power(V1, m)
            @test dim(V2) == binomial(dim(V1), m)

            V1 = exterior_power(standard_module(L), 2)
            V2 = exterior_power(V1, m)
            @test dim(V2) == binomial(dim(V1), m)
        end
    end

    @testset "conformance tests for $desc" for (desc, L, V, parentT, elemT) in [
        (
            "⋀^3 V of sl_4(QQ)",
            special_linear_liealgebra(QQ, 4),
            exterior_power(standard_module(special_linear_liealgebra(QQ, 4)), 3),
            LieAlgebraExteriorPowerModule{fmpq},
            LieAlgebraExteriorPowerModuleElem{fmpq},
        ),
        (
            "⋀^2 S^2 V of so_4(QQ)",
            special_orthogonal_liealgebra(QQ, 4),
            exterior_power(symmetric_power(standard_module(special_orthogonal_liealgebra(QQ, 4)), 2), 2),
            LieAlgebraExteriorPowerModule{fmpq},
            LieAlgebraExteriorPowerModuleElem{fmpq},
        ),
        (
            "⋀^2 V of so_4(CL(4))",
            special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4),
            exterior_power(standard_module(special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4)), 3),
            LieAlgebraExteriorPowerModule{nf_elem},
            LieAlgebraExteriorPowerModuleElem{nf_elem},
        ),
    ]
        liealgebra_module_conformance_test(L, V, parentT, elemT)
    end
end
