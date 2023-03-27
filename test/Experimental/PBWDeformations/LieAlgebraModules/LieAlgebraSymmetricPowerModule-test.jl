@testset ExtendedTestSet "All LieAlgebraSymmetricPowerModule.jl tests" begin

    @testset "constructors for R=$R, n=$n" for R in [QQ, cyclotomic_field(4)[1]],
        n in 1:5,
        L in [general_linear_liealgebra(R, n), special_linear_liealgebra(R, n), special_orthogonal_liealgebra(R, n)]

        for m in 1:min(3, n)
            V1 = standard_module(L)
            V2 = symmetric_power(V1, m)
            @test dim(V2) == binomial(dim(V1) + m - 1, m)

            V1 = symmetric_power(standard_module(L), 3)
            V2 = symmetric_power(V1, m)
            @test dim(V2) == binomial(dim(V1) + m - 1, m)

            V1 = exterior_power(standard_module(L), 2)
            V2 = symmetric_power(V1, m)
            @test dim(V2) == binomial(dim(V1) + m - 1, m)
        end
    end

    @testset "conformance tests for $desc" for (desc, L, V, parentT, elemT) in [
        (
            "S^3 V of sl_4(QQ)",
            special_linear_liealgebra(QQ, 4),
            symmetric_power(standard_module(special_linear_liealgebra(QQ, 4)), 3),
            LieAlgebraSymmetricPowerModule{fmpq},
            LieAlgebraSymmetricPowerModuleElem{fmpq},
        ),
        (
            "S^2 ⋀^2 V of so_4(QQ)",
            special_orthogonal_liealgebra(QQ, 4),
            symmetric_power(exterior_power(standard_module(special_orthogonal_liealgebra(QQ, 4)), 2), 2),
            LieAlgebraSymmetricPowerModule{fmpq},
            LieAlgebraSymmetricPowerModuleElem{fmpq},
        ),
        (
            "S^2 V of so_4(CL(4))",
            special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4),
            symmetric_power(standard_module(special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4)), 3),
            LieAlgebraSymmetricPowerModule{nf_elem},
            LieAlgebraSymmetricPowerModuleElem{nf_elem},
        ),
    ]
        liealgebra_module_conformance_test(L, V, parentT, elemT)
    end
end
