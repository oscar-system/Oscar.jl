@testset ExtendedTestSet "All LieAlgebraTensorPowerModule.jl tests" begin

    @testset "constructors for R=$R, n=$n" for R in [QQ, cyclotomic_field(4)[1]],
        n in 1:5,
        L in [general_linear_liealgebra(R, n), special_linear_liealgebra(R, n), special_orthogonal_liealgebra(R, n)]

        for m in 1:3
            V1 = standard_module(L)
            V2 = tensor_power(V1, m)
            @test dim(V2) == dim(V1)^m

            V1 = symmetric_power(standard_module(L), 3)
            V2 = tensor_power(V1, m)
            @test dim(V2) == dim(V1)^m

            V1 = tensor_power(standard_module(L), 2)
            V2 = tensor_power(V1, m)
            @test dim(V2) == dim(V1)^m
        end
    end

    @testset "conformance tests for $desc" for (desc, L, V, parentT, elemT) in [
        (
            "T^3 V of sl_4(QQ)",
            special_linear_liealgebra(QQ, 4),
            tensor_power(standard_module(special_linear_liealgebra(QQ, 4)), 3),
            LieAlgebraTensorPowerModule{fmpq},
            LieAlgebraTensorPowerModuleElem{fmpq},
        ),
        (
            "T^2 ⋀^2 V of so_4(QQ)",
            special_orthogonal_liealgebra(QQ, 4),
            tensor_power(exterior_power(standard_module(special_orthogonal_liealgebra(QQ, 4)), 2), 2),
            LieAlgebraTensorPowerModule{fmpq},
            LieAlgebraTensorPowerModuleElem{fmpq},
        ),
        (
            "T^2 V of so_4(CL(4))",
            special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4),
            tensor_power(standard_module(special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4)), 3),
            LieAlgebraTensorPowerModule{nf_elem},
            LieAlgebraTensorPowerModuleElem{nf_elem},
        ),
    ]
        liealgebra_module_conformance_test(L, V, parentT, elemT)
    end
end
