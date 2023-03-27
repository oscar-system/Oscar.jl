@testset ExtendedTestSet "All LieAlgebraStdModule.jl tests" begin

    @testset "constructors for R=$R, n=$n" for R in [QQ, cyclotomic_field(4)[1]], n in 1:5
        L = general_linear_liealgebra(R, n)
        V = standard_module(L)
        @test dim(V) == n

        L = special_linear_liealgebra(R, n)
        V = standard_module(L)
        @test dim(V) == n

        L = special_orthogonal_liealgebra(R, n)
        V = standard_module(L)
        @test dim(V) == n
    end

    @testset "conformance tests for $desc" for (desc, L, V, parentT, elemT) in [
        (
            "V of sl_4(QQ)",
            special_linear_liealgebra(QQ, 4),
            standard_module(special_linear_liealgebra(QQ, 4)),
            LieAlgebraStdModule{fmpq},
            LieAlgebraStdModuleElem{fmpq},
        ),
        (
            "V of so_4(QQ)",
            special_orthogonal_liealgebra(QQ, 4),
            standard_module(special_orthogonal_liealgebra(QQ, 4)),
            LieAlgebraStdModule{fmpq},
            LieAlgebraStdModuleElem{fmpq},
        ),
        (
            "V of sl_4(CF(4))",
            special_linear_liealgebra(cyclotomic_field(4)[1], 4),
            standard_module(special_linear_liealgebra(cyclotomic_field(4)[1], 4)),
            LieAlgebraStdModule{nf_elem},
            LieAlgebraStdModuleElem{nf_elem},
        ),
        (
            "V of so_4(CF(4))",
            special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4),
            standard_module(special_orthogonal_liealgebra(cyclotomic_field(4)[1], 4)),
            LieAlgebraStdModule{nf_elem},
            LieAlgebraStdModuleElem{nf_elem},
        ),
    ]
        liealgebra_module_conformance_test(L, V, parentT, elemT)
    end
end
