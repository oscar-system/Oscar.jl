@testset "Complex reflection groups" begin

    # We run some consistency tests on the models

    # Note: the complex_reflection_group constructor sets attributes like order and
    # is_complex_reflection_group etc. To perform an honest test, we create a copy of the
    # matrix group which then has no attributes set and run the tests on this copy.

    # CHEVIE models for exceptional groups
    for n=4:37
        W = complex_reflection_group(n ; model=:CHEVIE)
        W_copy = matrix_group(gens(W))
        @test order(W) == order(W_copy)
        @test is_complex_reflection_group(W_copy) == true
    end

    # Magma models for exceptional groups
    for n=4:37
        W = complex_reflection_group(n ; model=:Magma)
        W_copy = matrix_group(gens(W))
        @test order(W) == order(W_copy)
        @test is_complex_reflection_group(W_copy) == true
    end

    # LT models for exceptional groups (>22 missing here at the moment)
    for n=4:22
        W = complex_reflection_group(n ; model=:LT)
        W_copy = matrix_group(gens(W))
        @test order(W) == order(W_copy)
        @test is_complex_reflection_group(W_copy) == true
        @test is_unitary(W_copy) == true #LT models are unitary
    end

end