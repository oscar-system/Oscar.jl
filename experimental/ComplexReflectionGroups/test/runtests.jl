@testset "Complex reflection groups" begin

    # We run some consistency tests on the models

    # For example, we check the order of the group. This is set as an attribute already
    # upon creation, so we create a copy of the matrix group, determine its order and 
    # check whether the numbers match.

    # CHEVIE models for exceptional groups
    for n=4:37
        W = complex_reflection_group(n ; model=:CHEVIE)
        W_copy = matrix_group(gens(W))
        @test order(W) == order(W_copy)
    end

    # Magma models for exceptional groups
    for n=4:37
        W = complex_reflection_group(n ; model=:Magma)
        W_copy = matrix_group(gens(W))
        @test order(W) == order(W_copy)
    end

    # LT models for exceptional groups (>22 missing here at the moment)
    for n=4:22
        W = complex_reflection_group(n ; model=:LT)
        W_copy = matrix_group(gens(W))
        @test order(W) == order(W_copy)
        @test is_unitary(W) == true #LT models are unitary
    end

end