@testset "Complex reflection groups" begin

    #######################################################################################
    # Testing ComplexReflection
    #######################################################################################

    # Example 1
    V = vector_space(QQ,2)
    s = transpose(matrix(QQ,2,2,[-1 1; 0 1]))
    b,s_data = is_complex_reflection_with_data(s)
    @test b == true
    @test root(s_data) == V([1,0])
    @test coroot(s_data) == V([2,-1])
    @test hyperplane_basis(s_data) == [V([1,2])]
    @test matrix(complex_reflection(root(s_data), coroot(s_data))) == s
    @test is_unitary(s_data) == false

    # Example 2
    V = vector_space(QQ,2)
    s = transpose(matrix(QQ,2,2,[-1 1; 0 1]))
    b,s_data = is_complex_reflection_with_data(s)
    @test b == true
    @test root(s_data) == V([1,0])
    @test coroot(s_data) == V([2,-1])
    @test hyperplane_basis(s_data) == [V([1,2])]
    @test matrix(complex_reflection(root(s_data), coroot(s_data))) == s
    @test is_unitary(s_data) == false

    # Example 3 (orthogonal reflection in (1,1))
    V = vector_space(QQ,2)
    w = unitary_reflection(V([1,1]))
    @test root(w) == V([1,1])
    @test hyperplane_basis(w) == [V([1,-1])]
    @test matrix(w) == matrix(QQ,2,2,[0 -1; -1 0])
    @test is_unitary(w) == true
    @test order(w) == 2
    @test matrix(complex_reflection(root(w), coroot(w))) == matrix(w)

    # Example 4 (a transvection)
    t = transpose(matrix(QQ,2,2,[1 1; 0 1]))
    @test is_complex_reflection(t) == false

    # Example 5 (matrix fixing a hyperplane but not invertible)
    w = matrix(QQ,2,2,[1 0; 0 0])
    @test is_complex_reflection(w) == false

    # Example 6
    K,z = cyclotomic_field(3)
    V = vector_space(K,2)
    s = transpose(matrix(K, 2, 2, [z 0; z^-1 1]))
    b,s_data = is_complex_reflection_with_data(s)
    @test b == true
    @test order(s_data) == 3
    @test matrix(complex_reflection(root(s_data), coroot(s_data))) == s

    # Example 7
    K,z = cyclotomic_field(8)
    V = vector_space(K,2)
    s = transpose(matrix(K, 2, 2, [-1 0; -z+1 1]))
    b,s_data = is_complex_reflection_with_data(s)
    @test b == true
    @test order(s_data) == 2


    #######################################################################################
    # Exceptional complex reflection groups
    #######################################################################################

    # Note: the complex_reflection_group constructor sets attributes like order and
    # is_complex_reflection_group etc. To perform an honest test, we create a copy of the
    # matrix group which then has no attributes set and run the tests on this copy.

    # CHEVIE models for exceptional groups
    # for n=4:37
    #     W = complex_reflection_group(n ; model=:CHEVIE)
    #     W_copy = matrix_group(gens(W))
    #     @test order(W) == order(W_copy)
    #     @test is_complex_reflection_group(W_copy) == true
    # end

    # # Magma models for exceptional groups
    # for n=4:37
    #     W = complex_reflection_group(n ; model=:Magma)
    #     W_copy = matrix_group(gens(W))
    #     @test order(W) == order(W_copy)
    #     @test is_complex_reflection_group(W_copy) == true
    # end

    # # LT models for exceptional groups (>22 missing here at the moment)
    # for n=4:22
    #     W = complex_reflection_group(n ; model=:LT)
    #     W_copy = matrix_group(gens(W))
    #     @test order(W) == order(W_copy)
    #     @test is_complex_reflection_group(W_copy) == true
    #     @test is_unitary(W_copy) == true #LT models are unitary
    # end

end