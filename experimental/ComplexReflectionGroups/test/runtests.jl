@testset "Complex reflection groups" begin

    #######################################################################################
    # Complex reflections
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

    # Example 8 (matrix whose fixed space is a hyperplane and which is diagonalizble but
    # which is not of finite order)
    w = matrix(QQ,2,2,[2 0; 0 1])
    @test is_complex_reflection(t) == false

    #######################################################################################
    # is_complex_reflection_group
    #######################################################################################

    # The following is a fun example I found by brute force search: the following two
    # elements are contained in complex_reflection_group(4, :Magma) and they generate this
    # group but they are not reflections (ker(I-g) = 0). But nontheless the group they
    # generate is G4, so a complex reflection group. So, this is a nice test for 
    # is_complex_reflection_group.
    K,z = cyclotomic_field(3)
    g1 = matrix(K,2,2,[-1 -z-1; 0 -z])
    g2 = matrix(K,2,2,[0 -1; 1 0])
    @test is_complex_reflection(g1) == false
    @test is_complex_reflection(g2) == false
    G = matrix_group([g1,g2])
    @test is_complex_reflection_group(G) == true

    #######################################################################################
    # Exceptional complex reflection groups
    #######################################################################################

    # Note: the complex_reflection_group constructor sets attributes like order and
    # is_complex_reflection_group etc. To perform an honest test, we create a copy of the
    # matrix group which then has no attributes set and run the tests on this copy.

    for model in [:CHEVIE, :Magma, :LT]
        for n=4:37

            #LT not yet implemented for n > 23
            if model == :LT && n > 23
                continue
            end

            W = complex_reflection_group(n, model)
            W_copy = matrix_group(gens(W))
            @test order(W) == order(W_copy)
            @test is_complex_reflection_group(W_copy) == true

            # LT models are unitary
            if model == :LT
                @test is_unitary(W_copy) == true
            end

            # Test in a few cases if we find all reflections (would be too much work
            # for the big groups).
            if n in [4,12,23,25,28]
                @test length(complex_reflections(W_copy)) == number_of_reflections(complex_reflection_group_type(W))
            end
        end
    end

end
