# Tests the type catalogue, concrete models, and reflection geometry.
@testset "Complex reflection groups" begin

    @testset "Reflection group types" begin
        trivial_type = ComplexReflectionGroupType([])
        @test rank(trivial_type) == 0

        product_type = complex_reflection_group_type([4, (2, 1, 3)])
        @test rank(product_type) == 5
        @test product_type == complex_reflection_group_type([(2, 1, 3), 4])
        @test hash(product_type) == hash(complex_reflection_group_type([(2, 1, 3), 4]))

        @test complex_reflection_group_type(2, 2, 2) ==
              complex_reflection_group_type([(1, 1, 2), (1, 1, 2)])
        @test number_of_reflection_classes(complex_reflection_group_type(5, 1, 1)) == 4
        @test is_primitive(complex_reflection_group_type(5, 1, 1))
        @test is_imprimitive(complex_reflection_group_type(1, 1, 4))
        @test is_primitive(complex_reflection_group_type(1, 1, 5))
        @test_throws ArgumentError is_primitive(product_type)

        G4_type = complex_reflection_group_type(4)
        G4_models = [
            complex_reflection_group(4, model) for model in [:CHEVIE, :LT, :Magma]
        ]
        @test all(complex_reflection_group_type(W) == G4_type for W in G4_models)
        @test all(
            is_equivalent(complex_reflection_group_type(W), G4_type) for W in G4_models
        )
    end

    @testset "Type-level reflection data" begin
        T = complex_reflection_group_type(11)
        orbits = reflection_hyperplane_orbits(T)
        @test [(pointwise_stabilizer_order(O), number_of_hyperplanes(O)) for O in orbits] ==
              [(2, 12), (3, 8), (4, 6)]
        @test representative_word.(orbits) == [(1,), (2,), (3,)]
        @test reflection_exponent.(reflection_classes(T)) == [1, 1, 2, 1, 2, 3]
        @test orbits == reflection_hyperplane_orbits(complex_reflection_group_type(11))
        @test hash.(orbits) ==
              hash.(reflection_hyperplane_orbits(complex_reflection_group_type(11)))

        split_type = complex_reflection_group_type(6, 6, 2)
        split_orbits = reflection_hyperplane_orbits(split_type)
        @test [(pointwise_stabilizer_order(O), number_of_hyperplanes(O)) for O in split_orbits] ==
              [(2, 3), (2, 3)]
        @test length(reflection_classes(split_type)) == 2
    end

    @testset "Reducible matrix realizations" begin
        cyclotomic_group = complex_reflection_group([4, 8], :CHEVIE)
        @test cyclotomic_group isa MatGroup
        @test degree(cyclotomic_group) == 4
        @test [order(g) for g in gens(cyclotomic_group)] == [3, 3, 4, 4]
        is_cyclotomic, conductor =
            Hecke.is_cyclotomic_type(base_ring(cyclotomic_group))
        @test is_cyclotomic
        @test conductor == 12
        @test complex_reflection_group_model(cyclotomic_group) == [:CHEVIE, :CHEVIE]

        tower_group = complex_reflection_group([4, 8], :LT)
        factors = components(tower_group)
        embeddings = complex_reflection_group_component_embeddings(tower_group)
        @test degree(tower_group) == 4
        @test base_ring(tower_group) === base_ring(factors[1])
        @test base_field(base_ring(tower_group)) === base_ring(factors[2])
        @test all(codomain(f) === base_ring(tower_group) for f in embeddings)
        @test [order(g) for g in gens(tower_group)] == [3, 3, 4, 4]
        @test is_one(invariant_hermitian_form(tower_group))

        compositum_group = complex_reflection_group([4, 16], :LT)
        compositum_factors = components(compositum_group)
        compositum_embeddings =
            complex_reflection_group_component_embeddings(compositum_group)
        @test degree(compositum_group) == 4
        @test absolute_degree(base_ring(compositum_group)) == 16
        @test all(
            base_ring(compositum_group) !== base_ring(C) for
            C in compositum_factors
        )
        @test all(
            codomain(f) === base_ring(compositum_group) for
            f in compositum_embeddings
        )
        @test all(
            domain(f) === base_ring(C) for
            (f, C) in zip(compositum_embeddings, compositum_factors)
        )
        @test [order(g) for g in gens(compositum_group)] == [3, 3, 5, 5]
        @test is_unitary(
            compositum_group,
            invariant_hermitian_form(compositum_group),
        )
        @test length(complex_reflections(compositum_group)) ==
              number_of_reflections(complex_reflection_group_type(compositum_group))
    end

    @testset "Concrete reflection data" begin
        W = complex_reflection_group(11, :LT)
        K = base_ring(W)
        orbits = reflection_hyperplane_orbits(W)
        @test base_ring(W) === K
        @test [(pointwise_stabilizer_order(O), number_of_hyperplanes(O)) for O in orbits] ==
              [(2, 12), (3, 8), (4, 6)]
        classes = reflection_classes(W)
        @test length.(classes) == [12, 8, 8, 6, 6, 6]
        library = reflection_library(W)
        @test length.(library) == [12, 8, 6]
        class_offset = 0
        for (O, orbit_library) in zip(orbits, library)
            orbit_classes = classes[
                (class_offset + 1):(class_offset + pointwise_stabilizer_order(O) - 1)
            ]
            for (H, hyperplane_reflections) in
                zip(reflection_hyperplanes(O), orbit_library)
                @test all(
                    s -> hyperplane(complex_reflection(s)) == H,
                    hyperplane_reflections,
                )
                @test all(
                    s in elements(C) for
                    (s, C) in zip(hyperplane_reflections, orbit_classes)
                )
            end
            class_offset += pointwise_stabilizer_order(O) - 1
        end
        @test length(complex_reflections(W)) == 46

        J = invariant_hermitian_form(W)
        @test J !== nothing
        @test is_unitary(W, J)
        @test !is_unitary(W, zero_matrix(K, degree(W), degree(W)))
        nonhermitian_matrix = matrix(K, 2, 2, [1 1; 0 1])
        @test !is_unitary(W, nonhermitian_matrix)
        algebraic_closure_identity = identity_matrix(QQBar, 2)
        @test is_unitary(algebraic_closure_identity)
        @test is_unitary(algebraic_closure_identity, algebraic_closure_identity)
        @test is_unitary(matrix_group([algebraic_closure_identity]))

        cyclic_group = complex_reflection_group(5, 1, 1, :CHEVIE)
        cyclic_orbits = reflection_hyperplane_orbits(cyclic_group)
        @test length(cyclic_orbits) == 1
        @test number_of_hyperplanes(only(cyclic_orbits)) == 1
        @test pointwise_stabilizer_order(only(cyclic_orbits)) == 5
        @test length(reflection_classes(cyclic_group)) == 4
        @test length(complex_reflections(cyclic_group)) == 4
        @test isempty(hyperplane_basis(complex_reflection(first(gens(cyclic_group)))))

        split_group = complex_reflection_group(6, 6, 2, :CHEVIE)
        @test number_of_hyperplanes.(reflection_hyperplane_orbits(split_group)) == [3, 3]

        reordered_group = complex_reflection_group(10, :Magma)
        @test [(pointwise_stabilizer_order(O), number_of_hyperplanes(O)) for O in reflection_hyperplane_orbits(reordered_group)] ==
              [(3, 8), (4, 6)]

        # All three G(4, 2, 2) orbits have the same order and size, so their
        # reference marking cannot be recovered from those invariants alone.
        for (model, generator_indices) in
            [(:CHEVIE, [1, 2, 3]), (:LT, [2, 1, 3]), (:Magma, [3, 2, 1])]
            marked_group = complex_reflection_group(4, 2, 2, model)
            @test representative.(reflection_hyperplane_orbits(marked_group)) ==
                  gens(marked_group)[generator_indices]
        end

        custom_group = matrix_group(gens(complex_reflection_group(4, 2, 2, :Magma)))
        set_attribute!(
            custom_group,
            :complex_reflection_group_type,
            complex_reflection_group_type(4, 2, 2),
        )
        custom_marking = gens(custom_group)[[3, 2, 1]]
        set_reflection_hyperplane_orbit_marking!(custom_group, custom_marking)
        @test representative.(reflection_hyperplane_orbits(custom_group)) == custom_marking

        for model in [:CHEVIE, :Magma]
            G5 = complex_reflection_group(5, model)
            G7 = complex_reflection_group(7, model)
            @test representative.(reflection_hyperplane_orbits(G5)) == gens(G5)[[1, 2]]
            @test representative.(reflection_hyperplane_orbits(G7)) ==
                  gens(G7)[[1, 2, 3]]
        end
        G5 = complex_reflection_group(5, :LT)
        G7 = complex_reflection_group(7, :LT)
        @test representative.(reflection_hyperplane_orbits(G5)) ==
              [gens(G5)[1]^2, gens(G5)[2]^2]
        @test representative.(reflection_hyperplane_orbits(G7)) ==
              [gens(G7)[1], gens(G7)[2]^2, gens(G7)[3]^2]

        for model in [:CHEVIE, :Magma]
            G28 = complex_reflection_group(28, model)
            @test representative.(reflection_hyperplane_orbits(G28)) ==
                  gens(G28)[[1, 3]]
        end

        G4_LT = complex_reflection_group(4, :LT)
        @test representative(only(reflection_hyperplane_orbits(G4_LT))) ==
              gens(G4_LT)[1]^2
        G18_LT = complex_reflection_group(18, :LT)
        @test representative.(reflection_hyperplane_orbits(G18_LT)) ==
              [gens(G18_LT)[1], gens(G18_LT)[2]^4]
        G9_Magma = complex_reflection_group(9, :Magma)
        @test representative.(reflection_hyperplane_orbits(G9_Magma)) ==
              [gens(G9_Magma)[1], gens(G9_Magma)[2]^3]
        G10_Magma = complex_reflection_group(10, :Magma)
        @test representative.(reflection_hyperplane_orbits(G10_Magma)) ==
              [gens(G10_Magma)[2], gens(G10_Magma)[1]^3]
        G11_Magma = complex_reflection_group(11, :Magma)
        @test representative.(reflection_hyperplane_orbits(G11_Magma)) ==
              [gens(G11_Magma)[2], gens(G11_Magma)[3], gens(G11_Magma)[1]^3]

        custom_LT = matrix_group(gens(G4_LT))
        set_attribute!(
            custom_LT,
            :complex_reflection_group_type,
            complex_reflection_group_type(4),
        )
        set_reflection_hyperplane_orbit_marking!(
            custom_LT,
            [gens(custom_LT)[1]];
            class_exponents=[ZZ(2)],
        )
        @test representative(only(reflection_hyperplane_orbits(custom_LT))) ==
              gens(custom_LT)[1]^2

        unmarked_group = matrix_group(gens(G4_LT))
        set_attribute!(
            unmarked_group,
            :complex_reflection_group_type,
            complex_reflection_group_type(4),
        )
        @test_throws ArgumentError reflection_hyperplane_orbits(unmarked_group)

        cyclic_order_four = complex_reflection_group(4, 1, 1, :CHEVIE)
        @test_throws ArgumentError set_reflection_hyperplane_orbit_marking!(
            cyclic_order_four,
            [gens(cyclic_order_four)[1]];
            class_exponents=[2],
        )

        dual_group = complex_reflection_group_dual(G4_LT)
        @test complex_reflection_group_model(dual_group) == [:dual]
        @test complex_reflection_group_dual_source(dual_group) === G4_LT
        @test complex_reflection_group_dual_source_model(dual_group) == [:LT]
        expected_dual_reflections = [
            dual_group(transpose(matrix(representative(O)))) for
            O in reflection_hyperplane_orbits(G4_LT)
        ]
        @test representative.(reflection_hyperplane_orbits(dual_group)) ==
              expected_dual_reflections

        unmarked_dual_group = complex_reflection_group_dual(unmarked_group)
        @test complex_reflection_group_type(unmarked_dual_group) ==
              complex_reflection_group_type(unmarked_group)
        @test complex_reflection_group_model(unmarked_dual_group) == [:dual]
        @test_throws ArgumentError reflection_hyperplane_orbits(unmarked_dual_group)
    end

    @testset "Symplectic doubling" begin
        W = complex_reflection_group(4, :LT)
        doubled_group = symplectic_reflection_group(W)
        K = base_ring(W)
        n = degree(W)

        @test symplectic_doubling_source(doubled_group) === W
        @test symplectic_doubling_source_type(doubled_group) ==
              complex_reflection_group_type(W)
        @test symplectic_doubling_source_model(doubled_group) ==
              complex_reflection_group_model(W)
        @test complex_reflection_group_type(doubled_group) === nothing
        @test is_symplectic_reflection_group(doubled_group)
        @test base_ring(doubled_group) === K
        @test degree(doubled_group) == 2*n
        @test order(doubled_group) == order(W)

        primal_embedding, dual_embedding =
            symplectic_doubling_block_embeddings(doubled_group)
        omega = gram_matrix(symplectic_form(doubled_group))
        identity_block = identity_matrix(K, n)
        zero_block = zero_matrix(K, n, n)

        @test is_alternating(symplectic_form(doubled_group))
        @test primal_embedding*omega*transpose(primal_embedding) == zero_block
        @test dual_embedding*omega*transpose(dual_embedding) == zero_block
        @test primal_embedding*omega*transpose(dual_embedding) == identity_block

        for (g, doubled_generator) in zip(gens(W), gens(doubled_group))
            g_matrix = matrix(g)
            doubled_matrix = matrix(doubled_generator)
            @test doubled_matrix*omega*transpose(doubled_matrix) == omega
            @test primal_embedding*doubled_matrix == g_matrix*primal_embedding
            @test dual_embedding*doubled_matrix ==
                  transpose(inv(g_matrix))*dual_embedding
            @test rank(doubled_matrix - identity_matrix(K, 2*n)) == 2
        end
    end

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
    @test is_complex_reflection(w) == false

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

    typed_group = matrix_group([g1, g2])
    set_attribute!(
        typed_group,
        :complex_reflection_group_type,
        complex_reflection_group_type(4),
    )
    typed_dual_group = complex_reflection_group_dual(typed_group)
    @test complex_reflection_group_type(typed_dual_group) ==
          complex_reflection_group_type(typed_group)
    @test !has_attribute(typed_group, :is_complex_reflection_group)

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
