# Tests the type catalogue, concrete models, and reflection geometry.
@testset "Complex reflection groups" begin

    function fixed_space_dimension_from_generators(W)
        n = degree(W)
        isempty(gens(W)) && return n
        K = base_ring(W)
        identity = identity_matrix(K, n)
        equations = reduce(hcat, [matrix(g) - identity for g in gens(W)])
        return n - rank(equations)
    end

    @testset "Reflection group types" begin
        trivial_type = ComplexReflectionGroupType([])
        @test rank(trivial_type) == 0
        @test fixed_space_dimension(trivial_type) == 0
        @test degree(trivial_type) == 0

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

    @testset "Ambient reflection group types" begin
        trivial_zero = complex_reflection_group_type()
        trivial_one = complex_reflection_group_type(1, 1, 1)
        trivial_two = complex_reflection_group_type(; fixed_space_dimension=2)

        @test trivial_zero == complex_reflection_group_type([])
        @test trivial_one ==
              complex_reflection_group_type(; fixed_space_dimension=1)
        @test trivial_two ==
              complex_reflection_group_type([(1, 1, 1), (1, 1, 1)])
        @test_throws ArgumentError complex_reflection_group_type(
            4; fixed_space_dimension=-1
        )
        @test complex_reflection_group_type(7, 7, 1) == trivial_one
        @test complex_reflection_group_type(
            [(1, 1, 1), 4]; fixed_space_dimension=2
        ) == complex_reflection_group_type(4; fixed_space_dimension=3)
        @test is_trivial(trivial_zero)
        @test is_trivial(trivial_one)
        @test is_trivial(trivial_two)
        @test !is_irreducible(trivial_zero)
        @test is_irreducible(trivial_one)
        @test !is_irreducible(trivial_two)
        @test is_primitive(trivial_one)
        @test !is_imprimitive(trivial_one)
        @test_throws ArgumentError coxeter_number(trivial_one)
        @test length(Set([trivial_zero, trivial_one, trivial_two])) == 3
        @test rank.([trivial_zero, trivial_one, trivial_two]) == [0, 0, 0]
        @test fixed_space_dimension.([trivial_zero, trivial_one, trivial_two]) ==
              [0, 1, 2]
        @test degree.([trivial_zero, trivial_one, trivial_two]) == [0, 1, 2]
        @test isempty(degrees(trivial_zero))
        @test degrees(trivial_one) == [1]
        @test degrees(trivial_two) == [1, 1]
        @test isempty(codegrees(trivial_zero))
        @test codegrees(trivial_one) == [-1]
        @test codegrees(trivial_two) == [-1, -1]
        @test exponents(trivial_two) == [0, 0]
        @test coexponents(trivial_two) == [0, 0]
        @test number_of_hyperplanes(trivial_two) == 0
        @test number_of_reflections(trivial_two) == 0
        @test order(center(trivial_two)) == 1

        symmetric_reflection_type = complex_reflection_group_type(1, 1, 3)
        symmetric_permutation_type = complex_reflection_group_type(
            1, 1, 3; fixed_space_dimension=1
        )
        @test symmetric_permutation_type ==
              direct_product(symmetric_reflection_type, trivial_one)
        @test !is_equivalent(
            symmetric_reflection_type, symmetric_permutation_type
        )
        @test rank(symmetric_reflection_type) ==
              rank(symmetric_permutation_type) == 2
        @test fixed_space_dimension(symmetric_reflection_type) == 0
        @test fixed_space_dimension(symmetric_permutation_type) == 1
        @test degree(symmetric_reflection_type) == 2
        @test degree(symmetric_permutation_type) == 3
        @test order(symmetric_reflection_type) ==
              order(symmetric_permutation_type) == 6
        @test degrees(symmetric_reflection_type) == [2, 3]
        @test degrees(symmetric_permutation_type) == [1, 2, 3]
        @test codegrees(symmetric_reflection_type) == [0, 1]
        @test codegrees(symmetric_permutation_type) == [-1, 0, 1]

        irreducible_dihedral_type = complex_reflection_group_type(6, 6, 2)
        reducible_dihedral_type = complex_reflection_group_type(
            [(1, 1, 3), (2, 1, 1)]
        )
        @test order(irreducible_dihedral_type) ==
              order(reducible_dihedral_type) == 12
        @test rank(irreducible_dihedral_type) == 2
        @test rank(reducible_dihedral_type) == 3
        @test degree(irreducible_dihedral_type) == 2
        @test degree(reducible_dihedral_type) == 3
        @test is_irreducible(irreducible_dihedral_type)
        @test !is_irreducible(reducible_dihedral_type)
        @test !is_equivalent(
            irreducible_dihedral_type, reducible_dihedral_type
        )

        fixed_G4_type = complex_reflection_group_type(
            4; fixed_space_dimension=2
        )
        @test rank(fixed_G4_type) == 2
        @test fixed_space_dimension(fixed_G4_type) == 2
        @test degree(fixed_G4_type) == 4
        @test degrees(fixed_G4_type) == [1, 1, 4, 6]
        @test codegrees(fixed_G4_type) == [-1, -1, 0, 2]
        @test exponents(fixed_G4_type) == [0, 0, 3, 5]
        @test coexponents(fixed_G4_type) == [0, 0, 1, 3]
        @test length(degrees(fixed_G4_type)) == degree(fixed_G4_type)
        @test length(codegrees(fixed_G4_type)) == degree(fixed_G4_type)
        @test prod(degrees(fixed_G4_type); init=ZZ(1)) == order(fixed_G4_type)
        @test sum(exponents(fixed_G4_type); init=ZZ(0)) ==
              number_of_reflections(fixed_G4_type)
        @test sum(coexponents(fixed_G4_type); init=ZZ(0)) ==
              number_of_hyperplanes(fixed_G4_type)
        @test number_of_hyperplanes(fixed_G4_type) == 4
        @test number_of_reflections(fixed_G4_type) == 8
        @test order(center(fixed_G4_type)) == 2

        product_with_fixed_space = direct_product([
            complex_reflection_group_type(4; fixed_space_dimension=1),
            complex_reflection_group_type(8; fixed_space_dimension=2),
        ])
        @test product_with_fixed_space ==
              complex_reflection_group_type([4, 8]; fixed_space_dimension=3)
        @test number_of_components(product_with_fixed_space) == 2
        @test fixed_space_dimension(product_with_fixed_space) == 3
        @test rank(product_with_fixed_space) == 4
        @test degree(product_with_fixed_space) == 7

        @test is_well_generated(complex_reflection_group_type(6, 6, 2))
        @test is_well_generated(complex_reflection_group_type(8, 8, 3))
        @test is_well_generated(complex_reflection_group_type(6, 1, 3))
        @test !is_well_generated(complex_reflection_group_type(6, 3, 2))
        @test !is_well_generated(complex_reflection_group_type(4, 2, 3))
    end

    @testset "Ambient matrix realizations" begin
        trivial_types = [
            complex_reflection_group_type(),
            complex_reflection_group_type(1, 1, 1),
            complex_reflection_group_type(; fixed_space_dimension=2),
        ]
        trivial_groups = complex_reflection_group.(trivial_types)
        @test complex_reflection_group_type(complex_reflection_group(
            ; fixed_space_dimension=2
        )) == last(trivial_types)

        for (n, T, W) in zip(0:2, trivial_types, trivial_groups)
            @test base_ring(W) === QQ
            @test degree(W) == n
            @test isempty(gens(W))
            @test order(W) == 1
            @test complex_reflection_group_type(W) == T
            @test fixed_space_dimension(W) == n
            @test isempty(complex_reflection_group_model(W))
            @test isempty(components(W))
            @test isempty(complex_reflection_group_component_embeddings(W))
            @test fixed_space_dimension_from_generators(W) == n
            @test invariant_hermitian_form(W) == identity_matrix(QQ, n)
            @test isempty(reflection_hyperplane_orbits(W))
            @test isempty(reflection_classes(W))
            @test isempty(reflection_library(W))
            @test isempty(complex_reflections(W))

            WD = complex_reflection_group_dual(W)
            @test degree(WD) == n
            @test isempty(gens(WD))
            @test order(WD) == 1
            @test complex_reflection_group_type(WD) == T
            @test isempty(complex_reflection_group_model(WD))
            @test complex_reflection_group_dual_source(WD) === W
            @test isempty(components(WD))
            @test fixed_space_dimension_from_generators(WD) == n

            doubled_group = symplectic_reflection_group(W)
            @test degree(doubled_group) == 2*n
            @test isempty(gens(doubled_group))
            @test order(doubled_group) == 1
            @test symplectic_doubling_source(doubled_group) === W
            @test symplectic_doubling_source_type(doubled_group) == T
            @test symplectic_doubling_source_model(doubled_group) == Symbol[]
            @test complex_reflection_group_type(doubled_group) === nothing
            @test is_symplectic_reflection_group(doubled_group)
            primal_embedding, dual_embedding =
                symplectic_doubling_block_embeddings(doubled_group)
            @test size(primal_embedding) == (n, 2*n)
            @test size(dual_embedding) == (n, 2*n)
            @test size(gram_matrix(symplectic_form(doubled_group))) ==
                  (2*n, 2*n)
        end

        symmetric_reflection_type = complex_reflection_group_type(1, 1, 3)
        symmetric_permutation_type = complex_reflection_group_type(
            1, 1, 3; fixed_space_dimension=1
        )
        symmetric_reflection_group = complex_reflection_group(
            symmetric_reflection_type
        )
        symmetric_permutation_group = complex_reflection_group(
            symmetric_permutation_type
        )
        @test complex_reflection_group_type(complex_reflection_group(
            1, 1, 3; fixed_space_dimension=1
        )) == symmetric_permutation_type
        @test degree(symmetric_reflection_group) == 2
        @test degree(symmetric_permutation_group) == 3
        @test fixed_space_dimension_from_generators(
            symmetric_reflection_group
        ) == 0
        @test fixed_space_dimension_from_generators(
            symmetric_permutation_group
        ) == 1
        @test complex_reflection_group_type(symmetric_permutation_group) ==
              symmetric_permutation_type
        @test length(components(symmetric_permutation_group)) == 1
        @test degree(only(components(symmetric_permutation_group))) == 2
        @test length(complex_reflections(symmetric_reflection_group)) == 3
        @test length(complex_reflections(symmetric_permutation_group)) == 3

        fixed_G4_type = complex_reflection_group_type(
            4; fixed_space_dimension=2
        )
        fixed_G4 = complex_reflection_group(fixed_G4_type, :LT)
        @test degree(fixed_G4) == 4
        @test fixed_space_dimension_from_generators(fixed_G4) == 2
        @test complex_reflection_group_type(fixed_G4) == fixed_G4_type
        @test complex_reflection_group_model(fixed_G4) == [:LT]
        @test length(components(fixed_G4)) == 1
        @test degree(only(components(fixed_G4))) == 2
        @test length(complex_reflection_group_component_embeddings(fixed_G4)) == 1
        @test invariant_hermitian_form(fixed_G4) ==
              identity_matrix(base_ring(fixed_G4), 4)
        @test length(reflection_hyperplane_orbits(fixed_G4)) == 1
        @test dim(first(reflection_hyperplanes(
            only(reflection_hyperplane_orbits(fixed_G4))
        ))) == 3
        @test length(reflection_classes(fixed_G4)) == 2
        @test length(complex_reflections(fixed_G4)) == 8

        fixed_G4_dual = complex_reflection_group_dual(fixed_G4)
        @test degree(fixed_G4_dual) == 4
        @test complex_reflection_group_type(fixed_G4_dual) == fixed_G4_type
        @test fixed_space_dimension_from_generators(fixed_G4_dual) == 2
        @test length(components(fixed_G4_dual)) == 1
        @test degree(only(components(fixed_G4_dual))) == 2
        @test length(reflection_hyperplane_orbits(fixed_G4_dual)) == 1

        irreducible_dihedral_group = complex_reflection_group(6, 6, 2)
        reducible_dihedral_group = complex_reflection_group(
            [(1, 1, 3), (2, 1, 1)]
        )
        @test order(irreducible_dihedral_group) ==
              order(reducible_dihedral_group) == 12
        @test degree(irreducible_dihedral_group) == 2
        @test degree(reducible_dihedral_group) == 3
        @test fixed_space_dimension_from_generators(
            irreducible_dihedral_group
        ) == 0
        @test fixed_space_dimension_from_generators(
            reducible_dihedral_group
        ) == 0
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

        fixed_cyclotomic_group = complex_reflection_group(
            [4, 8], :CHEVIE; fixed_space_dimension=1
        )
        fixed_cyclotomic_type = complex_reflection_group_type(
            [4, 8]; fixed_space_dimension=1
        )
        @test degree(fixed_cyclotomic_group) == 5
        @test complex_reflection_group_type(fixed_cyclotomic_group) ==
              fixed_cyclotomic_type
        @test fixed_space_dimension_from_generators(fixed_cyclotomic_group) == 1
        @test degree.(components(fixed_cyclotomic_group)) == [2, 2]
        @test length(complex_reflection_group_component_embeddings(
            fixed_cyclotomic_group
        )) == 2
        @test length(reflection_hyperplane_orbits(fixed_cyclotomic_group)) ==
              length(reflection_hyperplane_orbits(fixed_cyclotomic_type))
        fixed_cyclotomic_dual = complex_reflection_group_dual(
            fixed_cyclotomic_group
        )
        @test degree(fixed_cyclotomic_dual) == 5
        @test fixed_space_dimension_from_generators(fixed_cyclotomic_dual) == 1
        @test degree.(components(fixed_cyclotomic_dual)) == [2, 2]
        @test length(reflection_hyperplane_orbits(fixed_cyclotomic_dual)) ==
              length(reflection_hyperplane_orbits(fixed_cyclotomic_type))

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

    @testset "Absolute coefficient-field presentations" begin
        # Field isomorphisms transport the LT matrices and their reference marking.
        flat = complex_reflection_group(4, :LT; field_presentation=:absolute)
        native = only(components(flat))
        embedding = only(complex_reflection_group_component_embeddings(flat))
        @test base_ring(native) isa Hecke.RelSimpleNumField
        @test base_ring(flat) isa AbsSimpleNumField
        @test absolute_degree(base_ring(flat)) == 4
        @test domain(embedding) === base_ring(native)
        @test codomain(embedding) === base_ring(flat)
        @test complex_reflection_group_model(flat) == [:LT]
        @test complex_reflection_group_type(flat) == complex_reflection_group_type(4)
        @test is_unitary(flat)
        @test is_one(invariant_hermitian_form(flat))
        conjugation_native = complex_conjugation(base_ring(native))
        conjugation_flat = complex_conjugation(base_ring(flat))
        for (g, h) in zip(gens(native), gens(flat))
            @test map_entries(embedding, matrix(g)) == matrix(h)
            @test map_entries(inv(embedding), matrix(h)) == matrix(g)
            @test all(
                conjugation_flat(embedding(x)) == embedding(conjugation_native(x))
                for x in matrix(g)
            )
        end
        for (C, D) in zip(reflection_classes(native), reflection_classes(flat))
            @test reflection_class_type(C) == reflection_class_type(D)
            @test Set(map_entries(embedding, matrix(s)) for s in elements(C)) ==
                  Set(matrix(s) for s in elements(D))
        end

        for alias in (:LT_flat, Symbol("LT-flat"))
            W = complex_reflection_group(4, alias)
            @test base_ring(W) isa AbsSimpleNumField
            @test complex_reflection_group_model(W) == [:LT]
            @test is_unitary(W)
        end
        for W in (
            complex_reflection_group(; fixed_space_dimension=2, field_presentation=:absolute),
            complex_reflection_group(4, 2, 2, :LT; field_presentation=:absolute),
            complex_reflection_group([4, 8], :LT; field_presentation=:absolute),
            complex_reflection_group(
                complex_reflection_group_type([4, 8]; fixed_space_dimension=1), :LT;
                field_presentation=:absolute,
            ),
        )
            @test base_ring(W) isa Union{QQField, AbsSimpleNumField}
            @test is_unitary(W)
            @test fixed_space_dimension_from_generators(W) == fixed_space_dimension(W)
            component_embeddings = complex_reflection_group_component_embeddings(W)
            @test all(codomain(f) === base_ring(W) for f in component_embeddings)
            @test all(
                domain(f) === base_ring(C) for
                (f, C) in zip(component_embeddings, components(W))
            )
            @test sum(number_of_hyperplanes, reflection_hyperplane_orbits(W); init=0) ==
                  number_of_hyperplanes(complex_reflection_group_type(W))
        end
        for (types, field_degree) in (([4, 16], 16), ([27, 33], 8))
            W = complex_reflection_group(types, :LT; field_presentation=:absolute)
            @test base_ring(W) isa AbsSimpleNumField
            @test absolute_degree(base_ring(W)) == field_degree
            @test is_unitary(W)
            conjugation = complex_conjugation(base_ring(W))
            for (C, f) in zip(components(W), complex_reflection_group_component_embeddings(W))
                source_conjugation = complex_conjugation(base_ring(C))
                @test all(
                    conjugation(f(x)) == f(source_conjugation(x))
                    for g in gens(C) for x in matrix(g)
                )
            end
        end
        @test_throws ArgumentError complex_reflection_group(
            4, :LT; field_presentation=:unknown
        )
    end

    @testset "Higher-rank Lehrer--Taylor models" begin
        # Independently tabulated line counts and stabilizer orders: LT09, Table D.2, p. 273.
        dimensions = [3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 6, 6, 7, 8]
        orbit_data = [
            [(2, 21)], [(3, 12)], [(2, 9), (3, 12)], [(2, 45)],
            [(2, 12), (2, 12)], [(2, 40)], [(2, 60)], [(2, 60)],
            [(3, 40)], [(2, 45)], [(2, 126)], [(2, 36)], [(2, 63)], [(2, 120)],
        ]
        for n in 24:37
            W = complex_reflection_group(n, :LT)
            @test degree(W) == dimensions[n - 23]
            @test is_unitary(W)
            @test is_one(invariant_hermitian_form(W))
            orbits = reflection_hyperplane_orbits(W)
            @test [(pointwise_stabilizer_order(O), length(O)) for O in orbits] ==
                  orbit_data[n - 23]
            @test sum(length, reflection_classes(W)) ==
                  sum((e - 1)*h for (e, h) in orbit_data[n - 23])
            if n in (25, 26, 32)
                K = base_ring(W)
                omega = gen(K)
                @test all(
                    eigenvalue(complex_reflection(representative(O))) == omega
                    for O in orbits if pointwise_stabilizer_order(O) == 3
                )
            end
        end

        # Explicit root-basis changes intertwine the LT and arithmetic Magma generators.
        for n in (25, 28, 32, 35, 36)
            W = complex_reflection_group(n, :LT)
            M = complex_reflection_group(n, :Magma)
            K = base_ring(W)
            if n == 28
                embedding = id_hom(QQ)
                basis_change = matrix(QQ, [
                    0 1 -1 0; 0 0 1 -1; 0 0 0 1; 1//2 -1//2 -1//2 -1//2
                ])
                indices = [1, 2, 3, 4]
            elseif n in (25, 32)
                omega = gen(K)
                theta = omega - omega^2
                embedding = hom(base_ring(M), K, omega)
                if n == 25
                    basis_change = matrix(K, [theta 0 0; 1 1 1; 0 theta 0])
                    indices = [1, 3, 2]
                else
                    basis_change = matrix(K, [
                        theta 0 0 0; -1 -1 -1 0; 0 theta 0 0; 0 -1 1 1
                    ])
                    indices = [1, 3, 2, 4]
                end
            else
                embedding = hom(QQ, K)
                radical = gen(K)
                if n == 35
                    basis_change = matrix(K, [
                        1 -1 0 0 0 0;
                        0 0 0 1 -1 0;
                        0 1 -1 0 0 0;
                        0 0 1 -1 0 0;
                        0 0 0 1 1 0;
                        -1//2 -1//2 -1//2 -1//2 -1//2 -radical//2
                    ])
                    indices = [1, 4, 2, 3, 5, 6]
                else
                    basis_change = matrix(K, [
                        -1//2 -1//2 -1//2 -1//2 -1//2 -1//2 -radical//2;
                        0 0 0 0 1 -1 0;
                        0 0 0 0 1 1 0;
                        0 0 0 1 -1 0 0;
                        0 0 1 -1 0 0 0;
                        0 1 -1 0 0 0 0;
                        1 -1 0 0 0 0 0
                    ])
                    indices = [8, 5, 6, 4, 3, 2, 1]
                end
            end
            @test !is_zero(det(basis_change))
            @test length(gens(M)) == length(indices)
            @test all(
                map_entries(embedding, matrix(g))*basis_change ==
                basis_change*matrix(gens(W)[j])
                for (g, j) in zip(gens(M), indices)
            )
        end

        for n in (2, 3, 4, 5, 6, 9)
            W = complex_reflection_group(1, 1, n, :LT)
            raw = matrix_group(gens(W))
            @test degree(W) == n - 1
            @test order(raw) == factorial(ZZ(n))
            @test is_unitary(W)
            @test (base_ring(W) isa QQField) == (n == 2 || is_square(ZZ(n)))
            if !(base_ring(W) isa QQField)
                @test absolute_degree(base_ring(W)) == 2
            end
            @test fixed_space_dimension_from_generators(W) == 0
            @test length(only(reflection_hyperplane_orbits(W))) == binomial(n, 2)
            generators = gens(W)
            @test all(is_one(g^2) for g in generators)
            @test all(is_one((generators[i]*generators[i + 1])^3) for i in 1:(n - 2))
            @test all(
                is_one((generators[i]*generators[j])^2)
                for i in 1:(n - 1) for j in (i + 2):(n - 1)
            )
        end
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
        algebraic_closure_identity = identity_matrix(algebraic_closure(QQ), 2)
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
