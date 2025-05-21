@testset "LieAlgebras.WeylGroup" begin
  function is_in_normal_form(x::WeylGroupElem)
    return word(parent(x)(word(x))) == word(x)
  end

  _isomorphic_group_on_gens = Oscar.LieAlgebras._isomorphic_group_on_gens

  @testset "apply_braid_move!" begin
    @test apply_braid_move!(UInt8[1, 3, 2, 1], (1, 2, 0)) == [3, 1, 2, 1]
    @test apply_braid_move!(UInt8[2, 3, 1, 4], (2, 2, 0)) == [2, 1, 3, 4]
    @test apply_braid_move!(UInt8[2, 3, 1, 2, 1, 4], (3, 3, -1)) == [2, 3, 2, 1, 2, 4]
    @test apply_braid_move!(UInt8[2, 3, 4, 3, 4, 2], (2, 4, -2)) == [2, 4, 3, 4, 3, 2]
    @test apply_braid_move!(UInt8[1, 2, 1, 2, 1, 2], (1, 6, -3)) == [2, 1, 2, 1, 2, 1]
  end

  @testset "braid_moves for $fam$rk" for (fam, rk) in [
    (:A, 3), (:B, 3), (:C, 3), (:D, 4), (:F, 4), (:G, 2)
  ]
    W = weyl_group(fam, rk)
    for _ in 1:10
      w = rand(W)
      w1 = word(w)
      w2 = copy(w1)

      i = 0
      n = rand(1:(10 * order(W)))
      for r in reduced_expressions(w)
        i += 1
        copy!(w2, r)
        if i == n
          break
        end
      end

      mvs = braid_moves(W, w1, w2)
      @test !isnothing(mvs)
      for mv in mvs
        apply_braid_move!(w2, mv)
      end
      @test w1 == w2
    end
  end

  # TODO: merge with conformance tests in test/LieTheory/WeylGroup.jl once this is moved to src
  @testset "WeylGroup Group isomorphism test for $(Wname)" for (Wname, W) in [
    ("A1", weyl_group(:A, 1)),
    (
      "A3 with non-canonical ordering of simple roots",
      weyl_group(root_system([2 -1 -1; -1 2 0; -1 0 2])),
    ),
    ("A5", weyl_group(:A, 5)),
    (
      "B4 with non-canonical ordering of simple roots",
      weyl_group(root_system([2 -1 -1 0; -1 2 0 -1; -2 0 2 0; 0 -1 0 2])),
    ),
    ("B4", weyl_group(root_system(:B, 4))),
    ("C3", weyl_group(:C, 3)),
    ("D5", weyl_group(cartan_matrix(:D, 5))),
    ("E6", weyl_group(:E, 6)),
    ("E7", weyl_group(:E, 7)),
    ("E8", weyl_group(:E, 8)),
    ("F4", weyl_group(:F, 4)),
    ("G2", weyl_group(:G, 2)),
    ("F4+G2", weyl_group((:F, 4), (:G, 2))),
    ("E6+C3", weyl_group([(:E, 6), (:C, 3)])),
    ("A_1^(1)", weyl_group(ZZ[2 -2; -2 2])), # TODO: replace with cartan_matrix(A_1^(1)), once functionality for affine type is added
    (
      "complicated case 1",
      begin
        cm = cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2))
        for _ in 1:50
          i, j = rand(1:nrows(cm), 2)
          if i != j
            swap_rows!(cm, i, j)
            swap_cols!(cm, i, j)
          end
        end
        weyl_group(cm)
      end,
    ),
    (
      "complicated case 2",
      begin
        cm = cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2))
        for _ in 1:50
          i, j = rand(1:nrows(cm), 2)
          if i != j
            swap_rows!(cm, i, j)
            swap_cols!(cm, i, j)
          end
        end
        weyl_group(root_system(cm))
      end,
    ),
  ]
    @testset "isomorphism(FPGroup, ::WeylGroup)" begin
      if is_finite(W) && ngens(W) < 6 #= for sane runtime =#
        G = _isomorphic_group_on_gens(FPGroup, W)
        @test is_finite(G) == is_finite(W)
        is_finite(W) && @test order(G) == order(W)
      end

      G = fp_group(W)
      @test is_finite(G) == is_finite(W)
      is_finite(W) && @test order(G) == order(W)

      iso = isomorphism(FPGroup, W)
      @test W === domain(iso)
      @test G === codomain(iso)
      @test iso === isomorphism(FPGroup, W) # test caching

      # test mapping of gens
      @test ngens(W) == ngens(G)
      @test iso.(gens(W)) == gens(G)
      @test inv(iso).(gens(G)) == gens(W)

      # test general mapping
      if ngens(W) < 10 #= for sane runtime =#
        for _ in 1:5
          if is_finite(W) # remove once rand(W) is implemented for infinite groups
            w = rand(W)
            @test w == inv(iso)(iso(w))
            v = rand(W)
            @test iso(v * w) == iso(v) * iso(w)
            @test v * w == inv(iso)(iso(v) * iso(w))
          end
          g = rand_pseudo(G)
          @test is_in_normal_form(inv(iso)(g))
          @test g == iso(inv(iso)(g))
          h = rand_pseudo(G)
          @test inv(iso)(h * g) == inv(iso)(h) * inv(iso)(g)
          @test h * g == iso(inv(iso)(h) * inv(iso)(g))
        end
      end
    end

    if has_root_system_type(root_system(W))
      type, ordering = root_system_type_with_ordering(root_system(W))
      @testset "isomorphism(PermGroup, ::WeylGroup)" begin
        if is_finite(W) && ngens(W) < 6 #= for sane runtime =#
          G = _isomorphic_group_on_gens(PermGroup, W)
          @test is_finite(G) == is_finite(W)
          is_finite(W) && @test order(G) == order(W)
        end

        G = permutation_group(W)
        @test is_finite(G) == is_finite(W)
        is_finite(W) && @test order(G) == order(W)

        iso = isomorphism(PermGroup, W)
        @test W === domain(iso)
        @test G === codomain(iso)
        @test iso === isomorphism(PermGroup, W) # test caching

        # test mapping of gens
        @test ngens(W) == ngens(G)
        @test iso.(gens(W)) == gens(G)
        @test inv(iso).(gens(G)) == gens(W)

        # test general mapping
        if ngens(W) < 10 #= for sane runtime =#
          for _ in 1:5
            if is_finite(W) # remove once rand(W) is implemented for infinite groups
              w = rand(W)
              @test w == inv(iso)(iso(w))
              v = rand(W)
              @test iso(v * w) == iso(v) * iso(w)
              @test v * w == inv(iso)(iso(v) * iso(w))
            end
            g = rand_pseudo(G)
            @test is_in_normal_form(inv(iso)(g))
            @test g == iso(inv(iso)(g))
            h = rand_pseudo(G)
            @test inv(iso)(h * g) == inv(iso)(h) * inv(iso)(g)
            @test h * g == iso(inv(iso)(h) * inv(iso)(g))
          end
        end
      end
    end
  end

  @testset "WeylGroup parabolic subgroup test for $(Wname)" for (
    Wname, W, vec, check_proj
  ) in [
    ("A1", weyl_group(:A, 1), [1], true),
    ("A5", weyl_group(:A, 5), [1, 5, 3], false),
    ("B2", weyl_group(:B, 2), [2, 1], false),
    ("F4", weyl_group(:F, 4), [2, 3], false),
    ("A5+E8+D4", weyl_group((:A, 5), (:E, 8), (:D, 4)), [6:13; 1:5], true),
    (
      "A3 with non-canonical ordering of simple roots",
      weyl_group(root_system([2 -1 -1; -1 2 0; -1 0 2])),
      [2, 3], false,
    ),
    (
      "B4 with non-canonical ordering of simple roots",
      weyl_group(root_system([2 -1 -1 0; -1 2 0 -1; -2 0 2 0; 0 -1 0 2])),
      [2, 4], false,
    ),
    (
      "complicated case 1",
      begin
        cm = cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2))
        for _ in 1:50
          i, j = rand(1:nrows(cm), 2)
          if i != j
            swap_rows!(cm, i, j)
            swap_cols!(cm, i, j)
          end
        end
        weyl_group(cm)
      end,
      unique(rand(1:14, 5)), false,
    ),
    (
      "complicated case 2",
      begin
        cm = cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2))
        for _ in 1:50
          i, j = rand(1:nrows(cm), 2)
          if i != j
            swap_rows!(cm, i, j)
            swap_cols!(cm, i, j)
          end
        end
        weyl_group(root_system(cm))
      end,
      unique(rand(1:15, 6)), false,
    ),
  ]
    for k in 1:4
      if k == 1
        # On the first run, test the standard parabolics and projections
        if check_proj
          para, emb, proj = parabolic_subgroup_with_projection(W, vec)
        else
          para, emb = parabolic_subgroup(W, vec)
        end
        genimgs = [gen(W, i) for i in vec] # Desired images of gens(para) in W
      else
        # On subsequent runs, conjugate by random elements
        r = rand(W)
        para, emb = parabolic_subgroup(W, vec, r)
        genimgs = [conj(W[i], r) for i in vec]
      end
      # Test that emb maps gens(para) to genimgs
      for i in 1:length(vec)
        @test emb(gen(para, i)) == genimgs[i]
      end
      # Test that emb is a homomorphism
      for _ in 1:5
        p1 = rand(para)
        p2 = rand(para)
        @test emb(p1) * emb(p2) == emb(p1 * p2)
      end
      # Test proj
      if k == 1 && check_proj
        # Test that proj maps gens(para) to gens(W)[vec]
        for i in 1:length(vec)
          @test proj(gen(W, vec[i])) == gen(para, i)
        end
        # Test that proj is a homomorphism
        for _ in 1:5
          w1 = rand(W)
          w2 = rand(W)
          @test proj(w1) * proj(w2) == proj(w1 * w2)
        end
        # Test that proj is the left-inverse of emb
        for _ in 1:5
          p = rand(para)
          @test proj(emb(p)) == p
        end
      end
    end
  end

  @testset "irreducible_factors and inner_direct_product for $(name)" for (
    name, factor_list
  ) in [
    ("Empty product", WeylGroup[]),
    ("Trivial Weyl group", [weyl_group(zero_matrix(ZZ, 0, 0))]),
    ("A1", [weyl_group(:A, 1)]),
    (
      "A3 with non-canonical ordering of simple roots",
      [weyl_group(root_system([2 -1 -1; -1 2 0; -1 0 2]))],
    ),
    (
      "B4 with non-canonical ordering of simple roots",
      [weyl_group(root_system([2 -1 -1 0; -1 2 0 -1; -2 0 2 0; 0 -1 0 2]))],
    ),
    ("B4", [weyl_group(root_system(:B, 4))]),
    ("F4", [weyl_group(:F, 4)]),
    ("G2+Trivial", [weyl_group(:G, 2), weyl_group(zero_matrix(ZZ, 0, 0))]),
    ("A2+E8", [weyl_group(:A, 2), weyl_group(:E, 8)]),
    ("C3+E6+G2", [weyl_group(:C, 3), weyl_group(:E, 6), weyl_group(:G, 2)]),
    ( # This case tests the inner direct product of reducible Weyl groups
      "complicated case",
      [
        begin
          cm = cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2))
          for _ in 1:50
            i, j = rand(1:nrows(cm), 2)
            if i != j
              swap_rows!(cm, i, j)
              swap_cols!(cm, i, j)
            end
          end
          weyl_group(cm)
        end,
        begin
          cm = cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2))
          for _ in 1:50
            i, j = rand(1:nrows(cm), 2)
            if i != j
              swap_rows!(cm, i, j)
              swap_cols!(cm, i, j)
            end
          end
          weyl_group(cm)
        end,
      ],
    ),
  ]
    for morphisms in (true, false)
      if morphisms
        W, emb_inner, proj_inner = inner_direct_product(factor_list; morphisms=morphisms)
        irred_factors, emb_irred, proj_irred = irreducible_factors(W; morphisms=morphisms)
      else
        W = inner_direct_product(factor_list; morphisms=morphisms)
        irred_factors = irreducible_factors(W; morphisms=morphisms)
      end

      types, ordering = root_system_type_with_ordering(root_system(W))

      # W has the same type as the inner direct product of its irreducible factors.
      # This is the only necessary test for inner_direct_product if morphisms=false.
      # Both type vectors must be in the same order because
      # cartan_type(cartan_matrix(type_vector)) == type_vector.
      @test root_system_type(root_system(inner_direct_product(irred_factors))) == types

      @test W isa WeylGroup

      @testset "irreducible_factors" begin
        # Test number of irreducible factors
        @test length(irred_factors) == length(types)
        if morphisms
          @test length(emb_irred) == length(proj_irred) == length(irred_factors)
        end
        # Test that each irreducible factor is correct
        start_index = 1
        for i in eachindex(types)
          irred_factor = irred_factors[i]
          type = types[i] # Type of irreducible factor
          rk = type[2] # Rank of irreducible factor
          # Test for correct type
          @test type == only(root_system_type(root_system(irred_factor)))
          # Test mapping of generators
          if morphisms
            irred_factor_gens = gens(irred_factor)
            image_gens = [W[ordering[j]] for j in start_index:(start_index + rk - 1)]
            @test emb_irred[i].(irred_factor_gens) == image_gens
            @test proj_irred[i].(image_gens) == irred_factor_gens
            start_index = start_index + rk
          end
          # We do not test that the maps are homomorphisms
          # because they come from parabolic_subgroup_with_projection,
          # for which this is tested separately
        end
      end

      if morphisms
        @testset "inner_direct_product" begin
          start_index = 1
          for (i, factor) in enumerate(factor_list)
            rs = root_system(factor)
            type = root_system_type(rs)
            rk = rank(rs)

            # Test mapping of generators
            factor_gens = gens(factor)
            image_gens = [W[j] for j in start_index:(start_index + rk - 1)]
            @test emb_inner[i].(factor_gens) == image_gens
            @test proj_inner[i].(image_gens) == factor_gens
            start_index = start_index + rk
            if morphisms
              # Test that emb_inner[i] is a homomorphism
              for _ in 1:5
                w1 = rand(factor)
                w2 = rand(factor)
                @test emb_inner[i](w1 * w2) == emb_inner[i](w1) * emb_inner[i](w2)
              end
              # Test that proj_inner[i] is a homomorphism
              for _ in 1:5
                w1 = rand(W)
                w2 = rand(W)
                @test proj_inner[i](w1 * w2) == proj_inner[i](w1) * proj_inner[i](w2)
              end
              # Test that proj_inner[i] is a left-inverse of emb_inner[i]
              for _ in 1:5
                w = rand(factor)
                @test proj_inner[i](emb_inner[i](w)) == w
              end
            end
          end
        end
      end
    end
  end

  @testset "G-sets of WeylGroups" begin
    # construction using RootSpaceElem elements from root system of W
    W = weyl_group(:A, 2)
    pts = roots(root_system(W))
    Omega = @inferred gset(W, pts)
    @test AbstractAlgebra.PrettyPrinting.repr_terse(Omega) == "G-set"
    @test isa(Omega, GSet)
    @test (@inferred length(Omega)) == 6
    @test (@inferred length(@inferred orbits(Omega))) == 1
    @test is_transitive(Omega)
    @test !is_primitive(Omega)
    @test is_regular(Omega)
    @test is_semiregular(Omega)

    W = weyl_group([(:A, 2), (:B, 3)])
    pts = roots(root_system(W))
    Omega = gset(W, pts)
    @test length(Omega) == 24
    @test length(orbits(Omega)) == 3
    @test !is_transitive(Omega)
    @test !is_primitive(Omega)
    @test !is_regular(Omega)
    @test !is_semiregular(Omega)

    # construction using WeightLatticeElem elements from root system
    R = root_system(:A, 2)
    w = WeightLatticeElem(R, [2, 2])
    W = weyl_group(R)
    Omega = @inferred gset(W, [w])
    @test isa(Omega, GSet)
    @test length(Omega) == 6
    @test length(orbits(Omega)) == 1
    @test is_transitive(Omega)
    @test !is_primitive(Omega)
    @test is_regular(Omega)
    @test is_semiregular(Omega)

    # orbit
    W = weyl_group([(:A, 2), (:B, 3)])
    pts = roots(root_system(W))
    Omega = gset(W, pts)
    orbs = orbits(Omega)
    @test length(orbs) == 3
    @test length(orbit(Omega, pts[2])) == 6
    @test sort(map(length, orbs)) == [6, 6, 12]

    # action homomorphism
    W = weyl_group(:A, 2)
    pts = roots(root_system(W))
    Omega = gset(W, pts)
    acthom = action_homomorphism(Omega)
    @test order(image(acthom)[1]) == order(W)

    # all_blocks
    bl = all_blocks(Omega)
    @test length(bl) == 4
    @test Set([pts[1], pts[5]]) in bl

    # blocks
    bl = blocks(Omega)
    @test length(bl) == 3
    @test elements(bl) == map(Set, [[pts[1], pts[4]], [pts[6], pts[3]], [pts[5], pts[2]]])
    @test length(orbits(bl)) == 1
  end
end
