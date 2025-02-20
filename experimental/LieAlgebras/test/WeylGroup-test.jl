@testset "LieAlgebras.WeylGroup" begin
  function is_in_normal_form(x::WeylGroupElem)
    return word(parent(x)(word(x))) == word(x)
  end

  _isomorphic_group_on_gens = Oscar.LieAlgebras._isomorphic_group_on_gens

  # TODO: merge with conformance tests in test/LieTheory/WeylGroup.jl once this is moved to src
  @testset "WeylGroup Group isomorphism test for $(Wname)" for (Wname, W) in [
    # Some test cases can be removed as soon as isomorphism has support for reducible types.
    # (Then nearly everything is covered by complicated case 1 and 2.)
    # Until then, all irreducible types and some explicit non-canonical orderings should be tested.
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
      if length(type) == 1
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
end
