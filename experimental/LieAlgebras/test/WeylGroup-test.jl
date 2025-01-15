@testset "LieAlgebras.WeylGroup" begin
  function is_in_normal_form(x::WeylGroupElem)
    return word(parent(x)(word(x))) == word(x)
  end

  # TODO: merge with conformance tests in test/LieTheory/WeylGroup.jl once this is moved to src
  @testset "WeylGroup Group isomorphism test for $(Wname)" for (Wname, W) in [
    ("A1", weyl_group(:A, 1)),
    ("A5", weyl_group(:A, 5)),
    ("B4", weyl_group(root_system(:B, 4))),
    ("D5", weyl_group(cartan_matrix(:D, 5))),
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
    @testset "isomorphism(FPGroup, ::WeylGroup; set_properties=$set_properties)" for set_properties in
                                                                                     [
      false, true
    ]
      G = fp_group(W; set_properties)
      if (is_finite(W) && ngens(W) < 6) || set_properties #= for sane runtime =#
        @test is_finite(G) == is_finite(W)
        is_finite(W) && @test order(G) == order(W)
      end

      iso = isomorphism(FPGroup, W; set_properties)
      @test W == domain(iso)
      G = codomain(iso)
      if (is_finite(W) && ngens(W) < 6) || set_properties #= for sane runtime =#
        @test is_finite(G) == is_finite(W)
        is_finite(W) && @test order(G) == order(W)
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

    if has_root_system_type(root_system(W))
      type, ordering = root_system_type_with_ordering(root_system(W))
      if length(type) == 1 && issorted(ordering) && only(type)[1] == :A # only implemented for A_n (yet)
        @testset "isomorphism(PermGroup, ::WeylGroup; set_properties=$set_properties)" for set_properties in
                                                                                           [
          false, true
        ]
          G = permutation_group(W; set_properties)
          if (is_finite(W) && ngens(W) < 6) || set_properties #= for sane runtime =#
            @test is_finite(G) == is_finite(W)
            is_finite(W) && @test order(G) == order(W)
          end

          iso = isomorphism(PermGroup, W; set_properties)
          @test W == domain(iso)
          G = codomain(iso)
          if (is_finite(W) && ngens(W) < 6) || set_properties #= for sane runtime =#
            @test is_finite(G) == is_finite(W)
            is_finite(W) && @test order(G) == order(W)
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
  end
end
