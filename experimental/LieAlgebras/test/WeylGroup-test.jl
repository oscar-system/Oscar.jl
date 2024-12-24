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

  @testset "Serialization" begin
    mktempdir() do path
      @testset "simple saving and loading" begin
        W = weyl_group((:A, 2), (:B, 4))

        test_save_load_roundtrip(path, W) do loaded
          # nothing, cause `W === loaded` anyway
        end

        x = rand(W)
        test_save_load_roundtrip(path, x) do loaded
          @test parent(loaded) === W
          @test word(loaded) == word(x)
        end

        test_save_load_roundtrip(path, gens(W)) do loaded
          @test length(loaded) == ngens(W)
          @test all(
            word(loaded[i]) == word(gen(W, i)) for i in 1:ngens(W)
          )
        end
      end

      @testset "cyclic reference between R and W survives" begin
        Oscar.reset_global_serializer_state()

        R_filename = joinpath(path, "R.mrdi")
        W_filename = joinpath(path, "W.mrdi")

        R = root_system(:D, 5)
        W = weyl_group(R)

        save(R_filename, R)
        save(W_filename, W)

        Oscar.reset_global_serializer_state()

        loaded_R = load(R_filename)
        loaded_W = load(W_filename)

        @test loaded_R === root_system(loaded_W)
        @test loaded_W === weyl_group(loaded_R)

        loaded_R = loaded_W = nothing # unset all references

        Oscar.reset_global_serializer_state()

        loaded_W = load(W_filename)
        loaded_R = load(R_filename)

        @test loaded_R === root_system(loaded_W)
        @test loaded_W === weyl_group(loaded_R)

        loaded_R = loaded_W = nothing # unset all references        
      end
    end
  end
end
