@testset "LieAlgebras.RootSystem" begin
  @testset "conjugate_dominant_weight_with_elem(w::WeightLatticeElem)" begin
    for (R, vec) in [
      (root_system(:A, 5), [1, -1, 2, 0, 2]),
      (root_system(:B, 3), [1, 1, 1]),
      (root_system(:C, 4), [2, 1, 0, 1]),
      (root_system(:D, 5), [-1, 2, 2, -1, -1]),
      (root_system(:E, 6), [1, 2, 0, 0, 2, 1]),
      (root_system(:F, 4), [1, 2, 3, 4]),
      (root_system(:G, 2), [-1, -1]),
    ]
      wt = WeightLatticeElem(R, vec)
      d, x = conjugate_dominant_weight_with_elem(wt)
      @test is_dominant(d)
      @test x * wt == d
    end
  end

  @testset "root_system(cartan_matrix::ZZMatrix)" begin
    R = root_system(:F, 4)
    @test n_positive_roots(R) == 24
    @test n_roots(R) == 48

    R = root_system(:G, 2)
    @test n_positive_roots(R) == 6
    @test n_roots(R) == 12

    @test coefficients(positive_root(R, 1)) == QQ[1 0]
    @test coefficients(positive_root(R, 2)) == QQ[0 1]
  end

  @testset "property tests" begin
    function root_system_property_tests(R::RootSystem, rk::Int, npositive_roots::Int)
      W = weyl_group(R)

      @testset "Property tests" begin
        @test rank(R) == rk
        @test n_simple_roots(R) == rk
        @test n_positive_roots(R) == npositive_roots
        @test n_roots(R) == 2 * npositive_roots

        @test length(simple_roots(R)) == n_simple_roots(R)
        @test length(positive_roots(R)) == n_positive_roots(R)
        @test length(negative_roots(R)) == n_positive_roots(R)
        @test length(roots(R)) == n_roots(R)
        @test all(i -> simple_root(R, i) == simple_roots(R)[i], 1:rk)
        @test all(i -> positive_root(R, i) == positive_roots(R)[i], 1:npositive_roots)
        @test all(i -> negative_root(R, i) == negative_roots(R)[i], 1:npositive_roots)
        @test simple_roots(R) == positive_roots(R)[1:rk]
        @test all(is_root, roots(R))
        n_roots(R) >= 1 && @test !is_root(root(R, 1) - root(R, 1))
        @test all(r -> !is_root(2 * r), roots(R))
        @test all(is_root_with_index(r) == (true, i) for (i, r) in enumerate(roots(R)))
        @test all(r -> is_positive_root(r) == is_positive_root_with_index(r)[1], roots(R))
        @test all(r -> is_negative_root(r) == is_negative_root_with_index(r)[1], roots(R))
        @test all(r -> is_simple_root(r) == is_simple_root_with_index(r)[1], roots(R))
        @test all(
          is_positive_root_with_index(r) == (true, i) for
          (i, r) in enumerate(positive_roots(R))
        )
        @test all(!is_negative_root, positive_roots(R))
        @test all(
          is_negative_root_with_index(r) == (true, i) for
          (i, r) in enumerate(negative_roots(R))
        )
        @test all(!is_positive_root, negative_roots(R))
        @test all(
          is_simple_root_with_index(r) == (true, i) for (i, r) in enumerate(simple_roots(R))
        )
        @test all(is_positive_root, simple_roots(R))
        @test all(!is_negative_root, simple_roots(R))
        @test all(iszero, positive_roots(R) + negative_roots(R))
        n_roots(R) >= 1 && for _ in 1:10
          r = root(R, rand(1:n_roots(R)))
          w = rand(W)
          @test is_root(w * r)
        end

        @test length(simple_coroots(R)) == n_simple_roots(R)
        @test length(positive_coroots(R)) == n_positive_roots(R)
        @test length(negative_coroots(R)) == n_positive_roots(R)
        @test length(coroots(R)) == n_roots(R)
        @test all(i -> simple_coroot(R, i) == simple_coroots(R)[i], 1:rk)
        @test all(i -> positive_coroot(R, i) == positive_coroots(R)[i], 1:npositive_roots)
        @test all(i -> negative_coroot(R, i) == negative_coroots(R)[i], 1:npositive_roots)
        @test simple_coroots(R) == positive_coroots(R)[1:rk]
        @test all(is_coroot, coroots(R))
        n_roots(R) >= 1 && @test !is_coroot(coroot(R, 1) - coroot(R, 1))
        @test all(r -> !is_coroot(2 * r), coroots(R))
        @test all(is_coroot_with_index(r) == (true, i) for (i, r) in enumerate(coroots(R)))
        @test all(
          r -> is_positive_coroot(r) == is_positive_coroot_with_index(r)[1], coroots(R)
        )
        @test all(
          r -> is_negative_coroot(r) == is_negative_coroot_with_index(r)[1], coroots(R)
        )
        @test all(r -> is_simple_coroot(r) == is_simple_coroot_with_index(r)[1], coroots(R))
        @test all(
          is_positive_coroot_with_index(r) == (true, i) for
          (i, r) in enumerate(positive_coroots(R))
        )
        @test all(!is_negative_coroot, positive_coroots(R))
        @test all(
          is_negative_coroot_with_index(r) == (true, i) for
          (i, r) in enumerate(negative_coroots(R))
        )
        @test all(!is_positive_coroot, negative_coroots(R))
        @test all(
          is_simple_coroot_with_index(r) == (true, i) for
          (i, r) in enumerate(simple_coroots(R))
        )
        @test all(is_positive_coroot, simple_coroots(R))
        @test all(!is_negative_coroot, simple_coroots(R))
        @test all(iszero, positive_coroots(R) + negative_coroots(R))

        @test issorted(height.(positive_roots(R))) # sorted by height

        @test all(
          i ->
            dot(coefficients(coroot(R, i)) * cartan_matrix(R), coefficients(root(R, i))) ==
            2,
          1:n_roots(R),
        )
      end

      @testset "Serialization" begin
        mktempdir() do path
          test_save_load_roundtrip(path, R) do loaded
            # nothing, cause `R === loaded` anyway
          end

          if n_roots(R) >= 1
            r = positive_root(R, n_positive_roots(R))
            test_save_load_roundtrip(path, r) do loaded
              @test root_system(loaded) === R
              @test coefficients(loaded) == coefficients(r)
            end
          end

          test_save_load_roundtrip(path, positive_roots(R)) do loaded
            @test length(loaded) == n_positive_roots(R)
            @test all(
              coefficients(loaded[i]) == coefficients(root(R, i)) for
              i in 1:n_positive_roots(R)
            )
          end
        end
      end
    end

    @testset "rk 0" begin
      R = root_system(zero_matrix(ZZ, 0, 0))
      root_system_property_tests(R, 0, 0)
    end

    @testset "A_$n" for n in [1, 2, 6]
      R = root_system(:A, n)
      root_system_property_tests(R, n, binomial(n + 1, 2))
    end

    @testset "B_$n" for n in [2, 3, 6]
      R = root_system(:B, n)
      root_system_property_tests(R, n, n^2)
    end

    @testset "C_$n" for n in [2, 3, 6]
      R = root_system(:C, n)
      root_system_property_tests(R, n, n^2)
    end

    @testset "D_$n" for n in [4, 6]
      R = root_system(:D, n)
      root_system_property_tests(R, n, n^2 - n)
    end

    @testset "E_6" begin
      R = root_system(:E, 6)
      root_system_property_tests(R, 6, 36)
    end

    @testset "E_7" begin
      R = root_system(:E, 7)
      root_system_property_tests(R, 7, 63)
    end

    @testset "E_8" begin
      R = root_system(:E, 8)
      root_system_property_tests(R, 8, 120)
    end

    @testset "F_4" begin
      R = root_system(:F, 4)
      root_system_property_tests(R, 4, 24)
    end

    @testset "G_2" begin
      R = root_system(:G, 2)
      root_system_property_tests(R, 2, 6)
    end

    @testset "something mixed 1" begin
      cm = cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2))
      for _ in 1:50
        i, j = rand(1:nrows(cm), 2)
        if i != j
          swap_rows!(cm, i, j)
          swap_cols!(cm, i, j)
        end
      end
      R = root_system(cm)
      root_system_property_tests(R, 3 + 3 + 6 + 2, binomial(3 + 1, 2) + 3^2 + 36 + 6)
    end

    @testset "something mixed 2" begin
      cm = cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2))
      for _ in 1:50
        i, j = rand(1:nrows(cm), 2)
        if i != j
          swap_rows!(cm, i, j)
          swap_cols!(cm, i, j)
        end
      end
      R = root_system(cm)
      root_system_property_tests(R, 4 + 2 + 7 + 2, 24 + 2^2 + 63 + 6)
    end
  end

  @testset "WeightLatticeElem" begin
    R = root_system(:A, 2)
    w = WeightLatticeElem(R, [2, 2])

    @test root_system(w) === R
  end
end
