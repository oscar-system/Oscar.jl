@testset "LieTheory.CartanMatrix" begin
  @testset "cartan_matrix" begin
    for fam in (:A, :B, :C, :D, :E, :F, :G, :X), i in -1:10
      if is_cartan_type(fam, i)
        cartan_matrix(fam, i)
      else
        @test_throws ArgumentError cartan_matrix(fam, i)
      end
    end

    @test cartan_matrix(:A, 1) == matrix(ZZ, 1, 1, [2])
    @test cartan_matrix(:A, 2) == ZZ[2 -1; -1 2]
    @test cartan_matrix(:B, 2) == ZZ[2 -1; -2 2]
    @test cartan_matrix(:B, 3) == ZZ[2 -1 0; -1 2 -1; 0 -2 2]
    @test cartan_matrix(:C, 2) == ZZ[2 -2; -1 2]
    @test cartan_matrix(:C, 3) == ZZ[2 -1 0; -1 2 -2; 0 -1 2]
    @test cartan_matrix(:D, 4) == ZZ[2 -1 0 0; -1 2 -1 -1; 0 -1 2 0; 0 -1 0 2]
  end

  @testset "cartan_matrix(type::Tuple{Symbol,Int}...)" begin
    @test cartan_matrix((:A, 2)) == cartan_matrix(:A, 2)
    @test cartan_matrix((:A, 2), (:B, 2)) ==
      block_diagonal_matrix([cartan_matrix(:A, 2), cartan_matrix(:B, 2)])

    @test_throws MethodError cartan_matrix()
  end

  @testset "is_cartan_matrix" begin
    # test non cartan
    @test is_cartan_matrix(ZZ[2 -1; -1 1]) == false
    @testset for i in 1:3
      is_cartan_matrix(ZZ[2 -i; 0 2]) == false
    end

    # test finite type
    function test_is_cartan_matrix_finite_type(cm::ZZMatrix)
      rk = nrows(cm)
      for _ in 1:50
        i, j = rand(1:rk, 2)
        if i != j
          swap_rows!(cm, i, j)
          swap_cols!(cm, i, j)
        end
      end

      @test is_cartan_matrix(cm) == true
      @test is_cartan_matrix(cm; generalized=false) == true
    end

    test_is_cartan_matrix_finite_type(cartan_matrix(:A, 1))
    test_is_cartan_matrix_finite_type(cartan_matrix(:A, 3))
    test_is_cartan_matrix_finite_type(cartan_matrix(:B, 4))
    test_is_cartan_matrix_finite_type(cartan_matrix(:C, 2))
    test_is_cartan_matrix_finite_type(cartan_matrix(:D, 5))
    test_is_cartan_matrix_finite_type(cartan_matrix(:E, 6))
    test_is_cartan_matrix_finite_type(cartan_matrix(:E, 7))
    test_is_cartan_matrix_finite_type(cartan_matrix(:E, 8))
    test_is_cartan_matrix_finite_type(cartan_matrix(:F, 4))
    test_is_cartan_matrix_finite_type(cartan_matrix(:G, 2))

    test_is_cartan_matrix_finite_type(cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2)))
    test_is_cartan_matrix_finite_type(cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2)))

    # test affine type
    @test is_cartan_matrix(ZZ[2 -2; -2 2]) == true
    @test is_cartan_matrix(ZZ[2 -2; -2 2]; generalized=false) == false

    @test is_cartan_matrix(ZZ[2 -4; -1 2]) == true
    @test is_cartan_matrix(ZZ[2 -4; -1 2]; generalized=false) == false

    @test is_cartan_matrix(ZZ[2 -1 -1; -1 2 -1; -1 -1 2]) == true
    @test is_cartan_matrix(ZZ[2 -1 -1; -1 2 -1; -1 -1 2]; generalized=false) == false
  end

  @testset "cartan_symmetrizer" begin
    # cartan_symmetrizer function is indirectly covered by the cartan_bilinear_form tests
    # here we simply guarantee that the choice of the symmetrizer is as expected

    # finite type
    @test cartan_symmetrizer(cartan_matrix(:A, 1)) == [1]
    @test cartan_symmetrizer(cartan_matrix(:A, 2)) == [1, 1]
    @test cartan_symmetrizer(cartan_matrix(:B, 2)) == [2, 1]
    @test cartan_symmetrizer(cartan_matrix(:B, 4)) == [2, 2, 2, 1]
    @test cartan_symmetrizer(cartan_matrix(:C, 2)) == [1, 2]
    @test cartan_symmetrizer(cartan_matrix(:C, 4)) == [1, 1, 1, 2]
    @test cartan_symmetrizer(cartan_matrix(:D, 4)) == [1, 1, 1, 1]
    @test cartan_symmetrizer(cartan_matrix(:E, 6)) == [1, 1, 1, 1, 1, 1]
    @test cartan_symmetrizer(cartan_matrix(:E, 7)) == [1, 1, 1, 1, 1, 1, 1]
    @test cartan_symmetrizer(cartan_matrix(:E, 8)) == [1, 1, 1, 1, 1, 1, 1, 1]
    @test cartan_symmetrizer(cartan_matrix(:F, 4)) == [2, 2, 1, 1]
    @test cartan_symmetrizer(cartan_matrix(:G, 2)) == [1, 3]

    # affine type
    #@test cartan_symmetrizer(ZZ[2 -2; -2 2]) == [1, 1] # A1~1
    #@test cartan_symmetrizer(ZZ[2 -2 0; -1 2 -1; 0 -2 2]) == [1, 2, 1] # D3~2
    #@test cartan_symmetrizer(ZZ[2 -4; -1 2]) == [1, 4] # A1~2

    # hyperbolic type
    #@test cartan_symmetrizer(ZZ[2 -2; -3 2]) == [3, 2]
  end

  @testset "cartan_bilinear_form" begin
    function test_cartan_bilinear_form(cm::ZZMatrix)
      rk = nrows(cm)
      for _ in 1:10
        for _ in 1:10
          i, j = rand(1:rk, 2)
          if i != j
            swap_rows!(cm, i, j)
            swap_cols!(cm, i, j)
          end
        end

        bil = cartan_bilinear_form(cm)
        @test is_symmetric(bil)
        @test all(cm[i, j] == div(2 * bil[i, j], bil[i, i]) for i in 1:rk, j in 1:rk)
      end
    end

    # irreducible matrices
    test_cartan_bilinear_form(cartan_matrix(:A, 1))
    test_cartan_bilinear_form(cartan_matrix(:A, 3))
    test_cartan_bilinear_form(cartan_matrix(:B, 4))
    test_cartan_bilinear_form(cartan_matrix(:C, 2))
    test_cartan_bilinear_form(cartan_matrix(:D, 5))
    test_cartan_bilinear_form(cartan_matrix(:E, 6))
    test_cartan_bilinear_form(cartan_matrix(:E, 7))
    test_cartan_bilinear_form(cartan_matrix(:E, 8))
    test_cartan_bilinear_form(cartan_matrix(:F, 4))
    test_cartan_bilinear_form(cartan_matrix(:G, 2))

    # reducible matrices
    test_cartan_bilinear_form(cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2)))
    test_cartan_bilinear_form(cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2)))

    # affine type
    #@test cartan_bilinear_form(ZZ[2 -2; -2 2]) == ZZ[2 -2; -2 2]
    #@test cartan_bilinear_form(ZZ[2 -2 0; -1 2 -1; 0 -2 2]) == ZZ[2 -2 0; -2 4 -2; 0 -2 2]
    #@test cartan_bilinear_form(ZZ[2 -4; -1 2]) == ZZ[2 -4; -4 8]
  end

  @testset "cartan_type_with_ordering" begin
    function test_cartan_type_with_ordering(
      fam::Symbol, n::Int; autos::Vector{PermGroupElem}=[one(symmetric_group(n))]
    )
      @req all(aut -> parent(aut) == symmetric_group(n), autos) "Incompatible permutation parent"
      cm = cartan_matrix(fam, n)
      for perm in symmetric_group(n)
        type, ord = cartan_type_with_ordering(
          permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm); check=false
        )
        @test type == [(fam, n)]
        @test ord in [Vector{Int}(aut * perm) for aut in autos]
        @test !is_one(perm) || ord == 1:n
      end
    end

    @testset "A_n" begin
      @testset "A_$n" for n in 1:8
        # automorphism group of Dynkin diagram ≅ S_2 (horizontal reflection)
        test_cartan_type_with_ordering(
          :A,
          n;
          autos=[
            cperm(symmetric_group(n), Int[]),
            cperm(symmetric_group(n), Vector{Int}[[i, n + 1 - i] for i in 1:div(n, 2)]),
          ],
        )
      end
    end

    @testset "B_n" begin
      @testset "B_2" begin
        n = 2
        # B_2 is isomorphic to C_2; we decide based on the ordering
        cm = cartan_matrix(:B, n)
        # identity permutation
        type, ord = cartan_type_with_ordering(cm; check=false)
        @test type == [(:B, n)]
        @test ord == [1, 2]
        # non-identity permutation
        type, ord = cartan_type_with_ordering(transpose(cm); check=false)
        @test type == [(:C, 2)]
        @test ord == [1, 2]
      end
      @testset "B_$n" for n in 3:8
        # no automorphisms of Dynkin diagram
        test_cartan_type_with_ordering(:B, n)
      end
    end

    @testset "C_n" begin
      @testset "C_2" begin
        n = 2
        # C_2 is isomorphic to B_2; we decide based on the ordering
        cm = cartan_matrix(:C, n)
        # identity permutation
        type, ord = cartan_type_with_ordering(cm; check=false)
        @test type == [(:C, n)]
        @test ord == [1, 2]
        # non-identity permutation
        type, ord = cartan_type_with_ordering(transpose(cm); check=false)
        @test type == [(:B, 2)]
        @test ord == [1, 2]
      end
      @testset "C_$n" for n in 3:8
        # no automorphisms of Dynkin diagram
        test_cartan_type_with_ordering(:C, n)
      end
    end

    @testset "D_n" begin
      @testset "D_$n" for n in 4:8
        if n == 4
          # automorphism group of Dynkin diagram ≅ S_3
          test_cartan_type_with_ordering(
            :D, n; autos=(@perm 4 [(), (1, 3), (1, 4), (3, 4), (1, 3, 4), (1, 4, 3)])
          )
        else
          # automorphism group of Dynkin diagram ≅ S_2 (vertical reflection)
          test_cartan_type_with_ordering(
            :D,
            n;
            autos=[
              cperm(symmetric_group(n), Int[]), cperm(symmetric_group(n), Int[n - 1, n])
            ],
          )
        end
      end
    end

    @testset "E_n" begin
      @testset "E_$n" for n in 6:8
        if n == 6
          # automorphism group of Dynkin diagram ≅ S_2 (horizontal reflection)
          test_cartan_type_with_ordering(:E, n; autos=(@perm 6 [(), (1, 6)(3, 5)]))
        else
          # no automorphisms of Dynkin diagram
          test_cartan_type_with_ordering(:E, n)
        end
      end
    end

    @testset "F_4" begin
      # no automorphisms of Dynkin diagram
      test_cartan_type_with_ordering(:F, 4)
    end

    @testset "G_2" begin
      # no automorphisms of Dynkin diagram
      test_cartan_type_with_ordering(:G, 2)
    end

    @testset "non-simple cases" begin
      type, ord = cartan_type_with_ordering(cartan_matrix((:A, 1), (:A, 2)); check=false)
      @test type == [(:A, 1), (:A, 2)]
      @test ord == 1:3

      type, ord = cartan_type_with_ordering(cartan_matrix((:A, 1), (:B, 2)); check=false)
      @test type == [(:A, 1), (:B, 2)]
      @test ord == 1:3

      type, ord = cartan_type_with_ordering(cartan_matrix((:C, 2), (:B, 2)); check=false)
      @test type == [(:C, 2), (:B, 2)]
      @test ord == 1:4

      type, ord = cartan_type_with_ordering(
        cartan_matrix((:E, 8), (:A, 5), (:D, 4), (:F, 4), (:B, 8)); check=false
      )
      @test type == [(:E, 8), (:A, 5), (:D, 4), (:F, 4), (:B, 8)]
      @test ord == 1:(8 + 5 + 4 + 4 + 8)

      type, ord = cartan_type_with_ordering(
        ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2]; check=false
      )
      @test type == [(:B, 2), (:C, 2)]
      @test ord == [1, 3, 2, 4]
    end
  end

  @testset "is_cartan_type" begin
    for i in -1:10
      @test is_cartan_type(:A, i) == (i >= 1)
      @test is_cartan_type(:B, i) == (i >= 2)
      @test is_cartan_type(:C, i) == (i >= 2)
      @test is_cartan_type(:D, i) == (i >= 4)
      @test is_cartan_type(:E, i) == (i in [6, 7, 8])
      @test is_cartan_type(:F, i) == (i == 4)
      @test is_cartan_type(:G, i) == (i == 2)
      @test is_cartan_type(:X, i) == false
    end
  end
end
