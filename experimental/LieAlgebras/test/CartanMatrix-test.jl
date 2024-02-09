@testset "LieAlgebras.CartanMatrix" begin
  @testset "cartan_matrix(fam::Symbol, rk::Int)" begin
    @test cartan_matrix(:A, 1) == matrix(ZZ, 1, 1, [2])
    @test cartan_matrix(:A, 2) == ZZ[2 -1; -1 2]

    @test cartan_matrix(:B, 2) == ZZ[2 -1; -2 2]
    @test cartan_matrix(:B, 3) == ZZ[2 -1 0; -1 2 -1; 0 -2 2]

    @test cartan_matrix(:C, 2) == ZZ[2 -2; -1 2]
    @test cartan_matrix(:C, 3) == ZZ[2 -1 0; -1 2 -2; 0 -1 2]

    @test cartan_matrix(:D, 4) == ZZ[2 -1 0 0; -1 2 -1 -1; 0 -1 2 0; 0 -1 0 2]

    @test_throws ArgumentError root_system(:X, 1)
    @test_throws ArgumentError root_system(:A, 0)
    @test_throws ArgumentError root_system(:B, 1)
    @test_throws ArgumentError root_system(:C, 1)
    @test_throws ArgumentError root_system(:D, 2)
    @test_throws ArgumentError root_system(:E, 5)
    @test_throws ArgumentError root_system(:E, 9)
    @test_throws ArgumentError root_system(:F, 3)
    @test_throws ArgumentError root_system(:F, 5)
    @test_throws ArgumentError root_system(:G, 1)
    @test_throws ArgumentError root_system(:G, 3)
  end

  @testset "cartan_matrix(type::Tuple{Symbol,Int}...)" begin
    @test cartan_matrix((:A, 2)) == cartan_matrix(:A, 2)
    @test cartan_matrix((:A, 2), (:B, 2)) ==
      block_diagonal_matrix([cartan_matrix(:A, 2), cartan_matrix(:B, 2)])

    @test_throws ArgumentError cartan_matrix()
  end

  @testset "is_cartan_matrix(mat::ZZMatrix; generalized::Bool=true)" begin
    # test finite type
    @test is_cartan_matrix(cartan_matrix(:A, 1)) == true
    @test is_cartan_matrix(cartan_matrix(:A, 1); generalized=false) == true

    @test is_cartan_matrix(cartan_matrix(:A, 2)) == true
    @test is_cartan_matrix(cartan_matrix(:A, 2); generalized=false) == true

    @test is_cartan_matrix(cartan_matrix(:B, 2)) == true
    @test is_cartan_matrix(cartan_matrix(:B, 2); generalized=false) == true

    @test is_cartan_matrix(cartan_matrix(:C, 2)) == true
    @test is_cartan_matrix(cartan_matrix(:C, 2); generalized=false) == true

    @test is_cartan_matrix(cartan_matrix(:D, 4)) == true
    @test is_cartan_matrix(cartan_matrix(:D, 4); generalized=false) == true

    @test is_cartan_matrix(cartan_matrix(:G, 2)) == true
    @test is_cartan_matrix(cartan_matrix(:G, 2); generalized=false) == true

    # test affine type
    @test is_cartan_matrix(ZZ[2 -2; -2 2]) == true
    @test is_cartan_matrix(ZZ[2 -2; -2 2]; generalized=false) == false

    @test is_cartan_matrix(ZZ[2 -4; -1 2]) == true
    @test is_cartan_matrix(ZZ[2 -4; -1 2]; generalized=false) == false

    @test is_cartan_matrix(ZZ[2 -1 -1; -1 2 -1; -1 -1 2]) == true
    @test is_cartan_matrix(ZZ[2 -1 -1; -1 2 -1; -1 -1 2]; generalized=false) == false
  end

  @testset "cartan_symmetrizer(gcm::ZZMatrix; check::Bool)" begin
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

    @test cartan_symmetrizer(cartan_matrix((:A, 3), (:B, 3))) == [[1, 1, 1]; [2, 2, 1]]
    @test cartan_symmetrizer(cartan_matrix((:F, 4), (:D, 4))) ==
      [[2, 2, 1, 1]; [1, 1, 1, 1]]
    @test cartan_symmetrizer(cartan_matrix((:C, 4), (:G, 2))) == [[1, 1, 1, 2]; [1, 3]]

    @test cartan_symmetrizer(ZZ[2 -2 0; -1 2 -1; 0 -1 2]) == [1, 2, 2]
    @test cartan_symmetrizer(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2]) == [2, 1, 1, 2]

    # affine type
    @test cartan_symmetrizer(ZZ[2 -2; -2 2]) == [1, 1] # A1~1
    @test cartan_symmetrizer(ZZ[2 -2 0; -1 2 -1; 0 -2 2]) == [1, 2, 1] # D3~2
    @test cartan_symmetrizer(ZZ[2 -4; -1 2]) == [1, 4] # A1~2

    # hyperbolic type
    @test cartan_symmetrizer(ZZ[2 -2; -3 2]) == [3, 2]
  end

  @testset "cartan_bilinear_form(gcm::ZZMatrix; check::Bool)" begin
    # finite type
    @test cartan_bilinear_form(cartan_matrix(:A, 2)) == ZZ[2 -1; -1 2]
    @test cartan_bilinear_form(cartan_matrix(:B, 4)) ==
      ZZ[4 -2 0 0; -2 4 -2 0; 0 -2 4 -2; 0 0 -2 2]
    @test cartan_bilinear_form(cartan_matrix(:C, 4)) ==
      ZZ[2 -1 0 0; -1 2 -1 0; 0 -1 2 -2; 0 0 -2 4]
    @test cartan_bilinear_form(cartan_matrix(:F, 4)) ==
      ZZ[4 -2 0 0; -2 4 -2 0; 0 -2 2 -1; 0 0 -1 2]
    @test cartan_bilinear_form(cartan_matrix(:G, 2)) == ZZ[2 -3; -3 6]

    @test cartan_bilinear_form(cartan_matrix((:B, 2), (:A, 1))) == ZZ[4 -2 0; -2 2 0; 0 0 2]
    @test cartan_bilinear_form(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2]) ==
      ZZ[4 0 -2 0; 0 2 0 -2; -2 0 2 0; 0 -2 0 4]

    # affine type
    @test cartan_bilinear_form(ZZ[2 -2; -2 2]) == ZZ[2 -2; -2 2]
    @test cartan_bilinear_form(ZZ[2 -2 0; -1 2 -1; 0 -2 2]) == ZZ[2 -2 0; -2 4 -2; 0 -2 2]
    @test cartan_bilinear_form(ZZ[2 -4; -1 2]) == ZZ[2 -4; -4 8]
  end

  @testset "cartan_type_with_ordering" begin
    @testset "A_n" begin
      @testset "A_$n" for n in 1:8
        cm = cartan_matrix(:A, n)
        for perm in symmetric_group(n)
          type, ord = cartan_type_with_ordering(
            permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm);
            check=false,
          )
          @test type == [(:A, n)]
          # automorphism group of Dynkin diagram is S_2
          @test ord in [Vector{Int}(perm), reverse(Vector{Int}(perm))]
          @test !is_one(perm) || ord == 1:n
        end
      end
    end

    @testset "B_n" begin
      @testset "B_$n" for n in 2:8
        cm = cartan_matrix(:B, n)
        for perm in symmetric_group(n)
          type, ord = cartan_type_with_ordering(
            permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm);
            check=false,
          )
          if n == 2
            # B_2 is isomorphic to C_2; we decide based on the ordering
            if is_one(perm)
              @test type == [(:B, n)]
              @test ord == [1, 2]
            else
              @test type == [(:C, 2)]
              @test ord == [1, 2]
            end
          else
            @test type == [(:B, n)]
            # no automorphisms of Dynkin diagram
            @test ord == Vector{Int}(perm)
          end
          @test !is_one(perm) || ord == 1:n
        end
      end
    end

    @testset "C_n" begin
      @testset "C_$n" for n in 2:8
        cm = cartan_matrix(:C, n)
        for perm in symmetric_group(n)
          type, ord = cartan_type_with_ordering(
            permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm);
            check=false,
          )
          if n == 2
            # C_2 is isomorphic to B_2; we decide based on the ordering
            if is_one(perm)
              @test type == [(:C, n)]
              @test ord == [1, 2]
            else
              @test type == [(:B, 2)]
              @test ord == [1, 2]
            end
          else
            @test type == [(:C, n)]
            # no automorphisms of Dynkin diagram
            @test ord == Vector{Int}(perm)
          end
        end
      end
    end

    @testset "D_n" begin
      @testset "D_$n" for n in 4:8
        cm = cartan_matrix(:D, n)
        for perm in symmetric_group(n)
          type, ord = cartan_type_with_ordering(
            permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm);
            check=false,
          )
          @test type == [(:D, n)]
          if n == 4
            # automorphism group of Dynkin diagram is S_3, only one fixpoint
            @test ord[2] == perm(2)
          else
            # automorphism group of Dynkin diagram is S_2
            aut = @perm (n - 1, n)
            @test ord in [Vector{Int}(perm), Vector{Int}(aut * perm)]
          end
          @test !is_one(perm) || ord == 1:n
        end
      end
    end

    @testset "E_n" begin
      @testset "E_$n" for n in 6:8
        cm = cartan_matrix(:E, n)
        for perm in symmetric_group(n)
          type, ord = cartan_type_with_ordering(
            permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm);
            check=false,
          )
          @test type == [(:E, n)]
          if n == 6
            # automorphism group of Dynkin diagram is S_2
            aut = @perm (1, 6)(3, 5)
            @test ord in [Vector{Int}(perm), Vector{Int}(aut * perm)]
          else
            # no automorphisms of Dynkin diagram
            @test ord == Vector{Int}(perm)
          end
          @test !is_one(perm) || ord == 1:n
        end
      end
    end

    @testset "F_$n" for n in 4:4
      cm = cartan_matrix(:F, n)
      for perm in symmetric_group(n)
        type, ord = cartan_type_with_ordering(
          permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm); check=false
        )
        @test type == [(:F, n)]
        # no automorphisms of Dynkin diagram
        @test ord == Vector{Int}(perm)
        @test !is_one(perm) || ord == 1:n
      end
    end

    @testset "G_$n" for n in 2:2
      cm = cartan_matrix(:G, n)
      for perm in symmetric_group(n)
        type, ord = cartan_type_with_ordering(
          permutation_matrix(ZZ, inv(perm)) * cm * permutation_matrix(ZZ, perm); check=false
        )
        @test type == [(:G, n)]
        # no automorphisms of Dynkin diagram
        @test ord == Vector{Int}(perm)
        @test !is_one(perm) || ord == 1:n
      end
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
