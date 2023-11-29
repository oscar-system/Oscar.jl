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

  @testset "cartan_type(gcm::ZZMatrix; check::Bool)" begin
    # test if we follow our own conventions
    @test cartan_type(cartan_matrix(:A, 3); check=false) == [(:A, 3)]
    @test cartan_type(cartan_matrix(:B, 2); check=false) == [(:B, 2)]
    @test cartan_type(cartan_matrix(:C, 2); check=false) == [(:C, 2)]

    # tests for irreducibles
    @test cartan_type(cartan_matrix(:A, 1); check=false) == [(:A, 1)]
    @test cartan_type(cartan_matrix(:A, 2); check=false) == [(:A, 2)]

    @test cartan_type(cartan_matrix(:B, 3); check=false) == [(:B, 3)]
    @test cartan_type(ZZ[2 -2 0; -1 2 -1; 0 -1 2]; check=false) == [(:B, 3)]

    # Cn
    @test cartan_type(cartan_matrix(:C, 3); check=false) == [(:C, 3)]
    @test cartan_type(ZZ[2 -1 0; -2 2 -1; 0 -1 2]; check=false) == [(:C, 3)]

    # Dn
    @test cartan_type(cartan_matrix(:D, 4); check=false) == [(:D, 4)]
    @test cartan_type(ZZ[2 -1 -1 -1; -1 2 0 0; -1 0 2 0; -1 0 0 2]) == [(:D, 4)]
    @test cartan_type(cartan_matrix(:D, 6); check=false) == [(:D, 6)]

    @test cartan_type(cartan_matrix(:E, 6); check=false) == [(:E, 6)]
    @test cartan_type(cartan_matrix(:E, 7); check=false) == [(:E, 7)]
    @test cartan_type(cartan_matrix(:E, 8); check=false) == [(:E, 8)]
    @test cartan_type(cartan_matrix(:F, 4); check=false) == [(:F, 4)]
    @test cartan_type(cartan_matrix(:G, 2); check=false) == [(:G, 2)]

    # for F4 and G2 we also allow the transposed Cartan matrix
    @test cartan_type(transpose(cartan_matrix(:F, 4))) == [(:F, 4)]
    @test cartan_type(transpose(cartan_matrix(:G, 2))) == [(:G, 2)]

    # test decomposable Cartan matrices
    @test cartan_type(cartan_matrix((:A, 1), (:A, 2)); check=false) == [(:A, 1), (:A, 2)]
    @test cartan_type(cartan_matrix((:A, 1), (:B, 2)); check=false) == [(:A, 1), (:B, 2)]
    @test cartan_type(cartan_matrix((:C, 2), (:B, 2)); check=false) == [(:C, 2), (:B, 2)]
    @test cartan_type(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2]; check=false) ==
      [(:B, 2), (:C, 2)]
  end

  @testset "cartan_type_with_ordering(gcm::ZZMatrix; check::Bool)" begin
    # we only test the ordering here the type detection is tested in "cartan_type(gcm::ZZMatrix; check::Bool)"

    # An
    _, ord = cartan_type_with_ordering(cartan_matrix(:A, 1))
    @test ord == [1]

    _, ord = cartan_type_with_ordering(cartan_matrix(:A, 2))
    @test ord == [1, 2]

    # Bn
    _, ord = cartan_type_with_ordering(cartan_matrix(:B, 2))
    @test ord == [1, 2]

    _, ord = cartan_type_with_ordering(cartan_matrix(:B, 3))
    @test ord == [1, 2, 3]

    # Cn
    _, ord = cartan_type_with_ordering(cartan_matrix(:C, 2))
    @test ord == [1, 2]

    _, ord = cartan_type_with_ordering(cartan_matrix(:C, 3))
    @test ord == [1, 2, 3]

    # Dn
    _, ord = cartan_type_with_ordering(cartan_matrix(:D, 4))
    @test ord == [1, 2, 3, 4]

    _, ord = cartan_type_with_ordering(ZZ[2 -1 -1 -1; -1 2 0 0; -1 0 2 0; -1 0 0 2])
    @test ord == [2, 1, 3, 4]

    # En
    _, ord = cartan_type_with_ordering(cartan_matrix(:E, 6))
    @test ord == [1, 3, 4, 2, 5, 6]

    _, ord = cartan_type_with_ordering(cartan_matrix(:E, 7))
    @test ord == [1, 3, 4, 2, 5, 6, 7]

    _, ord = cartan_type_with_ordering(cartan_matrix(:E, 8))
    @test ord == [1, 3, 4, 2, 5, 6, 7, 8]

    # F4
    _, ord = cartan_type_with_ordering(cartan_matrix(:F, 4))
    @test ord == [1, 2, 3, 4]

    _, ord = cartan_type_with_ordering(ZZ[2 -1 0 0; -1 2 -2 0; 0 -1 2 -1; 0 0 -1 2])
    @test ord == [4, 3, 2, 1]

    # G2
    _, ord = cartan_type_with_ordering(cartan_matrix(:G, 2))
    @test ord == [1, 2]

    _, ord = cartan_type_with_ordering(transpose(cartan_matrix(:G, 2)))
    @test ord == [2, 1]
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
