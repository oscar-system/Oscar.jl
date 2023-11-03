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

  @testset "cartan_to_coxeter_matrix(mat::ZZMatrix; generalized::Bool=true)" begin end

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

    # how test for loops ?
    @test is_cartan_matrix(ZZ[2 -1 -1; -1 2 -1; -1 -1 2]; generalized=false) == false broken =
      true
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

    @test cartan_type(cartan_matrix(:C, 3); check=false) == [(:C, 3)]
    @test cartan_type(ZZ[2 -1 0; -2 2 -1; 0 -1 2]; check=false) == [(:C, 3)]

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
end
