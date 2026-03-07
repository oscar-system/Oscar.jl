@testset "Finite group actions" begin
  E8 = root_lattice(:E, 8)
  A4 = lattice_in_same_ambient_space(E8, basis_matrix(E8)[4:7, :])
  f = matrix(QQ, 4, 4, [0 1 0 0; 0 0 1 0; 0 0 0 1; -1 -1 -1 -1])
  F, _ , _ = invariant_coinvariant_pair(A4, f; ambient_representation=false)
  @test rank(F) == 0

  G = matrix_group(f)
  @test Oscar.is_isometry_group(A4, G, false; is_special=true, is_stable=true)

  ok, Gsat = is_saturated_with_saturation(A4, matrix_group(f); ambient_representation=false, stable=true)
  @test !ok
  Gsat_amb = extend_to_ambient_space(A4, Gsat)
  _, C, H = invariant_coinvariant_pair(E8, Gsat_amb)
  @test C == A4

  H1 = pointwise_stabilizer_orthogonal_complement_in_orthogonal_group(E8, A4; special=true)
  @test describe(H1) == "A5"
  @test Oscar.is_isometry_group(ambient_space(A4), H1; is_special=true)

  H2 = pointwise_stabilizer_in_orthogonal_group(E8, A4)
  @test order(H2) == 120

  T = orthogonal_submodule(E8, A4)
  A2 = lattice_in_same_ambient_space(T, basis_matrix(T)[1:2, :])

  H3 = setwise_stabilizer_in_orthogonal_group(E8, A2; special=true)
  @test order(H3) == 311040

  H4 = restrict_to_lattice(A2, H3)
  @test order(H4) == 12

  H5 = setwise_stabilizer_in_orthogonal_group(E8, A2+A4)
  @test order(H5) == 2880
  
  AB = integer_lattice(gram=ZZ[4 -2 0 0; -2 4 0 0; 0 0 16 42; 0 0 42 112])
  BL = QQ[1 0 0 0; 0 1 0 0; 1//2 1//2 4 -3//2; 0 1//2 3//2 -1//2]
  A = lattice_in_same_ambient_space(AB, basis_matrix(AB)[1:2,:])
  B = lattice_in_same_ambient_space(AB, basis_matrix(AB)[3:4,:])
  GA = orthogonal_group(A)
  GB = orthogonal_group(B)
  L = lattice_in_same_ambient_space(AB, BL)
  GAB,_ = Oscar._isometry_group_via_decomposition(AB)
  G1 = Oscar._overlattice_stabilizer(G,AB,L)
  G2 = matrix_group(Oscar.stabilizer_in_diagonal_action(L,A,B,GA,GB))
  G3 = matrix_group(Oscar.stabilizer_in_diagonal_action(L,B,A,GB,GA))
  @test order(G1) == order(G2) == order(G3)
end
