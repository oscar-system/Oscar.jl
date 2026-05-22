@testset "Torsion quadratic modules with isometry" begin
  Tf = torsion_quadratic_module_with_isometry(QQ[-1//50;])
  T = underlying_module(Tf)
  @test order(T) == 50
  @test order_of_isometry(Tf) == 1
  @test order(underlying_module(primary_part(Tf, 5)[1])) == 25

  E8 = root_lattice(:E, 8)
  m = QQ[ 1  0  0  0  0  0  0  0;
          0  1  0  0  0  0  0  0;
         -2 -4 -5 -4 -3 -2 -1 -2;
          1  2  3  3  2  1  0  1;
         -1 -2 -3 -3 -2 -1  0 -2;
          1  2  3  2  1  0  0  2;
          0  0  0  0  0  1  0  0;
          2  4  6  5  4  3  2  3]
  Lf = integer_lattice_with_isometry(E8, m)
  F = kernel_lattice(Lf, 12)

  T, fT = discriminant_group(F)
  Tf1 = torsion_quadratic_module_with_isometry(T, hom(fT))
  @test length(submodules(Tf1)) == 2
  @test order(automorphism_group(Tf1)) == 3

  Tf2 = torsion_quadratic_module_with_isometry(T, fT)
  Tf3 = torsion_quadratic_module_with_isometry(T, matrix(fT))
  Tf4 = torsion_quadratic_module_with_isometry(T, abelian_group_homomorphism(hom(fT)))
  Tf5 = torsion_quadratic_module_with_isometry(T, first(gens(matrix_group([m]))))
  @test Tf1 == Tf2 == Tf3 == Tf4 == Tf5

  Tf6 = torsion_quadratic_module_with_isometry(gram_matrix_quadratic(T), matrix(fT))
  @test is_isomorphic_with_map(Tf1, Tf6)[1]

  Tf7 = torsion_quadratic_module_with_isometry(-gram_matrix_quadratic(T), matrix(fT))
  @test is_anti_isomorphic_with_map(Tf1, Tf7)[1]

  B = matrix(QQ, 5, 6, [1 2 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]);
  G = matrix(QQ, 6, 6, [2 -1 0 0 0 0; -1 2 -1 0 0 0; 0 -1 2 -1 0 -1; 0 0 -1 2 -1 0; 0 0 0 -1 2 0; 0 0 -1 0 0 2]);
  L = integer_lattice(B; gram=G);
  f = matrix(QQ, 6, 6, [1 0 0 0 0 0; -1 -1 -2 -2 -1 -1; 0 0 1 1 0 1; 0 0 -1 0 0 -1; 0 0 1 0 0 0; 0 0 0 0 1 0]);
  Lf = integer_lattice_with_isometry(L, f);
  Tf = discriminant_group(TorQuadModuleWithIsom, Lf)
  gene = [3*Tf.T[1]]

  S, j = sub(Tf, gene)
  @test is_injective(j)
  @test order(underlying_module(S)) == 2

  K, i = orthogonal_submodule(Tf, underlying_module(S))
  @test order(underlying_module(K)) == 3
  @test is_injective(i)
end
