@testset "Primitive extensions and embeddings" begin
  # Compute orbits of short vectors
  k = integer_lattice(; gram=matrix(QQ,1,1,[4]))
  E8 = root_lattice(:E, 8)
  ok, sv = primitive_embeddings(E8, k; classification=:sub)
  @test ok
  @test length(sv) == 1
  ok, sv = primitive_embeddings(rescale(E8, 2), rescale(k, QQ(1//2)); check=false)
  @test !ok
  @test is_empty(sv)
  @test isempty(primitive_embeddings(rescale(E8, -1), k; check=false)[2])

  k = integer_lattice(; gram=matrix(QQ,1,1,[6]))
  E7 = root_lattice(:E, 7)
  ok, sv = primitive_embeddings(E7, k; classification=:emb)
  @test ok
  @test length(sv) == 2
  q = discriminant_group(E7)
  p, z, n = signature_tuple(E7)
  ok, _ = primitive_embeddings(q, (p,n), E7; classification=:none)
  @test ok
  A5 = root_lattice(:A, 5)
  ok, sv = primitive_embeddings(A5, k)
  @test ok
  @test length(sv) == 2

  k = integer_lattice(; gram=matrix(QQ,1,1,[2]))
  ok, sv = primitive_embeddings(E7, k)
  @test ok
  @test length(sv) == 1

  ok, _ = primitive_embeddings(q, (p,n), E7; classification=:none)
  @test ok

  @test !primitive_embeddings(rescale(E7, 2), k; classification=:none, check=false)[1]

  B = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1]);
  G = matrix(QQ, 8, 8, [-4 2 0 0 0 0 0 0; 2 -4 2 0 0 0 0 0; 0 2 -4 2 0 0 0 2; 0 0 2 -4 2 0 0 0; 0 0 0 2 -4 2 0 0; 0 0 0 0 2 -4 2 0; 0 0 0 0 0 2 -4 0; 0 0 2 0 0 0 0 -4]);
  L = integer_lattice(B; gram=G);
  f = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; -2 -4 -6 -4 -3 -2 -1 -3; 2 4 6 5 4 3 2 3; -1 -2 -3 -3 -3 -2 -1 -1; 0 0 0 0 1 0 0 0; 1 2 3 3 2 1 0 2]);
  Lf = integer_lattice_with_isometry(L, f);
  F = invariant_lattice(Lf)
  C = coinvariant_lattice(Lf)
  reps = @inferred admissible_equivariant_primitive_extensions(F, C, Lf^0, 5)
  @test length(reps) == 1
  @test is_of_same_type(Lf, reps[1])
  _, reps = @inferred primitive_extensions(lattice(F), lattice(C); form_over=[discriminant_group(L)], classification=:embemb)
  @test length(reps) == 2
  @test reps[1] == (L, lattice(F), lattice(C))

  C = integer_lattice(; gram = QQ[-4 -2 -2  0  2 -2 -1 -3 -2 -1 -3 -2 -4  0 -2 -1;
                                  -2 -4 -1 -3  1 -1  1 -3 -1 -5 -3  2 -1  2 -2 -1;
                                  -2 -1 -4  0  2 -2 -1 -3 -2 -2 -3 -3 -4  0 -2  0;
                                   0 -3  0 -8  0 -2  4 -2  0 -7 -1  4  1  2  1  0;
                                   2  1  2  0 -4  0  0  2  2  1  1  2  3  0  1  0;
                                  -2 -1 -2 -2  0 -4  0 -4 -1 -2 -2 -2 -2  0  0 -1;
                                  -1  1 -1  4  0  0 -4  0 -1  3 -1 -3 -2 -1 -3 -1;
                                  -3 -3 -3 -2  2 -4  0 -8 -2 -4 -3 -2 -2  1 -1 -3;
                                  -2 -1 -2  0  2 -1 -1 -2 -4 -2 -3 -3 -4 -1 -3 -1;
                                  -1 -5 -2 -7  1 -2  3 -4 -2 -10 -4 3 -1  2 -2 -1;
                                  -3 -3 -3 -1  1 -2 -1 -3 -3 -4 -6 -1 -4  1 -5 -1;
                                  -2  2 -3  4  2 -2 -3 -2 -3  3 -1 -8 -5 -3 -1 -1;
                                  -4 -1 -4  1  3 -2 -2 -2 -4 -1 -4 -5 -8 -2 -4  0;
                                   0  2  0  2  0  0 -1  1 -1  2  1 -3 -2 -4  0 -1;
                                  -2 -2 -2  1  1  0 -3 -1 -3 -2 -5 -1 -4  0 -8 -2;
                                  -1 -1  0  0  0 -1 -1 -3 -1 -1 -1 -1  0 -1 -2 -4])

  L = integer_lattice(; gram = QQ[0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 1 -2 1 0 0 0 1 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 1 -2 0 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 1 0 0 0 0 -2 0 0 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2 1 0 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 1 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 1 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -2 0 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -2 0;
                                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4])

  ok, reps = primitive_embeddings(L, C; check=false)
  @test length(reps) == 1

  ## Odd case
  B = matrix(QQ, 5, 5 ,[1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]);
  G = matrix(QQ, 5, 5 ,[3, 1, 0, 0, 0, 1, 3, 1, 1, -1, 0, 1, 3, 0, 0, 0, 1, 0, 3, 0, 0, -1, 0, 0, 3]);
  L = integer_lattice(B, gram = G);
  k = lattice_in_same_ambient_space(L, B[2:2, :])
  N = orthogonal_submodule(L, k)
  ok, reps = primitive_extensions(k, N; glue_order=[3])
  @test ok
  @test length(reps) == 2

  GL = genus(torsion_quadratic_module(QQ[1;]), (3, 3))
  A2 = root_lattice(:A, 2)
  ok, reps = primitive_embeddings(GL, A2)
  @test ok
  @test !is_even(reps[1][3])

  k = integer_lattice(; gram=QQ[4;])
  I = direct_sum(hyperbolic_plane_lattice(), k)[1]
  @test primitive_embeddings(I, k; classification=:none)[1]
end

@testset "Equivariant primitive extensions" begin
  A2 = root_lattice(:A, 2)
  Fs = integer_lattice_with_isometry(root_lattice(:E, 6); neg=false)
  q = torsion_quadratic_module(QQ[0;])

  ok, reps = equivariant_primitive_extensions(A2, Fs; form_over=[q], classification=:embemb)
  @test ok
  @test length(reps) == 6

  ok, reps = equivariant_primitive_extensions(reps[1][2], reps[1][3]; form_over=[q], classification=:embemb)
  @test ok
  @test length(reps) == 2

  Fs2 = integer_lattice_with_isometry(integer_lattice(;gram = QQ[3;]); neg=false)
  ok, reps = equivariant_primitive_extensions(A2, Fs2; even=false)

  @test ok
  @test length(reps) == 9
end
