@testset "odd" begin
  R = residue_ring(ZZ, ZZ(3))[1]
  F = GF(3)
  G = zero_matrix(R, 0, 0)
  @test Oscar._orthogonal_grp_gens_odd(G, 3) == []
  G = diagonal_matrix([R(1)])
  @test Oscar._orthogonal_grp_gens_odd(G, 3) == [R[-1;]]
  Oscar._gens(G,2, 3)


  G = diagonal_matrix([R(1), R(1)])
  L = Oscar._orthogonal_grp_gens_odd(G, 3)
  L = [change_base_ring(F, lift(g)) for g in L]
  @test order(matrix_group(L)) == 8
  Oscar._gens(G,2, 3)


  G = diagonal_matrix([R(1),R(2)])
  L = Oscar._orthogonal_grp_gens_odd(G, 3)
  L = [change_base_ring(F, lift(g)) for g in L]
  @test order(matrix_group(L)) == 4
  Oscar._gens(G,2, 3)


  R = residue_ring(ZZ, 3^5)[1]
  F = GF(3)
  G = diagonal_matrix([R(x) for x in [3*1, 3*1]])
  L = [change_base_ring(F, lift(g)) for g in Oscar._gens_mod_p(G, 3)]
  @test order(matrix_group(L))==8
  # Oscar._gens(G,2, 3)



  G = diagonal_matrix([R(x) for x in [1, 3, 3, 9, 2*27]])
  L = [change_base_ring(F,lift(g)) for g in Oscar._gens_mod_p(G, 3)]
  @test order(matrix_group(L))==1259712
  #Oscar._gens(G,2, 3)

  R = residue_ring(ZZ, ZZ(3)^5)[1]  #det(::zzModMatrix) is broken and cannot be used here
  G = diagonal_matrix([R(x) for x in [2*27, 9, 3, 3, 1 ]])
  Oscar._gens(G,2, 3)

end
@testset "_gens_mod_2" begin
  R = residue_ring(ZZ, ZZ(2)^10)[1]
  F = GF(2)
  U = matrix(R, 2, 2, [0, 1, 1, 0])
  V = matrix(R, 2, 2, [2, 1, 1, 2])
  W0 = matrix(R, 2, 2, [1, 0, 0, 3])
  W1 = matrix(R, 2, 2, [1, 0, 0, 1])
  w1 = matrix(R, 1, 1, [1])
  w3 = matrix(R, 1, 1, [3])
  w5 = matrix(R, 1, 1, [5])
  w7 = matrix(R, 1, 1, [7])


  @test length(Oscar._orthogonal_gens_bilinear(zero_matrix(R,0,0))) == 0
  L = [[U], [V], [V, w1, w1], [U, w1, w3], [w1, w7], [U, V], [V, w3],
  [U, w1, w5]]
  GG = [diagonal_matrix(G) for G in L]
  GG = [Oscar._orthogonal_gens_bilinear(G) for G in GG]
  GG = [matrix_group([change_base_ring(F, lift(g)) for g in G]) for G in GG]
  GG = [order(G) for G in GG]
  @test GG == [6, 6, 48, 48, 2, 720, 6, 48]


  GG = [diagonal_matrix(G) for G in [[U,V], [U, W0],[V,W0],[U,U,W0],[U,V,W0], [W1]]]
  GG = [Oscar._orthogonal_grp_quadratic(G) for G in GG]
  GG = [matrix_group([change_base_ring(F, lift(g)) for g in G]) for G in GG]
  GG = [order(G) for G in GG]
  @test GG == [120, 8, 24, 1152,1920,2]


  _ker_gens = Oscar._ker_gens
  W0 = matrix(R, 2, 2, [0, 1, 1, 1])
  W1 = matrix(R, 1, 1, [1])
  W2 = matrix(R, 2, 2, [2, 1, 1, 1])
  G = diagonal_matrix([V*4, U*2, U])
  gen1 = _ker_gens(G, 2, 4, [0, 0, 0])
  gen1 = [Hecke.hensel_qf(G, g, 1, 2, 2) for g in gen1]
  @test length(gen1)==8

  G = diagonal_matrix([V*4, W1*2, W2])
  gen1 = _ker_gens(G, 2, 3, [0, 1, 1])
  gen1 = [Hecke.hensel_qf(G, g, 1, 2, 2) for g in gen1]
  @test length(gen1)==4

  G = diagonal_matrix([V*4, U*2, W1])
  gen1 = _ker_gens(G, 2, 4, [0, 0, 1])
  gen1 = [Hecke.hensel_qf(G, g, 1, 2, 2) for g in gen1]
  length(gen1)
  @test length(gen1) == 2

  G = diagonal_matrix([V*4, W2*2, U])
  gen1 = _ker_gens(G, 2, 4, [0, 1, 0])
  gen1 = [Hecke.hensel_qf(G, g, 1, 2, 2) for g in gen1]
  @test length(gen1) == 6

  G = diagonal_matrix([W1*4, U*2, U])
  gen1 = _ker_gens(G, 1, 3, [1,0,0])
  gen1 = [Hecke.hensel_qf(G, g, 1, 2, 2) for g in gen1]
  @test length(gen1)==6

  G = diagonal_matrix([W1*4, U*2, W0])
  gen1 = _ker_gens(G, 1, 3, [1,0,1])
  gen1 = [Hecke.hensel_qf(G, g, 1, 2, 2) for g in gen1]
  @test length(gen1)==5

  G = diagonal_matrix([W0*4, W0*2, U])
  gen1 = _ker_gens(G, 2, 4, [1,1,0])
  gen1 = [Hecke.hensel_qf(G, g, 1, 2, 2) for g in gen1]
  @test length(gen1)==6

  G = diagonal_matrix([W1*4, U*2, W2*2, W1])
  gen1 = _ker_gens(G, 2, 6, [1,1,1])
  @test length(gen1)==0



  L =[[2*U,2*V, w1], [2*U, w1,w7], [2*V,w5], [4*U,2*U,V,w5],
      [8*U,V,w1,w3], [4*U,2*U,w1,w1],
      [8*w1,4*V,2*U,2*w1,2*w5,U,V],
      #[8*U,4*V,2*U,2*w1,2*w5,U,V],
      [8*V,4*w1,2*U,2*w1,2*w5,U,V],
      [8*V,4*V,4*w7,2*U,2*w1,2*w5,w1]
    ]
  M = [diagonal_matrix(G) for G in L]
  GG = [Oscar._gens_mod_2(G) for G in M]
  GG = [matrix_group([change_base_ring(F, lift(g)) for g in G]) for G in GG]
  GG = [order(G) for G in GG]
  @test GG == [720,
      24,
      6,
      1179648,
      12288,
      24576,
      1781208836997120,
      #7295831396340203520,
      890604418498560,
      57982058496]
  for G in M
    Oscar._gens(G, 8,2)
  end
end
  
@testset "stabilisers of subspaces under isometries" begin  
  I = integer_lattice(gram=ZZ[1;])
  for i in 1:5
    @test Oscar._test_stabilizers(rescale(root_lattice(:A,i),2),2)
    @test Oscar._test_stabilizers(direct_sum(rescale(root_lattice(:A,i),2),I)[1],2)
  end
  for i in 4:6
    @test Oscar._test_stabilizers(rescale(root_lattice(:D,i),2),2)
    @test Oscar._test_stabilizers(direct_sum(rescale(root_lattice(:D,i),2),I)[1],2)
  end
  @test Oscar._test_stabilizers(rescale(root_lattice(:D,4),3),3)
  @test Oscar._test_stabilizers(rescale(root_lattice(:A,4),3),3)
  for i in 1:4
    @test Oscar._test_stabilizers(rescale(root_lattice(:A,i),7),7)
  end

  DD = TorQuadModule[]
  B = matrix(QQ, 1, 1 ,[1//2]);
  G = matrix(QQ, 1, 1 ,[2]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1;];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 1, 1 ,[1//2]);
  G = matrix(QQ, 1, 1 ,[6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1;];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 1, 1 ,[1//2]);
  G = matrix(QQ, 1, 1 ,[4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1;];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 1, 1 ,[1//2]);
  G = matrix(QQ, 1, 1 ,[8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1;];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[12, 0, 0, 6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1; -1//2 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[4, 2, 2, 4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1; -1 2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[0, 2, 2, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1; 1 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//4, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[80, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1; -1//4 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[14, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0; 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[2, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0; 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[2, 0, 0, 2]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0; 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//4, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[24, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1; -7//2 -3];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[4, 0, 0, 4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0; 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//2, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[8, 4, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1; -1 2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//8, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[160, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1; -1//4 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 2, 2 ,[1//8, 0, 0, 1//2]);
  G = matrix(QQ, 2, 2 ,[224, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1; 55//4 14];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//2, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[12, 0, 0, 0, 2, 0, 0, 0, 6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0; -1 0 1; 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//2, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[20, 0, 0, 0, 0, 2, 0, 2, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1; -1 1 0; 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//4, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[48, 0, 0, 0, 40, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 1; -1 -1//2 -3; 175//4 21 135];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//2, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[4, 2, 0, 2, 4, 0, 0, 0, 10]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0; -1 2 1; -5 10 6];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//2, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[0, 2, 0, 2, 0, 0, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0; 1 0 0; 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//2, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[4, 2, 0, 2, 4, 0, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0; -1 2 1; 7 -14 -6];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//2, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[0, 2, 0, 2, 0, 0, 0, 0, 2]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0; 1 0 0; 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[80, 0, 0, 0, 4, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0; -1 0 1; 1//4 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//4, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[0, 8, 0, 8, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 1; 1//2 0 0; 0 -7//2 -3];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//4, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[8, 0, 0, 0, 8, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1; 0 1//2 0; -7//2 0 -3];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[24, 0, 0, 0, 4, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0; -1 0 1; 7//2 0 -3];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//8, 0, 0, 0, 1//8, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[64, 32, 0, 32, 64, 0, 0, 0, 24]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 1; -1//4 1//2 0; 0 3//4 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//8, 0, 0, 0, 1//8, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[0, 32, 0, 32, 0, 0, 0, 0, 24]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 1; 1//4 0 0; 0 3//4 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[40, 0, 0, 0, 0, 4, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1; -1 1 0; 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[56, 0, 0, 0, 0, 4, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1; -1 1 0; 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//8, 0, 0, 0, 1//8, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[32, 0, 0, 0, 160, 0, 0, 0, 56]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1; 7//4 1//4 2; 105//4 7//2 30];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//8, 0, 0, 0, 1//8, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[96, 0, 0, 0, 224, 0, 0, 0, 56]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 1; 0 1//4 2; 7//4 3//2 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//8, 0, 0, 0, 1//8, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[32, 0, 0, 0, 96, 0, 0, 0, 56]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1; 7//4 1//4 2; -21//2 -7//4 -12];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//2, 0, 0, 0, 1//2, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[0, 4, 0, 4, 0, 0, 0, 0, 4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0; 1 0 0; 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 3, 3 ,[1//4, 0, 0, 0, 1//4, 0, 0, 0, 1//2]);
  G = matrix(QQ, 3, 3 ,[0, 16, 0, 16, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 1; 1//2 0 0; 0 -5//2 -2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[28, 0, 0, 0, 0, 4, 2, 0, 0, 2, 4, 0, 0, 0, 0, 6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0; -1 -1 2 0; 15 14 -28 1; -15//2 -7 14 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 14, 0, 0, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[4, 2, 0, 0, 2, 4, 0, 0, 0, 0, 2, 0, 0, 0, 0, 6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; -1 2 1 0; 1 -2 0 0; 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[16, 0, 0, 0, 0, 0, 8, 0, 0, 8, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 1; -1//2 1//2 0 0; -2 0 -5//2 -2; 15//4 0 5 4];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[112, 0, 0, 0, 0, 8, 0, 0, 0, 0, 24, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0 1; -1//2 0 1//2 0; 4 -7//2 -7 -3; -15//4 0 7//2 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[112, 0, 0, 0, 0, 56, 0, 0, 0, 0, 4, 0, 0, 0, 0, 12]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 1 0; -2 -1 0 1; 31//2 15//2 0 -7; -63//4 -15//2 0 7];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[16, 8, 0, 0, 8, 16, 0, 0, 0, 0, 8, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1; -1//2 1 1//2 0; 1//2 -9//2 0 -3; -7//2 28 0 18];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[16, 8, 0, 0, 8, 16, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1; -1 1 1 0; 1 -3//2 0 0; -1//2 1 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 8, 0, 0, 8, 0, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1; 1 0 1 0; 0 -1//2 0 0; -1//2 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[8, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1; -1 1 1 0; 1//2 -1 0 0; 0 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[8, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1; -1 1 1 0; 1//2 -1 0 0; 0 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[24, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0 1; -8 -7 1 0; 15//2 13//2 0 0; -7//2 -3 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[80, 0, 0, 0, 0, 0, 4, 0, 0, 4, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0; -1 1 0 0; 1 0 0 1; -1//4 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 32, 0, 0, 32, 0, 0, 0, 0, 0, 160, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1; 1//4 0 0 0; 0 0 1//4 0; 0 -5//4 0 -1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[64, 32, 0, 0, 32, 64, 0, 0, 0, 0, 160, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1; -1//4 1//2 1//4 0; -5//4 5//4 3//2 -1; 5//4 -5//2 -3//2 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 32, 0, 0, 32, 0, 0, 0, 0, 0, 224, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1; 1//4 0 0 0; 0 -5//4 1//4 -1; 0 35//2 -15//4 14];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[64, 32, 0, 0, 32, 64, 0, 0, 0, 0, 224, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1; -1//4 1//2 1//4 0; 7//4 -19//4 -3//2 -1; 35 -385//4 -30 -21];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[8, 0, 0, 0, 0, 24, 0, 0, 0, 0, 4, 0, 0, 0, 0, 12]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0; -1 1 0 1; 0 -1//2 0 0; 3//2 -3//2 0 -1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 8, 0, 0, 8, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 1 0; 1 0 0 1; 0 -1//2 0 0; -5//2 0 0 -2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//4, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[40, 0, 0, 0, 0, 8, 4, 0, 0, 4, 8, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0; -1 -1 2 0; -9 -10 20 1; 819//2 455 -910 -45];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 4, 0, 0, 4, 0, 0, 0, 0, 0, 28, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//8, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[224, 0, 0, 0, 0, 0, 16, 0, 0, 16, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 1; -1//2 1//2 0 0; -2 0 -5//2 -2; 55//4 0 35//2 14];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//8, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[32, 0, 0, 0, 0, 0, 16, 0, 0, 16, 0, 0, 0, 0, 0, 56]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 1; -1//2 1//2 0 0; -3 0 -7//2 -3; 35//4 0 21//2 9];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 4, 4 ,[1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 4, 4 ,[0, 4, 0, 0, 4, 0, 0, 0, 0, 0, 8, 4, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[28, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 2, 4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0; -1 1 0 0 0; 1 0 0 0 1; -1 0 0 -1 2; 15//2 0 0 7 -14];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[12, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0; -1 -1 2 1 0; 2 1 -2 0 1; -15 -7 14 0 -6; 91//2 21 -42 0 18];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 1 0 0; 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 2]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 -1 2 1; 0 0 1 -2 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 1 0 0; 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 -1 2 0; 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[112, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 1; -1//2 -1//2 1 0 0; 15//2 7 -14 1//2 0; -641 -595 2373//2 -42 -3; 1275//4 595//2 -595 21 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[80, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 1; -1 1 0 -1 2; 0 0 -1//2 0 0; 1 -3//2 0 1 -2; 21//4 -15//2 0 5 -10];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1//2 0 0 0 0; 0 0 0 1//2 0; 0 -5//2 -1//2 1 -2; 0 -15 -5//2 5 -12];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[48, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0 0 1; -1 0 1 1 0; 3 -1//2 -3 0 0; -3 0 5//2 0 0; 7//4 0 -3//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[16, 8, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; -1//2 1 1//2 0 0; 1//2 -9//2 0 0 -3; 0 0 0 1//2 0; -7//2 28 0 0 18];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1//4 0 0 0 0; 0 0 0 1//4 0; 0 0 1//4 0 0; 0 -1//4 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 64, 32, 0, 0, 0, 32, 64, 0, 0, 0, 0, 0, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1//4 0 0 0 0; 0 0 0 1//4 0; 0 -1//4 -1//4 1//2 0; 0 0 1//4 -1//2 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[80, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0 1 0; -22 -21 -1 2 0; 87 83 2 -4 1; 2730 5209//2 63 -126 -3; -6069//4 -2895//2 -35 70 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[16, 8, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; -1 1 0 1 0; 1//2 -1//2 1//2 0 0; 1//2 -1 0 0 0; 0 0 -1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1 0 0 1 0; 0 -1//2 1//2 0 0; -1//2 0 -1//2 0 0; 0 0 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[16, 8, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; -1 2 1 1 0; -5//2 9//2 3 0 0; 3 -6 -7//2 0 0; -5//2 5 3 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1 0 0 1 0; 0 -1//2 1//2 0 0; -1//2 0 -1//2 0 0; 0 0 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 12]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1 0; 1 0 0 0 1; 0 -1//2 1//2 0 0; 0 0 -1//2 0 0; -3//2 0 0 0 -1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 160, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1//4 0 0 0 0; 0 0 1//4 0 0; 0 0 0 1//4 0; 0 -5//4 0 0 -1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 96, 0, 0, 0, 0, 0, 224, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1//4 0 0 0 0; 0 0 1//4 1//4 0; 0 -5//4 -7//4 -3//2 -1; 0 -105//4 -35 -30 -21];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 96, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; 1//4 0 0 0 0; 0 0 1//4 0 0; 0 -5//4 0 1//4 -1; 0 -15//2 0 5//4 -6];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[64, 32, 0, 0, 0, 32, 64, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 224, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1; -1//4 1//2 1//4 0 0; 1//4 -1//2 0 1//4 0; -7//4 9//4 0 -3//2 -1; -35 175//4 0 -30 -21];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1 0; 1 0 -1 2 0; 0 -1 0 0 1; -3//2 5//2 1 -2 -2; -15//2 15 5 -10 -12];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[8, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1 0; -1 1 -1 2 0; 3 -4 2 -4 1; 0 1//2 0 0 -3; -21//2 27//2 -7 14 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[80, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0; -1 1 0 0 0; 1 0 0 1 0; -1 0 0 0 1; 1//4 0 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[56, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 12]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0; -1 -1 2 0 0; 15 14 -28 1 0; -15 -14 28 0 1; 15//2 7 -14 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[32, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1 1; 0 1//2 1//2 0 0; -3 0 0 -5//2 -2; 0 -1//4 0 0 0; 25//4 0 0 5 4];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[224, 0, 0, 0, 0, 0, 224, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1 1; 0 1//2 1//2 0 0; -3 0 0 -5//2 -2; 0 -1//4 0 0 0; 85//4 0 0 35//2 14];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[32, 0, 0, 0, 0, 0, 224, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1 1; -1//2 1//2 1//2 0 0; -2 -1 0 -5//2 -2; 15//4 11//4 0 5 4; -105//2 -155//4 0 -70 -56];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[40, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0; -1 1 0 0 0; 1 0 0 0 1; -1 0 0 -1 2; -9//2 0 0 -5 10];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//4, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[56, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0; -1 1 0 0 0; 1 0 0 0 1; -1 0 0 -1 2; 15//2 0 0 7 -14];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[64, 32, 0, 0, 0, 32, 64, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1 1; -1//2 1//2 1//2 0 0; 1 -4 0 -5//2 -2; -11//4 37//4 0 5 4; -65//4 55 0 30 24];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//8, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 32, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1 1; 1//2 0 1//2 0 0; 0 -3 0 -5//2 -2; -1//4 0 0 0 0; 0 25//4 0 5 4];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 5, 5 ,[1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 5, 5 ,[0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 1 0 0; 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[12, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 0, 10]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0 0; -1 1 0 0 0 0; 1 0 0 0 1 0; -1 0 0 -1 2 1; -4 0 0 -5 10 6; 15//2 0 0 10 -20 -12];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[16, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0 1; -1//2 1//2 0 0 0 0; 1//2 0 0 0 1//2 0; -1//2 0 0 -1//2 1 0; -1 0 -5//2 1 -2 -2; -15//4 0 -15//2 5//2 -5 -6];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 0 1 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 0, 2, 4]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 0 1 0 0 0; 0 0 0 0 0 1; 0 0 0 0 -1 2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[16, 0, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0 1; -1//2 -1//2 1 0 0 0; 3//2 1 -2 1//2 0 0; -5 -1 -3//2 0 1//2 -3; 33 8 5 0 -7//2 18; -399//4 -49//2 -14 0 21//2 -54];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 14, 0, 0, 0, 0, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 14]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 0, 2, 4, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 6]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 1 0 0; 0 0 -1 2 1 0; 0 0 1 -2 0 0; 0 0 0 0 0 1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[48, 0, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 12]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 1 0; -1 -1 2 1 0 1; -5//2 -5//2 9//2 3 0 0; 19//2 25//2 -25 -15 0 0; -6 -17//2 17 21//2 0 -1; 15//4 5 -10 -6 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1//2 0 0 0 0 0; 0 0 0 1//2 0 0; 0 -7//2 -1//2 1 0 -3; 0 0 0 0 1//2 0; 0 21 7//2 -7 0 18];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 12]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1 0; 1 0 0 0 0 1; 0 0 1//2 0 0 0; 0 -1//2 0 1//2 0 0; 0 0 0 -1//2 0 0; -3//2 0 0 0 0 -1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[80, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0 0 1 0; -1 0 1 -1 2 0; 20 -1 -20 0 0 1; -9 0 17//2 1 -2 0; -1379//2 7//2 655 70 -140 -3; 1239//4 0 -585//2 -35 70 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1 0 0 0 -1 2; 0 0 0 1//2 0 0; 0 0 1//2 0 0 0; 0 -1//2 0 0 0 0; -3//2 0 0 0 1 -2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1 0 0 0 -1 2; 0 0 0 1//2 0 0; 0 -1//2 -1//2 1 0 0; 0 0 1//2 -1 0 0; -3//2 0 0 0 1 -2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 160, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1//4 0 0 0 0 0; 0 0 0 1//4 0 0; 0 0 1//4 0 0 0; 0 0 0 0 1//4 0; 0 -5//4 0 0 0 -1];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 64, 32, 0, 0, 0, 0, 32, 64, 0, 0, 0, 0, 0, 0, 160, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1//4 0 0 0 0 0; 0 0 0 1//4 0 0; 0 0 -1//4 1//2 1//4 0; 0 -5//4 -5//4 5//2 3//2 -1; 0 0 5//4 -5//2 -3//2 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 224, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1//4 0 0 0 0 0; 0 0 0 1//4 0 0; 0 0 1//4 0 0 0; 0 -5//4 0 0 1//4 -1; 0 35//2 0 0 -15//4 14];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 64, 32, 0, 0, 0, 0, 32, 64, 0, 0, 0, 0, 0, 0, 224, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1//4 0 0 0 0 0; 0 0 0 1//4 0 0; 0 0 -1//4 1//2 1//4 0; 0 -5//4 7//4 -7//2 -3//2 -1; 0 -105//4 35 -70 -30 -21];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1 0; 1 0 0 0 0 1; 0 0 0 1//2 0 0; 0 -1//2 -1//2 1 0 0; -5//2 0 1//2 -1 0 -2; -15 0 5//2 -5 0 -12];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[16, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 1 0; -1 1 0 1 0 0; 0 0 -1 0 0 1; 0 -1//2 0 0 0 0; 1//2 0 5//2 0 0 -2; -5//4 0 -5 0 0 4];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1 0 0 0 1 0; 0 0 1//2 0 0 0; 0 -1//2 0 1//2 0 0; -1//2 0 0 -1//2 0 0; 0 0 0 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1 0 0 0 1 0; 0 0 1//2 1//2 0 0; 0 -1//2 -7//2 -3 0 0; -1//2 0 7//2 3 0 0; 0 0 -7//2 -3 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; 1 0 0 0 1 0; 0 0 1//2 0 0 0; 0 -1//2 0 1//2 0 0; -1//2 0 0 -1//2 0 0; 0 0 0 1//2 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 0 1; -1 2 1 0 1 0; 1//2 -1 0 1//2 0 0; -7//2 13//2 0 -3 0 0; 4 -8 -1//2 3 0 0; -7//2 7 0 -3 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1 0; 1 0 0 -1 2 0; 0 0 1 0 0 1; 0 -1//2 0 0 0 0; -3//2 0 -7//2 1 -2 -3; 21//2 0 21 -7 14 18];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[80, 0, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 28]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0 1 0 0; -22 -21 -1 2 0 0; 87 83 2 -4 1 0; 5982 5707 138 -276 0 1; -47943//2 -45739//2 -553 1106 0 -3; 6069//4 2895//2 35 -70 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[8, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1 0 0; -1 1 1 0 0 0; 1 -2 0 0 0 1; -1 3 0 0 -1 2; 3//2 -5 0 0 1 -2; -15//2 51//2 0 0 -5 10];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[24, 0, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 1 0 1 0 0; -8 -7 1 0 0 0; 15 13 0 0 0 1; -22 -19 0 0 -1 2; 73//2 63//2 0 0 1 -2; -7//2 -3 0 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[8, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1 0 0; -1 1 1 0 0 0; 1 -2 0 0 0 1; -1 3 0 0 -1 2; 0 -1//2 0 0 0 0; 3//2 -9//2 0 0 1 -2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[80, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 1 0 0 0; -1 1 0 0 0 0; 1 0 0 0 1 0; -1 0 0 1 0 0; 1 0 0 0 0 1; -1//4 0 0 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=1)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[8, 0, 0, 0, 0, 0, 0, 56, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 4, 8, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 12]);
  L = integer_lattice(B, gram = G);
  B2=QQ[1 0 0 1 0 0; -1 1 -1 2 0 0; 3 -4 2 -4 1 0; 18 -23 12 -24 0 1; -75//2 48 -25 50 0 -1; 21//2 -27//2 7 -14 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 20]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1 0 0; 1 0 1 0 0 0; 0 -1 0 0 1 0; -1 0 0 0 0 1; 0 1//2 0 0 0 0; 5//2 0 0 0 0 -2];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1 1; 1//2 0 0 1//2 0 0; 0 -5//2 1//2 0 -5//2 -2; 0 -1//4 0 0 0 0; -1//4 0 0 0 0 0; 0 5 -5//4 0 5 4];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[64, 32, 0, 0, 0, 0, 32, 64, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1 1; -1//2 1 1//2 1//2 0 0; 1//2 -4 0 0 -5//2 -2; 0 0 -1//4 0 0 0; -1 33//4 0 0 5 4; -25//4 50 0 0 30 24];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 32, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 96, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1 1; 1//2 0 0 1//2 0 0; 0 -3 1//2 0 -5//2 -2; -1//4 0 -1//4 0 0 0; 0 25//4 -1 0 5 4; 0 75//2 -25//4 0 30 24];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//8, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[64, 32, 0, 0, 0, 0, 32, 64, 0, 0, 0, 0, 0, 0, 96, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 40]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 0 1 1; -1//2 1//2 0 1//2 0 0; 1//2 -3 1//2 0 -5//2 -2; -1 23//4 -5//4 0 5 4; 0 0 -1//4 0 0 0; -25//4 35 -15//2 0 30 24];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[16, 8, 0, 0, 0, 0, 8, 16, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1 0 0; -1 1 1 0 0 0; 2 -3 0 0 0 1; -3 5 0 0 1 0; 2 -7//2 0 0 0 0; -1//2 1 0 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  B = matrix(QQ, 6, 6 ,[1//4, 0, 0, 0, 0, 0, 0, 1//4, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2, 0, 0, 0, 0, 0, 0, 1//2]);
  G = matrix(QQ, 6, 6 ,[0, 8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0]);
  L = integer_lattice(B, gram = G);
  B2=QQ[0 1 0 1 0 0; 1 0 1 0 0 0; 0 -1 0 0 0 1; -1 0 0 0 1 0; 0 1//2 0 0 0 0; 1//2 0 0 0 0 0];
  rels = lattice_in_same_ambient_space(L, B2);
  T = torsion_quadratic_module(L, rels; modulus=1, modulus_qf=2)
  push!(DD, T);
  
  for (k,D) in enumerate(DD)
    if k<=120
      @test Oscar._test_stabilizers(D)
    else 
      for i in 0:6
        @test Oscar._test_isotropic_stabilizer_orders(D, i)
      end
    end
  end
end 
