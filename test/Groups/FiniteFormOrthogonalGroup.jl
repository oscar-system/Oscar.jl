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
  
@test "stabilisers of subspaces under isometries" begin
  function _test(L::ZZLat, p)
    DL = discriminant_group(L)
    D,iD = kernel(hom(DL,DL, [p*i for i in gens(DL)]))
    
    G,iG = restrict_automorphism_group(image_in_Oq(L)[1],iD);
    OG = orthogonal_group(D);
    _,iG = is_subgroup(G, orthogonal_group(D));
    
    for k in 1:div(rank(L),2)
      @test Oscar._test_isotropic_stabilizer_orders(D, k, p)
      X = [i for (i,_) in subgroups(D; order=p^k) if is_totally_isotropic(i)]
      a = length(Oscar._isotropic_subspaces_representatives(D, iG, k))
      if length(X) > 0
        XG = gset(G, Oscar.on_subgroups_slow, X)
        b = length(orbits(XG))
      else 
        b = 0
      end
      @test a==b
    end
  end
  for i in 1:6
    _test(rescale(root_lattice(:A,i),2),2)
    _test(rescale(root_lattice(:A,i),3),3)
  end
  for i in 4:6
    _test(rescale(root_lattice(:D,i),2),2)
    _test(rescale(root_lattice(:D,i),3),3)
  end
  for i in 1:4
    _test(rescale(root_lattice(:A,i),7),7)
  end 
end
