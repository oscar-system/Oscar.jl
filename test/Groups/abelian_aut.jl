@testset "conversion" begin
  A = abelian_group([n for n in 1:6])
  Agap,to_gap,to_oscar = oscar._isomorphic_gap_group(A)
  @test all(to_oscar(to_gap(a))==a for a in A)
  @test all(to_gap(to_oscar(a))==a for a in Agap)
  @test all(to_oscar(a*b)==to_oscar(a)+to_oscar(b) for a in gens(Agap) for b in gens(Agap))
  @test all(to_gap(a+b)==to_gap(a)*to_gap(b) for a in gens(A) for b in gens(A))

  Agap,to_gap,to_oscar = oscar._isomorphic_gap_group(A;T=FPGroup)
  @test all(to_oscar(to_gap(a))==a for a in A)
  @test all(to_gap(to_oscar(a))==a for a in Agap)
  @test all(to_oscar(a*b)==to_oscar(a)+to_oscar(b) for a in gens(Agap) for b in gens(Agap))
  @test all(to_gap(a+b)==to_gap(a)*to_gap(b) for a in gens(A) for b in gens(A))

  Agap,to_gap,to_oscar = oscar._isomorphic_gap_group(A;T=PermGroup)
  @test all(to_oscar(to_gap(a))==a for a in A)
  @test all(to_gap(to_oscar(a))==a for a in Agap)
  @test all(to_oscar(a*b)==to_oscar(a)+to_oscar(b) for a in gens(Agap) for b in gens(Agap))
  @test all(to_gap(a+b)==to_gap(a)*to_gap(b) for a in gens(A) for b in gens(A))

  autA = automorphism_group(A)
  @test A[1]^(autA[2]*autA[3]) == (A[1]^autA[2])^autA[3]
  @test all(autA(hom(f)) == f for f in gens(autA))
  @test all(autA(matrix(f)) == f for f in gens(autA))
  @test all(defines_automorphism(domain(autA),matrix(f)) for f in gens(autA))

  A,_ = sub(A,[A[1],A[3],A[3]+A[2],A[2]-A[3]], false)
  Agap,to_gap,to_oscar = oscar._isomorphic_gap_group(A)
  @test all(to_oscar(to_gap(a))==a for a in A)
  @test all(to_gap(to_oscar(a))==a for a in Agap)
  @test all(to_oscar(a*b)==to_oscar(a)+to_oscar(b) for a in gens(Agap) for b in gens(Agap))
  @test all(to_gap(a+b)==to_gap(a)*to_gap(b) for a in gens(A) for b in gens(A))

  Agap,to_gap,to_oscar = oscar._isomorphic_gap_group(A;T=FPGroup)
  @test all(to_oscar(to_gap(a))==a for a in A)
  @test all(to_gap(to_oscar(a))==a for a in Agap)
  @test all(to_oscar(a*b)==to_oscar(a)+to_oscar(b) for a in gens(Agap) for b in gens(Agap))
  @test all(to_gap(a+b)==to_gap(a)*to_gap(b) for a in gens(A) for b in gens(A))

  Agap,to_gap,to_oscar = oscar._isomorphic_gap_group(A;T=PermGroup)
  @test all(to_oscar(to_gap(a))==a for a in A)
  @test all(to_gap(to_oscar(a))==a for a in Agap)
  @test all(to_oscar(a*b)==to_oscar(a)+to_oscar(b) for a in gens(Agap) for b in gens(Agap))
  @test all(to_gap(a+b)==to_gap(a)*to_gap(b) for a in gens(A) for b in gens(A))

  autA = automorphism_group(A)
  @test autA[1](A[1]) ==  oscar.apply_automorphism(autA[1],A[1]) ==  A[1]^autA[1]
  @test all(autA(hom(f)) == f for f in gens(autA))
  @test all(autA(matrix(f)) == f for f in gens(autA))
  @test all(defines_automorphism(domain(autA),matrix(f)) for f in gens(autA))

end

@testset "Orthogonal groups of torsion quadratic modules" begin

  L = Zlattice(gram=3*ZZ[2 1; 1 2])
  D = discriminant_group(L)
  G = orthogonal_group(D)
  d = gens(D)[1]
  f = gens(G)[1]
  f(d)
  @test G(matrix(f)) == f
  @test hom(f)(d) == f(d)
  @test G(hom(f)) == f


  L = Zlattice(gram=3*ZZ[2 1 0; 1 2 0; 0 0 1])
  @test order(orthogonal_group(discriminant_group(L))) == 72

  L = Zlattice(gram=ZZ[0 1; 1 2])
  # the trivial group
  D = discriminant_group(L)
  G = orthogonal_group(D)
  g = one(G)
  @test @inferred g ==G(matrix(g))

  L = root_lattice(:A, 2)
  q = discriminant_group(L)
  T = orthogonal_sum(q, q)[1]
  OT = orthogonal_group(T)
  f = matrix(ZZ, 2, 2, [1 1;0 1])
  fT = hom(T, T, f) # this works, we see it as a map of abelian group
  @test_throws ErrorException OT(fT) # this should not because fT does not preserve the bilinear form
  T = discriminant_group(root_lattice(:D, 13))
  Tsub, _ = sub(T, 4*gens(T))
  @test order(orthogonal_group(Tsub)) == 1
  L = Zlattice(gram=ZZ[1 0 0; 0 0 2; 0 2 0])
  @test order(orthogonal_group(discriminant_group(L)))==6
  # a test for odd lattices
end

@testset "Orthogonal groups of non-semiregular torquadmod" begin
  L = Zlattice(gram=matrix(ZZ, [[2, -1, 0, 0, 0, 0],[-1, 2, -1, -1, 0, 0],[0, -1, 2, 0, 0, 0],[0, -1, 0, 2, 0, 0],[0, 0, 0, 0, 6, 3],[0, 0, 0, 0, 3, 6]]))
  T = discriminant_group(L)
  Tsub, _ = sub(T, [2*T[1], 3*T[2]])
  TT = direct_sum(Tsub, Tsub)[1]
  r3 = radical_quadratic(primary_part(TT, 3)[1])[1]
  TT2 = primary_part(TT, 2)[1]
  @test order(orthogonal_group(Tsub)) == 12
  @test order(orthogonal_group(TT)) == 62208
  @test orthogonal_group(TT) === orthogonal_group(TT)
  @test order(orthogonal_group(TT2)) == 2
  @test order(orthogonal_group(r3)) == 48  # this is the order of GL_2(3)

  T = TorQuadMod(matrix(QQ, 1, 1, [1//27]))
  Tsub, _ = sub(T, 3*gens(T))
  @test order(orthogonal_group(Tsub)) == 6
  T2 = TorQuadMod(matrix(QQ, 1, 1, [21//25]))
  Tsub2, _ = sub(T2, 5*gens(T2))
  @test order(orthogonal_group(Tsub2)) == 4
  TT = direct_sum(Tsub, Tsub2)[1]
  @test order(orthogonal_group(TT)) == 24

  L = direct_sum(L, root_lattice(:A, 6))[1]
  T = discriminant_group(L)
  Tsub, _ = sub(T, [3*T[1], 3*T[2]])
  @assert !is_semi_regular(Tsub)
  @test order(orthogonal_group(Tsub)) == 24 # expected because for A_6, we have order 2
                                            # and the discriminant group is 7-elementary
end

@testset "Embedding of orthogonal groups" begin
  L = Zlattice(gram=matrix(ZZ, 6, 6, [ 2 -1  0  0  0  0;
                                      -1  2 -1 -1  0  0;
                                       0 -1  2  0  0  0;
                                       0 -1  0  2  0  0;
                                       0  0  0  0  6  3;
                                       0  0  0  0  3  6]))
  T = discriminant_group(L)
  i = id_hom(T)
  f = @inferred embedding_orthogonal_group(i)
  @test is_bijective(f)

  _, i = primary_part(T, 3)
  f = @inferred embedding_orthogonal_group(i)
  @test is_injective(f) && !is_surjective(f)
  @test order(domain(f)) == 12
  @test all(g -> order(f(g)) == order(g), domain(f))
end

