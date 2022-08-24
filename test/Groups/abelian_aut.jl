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

  A,_ = sub(A,[A[1],A[3],A[3]+A[2],A[2]-A[3]])
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
end
