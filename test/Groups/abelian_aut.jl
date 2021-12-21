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
<<<<<<< HEAD
  @test A[1]^(autA[2]*autA[3]) == (A[1]^autA[2])^autA[3]
  @test all(typeof(f) == elem_type(autA) for f in gens(autA))
  @test all(group_element(autA,f.X) == f for f in gens(autA))
  @test all(autA(hom(f)) == f for f in gens(autA))
  @test all(autA(f.map) == f for f in gens(autA))
  @test all(isautomorphism(autA.G,f.map) for f in gens(autA))
  @test all(inv(f) in autA for f in gens(autA))
  @test all(typeof(f*g) == elem_type(autA) for f in gens(autA) for g in gens(autA))
=======
  autA[1](A[1])
  @test A[1]^(autA[2]*autA[3]) == (A[1]^autA[2])^autA[3]
  A[1]^autA[1]
>>>>>>> 5db4abfde96a33f401986a368026321b7239735e

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

end
