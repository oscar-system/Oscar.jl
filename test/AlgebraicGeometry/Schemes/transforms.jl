@testset "different notions of transform" begin
  W = affine_space(QQ,3)
  R = OO(W)
  (x,y,z) = gens(R)

  I = ideal(R,[x*z,y*z,x^3+y^3])
  X = subscheme(W,I)

  C = ideal(R,[x,y,z])
  C2 = ideal(R,[x,y])
  bl = blow_up(W,C)
  E = exceptional_divisor(bl)
  IE = ideal_sheaf(E)
  W2 = domain(bl)

  Isheaf = ideal_sheaf(codomain(bl),W,gens(I))
  C2sheaf = ideal_sheaf(codomain(bl),W,gens(C2))

  Istrict = strict_transform(bl,Isheaf)
  @test is_smooth(subscheme(Istrict))
  
  Iweak = weak_transform(bl,Isheaf)
  Iweak2,mult = Oscar.weak_transform_with_multiplicity(bl, Isheaf)
  @test Iweak == Iweak2
  @test mult == 2

  Ictrl1 = controlled_transform(bl,Isheaf,1)
  @test IE*Istrict != Iweak
  @test IE*Iweak == Ictrl1

  C2strict = strict_transform(bl,C2sheaf)
  bl2 = blow_up(C2strict,covering=oscar.simplified_covering(W2),var_name="t")

  E2 = strict_transform(bl2,E)
end
