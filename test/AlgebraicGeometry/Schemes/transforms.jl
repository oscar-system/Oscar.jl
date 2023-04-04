@testset "different notions of transform" begin
  W = affine_space(QQ,3)
  R = OO(W)
  (x,y,z) = gens(R)

  I = ideal(R,[x*z,y*z,x^3+y^3])
  X = subscheme(W,I)

  C = ideal(R,[x,y,z])
  bl = blow_up(W,C)
  IE = ideal_sheaf(exceptional_divisor(bl))
  Isheaf = ideal_sheaf(codomain(bl),W,gens(I))

  Istrict = strict_transform(bl,Isheaf)
  @test is_smooth(subscheme(Istrict))
  
  Iweak = weak_transform(bl,Isheaf)
  Iweak2,mult = Oscar.weak_transform_with_multiplicity(bl, Isheaf)
  @test Iweak == Iweak2
  @test mult == 2

  Ictrl1 = controlled_transform(bl,Isheaf,1)
  @test IE*Istrict != Iweak
  @test IE*Iweak == Ictrl1

end