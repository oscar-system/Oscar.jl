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

  Isheaf = IdealSheaf(W, I; covered_scheme=codomain(bl))
  C2sheaf = IdealSheaf(W, C2; covered_scheme=codomain(bl))

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
  bl2 = blow_up(C2strict,covering=Oscar.simplified_covering(W2),var_name="t")

  E2 = strict_transform(bl2,E)
end

@testset "associated_points" begin
  # set up standard P2
  S, _ = graded_polynomial_ring(QQ,[:x,:y,:z])
  P2 = proj(S)
  X = covered_scheme(P2)

  C = standard_covering(P2)
  Ux = patches(C)[1]
  Uy = patches(C)[2]
  Uz = patches(C)[3]

  # define ideal on patches: one line, one isolated point, one embedded point
  Rx = OO(Ux)
  (y,z) = gens(Rx)
  Ix = ideal(Rx,[y^2,y*z])

  Ry = OO(Uy)
  (x,z) = gens(Ry)
  Iy = ideal(Ry,[x,z])

  Rz = OO(Uz)
  (x,y) = gens(Rz)
  Iz = ideal(Rz,[y])

  # glue together to an ideal sheaf on P2
  ID = IdDict{AbsAffineScheme, Oscar.Ideal}()
  ID[Ux]=Ix
  ID[Uy]=Iy
  ID[Uz]=Iz
  ISheaf=IdealSheaf(X, ID, check=true)

  # decompose and check that we obtain the expected results
  L1 = Oscar.associated_points(ISheaf)
  L2 = Oscar.maximal_associated_points(ISheaf)
  dimsL1 = [dim(a) for a in L1]
  dimsL2 = [dim(a) for a in L2]
  @test length(L1) == 3
  @test length(L2) == 2
  onesinL1 = findall(is_one, dimsL1)
  zerosinL1 = findall(is_zero, dimsL1)
  onesinL2 = findall(is_one, dimsL2)
  zerosinL2 = findall(is_zero, dimsL2)
  @test length(onesinL1) == 1
  @test length(zerosinL1) == 2
  @test length(zerosinL2) == 1
  @test length(onesinL2) == 1
end
