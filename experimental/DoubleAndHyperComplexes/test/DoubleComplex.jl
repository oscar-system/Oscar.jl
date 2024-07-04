@testset "double complexes" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  R1 = FreeMod(R, 1)
  Cx = ComplexOfMorphisms([hom(R1, R1, [x*R1[1]])])
  Cy = ComplexOfMorphisms([hom(R1, R1, [y*R1[1]])])
  Cz = ComplexOfMorphisms([hom(R1, R1, [z*R1[1]])])

  Cxy = tensor_product(Cx, Cy)
  Cxy[0, 0]
  Oscar.finalize!(Cxy)
  @test length(Cxy.chains) == 4
  @test length(Cxy.vertical_maps) == 2
  @test length(Cxy.horizontal_maps) == 2

  Cxy = tensor_product(Cx, Cy)
  Cxy[1, 0]
  Oscar.finalize!(Cxy)
  @test length(Cxy.chains) == 4
  @test length(Cxy.vertical_maps) == 2
  @test length(Cxy.horizontal_maps) == 2

  Cxy = tensor_product(Cx, Cy)
  Cxy[0, 1]
  Oscar.finalize!(Cxy)
  @test length(Cxy.chains) == 4
  @test length(Cxy.vertical_maps) == 2
  @test length(Cxy.horizontal_maps) == 2

  Cxy = tensor_product(Cx, Cy)
  Cxy[1, 1]
  Oscar.finalize!(Cxy)
  @test is_complete(Cxy)
  @test length(Cxy.chains) == 4
  @test length(Cxy.vertical_maps) == 2
  @test length(Cxy.horizontal_maps) == 2

  @test Oscar.right_bound(Cxy) == 1
  @test Oscar.left_bound(Cxy) == 0
  @test Oscar.upper_bound(Cxy) == 1
  @test Oscar.lower_bound(Cxy) == 0
  @test can_compute_index(Cxy, 0, 0) && Cxy[0, 0] isa FreeMod
  @test can_compute_index(Cxy, 1, 0) && Cxy[1, 0] isa FreeMod
  @test can_compute_index(Cxy, 0, 1) && Cxy[0, 1] isa FreeMod
  @test can_compute_index(Cxy, 1, 1) && Cxy[1, 1] isa FreeMod
  @test_throws ErrorException Cxy[0, 2]
  @test !can_compute_index(Cxy, 0, 2)
  @test_throws ErrorException Cxy[0, -1]
  @test !can_compute_index(Cxy, 0, -1)
  @test_throws ErrorException Cxy[-1, 0]
  @test !can_compute_index(Cxy, -1, 0)
  @test_throws ErrorException Cxy[2, 0]
  @test !can_compute_index(Cxy, 2, 0)

  @test !can_compute_horizontal_map(Cxy, 1, 100)
  @test can_compute_horizontal_map(Cxy, 1, 0)
  f = Oscar.horizontal_map(Cxy, 1, 0)
  @test domain(f) === Cxy[1, 0]
  @test codomain(f) === Cxy[0, 0]

  @test !can_compute_horizontal_map(Cxy, 1, 100)
  @test can_compute_horizontal_map(Cxy, 1, 1)
  f = Oscar.horizontal_map(Cxy, 1, 1)
  @test domain(f) === Cxy[1, 1]
  @test codomain(f) === Cxy[0, 1]

  @test !can_compute_vertical_map(Cxy, 0, 100)
  @test can_compute_vertical_map(Cxy, 0, 1)
  f = Oscar.vertical_map(Cxy, 0, 1)
  @test domain(f) === Cxy[0, 1]
  @test codomain(f) === Cxy[0, 0]

  @test !can_compute_vertical_map(Cxy, 1, 100)
  @test can_compute_vertical_map(Cxy, 1, 1)
  f = Oscar.vertical_map(Cxy, 1, 1)
  @test domain(f) === Cxy[1, 1]
  @test codomain(f) === Cxy[1, 0]

  tot_xy = Oscar.total_complex(Cxy)

  CC = tensor_product(tot_xy, Cz)
  @test Oscar.horizontal_range(CC) == 2:-1:0
  @test Oscar.vertical_range(CC) == 1:-1:0
  for I in Iterators.product(Oscar.horizontal_range(CC), Oscar.vertical_range(CC))
    @test CC[I] isa FreeMod
  end
  tot_xyz = Oscar.total_complex(CC)

  # Check that we really got the Koszul complex
  @test matrix(map(tot_xyz, 1)) == R[z; y; x]
  @test matrix(map(tot_xyz, 2)) == R[-y z 0; -x 0 z; 0 -x y]
  @test matrix(map(tot_xyz, 3)) == R[x -y z]

  # Test behavior w.r.t zero entries
  Z = FreeMod(R, 0)
  Kx = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [x*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Ky = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [y*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Kz = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [z*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Kxy = tensor_product(Kx, Ky)

  Kxy[0, 0]
  @test !is_complete(Kxy)
  Oscar.finalize!(Kxy)
  @test length(Kxy.chains) == 12
  @test length(Kxy.vertical_maps) == 2
  @test length(Kxy.horizontal_maps) == 2
  @test is_complete(Kxy)

  Kxy = tensor_product(Kx, Ky)
  Kxy[1, 0]
  @test !is_complete(Kxy)
  Oscar.finalize!(Kxy)
  @test length(Kxy.chains) == 12
  @test length(Kxy.vertical_maps) == 2
  @test length(Kxy.horizontal_maps) == 2
  @test is_complete(Kxy)

  Kxy = tensor_product(Kx, Ky)
  Kxy[0, 1]
  @test !is_complete(Kxy)
  Oscar.finalize!(Kxy)
  @test length(Kxy.chains) == 12
  @test length(Kxy.vertical_maps) == 2
  @test length(Kxy.horizontal_maps) == 2
  @test is_complete(Kxy)

  Kxy = tensor_product(Kx, Ky)
  Kxy[1, 1]
  @test !is_complete(Kxy)
  Oscar.finalize!(Kxy)
  @test length(Kxy.chains) == 12
  @test length(Kxy.vertical_maps) == 2
  @test length(Kxy.horizontal_maps) == 2
  @test is_complete(Kxy)

  tot_xy = Oscar.total_complex(Kxy)

  KK = tensor_product(tot_xy, Kz)
  @test !is_complete(KK)
  @test_throws ErrorException Oscar.finalize!(KK)
  KK[0, 0]
  @test !is_complete(KK)
  Oscar.finalize!(KK)
  @test is_complete(KK)
  tot_xyz_alt = Oscar.total_complex(KK)

  @test matrix(map(tot_xyz_alt, 1)) == matrix(map(tot_xyz, 1))
  @test matrix(map(tot_xyz_alt, 2)) == matrix(map(tot_xyz, 2))
  @test matrix(map(tot_xyz_alt, 3)) == matrix(map(tot_xyz, 3))
end

@testset "the same for quotient rings" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  R, _ = quo(R, ideal(R, [x + y + z]))
  R1 = FreeMod(R, 1)
  Cx = ComplexOfMorphisms([hom(R1, R1, [x*R1[1]])])
  Cy = ComplexOfMorphisms([hom(R1, R1, [y*R1[1]])])
  Cz = ComplexOfMorphisms([hom(R1, R1, [z*R1[1]])])

  Cxy = tensor_product(Cx, Cy)

  @test Oscar.right_bound(Cxy) == 1
  @test Oscar.left_bound(Cxy) == 0
  @test Oscar.upper_bound(Cxy) == 1
  @test Oscar.lower_bound(Cxy) == 0
  @test can_compute_index(Cxy, 0, 0) && Cxy[0, 0] isa FreeMod
  @test can_compute_index(Cxy, 1, 0) && Cxy[1, 0] isa FreeMod
  @test can_compute_index(Cxy, 0, 1) && Cxy[0, 1] isa FreeMod
  @test can_compute_index(Cxy, 1, 1) && Cxy[1, 1] isa FreeMod
  @test_throws ErrorException Cxy[0, 2]
  @test !can_compute_index(Cxy, 0, 2)
  @test_throws ErrorException Cxy[0, -1]
  @test !can_compute_index(Cxy, 0, -1)
  @test_throws ErrorException Cxy[-1, 0]
  @test !can_compute_index(Cxy, -1, 0)
  @test_throws ErrorException Cxy[2, 0]
  @test !can_compute_index(Cxy, 2, 0)

  f = Oscar.horizontal_map(Cxy, 1, 0)
  @test domain(f) === Cxy[1, 0]
  @test codomain(f) === Cxy[0, 0]

  f = Oscar.horizontal_map(Cxy, 1, 1)
  @test domain(f) === Cxy[1, 1]
  @test codomain(f) === Cxy[0, 1]

  f = Oscar.vertical_map(Cxy, 0, 1)
  @test domain(f) === Cxy[0, 1]
  @test codomain(f) === Cxy[0, 0]

  f = Oscar.vertical_map(Cxy, 1, 1)
  @test domain(f) === Cxy[1, 1]
  @test codomain(f) === Cxy[1, 0]

  tot_xy = Oscar.total_complex(Cxy)

  CC = tensor_product(tot_xy, Cz)
  @test Oscar.horizontal_range(CC) == 2:-1:0
  @test Oscar.vertical_range(CC) == 1:-1:0
  for I in Iterators.product(Oscar.horizontal_range(CC), Oscar.vertical_range(CC))
    @test CC[I] isa FreeMod
  end
  tot_xyz = Oscar.total_complex(CC)

  # Check that we really got the Koszul complex
  @test matrix(map(tot_xyz, 1)) == R[z; y; x]
  @test matrix(map(tot_xyz, 2)) == R[-y z 0; -x 0 z; 0 -x y]
  @test matrix(map(tot_xyz, 3)) == R[x -y z]

  # Test behavior w.r.t zero entries
  Z = FreeMod(R, 0)
  Kx = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [x*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Ky = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [y*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Kz = ComplexOfMorphisms([hom(Z, R1, elem_type(R1)[]), hom(R1, R1, [z*R1[1]]), hom(R1, Z, [zero(Z)])], seed = -1)
  Kxy = tensor_product(Kx, Ky)

  tot_xy = Oscar.total_complex(Kxy)

  KK = tensor_product(tot_xy, Kz)
  tot_xyz_alt = Oscar.total_complex(KK)

  @test matrix(map(tot_xyz_alt, 1)) == matrix(map(tot_xyz, 1))
  @test matrix(map(tot_xyz_alt, 2)) == matrix(map(tot_xyz, 2))
  @test matrix(map(tot_xyz_alt, 3)) == matrix(map(tot_xyz, 3))
end

@testset "extended functionality for tensor products" begin
  R, (x, y) = QQ[:x, :y]
  F = FreeMod(R, 1)
  phi1 = hom(F, F, [x*F[1]])

  G = FreeMod(R, 1)
  phi2 = hom(G, G, [y*G[1]])

  @test matrix(Oscar.tensor_product([phi1, phi2]))[1, 1] == x*y

  R4 = FreeMod(R, 4)
  R2 = FreeMod(R, 2)
  phi = hom(R2, R4, R[x y x y; 2*x y 3*x 7])

  R2x2 = tensor_product(R2, R2)
  R4x4 = tensor_product(R4, R4)

  psi = tensor_product(R2x2, R4x4, [phi, phi])
  @test domain(psi) === R2x2
  @test codomain(psi) === R4x4
  psi_alt = tensor_product(phi, phi)
  @test domain(psi_alt) !== R2x2
  @test codomain(psi_alt) !== R4x4

  @test matrix(psi) == matrix(tensor_product(phi, phi))

  R2x2x2 = tensor_product(R2, R2, R2)
  R4x4x4 = tensor_product(R4, R4, R4)

  psi = tensor_product(R2x2x2, R4x4x4, [phi, phi, phi])
end

@testset "the same for quotient rings again" begin
  R, (x, y) = QQ[:x, :y]
  R, _ = quo(R, ideal(R, [x^2 + y^2 - 1]))
  F = FreeMod(R, 1)
  phi1 = hom(F, F, [x*F[1]])

  G = FreeMod(R, 1)
  phi2 = hom(G, G, [y*G[1]])

  @test matrix(Oscar.tensor_product([phi1, phi2]))[1, 1] == x*y

  R4 = FreeMod(R, 4)
  R2 = FreeMod(R, 2)
  phi = hom(R2, R4, R[x y x y; 2*x y 3*x 7])

  R2x2 = tensor_product(R2, R2)
  R4x4 = tensor_product(R4, R4)

  psi = tensor_product(R2x2, R4x4, [phi, phi])
  @test domain(psi) === R2x2
  @test codomain(psi) === R4x4
  psi_alt = tensor_product(phi, phi)
  @test domain(psi_alt) !== R2x2
  @test codomain(psi_alt) !== R4x4

  @test matrix(psi) == matrix(tensor_product(phi, phi))

  R2x2x2 = tensor_product(R2, R2, R2)
  R4x4x4 = tensor_product(R4, R4, R4)

  psi = tensor_product(R2x2x2, R4x4x4, [phi, phi, phi])
end
