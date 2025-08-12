@testset "generic Ext" begin
  # We create several very simple Koszul complexes and their duals
  R, (x, y) = QQ[:x, :y]
  R1 = FreeMod(R, 1)
  Kx = Oscar.SimpleComplexWrapper(koszul_complex(x*R1[1]))
  Ky = Oscar.SimpleComplexWrapper(koszul_complex(y*R1[1]))

  Kx_dual = hom(Kx, R1) # Returns a cochain complex
  @test upper_bound(Kx_dual) == -lower_bound(Kx)
  @test lower_bound(Kx_dual) == -upper_bound(Kx)
  @test Oscar.direction(Kx_dual) == Oscar.direction(Kx)
  p = map(Kx_dual, 0)
  @test domain(p) === Kx_dual[0]
  @test codomain(p) === Kx_dual[-1]
  @test matrix(map(Kx, 1)) == matrix(map(Kx_dual, 0)) == R[x; ]

  KKx = hom(R1, Kx) # Returns a chain complex
  @test upper_bound(KKx) == upper_bound(Kx)
  @test lower_bound(KKx) == lower_bound(Kx)
  @test Oscar.direction(KKx) == Oscar.direction(Kx)
  @test matrix(map(Kx, 1)) == matrix(map(KKx, 1)) == R[x; ]

  Kxy = hom(Kx, Ky) # Returns a 2-dimensional hyper complex which 
  @test upper_bound(Kxy, 1) == -lower_bound(Kx)
  @test lower_bound(Kxy, 1) == -upper_bound(Kx)
  @test upper_bound(Kxy, 2) == upper_bound(Ky)
  @test lower_bound(Kxy, 2) == lower_bound(Ky)

  @test Oscar.direction(Kxy, 1) == :chain
  @test Oscar.direction(Kxy, 2) == :chain
  @test matrix(map(Kxy, 1, (0, 0))) == R[x;]
  @test matrix(map(Kxy, 2, (0, 1))) == R[y;]

  Kz = Oscar.SimpleComplexWrapper(koszul_complex((x+y)*R1[1]))
  C1 = hom(Kxy, Kz) # Returns a 3-dimensional hypercomplex 
  @test dim(C1) == 3
  @test upper_bound(C1, 1) == upper_bound(Kx)
  @test lower_bound(C1, 1) == lower_bound(Kx)
  @test upper_bound(C1, 2) == -lower_bound(Ky)
  @test lower_bound(C1, 2) == -upper_bound(Ky)
  @test upper_bound(C1, 3) == upper_bound(Kz)
  @test lower_bound(C1, 3) == lower_bound(Kz)

  @test Oscar.direction(C1, 1) == :chain
  @test Oscar.direction(C1, 2) == :chain
  @test Oscar.direction(C1, 3) == :chain
  @test isone(rank(C1[(0, 0, 0)]))
  @test isone(rank(C1[(1, 0, 0)]))
  p = map(C1, 1, (1, 0, 0))
  @test domain(p) === C1[1, 0, 0]
  @test codomain(p) === C1[0, 0, 0]
  @test matrix(p) == R[x;]
  q = map(C1, 2, (0, 0, 0))
  @test domain(q) === C1[(0,0,0)]
  @test codomain(q) === C1[(0, -1, 0)]
  @test matrix(q) == R[y;]
  r = map(C1, 3, (0, 0, 1))
  @test domain(r) === C1[(0,0,1)]
  @test codomain(r) === C1[(0, 0, 0)]
  @test matrix(r) == R[x+y;]

  C2 = hom(Kz, Kxy) # Returns a 3-dimensional hypercomplex
                    # reversing the direction in the first 
                    # dimension.
  @test dim(C2) == 3
  @test upper_bound(C2, 1) == -lower_bound(Kz)
  @test lower_bound(C2, 1) == -upper_bound(Kz)
  @test upper_bound(C2, 2) == -lower_bound(Kx)
  @test lower_bound(C2, 2) == -upper_bound(Kx)
  @test upper_bound(C2, 3) == upper_bound(Ky)
  @test lower_bound(C2, 3) == lower_bound(Ky)

  @test Oscar.direction(C2, 1) == :chain
  @test Oscar.direction(C2, 2) == :chain
  @test Oscar.direction(C2, 3) == :chain
  @test isone(rank(C2[(0, 0, 0)]))
  @test isone(rank(C2[(-1, 0, 0)]))
  p = map(C2, 1, (0, 0, 0))
  @test domain(p) === C2[(0, 0, 0)]
  @test codomain(p) === C2[(-1, 0, 0)]
  @test matrix(p) == R[x+y;]
  q = map(C2, 2, (0, 0, 0))
  @test domain(q) === C2[(0,0,0)]
  @test codomain(q) === C2[(0, -1, 0)]
  @test matrix(q) == R[x;]
  r = map(C2, 3, (0, 0, 1))
  @test domain(r) === C2[(0,0,1)]
  @test codomain(r) === C2[(0, 0, 0)]
  @test matrix(r) == R[y;]
  @test !is_complete(C2)
end
  

@testset "interpretation map for ext elements" begin
  S, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])

  # Hom of Koszul complexes
  R3 = FreeMod(S, 3)
  v = x*R3[1] + y*R3[2] + z*R3[3]
  K = koszul_complex(Oscar.KoszulComplex, v)
  sh = 0
  C = hom(shift(K, sh), shift(K, sh))
  @test C isa Oscar.HomComplex
  H_tot = total_complex(C);
  @test H_tot isa Oscar.TotalComplex{<:ModuleFP, <:ModuleFPHom, <:Oscar.HomComplex}
  for i in -10:10
    Oscar.can_compute_index(H_tot, i) || continue
    H0, _ = kernel(H_tot, i)
    iszero(H0) && continue
    v = first(gens(H0))
    interp = Oscar.element_to_homomorphism_map(H_tot, i)
    phi = interp(v)
    Oscar.can_compute_index(H_tot, i) && H_tot[i]
    Oscar.can_compute_map(H_tot, i) && map(H_tot, i)
    for j in -10:10
      if Oscar.can_compute_index(phi, j-1) && Oscar.can_compute_index(phi, j)
        @test matrix(compose(map(domain(phi), j), phi[j-1])) == matrix(compose(phi[j], map(codomain(phi), j + i)))
      end
    end
  end

  # Hom of higher complexes
  S1 = FreeMod(S, 1)
  vx = x*S1[1]
  Kx = koszul_complex(Oscar.KoszulComplex, vx)
  vy = y*S1[1]
  Ky = koszul_complex(Oscar.KoszulComplex, vy)
  vz = z*S1[1]
  Kz = koszul_complex(Oscar.KoszulComplex, vz)

  # Produce the homomorphisms for 
  #         ⋅x
  #   0 ← R ← R ← 0
  #       ↓   ↓
  #   0 ← R ← R ← 0
  #         ⋅y
  C = hom(Kx, Ky)
  tot = total_complex(C)
  H0, _ = kernel(tot, 0)
  v = first(gens(H0))
  interp = Oscar.element_to_homomorphism_map(tot, 0)
  phi = interp(v)
  @test Oscar.original_complex(domain(phi)) === Kx
  @test Oscar.original_complex(codomain(phi)) === Ky
  @test domain(phi[0]) === domain(phi)[0]
  for i in -10:10
      Oscar.can_compute_index(tot, i) && tot[i]
      Oscar.can_compute_map(tot, i) && map(tot, i)
  end
  @test compose(phi[1], map(codomain(phi), 1, (1,))) == compose(map(domain(phi), 1, (1,)), phi[0])
  @test !iszero(matrix(phi[0]))
  @test !iszero((phi[1]))

  # Produce the homomorphisms for 
  #             ⋅x
  #       0 ← R ← R ← 0
  #       ↓   ↓   ↓
  #   0 ← R ← R ← 0
  #         ⋅y
  H1, _ = kernel(tot, 1)
  @test iszero(H1)
  interp = Oscar.element_to_homomorphism_map(tot, 1)
  v = zero(H1)
  phi = interp(v)
  @test iszero(phi[0])

  # Produce the homomorphisms for 
  #         ⋅x
  #   0 ← R ← R ← 0
  #       ↓   ↓   ↓
  #       0 ← R ← R ← 0
  #             ⋅y
  Hm1, _ = kernel(tot, -1)
  @test iszero(H1)
  interp = Oscar.element_to_homomorphism_map(tot, -1)
  v = first(gens(Hm1))
  phi = interp(v)
  @test isone(matrix(phi[1]))

  # Produce something like a convoluted koszul complex
  Kxy = hom(Kx, Ky)
  Kxyz = hom(Kxy, Kz)

  tot = total_complex(Kxyz)

  # The morphisms
  #           -y    (x y)
  #            x
  #   0  ←  R  ←  R²  ←  R  ←  0
  #         ↓     ↓      ↓     ↓
  #         0  ←  R   ←  R  ←  0
  #                   z
  H0, _ = kernel(tot, 0)
  interp = Oscar.element_to_homomorphism_map(tot, 0)

  # The particular results below depend on the outcome of the kernel computation 
  # and might change!
  phi = interp(H0[1])
  dom_tot = domain(phi)
  cod_tot = codomain(phi)
  @test !iszero(phi[0])
  @test !iszero(phi[1])
  @test compose(map(dom_tot, 1), phi[0]) == compose(phi[1], map(cod_tot, 1))

  phi = interp(H0[2])
  dom_tot = domain(phi)
  cod_tot = codomain(phi)
  @test !iszero(phi[0])
  @test compose(map(dom_tot, 1), phi[0]) == compose(phi[1], map(cod_tot, 1))

  phi = interp(H0[3])
  dom_tot = domain(phi)
  cod_tot = codomain(phi)
  @test !iszero(phi[0])
  @test !iszero(phi[1])
  @test compose(map(dom_tot, 1), phi[0]) == compose(phi[1], map(cod_tot, 1))

  # The morphisms
  #           -y    (x y)
  #            x
  #   0  ←  R  ←  R² ←  R  ←  0
  #            ↘     ↘
  #         0  ←  R  ←  R  ←  0
  #                   z
  H1, _ = kernel(tot, 1)
  interp = Oscar.element_to_homomorphism_map(tot, 1)

  for v in gens(H1)
    phi = interp(v)
    dom_tot = domain(phi)
    cod_tot = codomain(phi)
    @test compose(map(dom_tot, 0), phi[-1]) == compose(phi[0], map(cod_tot, 1))
  end

  # Even more complicated hyper complexes
  for H in [hom(Kxyz, Kxyz), hom(Kxyz, Kz), hom(Kz, Kxyz), hom(C, Kxyz), hom(Kxyz, C), hom(C, C)]
    H = hom(Kxyz, Kxyz)
    H_tot = total_complex(H)

    for d in -3:3
      interp = Oscar.element_to_homomorphism_map(H_tot, d)

      Hd, _ = kernel(H_tot, d)
      v = iszero(ngens(Hd)) ? zero(Hd) : first(gens(Hd))

      phi_v = interp(v)

      @test Oscar.original_complex(domain(phi_v)) === Kxyz
      @test Oscar.original_complex(codomain(phi_v)) === Kxyz

      dom = domain(phi_v)
      cod = codomain(phi_v)

      for i in -10:10
        Oscar.can_compute_index(phi_v, (i,)) || continue
        Oscar.can_compute_index(phi_v, (i-1,)) || continue
        @test matrix(compose(map(dom, i), phi_v[i-1])) == matrix(compose(phi_v[i], map(cod, i+d)))
      end
    end
  end
end

@testset "induced morphisms on hom-complexes" begin
  R, (x, y, z) = QQ[:x, :y, :z]

  F = FreeMod(R, 3)
  v = x^2 * F[1] + y^2 * F[2] + z^2 * F[3]
  K = koszul_complex(Oscar.KoszulComplex, v)

  R1 = FreeMod(R, 1)
  K_dual = Oscar.hom(K, R1)

  R1_comp = codomain(K_dual)

  R2 = FreeMod(R, 2)
  w = (x^3 + y^3) * R2[1] + z^5 * R2[2]

  K2 = koszul_complex(Oscar.KoszulComplex, w)

  K2_dual = hom(K2, codomain(K_dual))

  phi = hom(K2[0], K[0], [K[0][1]]) # The identity map on a free module of rank 1

  f = Oscar.lift_map(K2, K, phi);

  ind_dom = Oscar.HomComplexMorphism(Any, K_dual, K2_dual, f, nothing)
  @test !Oscar.can_compute_index(ind_dom, (1,))
  @test !Oscar.can_compute_index(ind_dom, (-4,))
  @test compose(ind_dom[0], map(K2_dual, 1, (0,))) == compose(map(K_dual, 1, (0,)), ind_dom[(-1)])
  @test compose(ind_dom[-1], map(K2_dual, 1, (-1,))) == compose(map(K_dual, 1, (-1,)), ind_dom[(-2)])
  @test compose(ind_dom[-2], map(K2_dual, 1, (-2,))) == compose(map(K_dual, 1, (-2,)), ind_dom[(-3)])

  KK = hom(R1_comp, K)
  KK2 = hom(R1_comp, K2)

  ind_cod = Oscar.HomComplexMorphism(Any, KK2, KK, nothing, f)
  @test !Oscar.can_compute_index(ind_cod, (-1,))
  @test !Oscar.can_compute_index(ind_cod, (4,))
  @test compose(ind_cod[1], map(KK, 1, (1,))) == compose(map(KK2, 1, (1,)), ind_cod[(0)])
  @test compose(ind_cod[2], map(KK, 1, (2,))) == compose(map(KK2, 1, (2,)), ind_cod[(1)])
  @test compose(ind_cod[3], map(KK, 1, (3,))) == compose(map(KK2, 1, (3,)), ind_cod[(2)])


  K2_shift = shift(K2, -3)

  phi = hom(K2_shift[3], K[0], [K[0][1]]) # The identity map on a free module of rank 1

  g = Oscar.lift_map(K2_shift, K, phi, start_index=3, offset=-3);

  K2_shift_dual = hom(K2_shift, codomain(K_dual))

  ind_dom_shift = Oscar.HomComplexMorphism(Any, K_dual, K2_shift_dual, g, nothing)
  @test !Oscar.can_compute_index(ind_dom_shift, (1,))
  @test !Oscar.can_compute_index(ind_dom_shift, (-4,))
  @test compose(map(K_dual, 1, (0,)), ind_dom_shift[-1]) == compose(ind_dom_shift[0], map(K2_shift_dual, 1, (-3,)))

  K_shift = shift(K, -7)

  phi = hom(K2[0], K_shift[7], [K_shift[7][1]]) # The identity map on a free module of rank 1

  g = Oscar.lift_map(K2, K_shift, phi, start_index=0, offset=7);

  K_shift_dual = hom(K_shift, codomain(K2_dual))

  ind_dom_shift = Oscar.HomComplexMorphism(Any, K_shift_dual, K2_dual, g, nothing)

  @test !Oscar.can_compute_index(ind_dom_shift, (-6,))
  @test !Oscar.can_compute_index(ind_dom_shift, (-12,))
  @test compose(map(K_shift_dual, 1, (-7,)), ind_dom_shift[-8]) == compose(ind_dom_shift[-7], map(K2_dual, 1, (0,)))


  KK = hom(R1_comp, K_dual)
  KK2 = hom(R1_comp, K2_dual)

  ind_cod = Oscar.HomComplexMorphism(Any, KK, KK2, nothing, ind_dom)
  @test !Oscar.can_compute_index(ind_dom, (1,))
  @test !Oscar.can_compute_index(ind_dom, (-4,))
  @test compose(ind_cod[0], map(KK2, 1, (0,))) == compose(map(KK, 1, (0,)), ind_cod[(-1)])
  @test compose(ind_cod[-1], map(KK2, 1, (-1,))) == compose(map(KK, 1, (-1,)), ind_cod[(-2)])
  @test compose(ind_cod[-2], map(KK2, 1, (-2,))) == compose(map(KK, 1, (-2,)), ind_cod[(-3)])

  C1 = hom(K, K2)
  C2 = hom(K2, K)

  double_ind = Oscar.HomComplexMorphism(Any, C1, C2, f, f)
  @test !Oscar.can_compute_index(double_ind, (1, 0))
  @test !Oscar.can_compute_index(double_ind, (0, -1))
  @test double_ind[0, 0] isa ModuleFPHom
  @test double_ind[-1, 0] isa ModuleFPHom
  @test double_ind[-2, 0] isa ModuleFPHom
  @test double_ind[0, 1] isa ModuleFPHom
  @test double_ind[-1, 1] isa ModuleFPHom
  @test double_ind[-2, 1] isa ModuleFPHom

end

@testset "canonical bundle" begin
  R, z = polynomial_ring(GF(3), :z => 1:10)
  S, z = grade(R)
  V = [z[1]*z[7] + z[1]*z[9] + 2*z[1]*z[10] + z[2]*z[7] + 2*z[2]*z[8] + 2*z[2]*z[9] + z[2]*z[10] + z[3]*z[7] + 2*z[3]*z[8] + 2*z[3]*z[10] + 2*z[4]*z[7] + z[5]*z[6] + z[5]*z[7] + 2*z[5]*z[8] + z[5]*z[10] + 2*z[6]^2 + 2*z[6]*z[7] + z[7]^2 + 2*z[7]*z[8] + z[7]*z[9] + 2*z[7]*z[10] + z[8]^2 + 2*z[8]*z[9] + 2*z[8]*z[10] + z[9]*z[10] + z[10]^2, z[1]*z[7] + z[1]*z[8] + z[1]*z[9] + z[2]*z[6] + z[2]*z[7] + 2*z[2]*z[8] + z[2]*z[9] + z[2]*z[10] + 2*z[3]*z[6] + 2*z[3]*z[7] + z[3]*z[8] + z[3]*z[9] + 2*z[4]*z[6] + z[4]*z[7] + 2*z[4]*z[9] + 2*z[4]*z[10] + z[5]*z[9] + z[5]*z[10] + z[6]^2 + z[6]*z[7] + z[6]*z[9] + z[7]^2 + 2*z[7]*z[8] + z[7]*z[9] + 2*z[7]*z[10] + 2*z[8]^2 + z[8]*z[9] + 2*z[8]*z[10] + z[10]^2, z[1]*z[6] + z[1]*z[8] + z[1]*z[9] + z[1]*z[10] + 2*z[2]*z[7] + 2*z[2]*z[8] + 2*z[2]*z[9] + 2*z[3]*z[7] + 2*z[3]*z[10] + 2*z[4]*z[6] + z[4]*z[7] + 2*z[4]*z[9] + z[4]*z[10] + z[5]*z[7] + 2*z[5]*z[9] + 2*z[5]*z[10] + 2*z[6]^2 + z[6]*z[7] + z[6]*z[8] + 2*z[6]*z[9] + 2*z[7]^2 + z[7]*z[8] + z[7]*z[10] + z[8]*z[9] + z[8]*z[10] + z[9]^2 + 2*z[9]*z[10] + 2*z[10]^2, z[1]*z[7] + 2*z[1]*z[8] + z[1]*z[9] + z[1]*z[10] + z[2]*z[8] + z[2]*z[10] + 2*z[3]*z[8] + 2*z[3]*z[10] + z[4]*z[7] + 2*z[4]*z[10] + z[5]^2 + 2*z[5]*z[7] + 2*z[5]*z[8] + 2*z[5]*z[10] + 2*z[6]^2 + z[6]*z[8] + 2*z[6]*z[9] + 2*z[6]*z[10] + z[7]*z[8] + 2*z[7]*z[9] + 2*z[7]*z[10] + z[8]*z[9] + 2*z[9]^2 + 2*z[9]*z[10] + 2*z[10]^2, z[1]*z[8] + 2*z[2]*z[7] + 2*z[2]*z[8] + 2*z[2]*z[9] + 2*z[2]*z[10] + z[3]*z[7] + z[3]*z[8] + 2*z[3]*z[9] + z[4]*z[5] + 2*z[4]*z[6] + 2*z[4]*z[10] + 2*z[5]*z[7] + z[5]*z[9] + z[5]*z[10] + 2*z[6]^2 + z[6]*z[7] + z[6]*z[10] + 2*z[7]^2 + 2*z[7]*z[8] + z[7]*z[9] + 2*z[7]*z[10] + z[8]^2 + 2*z[8]*z[9] + 2*z[8]*z[10] + z[9]^2, z[1]*z[7] + z[1]*z[9] + 2*z[1]*z[10] + z[2]*z[8] + z[2]*z[9] + z[2]*z[10] + z[3]*z[5] + z[3]*z[6] + 2*z[3]*z[7] + z[3]*z[9] + 2*z[4]*z[8] + 2*z[4]*z[9] + z[5]*z[7] + 2*z[5]*z[8] + 2*z[5]*z[9] + 2*z[5]*z[10] + 2*z[6]^2 + z[6]*z[8] + 2*z[6]*z[9] + 2*z[6]*z[10] + 2*z[7]*z[8] + z[7]*z[9] + z[7]*z[10] + z[9]*z[10] + z[10]^2, z[1]*z[9] + z[1]*z[10] + z[2]*z[5] + z[2]*z[7] + z[2]*z[8] + 2*z[2]*z[10] + z[3]*z[6] + 2*z[3]*z[7] + z[3]*z[8] + z[4]*z[9] + z[5]*z[7] + 2*z[6]^2 + 2*z[6]*z[8] + z[6]*z[10] + z[7]^2 + 2*z[7]*z[8] + 2*z[7]*z[9] + z[7]*z[10] + 2*z[8]*z[10] + z[9]*z[10], z[1]*z[5] + 2*z[1]*z[10] + 2*z[2]*z[7] + 2*z[2]*z[8] + z[3]*z[8] + 2*z[4]*z[7] + 2*z[4]*z[8] + 2*z[4]*z[9] + z[4]*z[10] + z[5]*z[9] + 2*z[5]*z[10] + z[6]*z[7] + z[6]*z[9] + z[6]*z[10] + 2*z[7]^2 + z[7]*z[8] + z[7]*z[9] + z[7]*z[10] + 2*z[8]^2 + 2*z[8]*z[10] + 2*z[9]^2 + 2*z[9]*z[10] + 2*z[10]^2, 2*z[1]*z[8] + 2*z[1]*z[9] + 2*z[1]*z[10] + z[2]*z[7] + 2*z[2]*z[8] + 2*z[2]*z[9] + 2*z[2]*z[10] + z[3]*z[8] + 2*z[3]*z[9] + z[4]^2 + z[4]*z[6] + 2*z[4]*z[7] + z[4]*z[10] + z[5]*z[7] + 2*z[5]*z[8] + 2*z[5]*z[9] + 2*z[6]^2 + 2*z[6]*z[7] + z[6]*z[8] + 2*z[6]*z[9] + 2*z[6]*z[10] + z[7]^2 + 2*z[7]*z[8] + z[7]*z[9] + 2*z[7]*z[10] + z[8]^2 + z[8]*z[9] + 2*z[9]*z[10] + z[10]^2, z[1]*z[7] + z[1]*z[8] + z[2]*z[7] + 2*z[2]*z[8] + z[3]*z[4] + 2*z[3]*z[6] + 2*z[3]*z[7] + z[3]*z[8] + z[3]*z[9] + z[4]*z[6] + z[4]*z[7] + 2*z[4]*z[8] + z[4]*z[10] + z[5]*z[7] + 2*z[5]*z[8] + 2*z[5]*z[9] + z[5]*z[10] + 2*z[6]^2 + z[6]*z[7] + 2*z[6]*z[8] + 2*z[6]*z[9] + z[7]^2 + 2*z[7]*z[8] + z[7]*z[10] + z[8]^2 + z[8]*z[10] + 2*z[9]^2, 2*z[1]*z[8] + 2*z[1]*z[9] + z[1]*z[10] + z[2]*z[4] + z[2]*z[7] + z[2]*z[8] + z[2]*z[9] + 2*z[3]*z[6] + z[3]*z[8] + z[3]*z[10] + z[4]*z[6] + z[4]*z[7] + 2*z[4]*z[9] + z[4]*z[10] + z[5]*z[7] + z[5]*z[9] + z[6]^2 + 2*z[6]*z[7] + z[7]^2 + z[7]*z[10] + 2*z[8]^2 + 2*z[8]*z[9] + 2*z[8]*z[10] + 2*z[9]^2 + z[9]*z[10], z[1]*z[4] + z[1]*z[8] + z[2]*z[7] + z[2]*z[8] + z[2]*z[10] + z[3]*z[7] + 2*z[3]*z[8] + 2*z[3]*z[10] + 2*z[4]*z[6] + z[4]*z[7] + z[4]*z[8] + 2*z[4]*z[9] + 2*z[4]*z[10] + 2*z[5]*z[9] + z[5]*z[10] + 2*z[6]^2 + 2*z[6]*z[7] + 2*z[6]*z[9] + z[6]*z[10] + z[7]^2 + 2*z[7]*z[8] + z[7]*z[9] + z[8]^2 + z[8]*z[10] + z[9]^2 + 2*z[10]^2, 2*z[1]*z[7] + z[1]*z[8] + 2*z[1]*z[9] + z[1]*z[10] + z[2]*z[3] + z[2]*z[7] + 2*z[2]*z[9] + 2*z[2]*z[10] + 2*z[3]^2 + 2*z[3]*z[7] + 2*z[3]*z[8] + 2*z[3]*z[9] + z[3]*z[10] + 2*z[4]*z[7] + 2*z[4]*z[8] + 2*z[4]*z[9] + z[4]*z[10] + z[5]*z[7] + z[5]*z[8] + z[5]*z[9] + z[5]*z[10] + z[6]*z[7] + 2*z[6]*z[8] + z[6]*z[10] + z[7]^2 + 2*z[7]*z[9] + z[7]*z[10] + 2*z[9]^2 + z[9]*z[10], z[1]*z[3] + 2*z[1]*z[7] + z[1]*z[8] + 2*z[1]*z[9] + 2*z[1]*z[10] + z[2]*z[7] + z[2]*z[8] + z[2]*z[9] + 2*z[3]*z[6] + 2*z[3]*z[7] + 2*z[3]*z[8] + 2*z[3]*z[9] + 2*z[3]*z[10] + 2*z[4]*z[6] + z[4]*z[8] + z[4]*z[10] + 2*z[5]*z[7] + 2*z[5]*z[8] + z[5]*z[9] + z[6]*z[7] + 2*z[6]*z[10] + z[7]^2 + z[7]*z[8] + z[7]*z[9] + 2*z[7]*z[10] + 2*z[8]^2 + 2*z[8]*z[9] + 2*z[8]*z[10] + z[9]^2 + 2*z[9]*z[10] + z[10]^2, z[1]*z[7] + 2*z[1]*z[8] + z[1]*z[9] + z[2]^2 + 2*z[2]*z[7] + 2*z[2]*z[9] + 2*z[3]^2 + 2*z[3]*z[7] + 2*z[3]*z[8] + 2*z[3]*z[10] + 2*z[4]*z[6] + z[4]*z[8] + z[4]*z[10] + 2*z[5]*z[7] + 2*z[5]*z[9] + z[5]*z[10] + z[6]^2 + z[6]*z[7] + z[7]^2 + z[7]*z[9] + z[7]*z[10] + 2*z[8]^2 + z[8]*z[10] + 2*z[9]^2, z[1]*z[2] + 2*z[1]*z[7] + 2*z[1]*z[9] + z[1]*z[10] + z[2]*z[9] + z[2]*z[10] + 2*z[3]*z[6] + 2*z[3]*z[7] + z[3]*z[8] + z[3]*z[9] + z[3]*z[10] + z[4]*z[6] + z[4]*z[7] + z[4]*z[10] + 2*z[5]*z[8] + 2*z[5]*z[9] + 2*z[6]^2 + 2*z[6]*z[7] + z[6]*z[8] + z[7]*z[9] + 2*z[7]*z[10] + z[8]^2 + z[8]*z[9] + 2*z[9]^2 + z[9]*z[10] + 2*z[10]^2, z[1]^2 + 2*z[1]*z[7] + 2*z[1]*z[8] + 2*z[1]*z[9] + 2*z[1]*z[10] + z[2]*z[7] + 2*z[2]*z[8] + 2*z[2]*z[9] + z[3]*z[7] + 2*z[3]*z[8] + z[3]*z[9] + 2*z[3]*z[10] + 2*z[4]*z[6] + z[4]*z[7] + 2*z[4]*z[8] + z[5]*z[7] + z[5]*z[9] + z[5]*z[10] + z[6]^2 + 2*z[6]*z[7] + z[6]*z[8] + 2*z[6]*z[9] + 2*z[6]*z[10] + z[7]^2 + 2*z[7]*z[8] + z[7]*z[9] + 2*z[7]*z[10] + z[8]^2 + 2*z[8]*z[9] + 2*z[8]*z[10] + 2*z[9]*z[10] + z[10]^2, z[1]*z[7]*z[8] + 2*z[1]*z[7]*z[10] + 2*z[1]*z[8]*z[9] + z[1]*z[9]*z[10] + 2*z[1]*z[10]^2 + 2*z[2]*z[7]*z[9] + z[2]*z[8]^2 + z[2]*z[8]*z[10] + z[2]*z[9]^2 + z[2]*z[9]*z[10] + 2*z[2]*z[10]^2 + 2*z[3]*z[7]*z[8] + z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + 2*z[3]*z[8]*z[9] + z[3]*z[8]*z[10] + 2*z[3]*z[9]^2 + z[3]*z[9]*z[10] + z[4]*z[6]*z[9] + z[4]*z[6]*z[10] + 2*z[4]*z[7]*z[8] + z[4]*z[8]^2 + z[4]*z[8]*z[10] + z[4]*z[9]^2 + z[4]*z[9]*z[10] + z[4]*z[10]^2 + z[5]*z[7]*z[8] + z[5]*z[7]*z[9] + 2*z[5]*z[7]*z[10] + 2*z[5]*z[8]^2 + z[5]*z[8]*z[10] + 2*z[5]*z[9]^2 + z[5]*z[9]*z[10] + z[6]^2*z[8] + z[6]^2*z[9] + 2*z[6]^2*z[10] + 2*z[6]*z[7]*z[8] + z[6]*z[7]*z[9] + z[6]*z[8]*z[9] + 2*z[6]*z[9]^2 + z[6]*z[10]^2 + 2*z[7]^2*z[8] + z[7]^2*z[9] + z[7]^2*z[10] + z[7]*z[8]^2 + z[7]*z[8]*z[10] + 2*z[7]*z[9]^2 + z[7]*z[10]^2 + 2*z[8]^3 + z[8]^2*z[9] + 2*z[8]^2*z[10] + z[8]*z[9]^2 + z[8]*z[9]*z[10] + 2*z[8]*z[10]^2 + 2*z[9]^2*z[10] + 2*z[9]*z[10]^2, 2*z[1]*z[7]*z[8] + 2*z[1]*z[7]*z[9] + z[1]*z[7]*z[10] + 2*z[1]*z[8]^2 + 2*z[1]*z[10]^2 + z[2]*z[7]*z[8] + 2*z[2]*z[7]*z[9] + z[2]*z[7]*z[10] + 2*z[2]*z[8]^2 + 2*z[2]*z[8]*z[10] + z[2]*z[9]^2 + z[2]*z[10]^2 + z[3]*z[6]*z[9] + 2*z[3]*z[6]*z[10] + z[3]*z[7]*z[9] + 2*z[3]*z[8]^2 + 2*z[3]*z[8]*z[9] + z[3]*z[9]^2 + 2*z[3]*z[9]*z[10] + 2*z[3]*z[10]^2 + z[4]*z[6]*z[8] + 2*z[4]*z[6]*z[9] + z[4]*z[6]*z[10] + 2*z[4]*z[7]*z[8] + 2*z[4]*z[8]^2 + z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + z[5]*z[7]*z[8] + z[5]*z[7]*z[9] + z[5]*z[7]*z[10] + 2*z[5]*z[8]^2 + z[5]*z[8]*z[9] + 2*z[5]*z[8]*z[10] + 2*z[5]*z[9]^2 + z[5]*z[9]*z[10] + z[5]*z[10]^2 + z[6]^2*z[10] + z[6]*z[7]*z[9] + z[6]*z[7]*z[10] + 2*z[6]*z[8]^2 + z[6]*z[8]*z[9] + 2*z[6]*z[9]*z[10] + z[6]*z[10]^2 + 2*z[7]^2*z[8] + z[7]*z[8]^2 + z[7]*z[8]*z[9] + 2*z[7]*z[9]*z[10] + 2*z[7]*z[10]^2 + 2*z[8]^3 + z[8]^2*z[9] + z[8]*z[9]*z[10] + 2*z[8]*z[10]^2 + 2*z[9]^3 + 2*z[9]^2*z[10] + 2*z[10]^3, 2*z[1]*z[7]*z[8] + 2*z[1]*z[7]*z[10] + z[1]*z[8]^2 + z[1]*z[8]*z[9] + 2*z[1]*z[9]^2 + 2*z[1]*z[9]*z[10] + z[1]*z[10]^2 + 2*z[2]*z[7]*z[8] + 2*z[2]*z[7]*z[9] + 2*z[2]*z[7]*z[10] + 2*z[2]*z[8]*z[9] + 2*z[2]*z[8]*z[10] + 2*z[2]*z[9]^2 + z[2]*z[9]*z[10] + z[3]^2*z[9] + 2*z[3]^2*z[10] + z[3]*z[6]*z[8] + 2*z[3]*z[6]*z[10] + z[3]*z[7]*z[8] + 2*z[3]*z[7]*z[10] + z[3]*z[8]*z[10] + 2*z[3]*z[9]^2 + 2*z[3]*z[9]*z[10] + z[3]*z[10]^2 + z[4]*z[6]*z[9] + z[4]*z[6]*z[10] + 2*z[4]*z[7]*z[8] + 2*z[4]*z[7]*z[9] + z[4]*z[7]*z[10] + 2*z[4]*z[8]*z[10] + z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + 2*z[4]*z[10]^2 + z[5]*z[7]*z[9] + z[5]*z[8]*z[9] + 2*z[5]*z[8]*z[10] + 2*z[5]*z[9]^2 + 2*z[5]*z[9]*z[10] + z[5]*z[10]^2 + z[6]^2*z[9] + z[6]^2*z[10] + 2*z[6]*z[7]*z[8] + 2*z[6]*z[8]*z[9] + 2*z[6]*z[8]*z[10] + 2*z[6]*z[9]^2 + z[6]*z[9]*z[10] + 2*z[7]^2*z[8] + 2*z[7]^2*z[9] + 2*z[7]^2*z[10] + z[7]*z[8]^2 + 2*z[7]*z[8]*z[9] + z[7]*z[8]*z[10] + z[7]*z[9]*z[10] + z[7]*z[10]^2 + z[8]^3 + z[8]^2*z[10] + z[8]*z[9]*z[10] + 2*z[9]^3 + z[9]*z[10]^2 + 2*z[10]^3, z[1]*z[7]*z[9] + z[1]*z[7]*z[10] + 2*z[1]*z[8]^2 + 2*z[1]*z[8]*z[10] + z[2]*z[7]*z[8] + 2*z[2]*z[7]*z[9] + z[2]*z[7]*z[10] + z[2]*z[8]^2 + z[2]*z[8]*z[9] + 2*z[2]*z[9]*z[10] + z[2]*z[10]^2 + 2*z[3]*z[7]*z[8] + 2*z[3]*z[8]^2 + z[3]*z[8]*z[10] + 2*z[3]*z[9]*z[10] + 2*z[3]*z[10]^2 + z[4]*z[6]*z[9] + 2*z[4]*z[7]*z[8] + 2*z[4]*z[8]^2 + 2*z[4]*z[10]^2 + z[5]*z[7]*z[8] + z[5]*z[7]*z[10] + 2*z[5]*z[8]^2 + 2*z[5]*z[8]*z[9] + 2*z[5]*z[9]^2 + z[5]*z[9]*z[10] + z[5]*z[10]^2 + 2*z[6]^2*z[9] + z[6]^2*z[10] + z[6]*z[7]^2 + z[6]*z[7]*z[8] + z[6]*z[7]*z[9] + 2*z[6]*z[7]*z[10] + z[6]*z[8]^2 + 2*z[6]*z[8]*z[9] + z[6]*z[9]^2 + 2*z[6]*z[10]^2 + z[7]^2*z[9] + 2*z[7]^2*z[10] + z[7]*z[8]^2 + 2*z[7]*z[8]*z[9] + 2*z[7]*z[8]*z[10] + 2*z[7]*z[9]*z[10] + z[8]^3 + z[8]^2*z[9] + 2*z[8]^2*z[10] + z[8]*z[9]^2 + z[8]*z[10]^2 + z[9]^3 + z[9]^2*z[10] + 2*z[9]*z[10]^2 + z[10]^3, 2*z[1]*z[7]*z[8] + z[1]*z[8]^2 + z[1]*z[8]*z[9] + z[1]*z[8]*z[10] + z[1]*z[9]^2 + z[1]*z[10]^2 + z[2]*z[7]*z[8] + z[2]*z[7]*z[9] + 2*z[2]*z[7]*z[10] + 2*z[2]*z[8]^2 + 2*z[2]*z[8]*z[9] + 2*z[2]*z[9]^2 + 2*z[2]*z[9]*z[10] + z[2]*z[10]^2 + z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + z[3]*z[8]^2 + 2*z[3]*z[8]*z[10] + z[3]*z[9]^2 + 2*z[3]*z[10]^2 + z[4]*z[6]*z[8] + 2*z[4]*z[7]*z[8] + 2*z[4]*z[7]*z[9] + 2*z[4]*z[7]*z[10] + z[4]*z[8]^2 + 2*z[4]*z[8]*z[9] + z[4]*z[8]*z[10] + 2*z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + z[5]*z[7]^2 + 2*z[5]*z[7]*z[8] + 2*z[5]*z[7]*z[9] + 2*z[5]*z[8]*z[9] + z[5]*z[9]*z[10] + 2*z[5]*z[10]^2 + z[6]^2*z[9] + 2*z[6]*z[7]*z[8] + 2*z[6]*z[7]*z[9] + z[6]*z[8]^2 + z[6]*z[8]*z[9] + 2*z[6]*z[8]*z[10] + 2*z[6]*z[9]^2 + z[7]^2*z[8] + z[7]^2*z[9] + 2*z[7]^2*z[10] + 2*z[7]*z[8]*z[9] + 2*z[7]*z[8]*z[10] + 2*z[7]*z[10]^2 + z[8]^3 + 2*z[8]^2*z[10] + 2*z[8]*z[9]^2 + z[8]*z[9]*z[10] + z[8]*z[10]^2 + z[9]^3, 2*z[1]*z[7]*z[9] + 2*z[1]*z[7]*z[10] + 2*z[1]*z[8]^2 + 2*z[1]*z[8]*z[10] + z[1]*z[9]^2 + z[1]*z[9]*z[10] + 2*z[1]*z[10]^2 + 2*z[2]*z[7]*z[9] + 2*z[2]*z[7]*z[10] + z[2]*z[8]^2 + 2*z[2]*z[8]*z[9] + 2*z[2]*z[8]*z[10] + 2*z[2]*z[10]^2 + z[3]*z[7]*z[8] + 2*z[3]*z[7]*z[9] + 2*z[3]*z[7]*z[10] + z[3]*z[8]^2 + 2*z[3]*z[8]*z[9] + z[3]*z[8]*z[10] + 2*z[3]*z[9]^2 + 2*z[3]*z[9]*z[10] + 2*z[4]*z[6]*z[9] + z[4]*z[7]^2 + 2*z[4]*z[7]*z[9] + z[4]*z[7]*z[10] + 2*z[4]*z[8]^2 + z[4]*z[8]*z[9] + z[4]*z[8]*z[10] + z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + 2*z[5]*z[7]*z[8] + z[5]*z[7]*z[9] + z[5]*z[8]^2 + 2*z[5]*z[8]*z[9] + z[5]*z[9]^2 + 2*z[5]*z[9]*z[10] + 2*z[6]^2*z[8] + 2*z[6]^2*z[9] + z[6]^2*z[10] + 2*z[6]*z[7]^2 + 2*z[6]*z[7]*z[8] + z[6]*z[7]*z[9] + z[6]*z[8]^2 + 2*z[6]*z[8]*z[9] + z[6]*z[10]^2 + 2*z[7]^2*z[8] + 2*z[7]^2*z[9] + z[7]*z[8]^2 + 2*z[7]*z[8]*z[9] + z[7]*z[8]*z[10] + 2*z[7]*z[9]^2 + z[7]*z[10]^2 + z[8]^2*z[9] + 2*z[8]*z[9]^2 + z[8]*z[9]*z[10] + z[9]^3 + 2*z[9]^2*z[10] + 2*z[9]*z[10]^2 + z[10]^3, 2*z[1]*z[7]*z[9] + 2*z[1]*z[7]*z[10] + z[1]*z[8]^2 + z[1]*z[8]*z[10] + z[1]*z[10]^2 + 2*z[2]*z[7]*z[8] + z[2]*z[7]*z[9] + 2*z[2]*z[9]^2 + 2*z[2]*z[9]*z[10] + z[2]*z[10]^2 + 2*z[3]*z[6]*z[9] + z[3]*z[6]*z[10] + z[3]*z[7]^2 + 2*z[3]*z[7]*z[8] + 2*z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + 2*z[3]*z[8]^2 + 2*z[3]*z[10]^2 + z[4]*z[6]*z[8] + z[4]*z[6]*z[9] + 2*z[4]*z[6]*z[10] + 2*z[4]*z[7]*z[8] + z[4]*z[7]*z[9] + 2*z[4]*z[7]*z[10] + z[4]*z[8]^2 + z[4]*z[8]*z[10] + z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + 2*z[4]*z[10]^2 + 2*z[5]*z[7]^2 + z[5]*z[7]*z[8] + z[5]*z[7]*z[9] + 2*z[5]*z[7]*z[10] + z[5]*z[8]^2 + z[5]*z[8]*z[9] + 2*z[5]*z[8]*z[10] + z[5]*z[10]^2 + 2*z[6]^2*z[8] + 2*z[6]^2*z[9] + 2*z[6]^2*z[10] + 2*z[6]*z[7]^2 + z[6]*z[7]*z[8] + 2*z[6]*z[7]*z[9] + 2*z[6]*z[7]*z[10] + z[6]*z[8]^2 + z[6]*z[8]*z[9] + z[6]*z[9]^2 + z[7]^2*z[9] + 2*z[7]^2*z[10] + z[7]*z[8]^2 + z[7]*z[8]*z[9] + z[7]*z[9]^2 + 2*z[7]*z[9]*z[10] + z[7]*z[10]^2 + z[8]^3 + z[8]^2*z[9] + z[8]^2*z[10] + 2*z[8]*z[9]^2 + 2*z[8]*z[10]^2 + 2*z[9]^2*z[10] + 2*z[10]^3, 2*z[1]*z[7]*z[8] + 2*z[1]*z[7]*z[9] + 2*z[1]*z[9]*z[10] + z[2]*z[7]^2 + 2*z[2]*z[7]*z[8] + 2*z[2]*z[7]*z[10] + z[2]*z[8]^2 + z[2]*z[8]*z[9] + 2*z[2]*z[9]^2 + z[2]*z[9]*z[10] + z[3]*z[6]*z[9] + 2*z[3]*z[6]*z[10] + z[3]*z[7]^2 + 2*z[3]*z[7]*z[8] + z[3]*z[7]*z[9] + 2*z[3]*z[7]*z[10] + 2*z[3]*z[8]^2 + z[3]*z[8]*z[9] + z[3]*z[10]^2 + z[4]*z[6]*z[10] + 2*z[4]*z[7]^2 + 2*z[4]*z[7]*z[8] + 2*z[4]*z[7]*z[9] + z[4]*z[7]*z[10] + 2*z[4]*z[8]^2 + 2*z[4]*z[8]*z[10] + 2*z[4]*z[9]*z[10] + 2*z[5]*z[7]^2 + 2*z[5]*z[7]*z[8] + 2*z[5]*z[7]*z[9] + z[5]*z[7]*z[10] + z[5]*z[8]*z[9] + z[5]*z[9]^2 + 2*z[5]*z[9]*z[10] + 2*z[5]*z[10]^2 + z[6]^2*z[9] + 2*z[6]^2*z[10] + 2*z[6]*z[7]^2 + 2*z[6]*z[7]*z[8] + z[6]*z[7]*z[9] + 2*z[6]*z[8]^2 + z[6]*z[9]^2 + 2*z[6]*z[10]^2 + z[7]^3 + 2*z[7]^2*z[9] + z[7]^2*z[10] + z[7]*z[8]^2 + z[7]*z[9]^2 + 2*z[7]*z[9]*z[10] + z[7]*z[10]^2 + 2*z[8]^3 + 2*z[8]^2*z[9] + z[8]^2*z[10] + z[8]*z[10]^2 + 2*z[9]^3 + 2*z[9]^2*z[10] + z[9]*z[10]^2 + z[10]^3, z[1]*z[7]^2 + 2*z[1]*z[7]*z[8] + z[1]*z[7]*z[9] + z[1]*z[7]*z[10] + 2*z[1]*z[8]^2 + z[1]*z[8]*z[9] + z[1]*z[9]*z[10] + 2*z[1]*z[10]^2 + z[2]*z[7]^2 + z[2]*z[7]*z[9] + 2*z[2]*z[7]*z[10] + 2*z[2]*z[8]^2 + 2*z[2]*z[8]*z[9] + z[2]*z[8]*z[10] + z[2]*z[9]^2 + z[3]*z[6]*z[9] + 2*z[3]*z[6]*z[10] + z[3]*z[7]^2 + z[3]*z[7]*z[10] + z[3]*z[8]^2 + 2*z[3]*z[8]*z[10] + z[3]*z[9]*z[10] + 2*z[3]*z[10]^2 + z[4]*z[6]*z[8] + 2*z[4]*z[6]*z[9] + z[4]*z[6]*z[10] + 2*z[4]*z[7]^2 + z[4]*z[7]*z[8] + 2*z[4]*z[7]*z[9] + 2*z[4]*z[7]*z[10] + z[4]*z[8]*z[9] + 2*z[4]*z[8]*z[10] + z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + z[4]*z[10]^2 + 2*z[5]*z[7]^2 + 2*z[5]*z[7]*z[8] + 2*z[5]*z[7]*z[9] + z[5]*z[7]*z[10] + 2*z[5]*z[8]^2 + 2*z[5]*z[8]*z[9] + 2*z[5]*z[8]*z[10] + z[5]*z[9]*z[10] + 2*z[5]*z[10]^2 + z[6]^2*z[8] + z[6]^2*z[9] + 2*z[6]^2*z[10] + 2*z[6]*z[7]^2 + 2*z[6]*z[7]*z[9] + 2*z[6]*z[7]*z[10] + 2*z[6]*z[8]^2 + z[6]*z[8]*z[9] + 2*z[6]*z[8]*z[10] + z[6]*z[9]^2 + z[6]*z[9]*z[10] + z[6]*z[10]^2 + z[7]^3 + z[7]^2*z[8] + z[7]^2*z[9] + 2*z[7]^2*z[10] + 2*z[7]*z[8]^2 + z[7]*z[8]*z[9] + z[7]*z[9]^2 + z[8]^3 + 2*z[8]*z[9]^2 + z[9]^2*z[10] + z[9]*z[10]^2 + 2*z[10]^3, z[1]*z[7]*z[9] + z[1]*z[7]*z[10] + 2*z[1]*z[9]*z[10] + 2*z[1]*z[10]^2 + 2*z[2]*z[7]^2 + 2*z[2]*z[7]*z[8] + z[2]*z[8]*z[10] + 2*z[2]*z[9]^2 + z[2]*z[9]*z[10] + 2*z[2]*z[10]^2 + 2*z[3]*z[6]*z[9] + z[3]*z[6]*z[10] + z[3]*z[7]*z[8] + z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + 2*z[3]*z[8]^2 + z[3]*z[8]*z[9] + 2*z[4]*z[7]^2 + z[4]*z[7]*z[8] + z[4]*z[7]*z[10] + z[4]*z[9]^2 + z[4]*z[9]*z[10] + 2*z[4]*z[10]^2 + z[5]*z[7]*z[8] + 2*z[5]*z[7]*z[9] + z[5]*z[7]*z[10] + 2*z[5]*z[8]^2 + 2*z[5]*z[8]*z[9] + z[5]*z[9]^2 + 2*z[5]*z[10]^2 + z[6]^2*z[7] + 2*z[6]^2*z[8] + z[6]^2*z[9] + 2*z[6]^2*z[10] + 2*z[6]*z[7]^2 + z[6]*z[7]*z[8] + 2*z[6]*z[7]*z[9] + z[6]*z[8]*z[9] + 2*z[6]*z[9]*z[10] + 2*z[6]*z[10]^2 + 2*z[7]^3 + 2*z[7]^2*z[8] + 2*z[7]^2*z[9] + z[7]*z[8]*z[10] + 2*z[7]*z[9]*z[10] + 2*z[8]^2*z[9] + z[8]^2*z[10] + 2*z[8]*z[9]*z[10] + 2*z[9]^3 + z[9]*z[10]^2 + 2*z[10]^3, 2*z[1]*z[7]^2 + 2*z[1]*z[7]*z[8] + 2*z[1]*z[7]*z[9] + 2*z[1]*z[7]*z[10] + 2*z[1]*z[8]*z[9] + 2*z[1]*z[9]^2 + z[1]*z[9]*z[10] + z[2]*z[7]*z[9] + 2*z[2]*z[7]*z[10] + 2*z[2]*z[8]^2 + 2*z[2]*z[8]*z[9] + z[2]*z[8]*z[10] + z[2]*z[9]^2 + z[2]*z[9]*z[10] + 2*z[2]*z[10]^2 + z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + 2*z[3]*z[8]^2 + z[3]*z[8]*z[9] + z[3]*z[9]^2 + z[3]*z[9]*z[10] + z[3]*z[10]^2 + z[4]*z[6]*z[7] + 2*z[4]*z[6]*z[8] + 2*z[4]*z[6]*z[9] + 2*z[4]*z[6]*z[10] + z[4]*z[7]*z[8] + z[4]*z[7]*z[9] + 2*z[4]*z[7]*z[10] + 2*z[4]*z[8]^2 + z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + z[4]*z[10]^2 + z[5]*z[7]^2 + 2*z[5]*z[7]*z[8] + 2*z[5]*z[7]*z[9] + 2*z[5]*z[7]*z[10] + 2*z[5]*z[8]^2 + z[5]*z[8]*z[9] + 2*z[5]*z[8]*z[10] + z[5]*z[9]*z[10] + 2*z[5]*z[10]^2 + z[6]^2*z[7] + 2*z[6]^2*z[8] + z[6]^2*z[9] + 2*z[6]^2*z[10] + z[6]*z[7]^2 + z[6]*z[7]*z[8] + 2*z[6]*z[7]*z[9] + z[6]*z[8]^2 + z[6]*z[8]*z[9] + z[6]*z[8]*z[10] + z[6]*z[9]^2 + z[6]*z[9]*z[10] + 2*z[6]*z[10]^2 + 2*z[7]^2*z[8] + 2*z[7]^2*z[9] + 2*z[7]^2*z[10] + z[7]*z[8]^2 + 2*z[7]*z[8]*z[9] + 2*z[7]*z[8]*z[10] + 2*z[7]*z[9]^2 + 2*z[7]*z[10]^2 + z[8]^3 + 2*z[8]^2*z[9] + 2*z[8]^2*z[10] + 2*z[8]*z[9]^2 + z[8]*z[9]*z[10] + z[8]*z[10]^2 + z[9]^3 + z[9]^2*z[10] + z[9]*z[10]^2, z[1]*z[7]*z[9] + z[1]*z[7]*z[10] + 2*z[1]*z[8]^2 + z[1]*z[8]*z[10] + 2*z[1]*z[9]^2 + 2*z[1]*z[9]*z[10] + z[1]*z[10]^2 + z[2]*z[7]^2 + z[2]*z[7]*z[8] + 2*z[2]*z[7]*z[9] + z[2]*z[8]*z[9] + 2*z[2]*z[8]*z[10] + z[2]*z[9]^2 + z[2]*z[9]*z[10] + z[3]*z[6]*z[7] + z[3]*z[6]*z[8] + 2*z[3]*z[6]*z[10] + z[3]*z[7]^2 + 2*z[3]*z[7]*z[8] + z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + 2*z[3]*z[8]*z[10] + 2*z[3]*z[9]^2 + 2*z[4]*z[6]*z[8] + 2*z[4]*z[6]*z[9] + z[4]*z[7]*z[8] + z[4]*z[7]*z[9] + z[4]*z[7]*z[10] + z[4]*z[9]^2 + z[4]*z[9]*z[10] + z[4]*z[10]^2 + z[5]*z[7]^2 + z[5]*z[7]*z[9] + z[5]*z[8]^2 + z[5]*z[8]*z[9] + 2*z[5]*z[8]*z[10] + 2*z[5]*z[9]*z[10] + 2*z[5]*z[10]^2 + 2*z[6]^2*z[8] + 2*z[6]^2*z[9] + z[6]*z[7]*z[9] + z[6]*z[7]*z[10] + 2*z[6]*z[8]^2 + 2*z[6]*z[8]*z[9] + z[6]*z[8]*z[10] + 2*z[6]*z[9]*z[10] + z[7]^3 + z[7]^2*z[8] + z[7]^2*z[9] + 2*z[7]*z[8]^2 + z[7]*z[8]*z[9] + z[7]*z[8]*z[10] + 2*z[7]*z[9]^2 + z[7]*z[9]*z[10] + z[7]*z[10]^2 + 2*z[8]^2*z[9] + z[8]^2*z[10] + z[8]*z[9]^2 + 2*z[8]*z[9]*z[10] + z[8]*z[10]^2 + 2*z[9]^3 + 2*z[9]^2*z[10] + 2*z[9]*z[10]^2 + 2*z[10]^3, z[1]*z[7]*z[8] + 2*z[1]*z[7]*z[9] + z[1]*z[7]*z[10] + z[1]*z[8]^2 + z[1]*z[8]*z[9] + z[1]*z[8]*z[10] + z[1]*z[9]^2 + z[1]*z[9]*z[10] + 2*z[1]*z[10]^2 + 2*z[2]*z[7]*z[9] + 2*z[2]*z[8]*z[10] + z[2]*z[9]*z[10] + z[3]^2*z[7] + z[3]^2*z[8] + 2*z[3]^2*z[9] + 2*z[3]*z[6]*z[8] + 2*z[3]*z[6]*z[9] + 2*z[3]*z[6]*z[10] + 2*z[3]*z[7]*z[8] + z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + z[3]*z[8]^2 + z[3]*z[8]*z[9] + 2*z[3]*z[9]^2 + z[3]*z[9]*z[10] + 2*z[3]*z[10]^2 + 2*z[4]*z[6]*z[8] + z[4]*z[6]*z[10] + z[4]*z[7]*z[8] + z[4]*z[8]^2 + 2*z[4]*z[8]*z[10] + 2*z[4]*z[9]*z[10] + 2*z[5]*z[7]*z[8] + z[5]*z[7]*z[9] + 2*z[5]*z[8]^2 + 2*z[5]*z[8]*z[9] + z[5]*z[8]*z[10] + z[5]*z[9]*z[10] + 2*z[5]*z[10]^2 + 2*z[6]^2*z[10] + 2*z[6]*z[7]*z[9] + z[6]*z[7]*z[10] + 2*z[6]*z[8]^2 + z[6]*z[8]*z[9] + z[6]*z[9]^2 + 2*z[6]*z[10]^2 + 2*z[7]^2*z[8] + z[7]^2*z[9] + z[7]^2*z[10] + z[7]*z[8]^2 + z[7]*z[8]*z[10] + 2*z[7]*z[9]^2 + z[7]*z[9]*z[10] + z[7]*z[10]^2 + z[8]^3 + z[8]^2*z[9] + z[8]^2*z[10] + 2*z[8]*z[9]^2 + z[8]*z[9]*z[10] + z[9]^2*z[10], 2*z[1]*z[7]*z[8] + z[1]*z[7]*z[9] + z[1]*z[8]^2 + 2*z[1]*z[8]*z[9] + 2*z[1]*z[9]^2 + z[1]*z[9]*z[10] + z[1]*z[10]^2 + 2*z[2]*z[7]*z[9] + 2*z[2]*z[7]*z[10] + 2*z[2]*z[8]^2 + z[2]*z[8]*z[9] + z[2]*z[8]*z[10] + 2*z[2]*z[9]^2 + 2*z[2]*z[10]^2 + 2*z[3]*z[6]*z[7] + z[3]*z[6]*z[8] + z[3]*z[6]*z[9] + z[3]*z[6]*z[10] + 2*z[3]*z[7]^2 + 2*z[3]*z[7]*z[8] + 2*z[3]*z[7]*z[9] + 2*z[3]*z[7]*z[10] + 2*z[3]*z[8]^2 + 2*z[3]*z[8]*z[9] + 2*z[3]*z[8]*z[10] + 2*z[3]*z[9]^2 + 2*z[3]*z[9]*z[10] + z[4]*z[6]*z[7] + z[4]*z[6]*z[8] + 2*z[4]*z[6]*z[10] + z[4]*z[7]*z[8] + 2*z[4]*z[7]*z[9] + 2*z[4]*z[8]*z[10] + z[4]*z[9]^2 + 2*z[4]*z[9]*z[10] + 2*z[4]*z[10]^2 + 2*z[5]*z[7]^2 + 2*z[5]*z[8]^2 + z[5]*z[8]*z[9] + z[5]*z[8]*z[10] + z[5]*z[9]^2 + z[5]*z[10]^2 + z[6]^3 + z[6]^2*z[10] + z[6]*z[7]*z[8] + 2*z[6]*z[7]*z[10] + 2*z[6]*z[8]^2 + z[6]*z[8]*z[9] + z[6]*z[8]*z[10] + 2*z[6]*z[9]^2 + 2*z[6]*z[9]*z[10] + 2*z[7]^2*z[10] + 2*z[7]*z[8]^2 + z[7]*z[8]*z[9] + z[7]*z[8]*z[10] + 2*z[7]*z[9]^2 + z[7]*z[9]*z[10] + z[7]*z[10]^2 + 2*z[8]^3 + z[8]^2*z[9] + 2*z[8]^2*z[10] + z[8]*z[9]*z[10] + 2*z[8]*z[10]^2 + z[9]^3 + 2*z[9]^2*z[10] + z[9]*z[10]^2 + 2*z[10]^3, z[1]*z[7]*z[8] + z[1]*z[7]*z[9] + z[1]*z[7]*z[10] + 2*z[1]*z[8]^2 + z[1]*z[8]*z[9] + 2*z[1]*z[8]*z[10] + z[2]*z[7]^2 + 2*z[2]*z[7]*z[10] + 2*z[2]*z[8]^2 + 2*z[2]*z[8]*z[9] + 2*z[2]*z[8]*z[10] + 2*z[3]*z[6]*z[8] + z[3]*z[6]*z[10] + 2*z[3]*z[7]*z[8] + z[3]*z[7]*z[9] + 2*z[3]*z[7]*z[10] + 2*z[3]*z[8]*z[9] + 2*z[3]*z[8]*z[10] + 2*z[3]*z[9]^2 + 2*z[3]*z[9]*z[10] + z[3]*z[10]^2 + z[4]*z[6]^2 + 2*z[4]*z[6]*z[10] + 2*z[4]*z[7]^2 + 2*z[4]*z[7]*z[8] + z[4]*z[7]*z[9] + z[4]*z[8]^2 + z[4]*z[8]*z[9] + 2*z[4]*z[8]*z[10] + z[4]*z[9]*z[10] + 2*z[4]*z[10]^2 + z[5]*z[7]^2 + z[5]*z[7]*z[8] + 2*z[5]*z[8]^2 + z[5]*z[8]*z[10] + 2*z[5]*z[9]^2 + 2*z[6]^2*z[8] + 2*z[6]^2*z[9] + z[6]^2*z[10] + 2*z[6]*z[7]^2 + 2*z[6]*z[7]*z[8] + z[6]*z[7]*z[9] + 2*z[6]*z[7]*z[10] + z[6]*z[8]*z[9] + 2*z[6]*z[8]*z[10] + z[6]*z[9]*z[10] + z[7]^3 + 2*z[7]^2*z[8] + z[7]^2*z[9] + 2*z[7]*z[8]*z[9] + 2*z[7]*z[8]*z[10] + z[7]*z[9]^2 + z[7]*z[10]^2 + 2*z[8]^2*z[10] + z[8]*z[9]^2 + 2*z[9]^2*z[10] + z[9]*z[10]^2 + 2*z[10]^3, z[1]*z[7]*z[8] + 2*z[1]*z[7]*z[9] + z[1]*z[8]^2 + z[1]*z[8]*z[10] + 2*z[1]*z[9]^2 + z[1]*z[9]*z[10] + z[1]*z[10]^2 + 2*z[2]*z[7]^2 + 2*z[2]*z[7]*z[8] + 2*z[2]*z[8]^2 + 2*z[2]*z[8]*z[9] + z[2]*z[10]^2 + 2*z[3]^2*z[7] + z[3]^2*z[8] + 2*z[3]^2*z[9] + z[3]*z[6]^2 + 2*z[3]*z[6]*z[8] + z[3]*z[6]*z[9] + 2*z[3]*z[6]*z[10] + z[3]*z[7]^2 + z[3]*z[7]*z[8] + z[3]*z[7]*z[9] + z[3]*z[7]*z[10] + z[3]*z[8]*z[10] + 2*z[3]*z[9]^2 + z[3]*z[9]*z[10] + z[3]*z[10]^2 + 2*z[4]*z[6]*z[8] + 2*z[4]*z[7]^2 + z[4]*z[7]*z[9] + 2*z[4]*z[7]*z[10] + z[4]*z[8]^2 + z[4]*z[8]*z[9] + z[4]*z[8]*z[10] + 2*z[4]*z[9]^2 + 2*z[4]*z[10]^2 + 2*z[5]*z[7]*z[9] + 2*z[5]*z[7]*z[10] + 2*z[5]*z[8]^2 + 2*z[5]*z[8]*z[10] + 2*z[5]*z[9]^2 + z[6]^2*z[7] + 2*z[6]^2*z[9] + z[6]^2*z[10] + 2*z[6]*z[7]^2 + z[6]*z[7]*z[9] + 2*z[6]*z[7]*z[10] + z[6]*z[8]^2 + 2*z[6]*z[8]*z[9] + z[6]*z[8]*z[10] + z[6]*z[9]^2 + 2*z[6]*z[9]*z[10] + z[6]*z[10]^2 + 2*z[7]^3 + 2*z[7]^2*z[9] + 2*z[7]^2*z[10] + z[7]*z[8]*z[10] + 2*z[7]*z[9]^2 + 2*z[7]*z[9]*z[10] + 2*z[7]*z[10]^2 + 2*z[8]^3 + z[8]^2*z[10] + 2*z[8]*z[9]*z[10] + z[9]^3 + 2*z[9]^2*z[10] + 2*z[10]^3, 2*z[1]*z[7]*z[8]*z[9] + z[1]*z[7]*z[8]*z[10] + 2*z[1]*z[7]*z[9]^2 + 2*z[1]*z[7]*z[9]*z[10] + 2*z[1]*z[7]*z[10]^2 + 2*z[1]*z[8]^2*z[10] + 2*z[1]*z[8]*z[9]^2 + z[1]*z[8]*z[10]^2 + 2*z[1]*z[9]^3 + 2*z[1]*z[10]^3 + z[2]*z[7]*z[8]^2 + z[2]*z[7]*z[9]*z[10] + z[2]*z[8]^3 + z[2]*z[8]^2*z[9] + 2*z[2]*z[8]^2*z[10] + 2*z[2]*z[8]*z[9]^2 + 2*z[2]*z[8]*z[10]^2 + 2*z[2]*z[9]^3 + 2*z[2]*z[9]^2*z[10] + z[2]*z[10]^3 + z[3]^2*z[6]*z[9] + 2*z[3]^2*z[6]*z[10] + z[3]^2*z[8]^2 + z[3]^2*z[8]*z[10] + z[3]^2*z[9]^2 + z[3]*z[6]*z[9]*z[10] + 2*z[3]*z[7]*z[8]*z[9] + z[3]*z[7]*z[8]*z[10] + z[3]*z[7]*z[9]^2 + z[3]*z[7]*z[9]*z[10] + 2*z[3]*z[7]*z[10]^2 + z[3]*z[8]^2*z[9] + z[3]*z[8]^2*z[10] + z[3]*z[8]*z[9]^2 + 2*z[3]*z[8]*z[9]*z[10] + z[3]*z[8]*z[10]^2 + 2*z[3]*z[9]^3 + z[3]*z[9]^2*z[10] + z[3]*z[9]*z[10]^2 + z[4]*z[6]*z[9]^2 + z[4]*z[6]*z[9]*z[10] + z[4]*z[6]*z[10]^2 + z[4]*z[7]*z[8]*z[9] + 2*z[4]*z[7]*z[8]*z[10] + z[4]*z[8]^3 + 2*z[4]*z[8]^2*z[10] + z[4]*z[8]*z[9]^2 + z[4]*z[8]*z[10]^2 + 2*z[4]*z[9]*z[10]^2 + 2*z[5]*z[7]*z[8]^2 + 2*z[5]*z[7]*z[8]*z[10] + 2*z[5]*z[7]*z[9]^2 + z[5]*z[7]*z[9]*z[10] + 2*z[5]*z[8]^3 + z[5]*z[8]^2*z[9] + 2*z[5]*z[8]^2*z[10] + z[5]*z[8]*z[9]*z[10] + z[5]*z[9]^3 + 2*z[5]*z[9]^2*z[10] + z[5]*z[9]*z[10]^2 + 2*z[6]^2*z[9]*z[10] + 2*z[6]*z[7]*z[8]^2 + 2*z[6]*z[7]*z[8]*z[9] + z[6]*z[7]*z[8]*z[10] + 2*z[6]*z[7]*z[9]^2 + 2*z[6]*z[7]*z[9]*z[10] + z[6]*z[8]^2*z[9] + 2*z[6]*z[8]^2*z[10] + z[6]*z[8]*z[9]*z[10] + z[6]*z[9]^3 + 2*z[6]*z[9]*z[10]^2 + 2*z[6]*z[10]^3 + 2*z[7]^2*z[8]^2 + z[7]^2*z[8]*z[9] + 2*z[7]^2*z[9]*z[10] + 2*z[7]^2*z[10]^2 + z[7]*z[8]^3 + z[7]*z[8]^2*z[9] + 2*z[7]*z[8]^2*z[10] + z[7]*z[8]*z[10]^2 + 2*z[7]*z[9]^2*z[10] + z[7]*z[9]*z[10]^2 + z[7]*z[10]^3 + z[8]^4 + z[8]^3*z[9] + 2*z[8]^3*z[10] + z[8]^2*z[9]^2 + 2*z[8]^2*z[9]*z[10] + 2*z[8]*z[9]*z[10]^2 + z[9]^3*z[10] + z[9]^2*z[10]^2 + 2*z[10]^4, z[1]*z[7]*z[8]*z[9] + z[1]*z[7]*z[9]^2 + z[1]*z[8]^3 + 2*z[1]*z[8]^2*z[9] + z[1]*z[8]^2*z[10] + z[1]*z[8]*z[9]*z[10] + z[1]*z[8]*z[10]^2 + 2*z[1]*z[9]^3 + 2*z[1]*z[9]^2*z[10] + 2*z[1]*z[10]^3 + z[2]*z[7]*z[8]^2 + 2*z[2]*z[7]*z[8]*z[9] + z[2]*z[8]^3 + z[2]*z[8]^2*z[10] + z[2]*z[8]*z[9]^2 + z[2]*z[9]^2*z[10] + 2*z[2]*z[10]^3 + 2*z[3]*z[6]*z[9]*z[10] + z[3]*z[6]*z[10]^2 + 2*z[3]*z[7]*z[8]^2 + z[3]*z[7]*z[8]*z[10] + 2*z[3]*z[7]*z[9]^2 + 2*z[3]*z[7]*z[10]^2 + 2*z[3]*z[8]^3 + z[3]*z[8]^2*z[10] + 2*z[3]*z[8]*z[9]^2 + z[3]*z[8]*z[9]*z[10] + z[3]*z[8]*z[10]^2 + z[3]*z[9]^3 + z[3]*z[9]^2*z[10] + 2*z[3]*z[9]*z[10]^2 + z[3]*z[10]^3 + z[4]*z[7]*z[8]*z[9] + z[4]*z[7]*z[8]*z[10] + 2*z[4]*z[7]*z[9]*z[10] + z[4]*z[8]^3 + z[4]*z[8]*z[9]*z[10] + 2*z[4]*z[8]*z[10]^2 + z[4]*z[9]^3 + 2*z[4]*z[9]^2*z[10] + 2*z[4]*z[9]*z[10]^2 + z[4]*z[10]^3 + z[5]*z[7]*z[8]^2 + z[5]*z[7]*z[8]*z[9] + 2*z[5]*z[7]*z[9]*z[10] + 2*z[5]*z[7]*z[10]^2 + 2*z[5]*z[8]^2*z[9] + 2*z[5]*z[8]^2*z[10] + 2*z[5]*z[8]*z[9]^2 + z[5]*z[8]*z[10]^2 + z[5]*z[9]^3 + z[5]*z[9]*z[10]^2 + 2*z[6]^2*z[9]^2 + 2*z[6]^2*z[9]*z[10] + 2*z[6]*z[7]*z[8]^2 + z[6]*z[7]*z[8]*z[10] + z[6]*z[7]*z[9]^2 + 2*z[6]*z[7]*z[9]*z[10] + z[6]*z[7]*z[10]^2 + z[6]*z[8]^2*z[10] + 2*z[6]*z[8]*z[9]^2 + 2*z[6]*z[9]^3 + z[6]*z[9]^2*z[10] + z[6]*z[9]*z[10]^2 + 2*z[6]*z[10]^3 + z[7]^3*z[8] + z[7]^3*z[9] + 2*z[7]^3*z[10] + z[7]^2*z[8]^2 + 2*z[7]^2*z[8]*z[9] + 2*z[7]^2*z[9]*z[10] + z[7]^2*z[10]^2 + 2*z[7]*z[8]^3 + 2*z[7]*z[8]^2*z[9] + z[7]*z[8]*z[9]*z[10] + 2*z[7]*z[8]*z[10]^2 + 2*z[7]*z[9]^3 + z[7]*z[9]^2*z[10] + 2*z[7]*z[9]*z[10]^2 + 2*z[7]*z[10]^3 + 2*z[8]^4 + 2*z[8]^3*z[10] + z[8]^2*z[9]*z[10] + z[8]*z[9]^2*z[10] + 2*z[8]*z[9]*z[10]^2 + z[8]*z[10]^3 + 2*z[9]^3*z[10]];

  I = ideal(S, V);
  X = variety(I, check = false, is_radical = true);
  
  @test !is_zero(canonical_bundle(X))
end

