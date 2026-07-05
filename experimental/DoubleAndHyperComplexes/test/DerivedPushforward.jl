@testset "Variation of Hodge structures on elliptic curves" begin
  # We define a family of elliptic curves over a 1-dimensional base
  A, (t,) = polynomial_ring(QQ, [:t])
  P = projective_space(A, [:x, :y, :z])
  S = homogeneous_coordinate_ring(P)
  n = ngens(S)
  f1 = sum(x^n for x in gens(S))
  f0 = prod(gens(S))
  #f0 = S[1]^n
  # For t=0 this degenerates to three lines and we expect jumps in cohomology there
  f = (1-t)*f0 + t*f1

  # derived pushforward of the structure sheaf of the curve
  IPX, inc_X = sub(P, f)
  S1 = graded_free_module(S, [0])
  I = ideal(S, f)
  IS1, inc = I*S1
  M, _ = cokernel(inc)
  # Nothing is expected here to jump: H^0 is always free of rank one and 
  # so the remainder of this short complex must also be of constant rank.
  st1 = Oscar._derived_pushforward(M)

  Z0, _ = kernel(map(st1, 1, (0,)))
  @test !iszero(Z0)
  Z1, _ = kernel(map(st1, 1, (-1,)))
  B1, _ = image(map(st1, 1, (0,)))
  H1, _ = quo(Z1, B1)
  @test !iszero(H1)
  B2, _ = image(map(st1, 1, (-1,)))
  H2, _ = quo(st1[-2], B2)
  @test iszero(H2)
  

  # derived pushforward of the relative Ω¹
  M1 = Oscar.relative_cotangent_module(IPX)

  M1, _ = pushforward(inc_X, M1)

  #M1 = simplify(M1)[1] # Does not preserve gradings and can not be used

  st2 = Oscar._derived_pushforward(M1)

  Z0, _ = kernel(map(st2, 1, (0,)))
  @test !iszero(Z1)
  Z1, _ = kernel(map(st2, 1, (-1,)))
  B1, _ = image(map(st2, 1, (0,)))
  H1, _ = quo(Z1, B1)
  @test !iszero(H1)
  B2, _ = image(map(st2, 1, (-1,)))
  H2, _ = quo(st2[-2], B2)
  @test iszero(H2)
  

  # The relative R¹π_* Ω¹_C
  H11 = homology(st2.internal_complex, 1, (-1,))[1]
  @test !iszero(H11)

  # look at how the matrices jump
  B1 = matrix(map(st2, 1, (0,)))
  B2 = matrix(map(st2, 1, (-1,)))

  # H^{0, 1} jumps in the singular fiber
  B10 = map_entries(g->evaluate(g, [0]), B1)
  B11 = map_entries(g->evaluate(g, [1]), B1)
  @test rank(B10) != rank(B11)

  # H^{1, 1} stays constant throughout the family
  B20 = map_entries(g->evaluate(g, [0]), B2)
  B21 = map_entries(g->evaluate(g, [1]), B2)
  @test rank(B20) == rank(B21)
end

@testset "simplification of matrices" begin
  R, (x, y) = QQ[:x, :y]
  A = sparse_matrix(R, [x^3 x y 3*x y^2; 2*x 2*y 5 x^2+1 5*y; x x x x x])
  B = deepcopy(A)
  S, Sinv, T, Tinv, ind = Oscar._simplify_matrix!(B)
  @test isone(matrix(S)*matrix(Sinv))
  @test isone(matrix(T)*matrix(Tinv))
  @test inv(matrix(S))*matrix(B)*matrix(T) == matrix(A)

  A = sparse_matrix(R, [-1 x y; 2*x 2*y 0])
  B = deepcopy(A)
  S, Sinv, T, Tinv, ind = Oscar._simplify_matrix!(B)
  @test isone(matrix(S)*matrix(Sinv))
  @test isone(matrix(T)*matrix(Tinv))
  @test inv(matrix(S))*matrix(B)*matrix(T) == matrix(A)

  A = sparse_matrix(R, [2*x 2*y 0; -1 x y])
  B = deepcopy(A)
  S, Sinv, T, Tinv, ind = Oscar._simplify_matrix!(B)
  @test isone(matrix(S)*matrix(Sinv))
  @test isone(matrix(T)*matrix(Tinv))
  @test inv(matrix(S))*matrix(B)*matrix(T) == matrix(A)
end

@testset "simplification of complexes" begin
  R, (x, y) = QQ[:x, :y]
  R3 = FreeMod(R, 3)
  v = x*R3[1] + 5* R3[2] + y*R3[3]

  C = koszul_complex(v)
  g, phi, psi = simplify(C)
  @test all(iszero(g[i]) for i in range(g))
  @test all(iszero(map(g, i)) for i in map_range(g))

  U = powers_of_element(y)
  L, _ = localization(R, U)
  L3 = FreeMod(L, 3)
  v = x*L3[1] + y^2* L3[2] + y*L3[3]
  C = koszul_complex(v)
  g, phi, psi = simplify(C)
  @test all(iszero(g[i]) for i in range(g))
  @test all(iszero(map(g, i)) for i in map_range(g))
end

@testset "simplified cohomology of K3-surfaces" begin
  # long running test. Possibly to be disabled.
  A, (t,) = polynomial_ring(QQ, [:t])
  P = projective_space(A, [:x, :y, :z, :w])
  S = homogeneous_coordinate_ring(P)
  n = ngens(S)
  f1 = sum(x^n for x in gens(S))
  f0 = prod(gens(S))
  f = (1-t)*f0 + t*f1
  #f = f1

  IPX, inc_X = sub(P, f)
  Sm4 = graded_free_module(S, [4])

  # Test the routine for free modules
  st0 = Oscar._derived_pushforward(Sm4);

  Z0, _ = kernel(map(st0, 1, (0,)))
  @test iszero(Z0)

  B1, _ = image(map(st0, 1, (0,)))
  Z1, _ = kernel(map(st0, 1, (-1,)))

  H1, _ = quo(Z1, B1)
  @test iszero(H1) 

  B2, _ = image(map(st0, 1, (-1,)))
  Z2, _ = kernel(map(st0, 1, (-2,)))

  H2, _ = quo(Z2, B2)
  @test iszero(H2)

  B3, _ = image(map(st0, 1, (-2,)))
  H3, _ = quo(st0[-3], B3)
  @test !iszero(H3)

  # Compute the higher direct images of the structure sheaf
  S1 = graded_free_module(S, [0])
  I = ideal(S, f)
  IS1, inc = I*S1
  M, _ = cokernel(inc)
  st0 = Oscar._derived_pushforward(M);

  Z0, _ = kernel(map(st0, 1, (0,)))
  @test !iszero(Z0)
  @test is_projective(Z0)[1] # H^{0, 0} does not change throughout the family

  B1, _ = image(map(st0, 1, (0,)))
  Z1, _ = kernel(map(st0, 1, (-1,)))

  H1, _ = quo(Z1, B1)
  @test iszero(H1) # H^{1, 0} = 0 for K3-surfaces

  B2, _ = image(map(st0, 1, (-1,)))
  Z2, _ = kernel(map(st0, 1, (-2,)))

  H2, _ = quo(Z2, B2)
  @test !iszero(H2)
  @test is_projective(H2)[1] # H^{2, 0} is constant in this family

  B3, _ = image(map(st0, 1, (-2,)))
  H3, _ = quo(st0[-3], B3)
  @test iszero(H3)
  

  M1 = Oscar.relative_cotangent_module(IPX)
  M1, _ = pushforward(inc_X, M1)
  st = Oscar._derived_pushforward(M1);

  stsp = simplify(Oscar.SimpleComplexWrapper(st));
  A = matrix(map(stsp, 1, (0,)))
  B = matrix(map(stsp, 1, (-1,)))

  # Test for a smooth fiber: H¹(X₁, Ω¹) = ℚ ²⁰.
  A1 = map_entries(f->f(one(ZZ)), A)
  B1 = map_entries(f->f(one(ZZ)), B)
  @test rank(stsp[-1]) - rank(A1) - rank(B1) == 20

  # Test for the singular fiber: H¹(X₀, Ω¹) = ℚ ²³.
  A0 = map_entries(f->f(zero(ZZ)), A)
  B0 = map_entries(f->f(zero(ZZ)), B)
  @test rank(stsp[-1]) - rank(A0) - rank(B0) == 23
end

@testset "K3 surface over a field" begin
  A = QQ
  P = projective_space(A, [:x, :y, :z, :w])
  S = homogeneous_coordinate_ring(P)
  n = ngens(S)
  f1 = sum(x^n for x in gens(S))
  f0 = prod(gens(S))
  f = f1 + 2*f0

  IPX, inc_X = sub(P, f)
  Sm4 = graded_free_module(S, [4])

  # Compute the higher direct images of the structure sheaf
  S1 = graded_free_module(S, [0])
  I = ideal(S, f)
  IS1, inc = I*S1
  M, _ = cokernel(inc)
  st0 = Oscar._derived_pushforward(M);
  st0_simp = simplify(st0);
  
  @test rank(st0_simp[0]) == 1
  @test iszero(st0_simp[-1])
  @test rank(st0_simp[-2]) == 1
  @test iszero(st0_simp[-3])

  # Compute the higher direct images of \Omega^1
  M1 = Oscar.relative_cotangent_module(IPX)
  M1, _ = pushforward(inc_X, M1)
  st = Oscar._derived_pushforward(M1);
  stsp = simplify(st);

  @test iszero(stsp[0])
  @test rank(stsp[-1]) == 20
  @test iszero(stsp[-2])
end

@testset "mild log singularities" begin
  B, (t,) = polynomial_ring(QQ, [:t])

  IP5 = projective_space(B, [:u, :v, :w, :x, :y, :z])
  S = homogeneous_coordinate_ring(IP5)
  (u, v, w, x, y, z) = gens(S)

  f = x*y - t^2*v^2 - w^2
  g = z*w - t*u*v

  I = ideal(S, [f, g])

  X, inc_X = sub(IP5, I)

  # Testing the structure sheaf
  M0, _ = pushforward(inc_X, graded_free_module(homogeneous_coordinate_ring(X), [0]))

  C = simplify(Oscar._derived_pushforward(M0));

  @test rank(C[0]) == 1 # Free of rank one
  @test all(k->iszero(C[-k]), 1:5)

  # Testing the Kaehler differentials
  M1, _ = pushforward(inc_X, Oscar.relative_cotangent_module(X))

  C = simplify(Oscar._derived_pushforward(M1));

  @test iszero(C[0])
  @test rank(C[-1]) == 1
  @test rank(C[-2]) == 2
  @test iszero(C[-3])
  @test iszero(C[-4])
  @test iszero(C[-5])

  @test iszero(map(C, 1, (-1,)))
end

@testset "Weyman complexes backbone" begin
  X = abelian_d10_pi6() # Should run, too, but is too long for regular test suite.
  #X = cubic_scroll()
  #X = bordiga() # Yet another possibility which runs fast.
  Q = homogeneous_coordinate_ring(X)
  S = base_ring(Q);
  M = canonical_bundle(X);
  res, aug = free_resolution(Oscar.SimpleFreeResolution, M);
  res = Oscar.ReflectedComplex(res);

  ctx = Oscar.PushForwardCtx(S);
  wctx = Oscar.WeymanCtx(ctx, res);

  mb = Oscar.get_macro_block!(wctx, -1, 4, :cohomology);
  F = res[-1]
  alpha = -degree(F[85])
  e = Oscar._minimal_exponent_vector(ctx, alpha)
  h4 = Oscar.simplified_strand(ctx, e, alpha)[-4]
  v = Oscar.MicroVec(mb, 85, e, h4[1]);
  mv = Oscar.MacroVec(v)
  Imv = Oscar.apply_I(mv)
  w = addmul!(mv, Oscar.apply_P(Imv)[2], -1)
  @test is_zero(w)

  F = res[-5]
  for k in 1:5
    mb = Oscar.get_macro_block!(wctx, -5, 4, :cohomology);
    alpha = -degree(F[k])
    e = Oscar._minimal_exponent_vector(ctx, alpha)
    h4 = Oscar.simplified_strand(ctx, e, alpha)[-4]
    v = Oscar.MicroVec(mb, k, e, h4[1]);
    mv = Oscar.MacroVec(v)
    Imv = Oscar.apply_I(mv)
    PImv = Oscar.apply_P(Imv)
    w = mv + (-1)*PImv[end]
    @test all(is_zero, PImv[1:end-1])
    @test is_zero(w)

    CImv = Oscar.apply_cech_map(Imv)
    HCImv = Oscar.apply_homotopy_map(CImv)
    HCImv2 = Oscar.apply_homotopy_map(Oscar.apply_cech_map(Imv))
    @test is_zero(HCImv + (-1)*HCImv2)

    phiv = Oscar.apply_weyman_differential(mv)
    @test is_zero(first(phiv))
  end

# F = res[-5]
# for k in 1:ngens(F)
#   # k = 26
#   mv = Oscar.MacroVec(wctx, -5, 4, :cohomology);
#   alpha = -degree(F[k])
#   e = Oscar._minimal_exponent_vector(ctx, alpha)
#   h4 = Oscar.simplified_strand(ctx, e, alpha)[-4]
#   for w in gens(h4)
#     # w = h4[23]
#     v = Oscar.MicroVec(mv, k, e, w);
#     a = Oscar.apply_weyman_differential(Oscar.MacroVec(v))
#   end
# end

end

@testset "new weyman complexes" begin
  X = abelian_d10_pi6() # Should run, too, but is too long for regular test suite.
  #X = cubic_scroll()
  #X = bordiga() # Yet another possibility which runs fast.
  Q = homogeneous_coordinate_ring(X)
  S = base_ring(Q);
  M = canonical_bundle(X);
  res, aug = free_resolution(Oscar.SimpleFreeResolution, M);
  ctx = Oscar.PushForwardCtx(S);
  #ctx.fixed_exponent_vector = [5 for _ in 1:5]
  dir_im = Oscar.DirectImageComplex(ctx, res);
  wctx = Oscar.WeymanCtx(ctx, Oscar.ReflectedComplex(res));
  mb = Oscar.get_macro_block!(wctx, -5, 4, :cohomology);
  v = Oscar.MicroVec(mb, 6, 70)
  w = Oscar.apply_weyman_differential(Oscar.MacroVec(v))
  
  w0 = Oscar.MacroVec(v)
  w1 = Oscar.apply_I(w0)
  w10 = Oscar.apply_P(w1)
  w2_1 = (-1)*Oscar.apply_phi.(w1)
  w2_2 = Oscar.apply_cech_map(w1)
  w2 = w2_1[2:end] + w2_2
  pushfirst!(w2, first(w2_1))
  
  w3 = Oscar.apply_P(w2)
  w3_alt = Oscar.apply_P(w2_1)
  @assert w3 == w3_alt
  w3_aa = Oscar.apply_weyman_differential(w0)
  @assert w3 == w3_aa
  
  #[Oscar.project_to_weyman_complex!([Oscar.MacroVec(deepcopy(v))]) for (_, v) in Oscar.micro_vectors(w3[end])]
  
  w4 = Oscar.project_to_weyman_complex!(deepcopy(w3))
  w4_alt = Oscar.apply_weyman_differential(w0)
  @test w4 == w4_alt
end
  
@testset "NewToricCtx" begin
  X = hirzebruch_surface(NormalToricVariety, 3)

  ctx = Oscar.NewToricCtx(X)
  S = Oscar.graded_ring(ctx)
  G = grading_group(S)
  alpha = G[1]
  Oscar.all_monomials(ctx, alpha)
  Oscar.all_monomials_inv(ctx, alpha)
  p = [7, 2, 3, -5]
  str_inc = Oscar.fine_strand(ctx, p)
  str_simp = Oscar.simplified_fine_strand(ctx, p)
  @test Oscar.original_complex(str_simp) === str_inc
  q = [7, 2, -6, -5]
  str_inc = Oscar.fine_strand(ctx, q)
  str_simp = Oscar.simplified_fine_strand(ctx, q)
  @test Oscar.original_complex(str_simp) === str_inc

  @test Oscar.fine_strand(ctx, [1, 0, 0, -1]) === Oscar.fine_strand(ctx, [5, 0, 4, -70])
  trunc_str = ctx[10, G[1]]
  @test !is_zero(map(trunc_str, 0))

  # test the structure sheaf 
  ctx = Oscar.NewToricCtx(X);
  z = zero(G)
  str = ctx[0, z]
  [ngens(str[i]) for i in 0:-1:-2]

  str_simp = Oscar.simplified_strand(ctx, 0, z)
  betti_nums = [ngens(str_simp[i]) for i in 0:-1:-2]
  @test isone(first(betti_nums)) && all(is_zero, betti_nums[2:end])


  simp_trunc_str = Oscar.simplified_strand(ctx, 10, G[1])
  Oscar.induced_cohomology_map(ctx, 1, 2, G[1], 0)
  Oscar.simplified_strand_homotopy(ctx, 3, G[1], -1)
  Oscar.simplified_strand_inclusion(ctx, 3, G[1], -1)

  # test compatibility with the already established code 
  ctx[1, 2, G[1]]
  Oscar.multiplication_map(ctx, S[1], 1, G[1], 0)

  ctx_old = Oscar.ToricCtx(X);
  for _ in 1:100
    g = sum(rand(-5:5)*g for g in gens(G); init=zero(G))
    e = rand(0:10)
    ee = [e for _ in 1:ngens(irrelevant_ideal(X))]
    s1 = Oscar.simplified_strand(ctx, e, g)
    s2 = Oscar.simplified_strand(ctx_old, ee, g)
    @test all(ngens(s1[k]) == ngens(s2[k]) for k in 0:-1:-dim(X))
  end


  # Test computation of support sets in easy cases. 
  IP = projective_space(NormalToricVariety, 2)
  S = cox_ring(IP)
  G = grading_group(S)
  g = G[1]
  ctx = Oscar.NewToricCtx(IP);
  supps = Dict{Int, Vector{Int}}(i => Oscar.support_set(ctx, i) for i in 0:dim(IP))

  @test Oscar.optimal_k(ctx, 0, zero(G)) == 0
  @test Oscar.optimal_k(ctx, 1, -g) == 0
  @test Oscar.optimal_k(ctx, 2, -3*g) == 1
  @test Oscar.optimal_k(ctx, 2, -4*g) == 2
end
