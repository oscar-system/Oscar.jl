@testset "Graded Modules 1" begin
  Qx, g = polynomial_ring(QQ, 3)
  R = grade(Qx, [1,2,3])[1]
  F = FreeModule(R, [i^2*R.D[1] for i=1:4])
  s = sub(F, [R(g[1])^2*gen(F, 1), R(g[2])*gen(F, 1)])
  q = quo(s, [R(g[1])^4*gen(F, 1), R(g[2])^2*gen(F, 1)]) 
  @test length(presentation(q)) == 3

  h = hom(F,F, [R(g[1]), R(g[2]), R(g[2]), R(g[3])] .* basis(F)) 

  @test length([h(g) for g = basis(F)]) == dim(F)

  @test !Oscar.is_homogeneous(h)

  c = homogeneous_components(h)
  @test length(c) == 3

  for g = keys(c)
    @test degree(c[g]) == g
  end

  H, phi = hom(F, F)
  
  phi(rand(gens(H)))

  H, mH = hom(q,q )

  free_resolution(H)

  homogeneous_component(H, grading_group(F)[0])
end

@testset "Graded Modules 2" begin

  D = free_abelian_group(2)
  Qx, x, y = polynomial_ring(QQ, :x=>1:2, :y=>1:2)
  R = grade(Qx, [D[1], D[1], D[2], D[2]])[1]

  f = x[1]^3-x[2]^3
  g = y[1]^5+y[2]^5
  h = f+g

  F = FreeModule(R, 1)

  s = sub(F, [derivative(h, i)*F[1] for i=1:4])

  free_resolution(s)

  H, mH = hom(s, quo(F, s))

  free_resolution(H)

  tensor_product(s, s, task = :none)

end

@testset "Graded modules kernel" begin
  R, (x,y) = polynomial_ring(QQ, ["x", "y"])
  R = grade(R, [0,0])

  F1 = FreeModule(R, 1)
  M1 = quo(F1, [x*F1[1], y*F1[1]])

  F2 = FreeModule(R, 2)
  M2 = quo(F2, [x*F2[1], y*F2[1], x*F2[2], y*F2[2]])

  phi = hom(M2,M1, [y*M1[1], x*M1[1]]) # zero-morphism
  K = kernel(phi)[1]
  @test (iszero(quo(M2,K)) && iszero(quo(K,M2))) # mathematical equality
end

@testset "all_monomials for graded free modules" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]
  S, (x, y, u, v, w) = grade(R)

  F = graded_free_module(S, [-1, 2])

  for d in -4:4
    amm = Oscar.all_monomials(F, d)
    length(amm)
    eltype(amm)
    @test d < -1 || !isempty(amm)
    @test all(x->degree(x) == grading_group(F)([d]), amm)
    ae = Oscar.all_exponents(F, d)
    @test d < -1 || !isempty(ae)
    @test all(x->sum(x[1]) + Int(degree(F[x[2]])[1]) == d, ae)
  end
end

#= Disabled for the moment, but continued soon.
@testset "monomials of subquos" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])

  S1 = graded_free_module(S, [0])
  I = ideal(S, [u^2 for u in gens(S)])
  IS1, inc = I*S1
  M = cokernel(inc)

  a = Oscar.AllSubquoMonomials(M, 3)
  b = Oscar.all_exponents(M, 3)

  v = collect(a)
  @test length(v) == length(a) == 1 == length(collect(b))
  @test !(M(x^3 * S1[1]) in v)
  @test M(x*y*z * S1[1]) in v

  I = ideal(S, [u^3 for u in gens(S)])
  IS1, inc = I*S1
  M = cokernel(inc)

  a = Oscar.AllSubquoMonomials(M, 3)
  b = Oscar.all_exponents(M, 3)

  v = collect(a)
  @test length(v) == length(a) == length(collect(b))
  @test !(M(x^3 * S1[1]) in v)
  @test M(x*y*z * S1[1]) in v
  @test M(x^2*y * S1[1]) in v

  J, _ = sub(S1, [x*y*z*S1[1]])
  I = ideal(S, [u^4 for u in gens(S)])
  IS1, inc = I*S1
  M, _ = quo(J, IS1)
  a = Oscar.AllSubquoMonomials(M, 4)
  b = Oscar.all_exponents(M, 4)
  @test length(collect(a)) == length(a) == 3 == length(collect(b))

  a = Oscar.AllSubquoMonomials(M, 6)
  b = Oscar.all_exponents(M, 6)
  @test length(collect(a)) == length(a) == 7 == length(collect(b))
end
=#

@testset "simplification of graded subquos, issue #3108" begin
  X = rational_d10_pi9_quart_1();
  I = defining_ideal(X);
  Pn = base_ring(I)
  n = ngens(Pn)-1
  c = codim(I)
  FI = free_resolution(I)
  FIC = FI.C;
  r = range(FIC)
  C = shift(FIC[first(r):-1:1], -c)
  F = free_module(Pn, 1)
  OmegaPn = grade(F, [n+1])
  D = hom(C, OmegaPn)
  Omega = homology(D, 0);
  is_graded(Omega)
  SOmega, b = simplify(Omega)
  a = get_attribute(b, :inverse)
  @test is_graded(SOmega)
  @test is_isomorphism(a)
  @test is_isomorphism(b)
  M, iso = forget_grading(Omega)
  @test is_isomorphism(iso)
  @test is_isomorphism(inverse(iso))
  M, _ = forget_grading(Omega)
  prune_with_map(M)
end
