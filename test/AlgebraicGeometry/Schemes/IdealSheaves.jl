@testset "Ideal sheaves" begin
  kk = GF(29)

  # Set up the base ℙ¹ with coordinates s and t
  S, _ = graded_polynomial_ring(kk, [:s, :t])

  base_P1 = proj(S)

  # split this into the standard covering
  base_covering = standard_covering(base_P1)

  A1s = patches(base_covering)[1]
  A1t = patches(base_covering)[2]

  # Set up relative projective space of relative dimension 2
  # over both base patches
  P2_s = projective_space(OO(A1s), ["xs", "ys", "zs"])

  Cs = standard_covering(P2_s)

  P2_t = projective_space(OO(A1t), ["xt", "yt", "zt"])

  Ct = standard_covering(P2_t)

  # Join the resulting schemes in a disjoint union with two
  # components
  C = disjoint_union(Cs, Ct)

  # Manually glue two (dense) patches of the two components
  X = Cs[3]
  Y = Ct[3]
  x = gens(OO(X))
  y = gens(OO(Y))
  f = maximal_extension(X, Y, [x[1]//(x[3])^4, x[2]//(x[3])^6, 1//x[3]])
  g = maximal_extension(Y, X, [y[1]//(y[3])^4, y[2]//(y[3])^6, 1//y[3]])
  add_gluing!(C, Gluing(X, Y, restrict(f, domain(f), domain(g)), restrict(g, domain(g), domain(f))))

  # Extend the gluing to the whole covered scheme
  Oscar.fill_transitions!(C)

  X = CoveredScheme(C)

  U = C[1]
  x = gens(ambient_coordinate_ring(U))
  I = IdealSheaf(X, U, OO(U).([x[1]-1, x[2]-2, x[3]-3]))
  J = IdealSheaf(X, U, OO(U).([x[1]-5, x[2]-1, x[3]]))

  @test I+J == Oscar.unit_ideal_sheaf(X)
  K = I*J
end

@testset "ideal sheaves II" begin
  IP2 = projective_space(QQ, 2)
  X = covered_scheme(IP2)
  S = homogeneous_coordinate_ring(IP2)
  (u,v,w) = gens(S)
  Ihom = ideal(S, u^2 - v*w)
  I = IdealSheaf(IP2, Ihom)
  @test is_prime(I)
  I2 = IdealSheaf(IP2, gens(Ihom))
  I3 = IdealSheaf(IP2, gen(Ihom, 1))
  @test I == I2 == I3
  U = patches(default_covering(X))
  @test I(U[1]) isa Oscar.Ideal
  @test I(U[2]) isa Oscar.Ideal
  @test I(U[3]) isa Oscar.Ideal
  V = PrincipalOpenSubset(U[1], gens(OO(U[1]))[1])
  rho = I(U[1], V)
  @test I(V) == ideal(OO(V), rho.(gens(I(U[1]))))
  Y = subscheme(I)
  simplify!(I)

  # run the check block in extend!
  ID = IdDict{AbsAffineScheme, Oscar.Ideal}()
  ID[X[1][1]] = I(X[1][1])
  J = IdealSheaf(X, extend!(default_covering(X), ID), check=true)
  ID[X[1][1]] = I(X[1][1])
  J2 = IdealSheaf(X, ID)
  @test J == J2
  @test sprint(show, J) isa String

  @test issubset(IdealSheaf(X), I) # Whether the zero ideal sheaf is a subset of I
end

@testset "pullbacks of ideal sheaves" begin
  P = projective_space(QQ, 2)
  SP = homogeneous_coordinate_ring(P)
  (x, y, z) = gens(SP)
  m = ideal(SP, gens(SP))
  m3 = m^3
  n = ngens(m3)
  Q = projective_space(QQ, n-1)
  SQ = homogeneous_coordinate_ring(Q)
  phi = hom(SQ, SP, gens(m3))
  f = ProjectiveSchemeMor(P, Q, phi)
  f_cov = covered_scheme_morphism(f)

  A = [1 2 4; 1 3 9; 1 5 25]
  v = A*[x, y, z]
  psi = hom(SP, SP, v)
  g = ProjectiveSchemeMor(P, P, psi)
  g_cov = covered_scheme_morphism(g)

  II = IdealSheaf(Q, [gen(SQ, 1)+gen(SQ, 2)])
  pbII = pullback(f_cov)(II)
  X = covered_scheme(P)
  U = affine_charts(X)
  @test pbII(U[1]) isa Ideal
  @test !haskey(Oscar.object_cache(pbII), U[2])
  @test pbII(U[2]) isa Ideal
  @test haskey(Oscar.object_cache(pbII), U[2])
end

@testset "colength of ideal sheaves" begin
  P3 = projective_space(QQ, 3)
  X = covered_scheme(P3)
  S = homogeneous_coordinate_ring(P3)
  (x,y, z, w) = gens(S)
  I = ideal(S, [x^3+y^3+z^3+w^3, x+y+z+w, 37*x^2-x*z+4*w^2-5*w*z])
  II = ideal_sheaf(P3, I)
  @test Oscar.colength(II) == 6
  @test Oscar.colength(II, covering=Oscar.simplified_covering(X)) == 6
  J = ideal(S, [x, y, z^4])
  JJ = ideal_sheaf(P3, J)
  @test Oscar.colength(JJ) == 4
  @test Oscar.colength(JJ, covering=Oscar.simplified_covering(X)) == 4
end

@testset "separation of ideal sheaves" begin
  P2 = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(P2)
  (x, y, z) = gens(S)
  I = Vector{Ideal}()
  push!(I, ideal(S, [x, y]))
  push!(I, ideal(S, [y, z]))
  push!(I, ideal(S, [x, z]))
  push!(I, ideal(S, [x+y, z]))
  push!(I, ideal(S, [x+z, y]))
  push!(I, ideal(S, [y+z, x]))
  push!(I, ideal(S, [y+2*z, x]))

  II = [IdealSheaf(P2, a) for a in I]
  X = covered_scheme(P2)
  C = Oscar._separate_disjoint_components(II, covering=Oscar.simplified_covering(X))
  for U in patches(C)
    @test sum(isone(I(U)) for I in II) == 6
  end
  C = Oscar._separate_disjoint_components(II)
  for U in patches(C)
    @test sum(isone(I(U)) for I in II) == 6
  end

  P3 = projective_space(QQ, 3)

  S = homogeneous_coordinate_ring(P3)
  (x, y, z, w) = gens(S)
  I = ideal(S, x^4 + y^4 + z^4 + w^4)
  X = subscheme(P3, I)
  S = homogeneous_coordinate_ring(X)
  (x, y, z, w) = gens(S)
  J = ideal(S, [x^2 - 3*y^2 + 4*y*z - 5*w^2 + 3*x*w, x+y+z+w])
  J = J*ideal(S, [5*x + 9*y - 5*z + 3*w, x+8*y+z+w])
  J = J*ideal(S, [-9*x + 3*y - 5*z + w, x+8*y+15*z+w])
  #J = ideal(S, [x^2 - 3*y^2 + 4*y*z - 5*w^2 + 3*x*w, 25*x^2*y + y^2*z + z^2*w + w^2*x])
  JJ = Oscar.maximal_associated_points(ideal_sheaf(X, J))

  X = covered_scheme(X)
  C = Oscar._separate_disjoint_components(JJ, covering=Oscar.simplified_covering(X))
  for U in patches(C)
    @test sum(!isone(I(U)) for I in JJ) == 1
  end

  CC = Oscar._one_patch_per_component(C, JJ)
  for P in JJ
    @test isone(sum(!isone(P(U)) for U in patches(CC)))
  end
end

@testset "saturation of ideal sheaves" begin
  IP2 = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(IP2)
  (x, y, z) = gens(S)

  I = IdealSheaf(IP2, [x+y])
  J = IdealSheaf(IP2, [z^2])

  @test J == saturation(I*J, I)
end

@testset "pushforward of ideal sheaves" begin
  IP2 = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(IP2)
  (x, y, z) = gens(S)
  II = IdealSheaf(IP2, ideal(S, [x, y]))
  bl = blow_up(II)
  E = ideal_sheaf(exceptional_divisor(bl))
  JJ = pushforward(bl, E)
  @test JJ == II
end

@testset "subschemes of ideal sheaves" begin
  IP2 = projective_space(QQ, [:x, :y, :z])

  X = covered_scheme(IP2)
  S = homogeneous_coordinate_ring(IP2)
  x, y, z = gens(S)
  I = ideal(S, [x*y, x*z, y*z])
  II = IdealSheaf(IP2, I)
  H = IdealSheaf(IP2, ideal(S, z))
  Y = subscheme(II + H)
  U = affine_charts(Y)
  @test length(U) == 2
  f, g = gluing_morphisms(Y[1][U[1], U[2]])
  A = domain(f)
  B = domain(g)
  @test isempty(A)
  @test isempty(B)
end

@testset "maximal associated points" begin
  IA3 = affine_space(QQ, [:x, :y, :z])
  R = OO(IA3)
  (x, y, z) = gens(R)

  I = ideal(R, z*(z-1))
  X, inc_X = sub(IA3, I)

  A = OO(X)
  J = ideal(A, [x, y])*ideal(A, [y-1, z-1])

  JJ = IdealSheaf(X, J)

  @test dim(JJ) == 1
  comp = Oscar.maximal_associated_points(JJ)
  @test length(comp) == 3
  @test 1 in dim.(comp)
  @test 0 in dim.(comp)

  Z = ideal(R, [x, y, z])

  bl = blow_up(IA3, Z)

  str_JJ = strict_transform(bl, IdealSheaf(IA3, pushforward(inc_X, J), covered_scheme=codomain(bl)))

  comp2 = Oscar.maximal_associated_points(str_JJ)
  @assert length(comp2) == 2
  @test 1 in dim.(comp2)
  @test 0 in dim.(comp2)
end

@testset "production of ideals" begin
  IP2 = projective_space(QQ, [:x, :y, :z])
  P2 = covered_scheme(IP2)
  U = first(affine_charts(P2))
  y, z = gens(OO(U))
  I = ideal(OO(U), y-z)
  II = IdealSheaf(P2, U, I)
  Oscar.produce_object(II, affine_charts(P2)[2]; algorithm=:pullback)
  Oscar.produce_object(II, affine_charts(P2)[3]; algorithm=:pushforward)
end

