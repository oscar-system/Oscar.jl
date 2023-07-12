@testset "Ideal sheaves" begin
  kk = GF(29)

  # Set up the base ℙ¹ with coordinates s and t
  R, (s,t) = polynomial_ring(kk, ["s", "t"])
  S, _ = grade(R, [1, 1])

  base_P1 = ProjectiveScheme(S)

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
  add_glueing!(C, Glueing(X, Y, restrict(f, domain(f), domain(g)), restrict(g, domain(g), domain(f))))

  # Extend the glueing to the whole covered scheme
  fill_transitions!(C)

  X = CoveredScheme(C)

  U = C[1]
  x = gens(ambient_coordinate_ring(U))
  I = IdealSheaf(X, U, OO(U).([x[1]-1, x[2]-2, x[3]-3]))
  J = IdealSheaf(X, U, OO(U).([x[1]-5, x[2]-1, x[3]]))

  @test I+J == IdealSheaf(X, U, [one(OO(U))])
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
  ID = IdDict{AbsSpec, Oscar.Ideal}()
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
  @test !haskey(pbII.I.obj_cache, U[2])
  @test pbII(U[2]) isa Ideal
  @test haskey(pbII.I.obj_cache, U[2])
end

@testset "colength of ideal sheaves" begin
  P3 = projective_space(QQ, 3)
  X = covered_scheme(P3)
  S = homogeneous_coordinate_ring(P3)
  (x,y, z, w) = gens(S)
  I = ideal(S, [x^3+y^3+z^3+w^3, x+y+z+w, 37*x^2-x*z+4*w^2-5*w*z])
  II = ideal_sheaf(P3, I)
  @test oscar.colength(II) == 6
  @test oscar.colength(II, covering=oscar.simplified_covering(X)) == 6
  J = ideal(S, [x, y, z^4])
  JJ = ideal_sheaf(P3, J)
  @test oscar.colength(JJ) == 4
  @test oscar.colength(JJ, covering=oscar.simplified_covering(X)) == 4
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
  C = oscar._separate_disjoint_components(II, covering=oscar.simplified_covering(X))
  for U in patches(C)
    @test sum(isone(I(U)) for I in II) == 6
  end
  C = oscar._separate_disjoint_components(II)
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
  JJ = oscar.maximal_associated_points(ideal_sheaf(X, J))
  X = covered_scheme(X)
  C = oscar._separate_disjoint_components(JJ, covering=oscar.simplified_covering(X))
  for U in patches(C)
    @test sum(!isone(I(U)) for I in JJ) == 1
  end

  CC = oscar._one_patch_per_component(C, JJ)
  for P in JJ
    @test isone(sum(!isone(P(U)) for U in patches(CC)))
    @show sum(!isone(P(U)) for U in patches(CC))
  end
end
