@testset "kernels via flattenings" begin
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z], cached=false)
  S, _ = polynomial_ring(R, [:s, :t], cached=false)
  phi = Oscar.flatten(S)
  @test vcat(phi.(gens(S)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))

  f = hom(S, S, gens(S))
  K = kernel(f)
  @test zero(S) in K

  f = hom(S, S, hom(R, R, [first(gens(R)), first(gens(R)), last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gen(R, 1) - gen(R, 2)) in K

  R, _ = quo(R, ideal(R, sum(gens(R))))
  S, _ = polynomial_ring(R, [:s, :t])
  phi = Oscar.flatten(S; cached=true)
  @test vcat(phi.(gens(S)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))

  f = hom(S, S, gens(S))
  K = kernel(f)
  @test zero(S) in K

  f = hom(S, S, hom(R, R, [first(gens(R)) + gen(R, 2), zero(R), last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gen(R, 2)) in K

  R = base_ring(R)
  U = powers_of_element(sum(gens(R)))
  R, _ = localization(R, U)
  S, _ = polynomial_ring(R, [:s, :t])
  phi = Oscar.flatten(S)
  @test vcat(phi.(gens(S)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))

  f = hom(S, S, gens(S))
  @test Oscar.flatten(f) isa Map
  K = kernel(f) 
  @test zero(S) in K 
 
  f = hom(S, S, identity, [S[1], S[1]])
  K = kernel(f)
  @test S[1] - S[2] in K
 
  I = ideal(R, first(gens(R))^2)
  R, _ = quo(R, I)
  S, _ = polynomial_ring(R, [:s, :t])
  phi = Oscar.flatten(S)
  @test vcat(phi.(gens(S)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))
 
  f = hom(S, S, [S[1], S[1]])
  K = kernel(f)
  @test S[1] - S[2] in K
 
  f = hom(S, S, hom(R, R, [zero(R), first(gens(R)) + gen(R, 2), last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gen(R, 1)) in K
end

@testset "cross kernel computations" begin
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z], cached=false)
  S, (s, t) = polynomial_ring(R, [:s, :t], cached=false)
  
  f = hom(R, S, [s, s, t])
  @test x-y in kernel(f)

  g = hom(S, R, identity, [x, y])
  @test s-x in kernel(g)

  h = hom(S, R, hom(R, R, [R[2], R[1], R[2]]), [x, y])
  @test s-y in kernel(h)
  @test S(x-z) in kernel(h)
end

@testset "flattenings of quotient rings" begin
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z], cached=false)

  S, _ = polynomial_ring(R, [:s, :t], cached=false)
  Q, _ = quo(S, ideal(S, S[1]))
  phi = Oscar.flatten(Q)
  @test vcat(phi.(gens(Q)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(Q)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(Q.(gens(R)))

  U = powers_of_element(x)
  L, _ = localization(R, U)
  S, _ = polynomial_ring(L, [:s, :t])
  Q, _ = quo(S, ideal(S, S[1]))
  phi = Oscar.flatten(Q);
  @test vcat(phi.(gens(Q)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(Q)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(Q.(S.(gens(R))))
  
  A, _ = quo(R, ideal(R, y))
  S, _ = polynomial_ring(A, [:s, :t])
  Q, _ = quo(S, ideal(S, S[1]))
  phi = Oscar.flatten(Q);
  @test vcat(phi.(gens(Q)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(Q)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(Q.(S.(gens(R))))

  LA, _ = localization(A, U)
  S, _ = polynomial_ring(LA, [:s, :t])
  Q, _ = quo(S, ideal(S, S[1]))
  phi = Oscar.flatten(Q);
  @test vcat(phi.(gens(Q)), Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(Q)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test Oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(Q.(S.(gens(R))))
end

@testset "flattenings of graded rings" begin
  # towers of polynomial rings
  R, (x, y) = polynomial_ring(QQ, [:x, :y], cached=false)
  S, (u, v) = grade(polynomial_ring(R, [:u, :v], cached=false)[1])

  phi = hom(R, S, [u, v])
  @test phi(x) == u
  
  flat = Oscar.flatten(S)
  S_flat = codomain(flat)

  @test is_graded(S_flat)
  G = grading_group(S_flat)

  @test degree.(gens(S_flat)) == [G[1], G[1], zero(G), zero(G)]
  @test default_ordering(S_flat) == degrevlex(gens(S_flat)[1:2])*degrevlex(gens(S_flat)[3:4])

  I = ideal(S, [u-v])
  I_flat = flat(I)
  @test is_graded(I_flat)
  @test !(u in I)
  @test u-v in I

  # graded rings over quotient rings
  A, _ = quo(R, [x-y])
  S, (u, v) = grade(A[:u, :v][1])

  flat = Oscar.flatten(S)

  S_flat = codomain(flat)
  S_poly = base_ring(S_flat)

  @test default_ordering(S_flat) == degrevlex(gens(S_poly)[1:2])*degrevlex(gens(S_poly)[3:4])

  I = ideal(S, [u*x])
  I_flat = flat(I)
  @test !(u in I)
  @test u*y in I

  # graded rings over localizations of quotient rings
  U = powers_of_element(x^2 + y^2)
  W, _ = localization(A, U)
  S, (u, v) = graded_polynomial_ring(W, [:u, :v])

  flat = Oscar.flatten(S)

  S_flat = codomain(flat)
  S_poly = base_ring(S_flat)

  @test default_ordering(S_poly) == degrevlex(gens(S_poly)[1:2])*degrevlex(gens(S_poly)[3:4])

  I = ideal(S, [u*x])
  I_flat = flat(I)
  @test u in I
  @test !(v in I)

  # graded rings over localizations of polynomial rings 
  L, _ = localization(R, U)
  S, (u, v) = graded_polynomial_ring(L, [:u, :v])

  flat = Oscar.flatten(S)

  S_flat = codomain(flat)
  S_poly = base_ring(S_flat)

  @test default_ordering(S_poly) == degrevlex(gens(S_poly)[1:2])*degrevlex(gens(S_poly)[3:4])

  I = ideal(S, [u*(x^2 + y^2)])
  I_flat = flat(I)
  @test u in I
  @test !(v in I)
end

@testset "flattenings of modules" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y], cached=false)
  S, (u, v) = grade(polynomial_ring(R, [:u, :v], cached=false)[1])

  flat = Oscar.flatten(S)

  S_flat = codomain(flat)

  F = graded_free_module(S, [-2, 1])
  F_flat, iso, iso_inv = Oscar.flatten(F)

  M, _ = sub(F, [F[1]])
  M_flat, iso_M, iso_M_inv = Oscar.flatten(M)
  @test is_graded(M_flat)
  @test ambient_free_module(M_flat) === F_flat
  @test F[1] in M
  @test coordinates(u*F[1], M)[1] == u

  @test is_graded(F_flat)

  @test iso_inv(iso(F[1])) == F[1]
  @test F_flat === Oscar.flatten(F)[1]

  A, _ = quo(R, [x-y])

  S, (u, v) = grade(A[:u, :v][1])

  flat = Oscar.flatten(S)

  S_flat = codomain(flat)

  F = graded_free_module(S, [-2, 1])
  F_flat, iso, iso_inv = Oscar.flatten(F)

  M, _ = sub(F, [F[1]])
  M_flat, iso_M, iso_M_inv = Oscar.flatten(M)
  @test is_graded(M_flat)
  @test ambient_free_module(M_flat) === F_flat
  @test F[1] in M
  @test coordinates(u*F[1], M)[1] == u

  @test iso_inv(iso(F[1])) == F[1]
  @test F_flat === Oscar.flatten(F)[1]
  @test is_graded(F_flat)

  U = powers_of_element(y)
  L, _ = localization(R, U)
  S, (u, v) = grade(L[:u, :v][1])

  flat = Oscar.flatten(S)
  S_flat = codomain(flat)
  U, V, X, Y = gens(S_flat)

  @test degree(U) == degree(u)
  @test iszero(degree(X))

  F = graded_free_module(S, [-2, 1])
  F_flat, iso, iso_inv = Oscar.flatten(F)

  M, _ = sub(F, [F[1]])
  M_flat, iso_M, iso_M_inv = Oscar.flatten(M)
  @test is_graded(M_flat)
  @test ambient_free_module(M_flat) === F_flat
  @test F[1] in M
  @test coordinates(u*F[1], M)[1] == u

  @test iso_inv(iso(F[1])) == F[1]
  @test F_flat === Oscar.flatten(F)[1]
  @test is_graded(F_flat)

  U = powers_of_element(y)
  W, _ = localization(A, U)
  S, (u, v) = grade(W[:u, :v][1])

  flat = Oscar.flatten(S)
  S_flat = codomain(flat)
  U, V, X, Y = gens(S_flat)

  @test degree(U) == degree(u)
  @test iszero(degree(X))

  F = graded_free_module(S, [-2, 1])
  F_flat, iso, iso_inv = Oscar.flatten(F)

  M, _ = sub(F, [F[1]])
  M_flat, iso_M, iso_M_inv = Oscar.flatten(M)
  @test is_graded(M_flat)
  @test ambient_free_module(M_flat) === F_flat
  @test F[1] in M
  @test coordinates(u*F[1], M)[1] == u

  @test iso_inv(iso(F[1])) == F[1]
  @test F_flat === Oscar.flatten(F)[1]
  @test is_graded(F_flat)
end

@testset "free (graded) resolutions" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y], cached=false)
  S, (u, v) = grade(polynomial_ring(R, [:u, :v], cached=false)[1])

  flat = Oscar.flatten(S)

  S_flat = codomain(flat)

  F = graded_free_module(S, [-2, 1])
  F_flat, iso, iso_inv = Oscar.flatten(F)

  M, _ = quo(F, [u^5*F[1], v*F[1]])
  M_flat, iso, iso_inv = Oscar.flatten(M)
  res = free_resolution(M)
  @test Oscar._regularity_bound(M) == 4
end

@testset "garbage collection" begin
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z], cached=false)
  S, _ = polynomial_ring(R, [:s, :t], cached=false)
  S, _ = quo(S, ideal(S, elem_type(S)[]))
  phi = Oscar.flatten(S; cached=true)
  # julia might not consider `x^2+y^2` as being unused until the end of the scope
  # see https://github.com/JuliaLang/julia/issues/51818
  # so we put this into a function
  function use_it(phi)
    a = x^2 + y^2
    b = phi(a)
    @test phi(a) === b
    a = nothing
  end
  use_it(phi)
  GC.gc()
  @test isempty(keys(Oscar.flat_counterparts(phi)))
end
  
