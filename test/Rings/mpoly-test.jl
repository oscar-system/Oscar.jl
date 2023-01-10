@testset "Polynomial ring constructor" begin
  R, x = @inferred PolynomialRing(QQ, "x" => 1:3)
  @test length(x) == 3
  @test x isa Vector
  @test x == gens(R)
  for i in 1:3
    @test sprint(show, "text/plain", x[i]) == "x[$i]"
  end

  R, x = @inferred PolynomialRing(QQ, "x" => (1:3, 1:4))
  @test length(x) == 12
  @test x isa Matrix{<: Any}
  @test size(x) == (3, 4)
  @test length(unique(x)) == 12
  @test Set(gens(R)) == Set(x)
  for i in 1:3
    for j in 1:4
      @test sprint(show, "text/plain", x[i, j]) == "x[$i, $j]"
    end
  end

  R, x, y = @inferred PolynomialRing(QQ, "x" => (1:3, 1:4), "y" => 1:2)
  @test length(x) == 12
  @test x isa Matrix{<: Any}
  @test size(x) == (3, 4)
  @test length(unique(x)) == 12
  for i in 1:3
    for j in 1:4
      @test sprint(show, "text/plain", x[i, j]) == "x[$i, $j]"
    end
  end
  @test y isa Vector{<: Any}
  @test length(y) == 2
  @test length(unique(y)) == 2
  for i in 1:2
    @test sprint(show, "text/plain", y[i]) == "y[$i]"
  end
  @test Set(gens(R)) == union(Set(x), Set(y))

  R, x, y, z = @inferred PolynomialRing(QQ, "x" => (1:3, 1:4),
                                            "y" => 1:2,
                                            "z" => (1:1, 1:1, 1:1))
  @test length(x) == 12
  @test x isa Matrix{<: Any}
  @test size(x) == (3, 4)
  @test length(unique(x)) == 12
  for i in 1:3
    for j in 1:4
      @test sprint(show, "text/plain", x[i, j]) == "x[$i, $j]"
    end
  end
  @test y isa Vector{<: Any}
  @test length(y) == 2
  @test length(unique(y)) == 2
  for i in 1:2
    @test sprint(show, "text/plain", y[i]) == "y[$i]"
  end
  @test z isa Array{<: Any, 3}
  @test sprint(show, "text/plain", z[1, 1, 1]) == "z[1, 1, 1]"
  @test Set(gens(R)) == union(Set(x), Set(y), Set(z))
end

@testset "Polynomial homs" begin
  R, (x, y) = PolynomialRing(QQ, ["x", "y"])
  I1 = x^2 + y^2
  I2 = x^2 * y^2
  I3 = x^3*y - x*y^3
  S, (a,b,c) = PolynomialRing(QQ, ["a", "b", "c"])
  h = hom(S, R, [I1, I2, I3])
  @test kernel(h) == ideal(S, [a^2*b - 4*b^2 - c^2])
  @test h(gen(S, 1)) == I1
  @test h.([a,b]) == [I1, I2]
  @test preimage(h, ideal(R, [I2, I3])) == ideal(S, [b, c])
end

@testset "Ideal operations" begin
  R, (x, y) = PolynomialRing(QQ, ["x", "y"])
  I = ideal(R, [x^2, y^2])
  J = ideal(R, [x*y^2, x^2])
  S = ideal(R, [x^2, y^2, x*y^2])
  P = ideal(R, [x^3*y^2, x^4, x*y^4, x^2*y^2])
  @test Oscar.check_base_rings(I, J) == nothing
  @test I+J == I-J == S
  @test I*J == P
  f = x^2 + y^2
  g = x^4*y - x*y^3
  I = [f, g]
  S, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"])
  J = ideal(S, [(c^2+1)*(c^3+2)^2, b-c^2])
  @test_throws ErrorException Oscar.check_base_rings(P, J)
  r1 = c^2-b
  r2 = b^2*c+c^3+2*c^2+2
  L = gens(radical(J))

  @test jacobi_ideal(f) == ideal(R, [2*x, 2*y])
  @test jacobi_matrix(f) == matrix(R, 2, 1, [2*x, 2*y])
  @test jacobi_matrix(I) == matrix(R, 2, 2, [2*x, 4*x^3*y-y^3, 2*y, x^4-3*x*y^2])
  @test length(L) == 2
  @test length(findall(x->x==r1, L)) == 1
  @test length(findall(x->x==r2, L)) == 1

  @test issubset(ideal(S, [a]), ideal(S, [a]))
  @test issubset(ideal(S, [a]), ideal(S, [a, b]))
  @test !issubset(ideal(S, [c]), ideal(S, [b]))
  @test !issubset(ideal(S, [a, b, c]), ideal(S, [a*b*c]))
end

@testset "Primary decomposition" begin

  # primary_decomposition
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [x, y*z^2])
  for method in (:GTZ, :SY)
    j = ideal(R, [R(1)])
    for (q, p) in primary_decomposition(i, alg=method)
      j = intersect(j, q)
      @test is_primary(q)
      @test is_prime(p)
      @test p == radical(q)
    end
    @test j == i
  end

  R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
  i = ideal(R, [9, (a+3)*(b+3)])
  l = primary_decomposition(i)
  @test length(l) == 2

  # minimal_primes
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  i1 = ideal(R, [z^3+2, -z^2+y])
  i2 = ideal(R, [z^2+1, -z^2+y])
  l = minimal_primes(i)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  l = minimal_primes(i, alg=:charSets)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
  i = ideal(R, [R(9), (a+3)*(b+3)])
  i1 = ideal(R, [R(3), a])
  i2 = ideal(R, [R(3), b])
  l = minimal_primes(i)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  # equidimensional_decomposition_weak
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]),
                ideal(R, [x^2, z^2]), ideal(R, [x^5, y^5, z^5]))
  l = equidimensional_decomposition_weak(i)
  @test length(l) == 3
  @test l[1] == ideal(R, [z^4, y^5, x^5, x^3*z^3, x^4*y^4])
  @test l[2] == ideal(R, [y*z, x*z, x^2])
  @test l[3] == ideal(R, [z])

 # equidimensional_decomposition_radical
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  l = equidimensional_decomposition_radical(i)
  @test length(l) == 1
  @test l[1] == ideal(R, [z^2-y, y^2*z+z^3+2*z^2+2])
  
  # equidimensional_hull
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]),
                ideal(R, [x^2, z^2]), ideal(R, [x^5, y^5, z^5]))
  @test equidimensional_hull(i) == ideal(R, [z])

  R, (a, b, c, d) = PolynomialRing(ZZ, ["a", "b", "c", "d"])
  i = intersect(ideal(R, [R(9), a, b]),
                ideal(R, [R(3), c]),
                ideal(R, [R(11), 2*a, 7*b]),
                ideal(R, [13*a^2, 17*b^4]),
                ideal(R, [9*c^5, 6*d^5]),
                ideal(R, [R(17), a^15, b^15, c^15, d^15]))
  @test equidimensional_hull(i) == ideal(R, [R(3)])

  # equidimensional_hull_radical
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  @test equidimensional_hull_radical(i) == ideal(R, [z^2-y, y^2*z+z^3+2*z^2+2])

  # absolute_primary_decomposition
  R,(x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
  I = ideal(R, [(z+1)*(z^2+1)*(z^3+2)^2, x-y*z^2])
  d = absolute_primary_decomposition(I)
  @test length(d) == 3

  R,(x,y,z) = GradedPolynomialRing(QQ, ["x", "y", "z"])
  I = ideal(R, [(z+y)*(z^2+y^2)*(z^3+2*y^3)^2, x^3-y*z^2])
  d = absolute_primary_decomposition(I)
  @test length(d) == 5
end

@testset "Groebner" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  I = ideal([x+y^2+4*z-5,x+y*z+5*z-2])
  @test gens(groebner_basis(I, ordering=lex(gens(R)))) == [y^2 - y*z - z - 3, x + y*z + 5*z - 2]
  @test gens(groebner_basis(I, ordering=degrevlex(gens(R)), complete_reduction = true)) == [x + y*z + 5*z - 2, x + y^2 + 4*z - 5, x*y - x*z - 5*x - 2*y - 4*z^2 - 20*z + 10, x^2 + x*z^2 + 10*x*z - 4*x + 4*z^3 + 20*z^2 - 20*z + 4]

  # Test coefficient rings that are actually fields for safety. The first three
  # are native to singular while gfp_fmpz_elem now has a proper wrapper
  for Zn in [ResidueRing(ZZ, 11), ResidueRing(ZZ, fmpz(10)^50+151), GF(11),
             GF(fmpz(10)^50+151)]
    R, (x, y) = PolynomialRing(Zn, ["x", "y"], ordering = :degrevlex)
    l = [x*y+x^3+1, x*y^2+x^2+1]
    g = gens(groebner_basis(ideal(R, l); ordering = degrevlex(gens(R))))
    @test iszero(divrem(l[1] + l[2], g)[2])
  end

  F, a = FiniteField(11, 2, "a")
  R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"], ordering = :degrevlex)
  l = [3*x^5 + a*x*y^2 + a^2*z^2, z^3*x^2 + 7*y^3 + z]
  gb = gens(groebner_basis(ideal(R, l); ordering = degrevlex(gens(R))))
  @test iszero(divrem(l[1] + l[2], gb)[2])
end

@testset "Primary decomposition" begin
  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  I = ideal(R, [x, y*z^2])
  J = ideal(R, [x, y^2])
  L = primary_decomposition(I)
  Q = ideal(R, [R(1)])
  @test is_prime(I) == false
  @test is_primary(I) == false
  @test is_prime(J) == false
  @test is_primary(J) == true
  for (q, p) in L
    Q = intersect(Q, q)
    @test is_primary(q)
    @test is_prime(p)
    @test p == radical(q)
  end
  @test Q == I
end

@testset "#795" begin
  R, = QQ["x", "y"]
  I = ideal(R, zero(R))
  @test issubset(I, I)
  @test I == I
end

@testset "#975" begin
  A, t = PolynomialRing(QQ, ["t"])
  R, (x,y,z) = PolynomialRing(A, ["x", "y", "z"])
  I = ideal(R, [x])
  @test x in I
end

@testset "#1668" begin
  R, (x,) = PolynomialRing(QQ, ["x"])
  @test is_one(evaluate(x//x, [QQ(1)]))
end

@testset "Fraction fields over polynomial rings" begin
  R, x = PolynomialRing(QQ, "x")
  F = FractionField(R)
  @test gen(F) == F(x)
  @test gens(F) == elem_type(F)[ F(x) ]

  R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
  F = FractionField(R)
  @test ngens(F) == 3
  @test gen(F, 2) == F(y)
  @test gens(F) == elem_type(F)[ F(x), F(y), F(z) ]
  @test F[1] == F(x)
  @test F[0] == zero(F)
end

@testset "Grassmann Plücker Relations" begin
    R, x = PolynomialRing(ResidueRing(ZZ, 7), "x" => (1:2, 1:3), ordering=:degrevlex)
    test_ideal = ideal([x[1, 2]*x[2, 2] + 6*x[2, 1]*x[1, 3] + x[1, 1]*x[2, 3]])
    @test grassmann_pluecker_ideal(R, 2, 4) == test_ideal
end
