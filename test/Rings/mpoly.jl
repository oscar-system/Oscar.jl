@testset "Polynomial ring constructor" begin
  R, x = @inferred polynomial_ring(QQ, :x => 1:3)
  @test length(x) == 3
  @test x isa Vector
  @test x == gens(R)
  for i in 1:3
    @test sprint(show, "text/plain", x[i]) == "x[$i]"
  end

  R, x = @inferred polynomial_ring(QQ, :x => (1:3, 1:4))
  @test length(x) == 12
  @test x isa Matrix{<: Any}
  @test size(x) == (3, 4)
  @test allunique(x)
  @test Set(gens(R)) == Set(x)
  for i in 1:3
    for j in 1:4
      @test sprint(show, "text/plain", x[i, j]) == "x[$i, $j]"
    end
  end

  R, x, y = @inferred polynomial_ring(QQ, :x => (1:3, 1:4), :y => 1:2)
  @test length(x) == 12
  @test x isa Matrix{<: Any}
  @test size(x) == (3, 4)
  @test allunique(x)
  for i in 1:3
    for j in 1:4
      @test sprint(show, "text/plain", x[i, j]) == "x[$i, $j]"
    end
  end
  @test y isa Vector{<: Any}
  @test length(y) == 2
  @test allunique(y)
  for i in 1:2
    @test sprint(show, "text/plain", y[i]) == "y[$i]"
  end
  @test Set(gens(R)) == union(Set(x), Set(y))

  R, x, y, z = @inferred polynomial_ring(QQ, :x => (1:3, 1:4),
                                             :y => 1:2,
                                             :z => (1:1, 1:1, 1:1))
  @test length(x) == 12
  @test x isa Matrix{<: Any}
  @test size(x) == (3, 4)
  @test allunique(x)
  for i in 1:3
    for j in 1:4
      @test sprint(show, "text/plain", x[i, j]) == "x[$i, $j]"
    end
  end
  @test y isa Vector{<: Any}
  @test length(y) == 2
  @test allunique(y)
  for i in 1:2
    @test sprint(show, "text/plain", y[i]) == "y[$i]"
  end
  @test z isa Array{<: Any, 3}
  @test sprint(show, "text/plain", z[1, 1, 1]) == "z[1, 1, 1]"
  @test Set(gens(R)) == union(Set(x), Set(y), Set(z))
end

@testset "Polynomial homs" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I1 = x^2 + y^2
  I2 = x^2 * y^2
  I3 = x^3*y - x*y^3
  S, (a,b,c) = polynomial_ring(QQ, [:a, :b, :c])
  h = hom(S, R, [I1, I2, I3])
  @test kernel(h) == ideal(S, [a^2*b - 4*b^2 - c^2])
  @test h(gen(S, 1)) == I1
  @test h.([a,b]) == [I1, I2]
  @test preimage(h, ideal(R, [I2, I3])) == ideal(S, [b, c])
end

@testset "Ideal operations" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [x^2, y^2])
  J = ideal(R, [x*y^2, x^2])
  S = ideal(R, [x^2, y^2, x*y^2])
  P = ideal(R, [x^3*y^2, x^4, x*y^4, x^2*y^2])
  @test Oscar.check_base_rings(I, J) === nothing
  @test I+J == I-J == S
  @test I*J == P
  @test intersect(I,J,P) == ideal(R,[x^2*y^2, x^4, x*y^4])
  @test intersect(I,J,P) == intersect([I,J,P])

  I = ideal(R, [x^3, x*y^2, x^2*y])
  J = ideal(R, [x,y])
  @test saturation(I) == ideal(R, [x])
  @test saturation(I) == saturation(I, J)
  @test saturation_with_index(I) == (ideal(R, [x]), 2)
  @test saturation_with_index(I) == saturation_with_index(I, J)

  @test I != J
  RR, (xx, yy) = grade(R, [1, 1])
  @test_throws ErrorException ideal(R, [x]) == ideal(RR, [xx])
  @test is_subset(I, I)
  RR, (xx, yy) = polynomial_ring(QQ, [:xx, :yy])
  @test_throws ErrorException is_subset(ideal(R, [x]), ideal(RR, [xx]))

  f = x^2 + y^2
  g = x^4*y - x*y^3
  I = [f, g]
  S, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
  J = ideal(S, [(c^2+1)*(c^3+2)^2, b-c^2])
  @test_throws ErrorException Oscar.check_base_rings(P, J)
  r1 = c^2-b
  r2 = b^2*c+c^3+2*c^2+2
  L = gens(radical(J))

  @test jacobian_ideal(f) == ideal(R, [2*x, 2*y])
  @test jacobian_matrix(f) == matrix(R, 2, 1, [2*x, 2*y])
  @test jacobian_matrix(I) == matrix(R, 2, 2, [2*x, 4*x^3*y-y^3, 2*y, x^4-3*x*y^2])
  @test length(L) == 2

  # Test disabled because it could not be reliably reproduced and also 
  # it is mathematical not rigorous
  # @test length(findall(==(r1), L)) == 1
  # @test length(findall(==(r2), L)) == 1
  @test ideal(parent(r1), L) == ideal(parent(r1), [r1, r2])

  @test issubset(ideal(S, [a]), ideal(S, [a]))
  @test issubset(ideal(S, [a]), ideal(S, [a, b]))
  @test !issubset(ideal(S, [c]), ideal(S, [b]))
  @test !issubset(ideal(S, [a, b, c]), ideal(S, [a*b*c]))
end

@testset "Primary decomposition" begin

  # primary_decomposition
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  i = ideal(R, [x, y*z^2])
  for method in (:GTZ, :SY)
    j = ideal(R, [R(1)])
    for (q, p) in primary_decomposition(i, algorithm=method)
      j = intersect(j, q)
      @test is_primary(q)
      @test is_prime(p)
      @test p == radical(q)
    end
    @test j == i
  end

  R, (a, b, c, d) = polynomial_ring(ZZ, [:a, :b, :c, :d])
  i = ideal(R, [9, (a+3)*(b+3)])
  l = primary_decomposition(i)
  @test length(l) == 2

  # minimal_primes
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  i1 = ideal(R, [z^3+2, -z^2+y])
  i2 = ideal(R, [z^2+1, -z^2+y])
  l = minimal_primes(i)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  l = minimal_primes(i, algorithm=:charSets)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1
  R, (a, b, c, d) = polynomial_ring(ZZ, [:a, :b, :c, :d])
  i = ideal(R, [R(9), (a+3)*(b+3)])
  i1 = ideal(R, [R(3), a])
  i2 = ideal(R, [R(3), b])
  l = minimal_primes(i)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  # equidimensional_decomposition_weak
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]),
                ideal(R, [x^2, z^2]), ideal(R, [x^5, y^5, z^5]))
  l = equidimensional_decomposition_weak(i)
  @test length(l) == 3
  @test l[1] == ideal(R, [z^4, y^5, x^5, x^3*z^3, x^4*y^4])
  @test l[2] == ideal(R, [y*z, x*z, x^2])
  @test l[3] == ideal(R, [z])

 # equidimensional_decomposition_radical
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  l = equidimensional_decomposition_radical(i)
  @test length(l) == 1
  @test l[1] == ideal(R, [z^2-y, y^2*z+z^3+2*z^2+2])
  
  # equidimensional_hull
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]),
                ideal(R, [x^2, z^2]), ideal(R, [x^5, y^5, z^5]))
  @test equidimensional_hull(i) == ideal(R, [z])

  R, (a, b, c, d) = polynomial_ring(ZZ, [:a, :b, :c, :d])
  i = intersect(ideal(R, [R(9), a, b]),
                ideal(R, [R(3), c]),
                ideal(R, [R(11), 2*a, 7*b]),
                ideal(R, [13*a^2, 17*b^4]),
                ideal(R, [9*c^5, 6*d^5]),
                ideal(R, [R(17), a^15, b^15, c^15, d^15]))
  @test equidimensional_hull(i) == ideal(R, [R(3)])

  # equidimensional_hull_radical
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  @test equidimensional_hull_radical(i) == ideal(R, [z^2-y, y^2*z+z^3+2*z^2+2])

  # absolute_primary_decomposition
  R,(x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [(z+1)*(z^2+1)*(z^3+2)^2, x-y*z^2])
  d = @inferred absolute_primary_decomposition(I)
  @test length(d) == 3

  R,(x,y,z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  I = ideal(R, [(z+y)*(z^2+y^2)*(z^3+2*y^3)^2, x^3-y*z^2])
  d = @inferred absolute_primary_decomposition(I)
  @test length(d) == 5

  d = @inferred absolute_primary_decomposition(ideal(R()))
  @test length(d) == 1

  d = @inferred absolute_primary_decomposition(ideal(R(1)))
  @test isempty(d)

  # Issue 4039
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [x + 1, y + 1, y])
  d = @inferred absolute_primary_decomposition(I)
  @test isempty(d)

  # is_prime
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [one(R)])
  @test is_prime(I) == false
  
  J = ideal(R, [x*(x-1), y*(y-1), x*y])
  l = minimal_primes(J)
  @test length(l) == 3
  
  QQt, t = QQ[:t]
  kk, a = extension_field(t^2 + 1)
  
  R, (x, y) = kk[:x, :y]
  J = ideal(R, [x^2 + 1, y^2 + 1, (x - a)*(y - a)])
  l = minimal_primes(J)
  @test length(l) == 3

end

@testset "Groebner" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  I = ideal([x+y^2+4*z-5,x+y*z+5*z-2])
  @test gens(groebner_basis(I, ordering=lex(gens(R)))) == [y^2 - y*z - z - 3, x + y*z + 5*z - 2]
  @test gens(groebner_basis(I, ordering=degrevlex(gens(R)), complete_reduction = true)) == [x + y*z + 5*z - 2, x + y^2 + 4*z - 5, x*y - x*z - 5*x - 2*y - 4*z^2 - 20*z + 10, x^2 + x*z^2 + 10*x*z - 4*x + 4*z^3 + 20*z^2 - 20*z + 4]

  # Test coefficient rings that are actually fields for safety. The first three
  # are native to singular while FpFieldElem now has a proper wrapper
  for Zn in [residue_ring(ZZ, 11)[1], residue_ring(ZZ, ZZRingElem(10)^50+151)[1], GF(11),
             GF(ZZRingElem(10)^50+151)]
    # Setting the internal_ordering is necessary for divrem to use the correct ordering
    R, (x, y) = polynomial_ring(Zn, [:x, :y], internal_ordering = :degrevlex)
    l = [x*y+x^3+1, x*y^2+x^2+1]
    g = gens(groebner_basis(ideal(R, l); ordering = degrevlex(gens(R))))
    @test iszero(divrem(l[1] + l[2], g)[2])
  end

  F, a = finite_field(11, 2, "a")
  # Setting the internal_ordering is necessary for divrem to use the correct ordering
  R, (x, y, z) = polynomial_ring(F, [:x, :y, :z], internal_ordering = :degrevlex)
  l = [3*x^5 + a*x*y^2 + a^2*z^2, z^3*x^2 + 7*y^3 + z]
  gb = gens(groebner_basis(ideal(R, l); ordering = degrevlex(gens(R))))
  @test iszero(divrem(l[1] + l[2], gb)[2])
end

@testset "Primary decomposition" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
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
  R, = QQ[:x, :y]
  I = ideal(R, zero(R))
  @test issubset(I, I)
  @test I == I
end

@testset "#975" begin
  A, t = polynomial_ring(QQ, [:t])
  R, (x,y,z) = polynomial_ring(A, [:x, :y, :z])
  I = ideal(R, [x])
  @test x in I
end

@testset "#1668" begin
  R, (x,) = polynomial_ring(QQ, [:x])
  @test is_one(evaluate(x//x, [QQ(1)]))
end

@testset "Fraction fields over polynomial rings" begin
  R, x = polynomial_ring(QQ, :x)
  F = fraction_field(R)
  @test gen(F) == F(x)
  @test gens(F) == elem_type(F)[ F(x) ]

  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  F = fraction_field(R)
  @test ngens(F) == 3
  @test gen(F, 2) == F(y)
  @test gens(F) == elem_type(F)[ F(x), F(y), F(z) ]
  @test F[1] == F(x)
end

@testset "Grassmann Plücker Relations" begin
  R, x = graded_polynomial_ring(residue_ring(ZZ, 7)[1], :x => (1:2, 1:3))
  test_ideal =ideal([x[1, 2]*x[2, 2] + 6*x[2, 1]*x[1, 3] + x[1, 1]*x[2, 3]])
  @test grassmann_pluecker_ideal(R, 2, 4) == test_ideal

  # Issue #4018
  I = grassmann_pluecker_ideal(2, 5)
  @test degree(I) == 5
  @test dim(I) == 7
  R, x = graded_polynomial_ring(GF(7), 10, :x)
  I = grassmann_pluecker_ideal(R, 2, 5)
  @test degree(I) == 5
  @test dim(I) == 7
end

@testset "IdealGens" begin
   R, (x0, x1, x2) = polynomial_ring(QQ, [:x0, :x1, :x2])
   I = ideal([x0*x1,x2])
   g = generating_system(I)
   @test elements(g) == [x0*x1, x2]
   @test g.isGB == false
   @test isdefined(g, :ord) == false
   h = set_ordering(g, degrevlex(gens(R)))
   @test h != g
   @test ordering(h) == degrevlex(gens(R))
   h1 = Oscar.IdealGens(base_ring(g), gens(g))
   Oscar.set_ordering!(h1, lex(R))
   @test ordering(h1) == lex(gens(R))
end

@testset "NonSimpleField" begin
  Qt, t = QQ[:t];
  K, (a1,a2) = number_field([t^2 - 2, t^2 - 3], "a");
  R, (x,y) = polynomial_ring(K,[:x,:y]);
  I = ideal(R, [x^2-a1])

  #just to see it working...
  radical(I)
  primary_decomposition(I)
  minimal_primes(I)
  equidimensional_hull(I)
  equidimensional_hull_radical(I)
  equidimensional_decomposition_weak(I)
end

@testset "absolute primary decomposition over number fields" begin
  P1, t1 = QQ[:t1];
  kk1, a = extension_field(t1^4 + 1);
  #P2, t2 = kk1[:t2]
  #kk2, b = extension_field(t2^3 - 7); # working over this field is too expensive.
  R, (x,) = polynomial_ring(kk1, [:x]);

  I = ideal(R, x^8 + 1)
  dec = absolute_primary_decomposition(I^2)
  @test length(dec) == 4
  for (Q, P, PP, d) in dec
    @test d == 2
    R_ext = base_ring(PP)
    L = coefficient_ring(R_ext)
    @test all(map_coefficients(L, g; parent=R_ext) in PP for g in gens(Q))
  end
end

@testset "Hessian matrix" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  f = x^2 + x*y^2 - z^3
  H = hessian_matrix(f)
  @test H == R[2 2*y 0; 2*y 2*x 0; 0 0 -6*z]
  @test hessian(f) == -24*x*z + 24*y^2*z
end

@testset "rerouting of primary decomposition over number fields" begin
  P, t = QQ[:t]
  kk = splitting_field(t^3 - 1)
  kks, s = kk[:s]
  kk2 = splitting_field(s^2 + 1)

  R, (x, y) = kk[:x, :y]
  alpha = first(gens(kk))
  R_flat, iso, iso_inv = Oscar._expand_coefficient_field(R)
  theta = first(gens(R_flat))
  @test iso_inv(R(alpha)) == theta
  @test iso(theta) == alpha
  id1 = compose(iso, iso_inv)
  @test all(x->x==id1(x), gens(R_flat))
  id2 = compose(iso_inv, iso)
  @test all(x->x==id2(x), gens(R))

  S, (u, v) = kk2[:u, :v]
  S_flat, iso, iso_inv = Oscar._expand_coefficient_field(S)
  S_ff, iso2, iso_inv2 = Oscar._expand_coefficient_field(S_flat)

  iso_tot = compose(iso2, iso)
  iso_inv_tot = compose(iso_inv, iso_inv2)

  beta = first(gens(kk2))
  @test iso_inv_tot(S(beta)) == S_ff[2]
  @test iso_inv_tot(S(alpha)) == S_ff[1]

  Sff, iso2, iso_inv2 = Oscar._expand_coefficient_field_to_QQ(S)
  @test iso_inv_tot(S(beta)) == S_ff[2]
  @test iso_inv_tot(S(alpha)) == S_ff[1]

  I = ideal(S, [(u^2 + 1)^3])
  dec = primary_decomposition(I)
  @test length(dec) == 2
end

@testset "primary decomposition in graded rings" begin
  Pt, t = QQ[:t]
  f = t^2 + 1
  kk, i = number_field(f)
  R, (x, y) = kk[:x, :y]

  S, _ = grade(R)

  S_exp, iso, iso_inv = Oscar._expand_coefficient_field(S)
  @test is_graded(S_exp)
  @test iszero(degree(S_exp[1]))
  @test domain(iso) === S_exp
  @test codomain(iso) === S
  @test domain(iso_inv) === S
  @test codomain(iso_inv) === S_exp
  @test compose(iso, iso_inv).(gens(S_exp)) == gens(S_exp)

  I = ideal(S, [x^2 + y^2])
  dec = primary_decomposition(I)
  @test length(dec) == 2

  Q, _ = quo(S, I)
  dec = primary_decomposition(ideal(Q, zero(Q)))
  @test length(dec) == 2
end
  
@testset "primary decomposition over number fields" begin
  Pt, t = QQ[:t]
  f = t^3 - 7
  kk, _ = number_field(f)

  # primary_decomposition
  R, (x, y, z) = polynomial_ring(kk, [:x, :y, :z])
  i = ideal(R, [x, y*z^2])
  for method in (:GTZ, :SY)
    j = ideal(R, [R(1)])
    for (q, p) in primary_decomposition(i, algorithm=method)
      j = intersect(j, q)
      @test is_primary(q)
      @test is_prime(p)
      @test p == radical(q)
    end
    @test j == i
  end

  R, (a, b, c, d) = polynomial_ring(ZZ, [:a, :b, :c, :d])
  i = ideal(R, [9, (a+3)*(b+3)])
  l = primary_decomposition(i)
  @test length(l) == 2

  # minimal_primes
  R, (x, y, z) = polynomial_ring(kk, [:x, :y, :z])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  l = minimal_primes(i)
  @test length(l) == 2

  l = minimal_primes(i, algorithm=:charSets)
  @test length(l) == 2

  R, (a, b, c, d) = polynomial_ring(ZZ, [:a, :b, :c, :d])
  i = ideal(R, [R(9), (a+3)*(b+3)])
  i1 = ideal(R, [R(3), a])
  i2 = ideal(R, [R(3), b])
  l = minimal_primes(i)
  @test length(l) == 2
  @test l[1] == i1 && l[2] == i2 || l[1] == i2 && l[2] == i1

  # equidimensional_decomposition_weak
  R, (x, y, z) = polynomial_ring(kk, [:x, :y, :z])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]),
                ideal(R, [x^2, z^2]), ideal(R, [x^5, y^5, z^5]))
  l = equidimensional_decomposition_weak(i)
  @test length(l) == 3
  @test l[1] == ideal(R, [z^4, y^5, x^5, x^3*z^3, x^4*y^4])
  @test l[2] == ideal(R, [y*z, x*z, x^2])
  @test l[3] == ideal(R, [z])

  # equidimensional_decomposition_radical
  R, (x, y, z) = polynomial_ring(kk, [:x, :y, :z])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  l = equidimensional_decomposition_radical(i)
  @test length(l) == 1

  # equidimensional_hull
  R, (x, y, z) = polynomial_ring(kk, [:x, :y, :z])
  i = intersect(ideal(R, [z]), ideal(R, [x, y]),
                ideal(R, [x^2, z^2]), ideal(R, [x^5, y^5, z^5]))
  @test equidimensional_hull(i) == ideal(R, [z])

  R, (a, b, c, d) = polynomial_ring(ZZ, [:a, :b, :c, :d])
  i = intersect(ideal(R, [R(9), a, b]),
                ideal(R, [R(3), c]),
                ideal(R, [R(11), 2*a, 7*b]),
                ideal(R, [13*a^2, 17*b^4]),
                ideal(R, [9*c^5, 6*d^5]),
                ideal(R, [R(17), a^15, b^15, c^15, d^15]))
  @test equidimensional_hull(i) == ideal(R, [R(3)])

  # equidimensional_hull_radical
  R, (x, y, z) = polynomial_ring(kk, [:x, :y, :z])
  i = ideal(R, [(z^2+1)*(z^3+2)^2, y-z^2])
  @test equidimensional_hull_radical(i) == ideal(R, [z^2-y, y^2*z+z^3+2*z^2+2])

  # absolute_primary_decomposition
  R,(x,y,z) = polynomial_ring(kk, [:x, :y, :z])
  I = ideal(R, [x^2 + 1])
  d = absolute_primary_decomposition(I)
  @test length(d) == 1

  #= Tests disabled because of too long runtime
  R,(x,y,z) = graded_polynomial_ring(kk, [:x, :y, :z])
  I = ideal(R, [(z+y)*(z^2+y^2)*(z^3+2*y^3)^2, x^3-y*z^2])
  d = absolute_primary_decomposition(I)
  @test length(d) == 5
  =#

  # is_prime
  R, (x, y) = polynomial_ring(kk, [:x, :y])
  I = ideal(R, [one(R)])
  @test is_prime(I) == false
end

@testset "primary decomposition over AbsNonSimpleNumField" begin
  _, x = QQ[:x]
  K, a = number_field([x - 1, x - 2]);
  Kt, t = K[:t];
  L, b = number_field(t - 1, "b");
  M, = number_field(t - 1, "b");
  Mx, = polynomial_ring(M, 2);
  primary_decomposition(ideal(Mx, [gen(Mx, 1)]))
  S, _ = grade(Mx)
  @test is_graded(Oscar._expand_coefficient_field_to_QQ(S)[1])
  primary_decomposition(ideal(S, [gen(S, 1)]))
end

@testset "primary decomposition over RelNonSimpleNumField" begin
  Pt, t = QQ[:t]
  f = t^2 + 1
  kk, i = number_field(f)
  _, T = kk[:T]
  g = T^3 - 5
  K, zeta = number_field([g], "zeta")
  @test K isa Hecke.RelNonSimpleNumField
  _, s = K[:s]
  h = s-1
  L, xi = number_field(h)
  R, (x, y) = L[:x, :y]
  I = ideal(R, x^2 + y^2)
  @test length(primary_decomposition(I)) == 2
  S, _ = grade(R)
  @test is_graded(Oscar._expand_coefficient_field_to_QQ(S)[1])
  IS = ideal(S, x^2 + y^2)
  @test length(primary_decomposition(IS)) == 2
end

@testset "flag pluecker ideal" begin
  dimension_vector = [2]
  ambient_dimension = 4
  I = flag_pluecker_ideal(dimension_vector, ambient_dimension)
  R = base_ring(I)
  @test dim(R) == 6
  x = gens(R)
  f1 = -x[1]*x[5]+x[2]*x[4]-x[3]*x[6] 
  @test [f1] == gens(I)
  I2 = flag_pluecker_ideal(GF(3),[1,2,3], 4);
  @test length(gens(I2)) == 10

  R, _ = polynomial_ring(QQ, 6)
  I = flag_pluecker_ideal(R,dimension_vector, ambient_dimension, minimal=false)
  @test isa(I, Ideal)
end

@testset "default ordering" begin
  R, _ = QQ[:x, :y, :z]
  S, _ = grade(R, [2, 1, 2])
  for T in [R, S]
    x, y, z = gens(T)
    f = x + y^2
    I = ideal(T, [y^2 - z, x - z])
    old_default = @inferred default_ordering(T)
    f = with_ordering(T, invlex(T)) do
          normal_form(f, I)
        end
    @test f == normal_form(f, I, ordering = invlex(T))
    @test default_ordering(T) == old_default

    # Make sure the ordering is reset in case of an error
    # The `@test_throws` is just here to catch the error
    @test_throws ErrorException with_ordering(T, invlex(T)) do
      error()
    end
    @test default_ordering(T) == old_default
  end
end

@testset "Issue 3992" begin
  P, (x, y) = QQ[:x, :y]
  I = ideal(P, elem_type(P)[])
  @test !radical_membership(x, I)
end

@testset "preprocessing for radical computations" begin
  kk7 = GF(7^3)
  P0, t0 = QQ[:t0]
  mipo1 = t0^2 + 1
  kk1, alpha_1 = extension_field(mipo1)

  P1, t1 = kk1[:t1]
  mipo2 = t1^2 - 2
  kk2, alpha2 = extension_field(mipo2)

  for kk in [QQ, GF(23), kk1, kk2, kk7]
    R0, (x0, y0) = kk[:x0, :y0]
    I0 = ideal(R0, [x0^4*(y0 + 5)^8])
    @test x0*(y0+5) in radical(I0^2)
    @test x0*(y0+5) in radical(I0^2; eliminate_variables=false)
    @test x0*(y0+5) in radical(I0^2; eliminate_variables=false, factor_generators=false)
    @test x0*(y0+5) in radical(I0^2; eliminate_variables=true, factor_generators=false)
  end
end

@testset "dimensions" begin
  # to address issue #2721
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, x)
  @test I.dim === nothing
  @test dim(I) == 1
  @test I.dim !== nothing
  
  I2 = ideal(R, [x, x+1])
  @test I2.dim === nothing
  @test dim(I2) == -inf
  @test I2.dim !== nothing
end

@testset "dimensions over number fields" begin
  P, t = QQ[:t]
  kk, i = extension_field(t^2 + 1)
  R, (x, y) = kk[:x, :y]
  @test dim(R) == 2
end

