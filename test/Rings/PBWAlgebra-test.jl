@testset "PBWAlgebra.constructor" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(r))

  @test elem_type(R) == PBWAlgElem{fmpq, Singular.n_Q}
  @test parent_type(x) == PBWAlgRing{fmpq, Singular.n_Q}
  @test coefficient_ring(R) == QQ
  @test coefficient_ring(x) == QQ
  @test base_ring(R) == r
  @test base_ring(x) == r
  @test symbols(R) == [:x, :y, :z]
  @test gens(R) == [x, y, z]
  @test ngens(R) == 3
  @test gen(R, 2) == y
  @test R[2] == y

  3*x^2 + y*z == R([3, 1], [[2, 0, 0], [0, 1, 1]])

  r, (x, y, z) = QQ["x", "y", "z"]
  @test_throws Exception pbw_algebra(r, [0 x*y+y x*z+z; 0 0 y*z+1; 0 0 0], lex(r))
  R, (x, y, z) = pbw_algebra(r, [0 x*y+y x*z+z; 0 0 y*z+1; 0 0 0], deglex(r); check = false)
  @test (z*y)*x != z*(y*x)
end

@testset "PBWAlgebra.printing" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(r))

  @test length(string(R)) > 2
  @test length(string(x + y)) > 2
  @test string(one(R)) == "1"
end

@testset "PBWAlgebra.iteration" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y x*z; 0 0 y*z + 1; 0 0 0], deglex(r))

  p = -((x*z*y)^6 - 1)

  @test iszero(constant_coefficient(p - constant_coefficient(p)))
  @test length(p) == length(collect(terms(p)))
  @test length(p) == length(collect(monomials(p)))

  s = zero(R)
  for t in terms(p)
    s += t
  end
  @test p == s
  @test p == leading_term(p) + tail(p)

  s = zero(R)
  for (c, m) in zip(coefficients(p), monomials(p))
    s += c*m
  end
  @test p == s
  @test p == leading_coefficient(p)*leading_monomial(p) + tail(p)

  s = build_ctx(R)
  for (c, e) in zip(coefficients(p), exponent_vectors(p))
    push_term!(s, c, e)
  end
  @test p == finish(s)

  @test p == R(collect(coefficients(p)), collect(exponent_vectors(p)))
end

@testset "PBWAlgebra.weyl_algebra" begin
  R, (x, dx) = weyl_algebra(QQ, ["x"])
  @test dx*x == 1 + x*dx

  R, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
  @test dx*x == 1 + x*dx
  @test dy*y == 1 + y*dy
  @test dx*y == y*dx
  @test dy*x == x*dy
  @test x*y == y*x
end

@testset "PBWAlgebra.opposite_algebra" begin
  R, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])
  opR, M = opposite_algebra(R)
  @test M(dy*dx*x*y) == M(y)*M(x)*M(dx)*M(dy)
  @test inv(M)(M(x)) == x
  @test M.([x, y, dx, dy]) == [M(x), M(y), M(dx), M(dy)]

  g = [1-dx, dy]
  @test base_ring(M.(right_ideal(g))) === opR
  @test M.(left_ideal(g)) == right_ideal(M.(g))
  @test M.(right_ideal(g)) == left_ideal(M.(g))
  @test M.(two_sided_ideal(g)) == two_sided_ideal(M.(g))
end

@testset "PBWAlgebra.ideals" begin
  R, (x, y, dx, dy) = weyl_algebra(QQ, ["x", "y"])

  I = left_ideal([x^2, y^2])
  @test length(string(I)) > 2
  @test ngens(I) > 1
  @test !iszero(I)
  @test !isone(I)
  @test x^2 - y^2 in I
  @test !(x + 1 in I)
  @test isone(I + left_ideal([dy^2]))
  @test gens(I) == [I[k] for k in 1:ngens(I)]

  @test y*dy in left_ideal(R, [dy])
  @test !(dy*y in left_ideal(R, [dy]))
  @test dy*y in right_ideal(R, [dy])
  @test !(y*dy in right_ideal(R, [dy]))

  I = two_sided_ideal(R, [dy])
  @test y*dy*y in I
  @test isone(I)

  I = intersect(left_ideal(R, [dy]), left_ideal(R, [dx]))
  @test x*dy*dx in I
  @test !(dy*dx*x in I)

  I = intersect(right_ideal(R, [dy]), right_ideal(R, [dx]))
  @test dy*dx*x in I
  @test !(x*dy*dx in I)

  @test is_one(intersect(two_sided_ideal(R, [dy]), two_sided_ideal(R, [x])))

  I = intersect(left_ideal([dx]), left_ideal([dy]), left_ideal([x]))
  @test x^2*dx == (x*dx-1)*x
  @test x^2*dx*dy in I
  @test dx^2*x == (x*dx+2)*dx
  @test dx^2*x*dy in I
  @test !(x in I)
  @test !(y in I)
  @test !(dx in I)
  @test !(dy in I)

  @test intersect(left_ideal([dx])) == left_ideal([dx])

  I = right_ideal(R, [dx^2])
  J = right_ideal(R, [dx^4*x, dx^2*y])
  @test intersect(I, J) != I
  @test intersect(I, J) == J
end

@testset "PBWAlgebra.ideals.multiplication" begin
  r, (x, y, z) = QQ["x", "y", "z"]
  R, (x, y, z) = pbw_algebra(r, [0 x*y+y x*z+z+y; 0 0 y*z; 0 0 0], lex(r))

  e1 = y*z
  e2 = x*y*z

  @test !is_one(two_sided_ideal([e1]))
  @test !is_one(two_sided_ideal([e2]))

  I = two_sided_ideal([e1*e2])

  @test I == left_ideal([e1])*right_ideal([e2])
  @test I != two_sided_ideal([e1])*right_ideal([e2])
  @test issubset(I, two_sided_ideal([e1])*right_ideal([e2]))
  @test issubset(I, two_sided_ideal([e1])*two_sided_ideal([e2]))

  @test isone(I^0)
  @test I^1 == I
  @test I^4 == (I^2)^2
end

@testset "PBWAlgebra.ideals.eliminate" begin
  r, (e, f, h, a) = QQ["e", "f", "h", "a"]
  rel = @pbw_relations(f*e == e*f-h, h*e == e*h+2*e, h*f == f*h-2*f)
  for o in [lex(r), deglex(r)]
    R, (e, f, h, a) = pbw_algebra(r, rel, o)
    I = left_ideal([e^3, f^3, h^3-4*h, 4*e*f+h^2-2*h-a])
    @test eliminate(I, [e, f, h]) == left_ideal([a^3 - 32*a^2 + 192*a])

    oo = weight_ordering([1,1,1,0], deglex(r))
    @test eliminate(I, [e, f, h]; ordering = oo) == left_ideal([a^3 - 32*a^2 + 192*a])

    @test_throws ErrorException eliminate(I, [h])
  end

  r, (a, h, f, e) = QQ["a", "h", "f", "e"]
  rel = @pbw_relations(e*f == f*e-h, e*h == h*e+2*e, f*h == h*f-2*f)
  for o in [revlex(r), deglex(r)]
    R, (a, h, f, e) = pbw_algebra(r, rel, o)
    I = right_ideal([e^3, f^3, h^3-4*h, 4*f*e+h^2-2*h-a])
    @test eliminate(I, [e, f, h]) == right_ideal([a^3 - 32*a^2 + 192*a])

    oo = weight_ordering([0,1,1,1], deglex(r))
    @test eliminate(I, [e, f, h]; ordering = oo) == right_ideal([a^3 - 32*a^2 + 192*a])

    @test_throws ErrorException eliminate(I, [h])
  end

  # forces the discovery of the weight [0,0,1,2]
  r, (a, b, x, d) = QQ["a", "b", "x", "d"]
  rel = @pbw_relations(b*a == a*b+3*a, d*a == a*d+3*x^2, x*b == b*x-x,
                       d*b == b*d+d, d*x == x*d+1)
  for o in [lex(r), deglex(r)]
    R, (a, b, x, d) = pbw_algebra(r, rel, o)
    I = left_ideal([a, x])
    @test eliminate(I, [x, d]) == left_ideal([a])

    oo = weight_ordering([0,0,1,2], deglex(r))
    @test eliminate(I, [x, d]; ordering = oo) == left_ideal([a])

    Rop, M = opposite_algebra(R)
    @test eliminate(M.(I), M.([x, d])) == M.(left_ideal([a]))
  end

  R, (x, dx) = weyl_algebra(QQ, ["x"])
  @test is_zero(eliminate(left_ideal([x*dx]), [x, dx]))
  @test is_one(eliminate(left_ideal([x, 1-x]), [x, dx]))

  # two-sided example
  R, (x, y, z) = QQ["x", "y", "z"];
  rel = @pbw_relations(y*x == -x*y, z*x == -x*z, z*y == -y*z)
  A, (x,y,z) = pbw_algebra(R, rel, lex(R))
  I = two_sided_ideal(A, [x^2+1, y^2+1, z^2+1, y*x+z, z*y+x, z*x-y])
  @test eliminate(I, [x, y]) == two_sided_ideal([z^2+1])
end
