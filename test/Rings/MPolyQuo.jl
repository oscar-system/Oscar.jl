@testset "MPolyQuoRing" begin
  R, (x,y) = polynomial_ring(QQ, [:x, :y])
  f = y^2+y+x^2
  C = ideal(R, [f])

  @test_throws ArgumentError P, = quo(R, C; ordering=neglex(R))
  P, = quo(R, C; ordering=lex(R))
  @test P.ordering == lex(R)

  Q, = quo(R, C)
  @test Q.ordering == degrevlex(R)
  @test one(Q) == 1
  @test zero(Q) == 0

  I = ideal([x^3 + y^3 - 3, x^5 + y^5 - 5])
  Q, = quo(R, I)
  @test length(Oscar._kbase(Q)) == 12
  b = inv(Q(x))
  @test isone(b*Q(x))

  xx, yy = gens(Q)

  @test xx^2 == xx * xx
  @test xx^0 == one(Q)

  J = ideal(R, [x-y^2, x^2-y^3+1])
  A, p = quo(R, J)
  f = A(y^2-y^4)
  @test simplify(f) == A(-x*y + x + 1)
  A, p = quo(R, J; ordering=lex(R))
  f = A(y^2-y^4)
  @test simplify(f) == A(-y^3 + y^2 + 1)

  A,p = quo(R, [x^3+y,y^2]; ordering=lex(R))
  @test A.I == ideal(R, [x^3+y, y^2])

  B,q = quo(R, x^3+y,y^2; ordering=lex(R))
  @test A.I == B.I

  @test characteristic(quo(R, ideal(x))[1]) == 0
  @test characteristic(quo(R, ideal(R, 1))[1]) == 1

  S, (s, t) = GF(5)[:s, :t]
  @test characteristic(quo(S, ideal([s, t]))[1]) == 5
end

@testset "MpolyQuo.manipulation" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [zero(R)])
  Q, q = quo(R,I)
  f = q(x*y)
  @test divides(one(Q), f) == (false, one(Q))
  @test_throws ErrorException divexact(one(Q), f)

  A, _ = quo(R, 2*x^2-5*y^3)
  (x, y) = (A(x), A(y))

  @test iszero(x-x)
  @test x*deepcopy(x) == x^2
  @test iszero(det(matrix(A, [x 5; y^3 2*x])))

  @test !divides(x, y)[1]
  @test divides(x, x) == (true, one(A))
  @test divexact(x, x) == one(A)
  @test divides(zero(A), x) == (true, zero(A))
  @test divexact(zero(A), x) == zero(A)

  # promote rule
  K = GF(2)
  Kx, (x, y) = K[:x, :y]
  I = ideal(Kx, elem_type(Kx)[])
  Kx, = quo(Kx, I)
  x = Kx(x)
  @test K(2) * x == Kx(2) * x

  # simplify
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [x^2])
  Q, q = quo(R,I)

  z = one(Q)
  simplify(z)
  mul!(z, z, Q(x))
  mul!(z, z, Q(x))
  @test iszero(z)

  R, (x,) = polynomial_ring(QQ, [:x])
  I = ideal(R, [x])
  Q, q = quo(R,I)

  z = zero(Q)
  simplify(z)
  z = add!(z, Q(x))
  @test iszero(z)
end

@testset "MPolyQuoRing.ideals" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  Q, _ = quo(R, ideal(R, [x*y, x*z]))
  (x, y, z) = map(Q, (x, y, z))

  I = ideal(Q, [x^2, y^2, z^2])
  @test I isa Oscar.Ideal
  @test y^2 in I
  @test x*(x-y) in I
  @test !(x in I)
  @test !(y*z in I)

  @test !iszero(ideal(Q, [x, y]))
  @test !iszero(ideal(Q, [y*z]))
  @test iszero(ideal(Q, [x*y, x*z]))
  @test iszero(ideal(Q, [x*y*z]))
  @test ideal(Q, [2*x]) + ideal(Q, [x*(y+z)]) == ideal(Q, [x])
  @test iszero(ideal(Q, [y*z])*ideal(Q, [x]))

  a = quotient(ideal(Q, [zero(Q)]), ideal(Q, [y*z]))
  @test a == ideal(Q, [x])
  @test a == ideal(Q, gens(a))

  b = a + ideal(Q, [z])
  @test b == ideal(Q, [z, x])
  @test b == ideal(Q, gens(b))

  b = a*ideal(Q, [z])
  @test b == ideal(Q, [z*x])
  @test b == ideal(Q, gens(b))

  I = ideal(Q, [x^2*y-x+y,y+1])
  simplify(I)
  SQ = Oscar.singular_poly_ring(Q)
  SI = I.gens.gens.S
  @test SI[1] == SQ(-x+y) && SI[2] == SQ(y+1)
  J = ideal(Q, [x+y+1,y+1])
  @test issubset(J, I) == true
  @test issubset(I, J) == false
  @test (I == J) == false
  @test dim(J)  == 1
  @test dim(J)  == J.dim  # test case if dim(J) is already set
  K = ideal(Q, [ x*y+1 ])
  @test intersect(I,J,K) == ideal(Q, [y+1, x])
  @test intersect(I,J,K) == intersect([I,J,K])

  R, (x, y) = graded_polynomial_ring(QQ, [ :x, :y ], [ 1, 2 ])
  I = ideal(R, [ x*y ])
  Q, RtoQ = quo(R, I)
  J = ideal(Q, [ x^3 + x*y, y, x^2 + y ])
  @test minimal_generating_set(J) == [ Q(y), Q(x^2) ]
  @test isdefined(J, :gb)
  @test minimal_generating_set(ideal(Q, [ Q() ])) == elem_type(Q)[]

  # 1530
  R, (x,y,z) = QQ[:x, :y, :z]
  I = ideal(R, zero(R))
  Q, _ = quo(R, I)
  J = ideal(Q, [x^7, y^2])
  @test !(one(Q) in J)

  # allowing empty set of generators 
  R, (x,y,z) = QQ[:x, :y, :z]
  I = ideal(R, x*y)
  Q, _ = quo(R, I)
  J = ideal(Q, elem_type(Q)[])
  @test Q(x*y) in J
  @test !(Q(x) in J)
  J2 = ideal(Q, elem_type(R)[])
  @test Q(x*y) in J2
  @test !(Q(x) in J2)
end

@testset "prime ideals" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  I = ideal(R, [x])
  A, _ = quo(R, I)
  J = ideal(A, [A(y)])
  @test is_prime(J)
  J2 = ideal(A, [A(z*y)])
  @test !is_prime(J2)
end

@testset "saturated ideal compatibility" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  I = ideal(R, [x+y+z])
  A, _ = quo(R,I)
  J = ideal(A, [x, y])
  @test z in saturated_ideal(J)
end

@testset "issue 1794" begin
  kk, _ = cyclotomic_field(5)
  R, (x, y) = kk[:x, :y]
  I = ideal(R, [x^2, y^2])
  A, _ = quo(R, I)
  a = inv(A(1-x*y))
  @test isone(a*(1-x*y))
end

@testset "issue #1901" begin
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  L, _ = localization(R, powers_of_element(R[1]))
  S, (s0, s1, s2) = polynomial_ring(L, [:s0, :s1, :s2])
  I = ideal(S, [x*s0 - y*s1^2, y*s0 - z*s2^7])
  Q, _ = quo(S, I)
  @test Q isa MPolyQuoRing
end

@testset "#2119/wrong promotion" begin
  R, (x,y) = QQ[:x, :y]
  I = ideal(R, x)
  Q, _ = quo(R, I)
  P, (u, v) = Q[:u, :v]
  J = ideal(P, u)
  A, _ = quo(P, J)
  @test iszero(x * u)
  @test iszero(u * x)
  @test iszero(Q(x) * u)
  @test iszero(u * Q(x))
  @test iszero(Q(x) * A(u))
  @test iszero(A(u) * A(u))
  @test iszero(A(x)*u)
end

@testset "issue #2292" begin
  R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  A, p = quo(R, ideal(R, [x-y]))
  V = [x, z^2, x^3+y^3, y^4, y*z^5]
  a = ideal(A, V)
  dim(a) # cashes a.gb
  gens(a.gb)
  @test Oscar.oscar_generators(a.gb) == MPolyDecRingElem[y, z^2]
end

@testset "quotients as vector spaces" begin
  R, (x,y) = QQ[:x, :y]
  I = ideal(R, [x^3- 2, y^2+3*x])
  A, pr = quo(R, I)
  V, id = vector_space(QQ, A)
  @test dim(V) == 6
  @test id.(gens(V)) == A.([x^2*y, x*y, y, x^2, x, one(x)])
  f = (x*3*y-4)^5
  f = A(f)
  @test id(preimage(id, f)) == f
  v = V[3] - 4*V[5]
  @test preimage(id, id(v)) == v
end

@testset "divides hack" begin
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, 1-x*y)
  o = invlex([x, y])
  Q = MPolyQuoRing(R, I, o)
  @test Oscar._divides_hack(one(Q), Q(y))[2] == Q(x)
end

@testset "representatives and hashing" begin
  R1, (x1, y1) = QQ[:x, :y]
  R2, (x2, y2) = ZZ[:x, :y]
  R3, (x3, y3) = GF(7)[:x, :y]

  Q1, p1 = quo(R1, ideal(R1, [x1-3]))
  Q2, p2 = quo(R2, ideal(R2, [x2-3]))
  Q3, p3 = quo(R3, ideal(R3, [x3-3]))


  @test hash(p1(x1)) == hash(Q1(3))
  @test hash(p2(x2)) == hash(Q2(3))
  @test hash(p3(x3)) == hash(Q3(3))

  @test simplify(p1(x1)).f == Q1(3).f
  @test simplify(p2(x2)).f == Q2(3).f
  @test simplify(p3(x3)).f == Q3(3).f

  # A case that uses the RingFlattening and does not have a singular backend:
  P, (u, v) = Q3[:u, :v]
  QP, pP = quo(P, ideal(P, [u]))
  @test_throws ErrorException hash(pP(v)) # Hashing is forbidden
  @test simplify(pP(v-3*u)).f != simplify(pP(v)).f # Simplification does not bring 
  # representatives to normal form
  @test pP(v-3*u) == pP(v) # Equality check works via ideal membership
  a = pP(v) # Simplification only checks for being zero
  @test !a.simplified
  simplify(a)
  @test a.simplified
  @test !iszero(a.f)
  b = pP(u) # Only in case the element is zero, the representative is changed
  simplify(b)
  @test iszero(b.f)

  P, (u, v) = R3[:u, :v]
  QP, pP = quo(P, ideal(P, [u]))
  @test_throws ErrorException hash(pP(v)) # Hashing is forbidden
  @test simplify(pP(v-3*u)).f != simplify(pP(v)).f # Simplification does not bring 
  # representatives to normal form
  @test pP(v-3*u) == pP(v) # Equality check works via ideal membership
  a = pP(v) # Simplification only checks for being zero
  @test !a.simplified
  simplify(a)
  @test a.simplified
  @test !iszero(a.f)
  b = pP(u) # Only in case the element is zero, the representative is changed
  simplify(b)
  @test iszero(b.f)

  L3, _ = localization(R3, powers_of_element(y3))
  P, (u, v) = L3[:u, :v]
  QP, pP = quo(P, ideal(P, [u]))
  @test_throws ErrorException hash(pP(v)) # Hashing is forbidden
  @test simplify(pP(v-3*u)).f != simplify(pP(v)).f # Simplification does not bring 
  # representatives to normal form
  @test pP(v-3*u) == pP(v) # Equality check works via ideal membership
  a = pP(v) # Simplification only checks for being zero
  @test !a.simplified
  simplify(a)
  @test a.simplified
  @test !iszero(a.f)
  b = pP(u) # Only in case the element is zero, the representative is changed
  simplify(b)
  @test iszero(b.f)

  W3, _ = localization(Q3, powers_of_element(y3))
  P, (u, v) = W3[:u, :v]
  QP, pP = quo(P, ideal(P, [u]))
  @test_throws ErrorException hash(pP(v)) # Hashing is forbidden
  @test simplify(pP(v-3*u)).f != simplify(pP(v)).f # Simplification does not bring 
  # representatives to normal form
  @test pP(v-3*u) == pP(v) # Equality check works via ideal membership
  a = pP(v) # Simplification only checks for being zero
  @test !a.simplified
  simplify(a)
  @test a.simplified
  @test !iszero(a.f)
  b = pP(u) # Only in case the element is zero, the representative is changed
  simplify(b)
  @test iszero(b.f)
  
  @test Oscar.HasNormalFormTrait(ZZ) isa Oscar.HasNoNormalForm
  @test Oscar.HasNormalFormTrait(zero(ZZ)) isa Oscar.HasNoNormalForm
  @test Oscar.HasNormalFormTrait(QQ) isa Oscar.HasSingularNormalForm
  @test Oscar.HasNormalFormTrait(zero(QQ)) isa Oscar.HasSingularNormalForm
  @test Oscar.HasGroebnerAlgorithmTrait(ZZ) isa Oscar.HasSingularGroebnerAlgorithm
  @test Oscar.HasGroebnerAlgorithmTrait(one(ZZ)) isa Oscar.HasSingularGroebnerAlgorithm
  @test Oscar.HasGroebnerAlgorithmTrait(QQ) isa Oscar.HasSingularGroebnerAlgorithm
  @test Oscar.HasGroebnerAlgorithmTrait(one(QQ)) isa Oscar.HasSingularGroebnerAlgorithm
  @test Oscar.HasNormalFormTrait(R1) isa Oscar.HasNoNormalForm
  @test Oscar.HasNormalFormTrait(one(R1)) isa Oscar.HasNoNormalForm
  @test Oscar.HasGroebnerAlgorithmTrait(R1) isa Oscar.HasRingFlattening
  @test Oscar.HasGroebnerAlgorithmTrait(one(R1)) isa Oscar.HasRingFlattening
end

@testset "Tensor product" begin
  R, x = polynomial_ring(QQ, 3, :x)
  S, y = polynomial_ring(QQ, 4, :y)

  @test_throws AssertionError tensor_product(R, GF(3)[:x, :y][1])

  I = ideal(R, [x[1]])
  J = ideal(S, [y[1]])
  for A in [R, quo(R, I)[1]]
    for B in [S, quo(S, J)[1]]
      T, AtoT, BtoT = tensor_product(A, B)

      @test ngens(T) == 7
      @test is_zero(kernel(AtoT))
      @test is_zero(kernel(BtoT))

      v = append!(AtoT.(gens(A)), BtoT.(gens(B)))
      @test v == gens(T)

      T, _, _ = tensor_product(A, B, use_product_ordering = true)
      @test default_ordering(T).o isa Oscar.Orderings.ProdOrdering
      @test cmp(default_ordering(T), lift(gen(T, ngens(A))), lift(gen(T, ngens(A) + 1))) == 1
    end
  end
end

@testset "dimensions" begin
  # to address issue #2721
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, x)
  Q, pr = quo(R, I)
  J = ideal(Q, y)
  @test J.dim === nothing
  @test dim(J) == 0
  @test J.dim !== nothing
  J2 = ideal(Q, [y, y+1])
  @test J2.dim === nothing
  @test dim(J2) == -inf
  @test J2.dim !== nothing
end

