
@testset "alltest" begin
  _set_even_odd_coeff! = Oscar._set_even_odd_coeff!
  @testset "rank zero corner cases" begin
    @testset "CliffordOrder" begin
      K, a = quadratic_field(-5)
      OK = ring_of_integers(K)
      Kzer = [K(0)]
      idzero = identity_matrix(K, 0)
      qs = quadratic_space(K, idzero)
      ls = lattice(qs, idzero)
      C = clifford_order(ls)
      @test is_commutative(C)
      @testset "Construction" begin
        @test typeof(C) == CliffordOrder{elem_type(base_ring_type(C)), typeof(ambient_algebra(C))}
        @test elem_type(C) == CliffordOrderElem{elem_type(base_ring_type(C)), typeof(ambient_algebra(C)), typeof(C)}
        @test elem_type(C) == typeof(C())
        @test base_ring_type(C) == typeof(OK)
        @test base_ring_type(C) == typeof(base_ring(C))

        @test (base_ring(C), rank(C), lattice(C), gram_matrix(C), coefficient_ideals(C)) == (OK, 1, ls, idzero, [fractional_ideal(OK, OK(1))])
        @test C() == C(0) && C() == zero(C)
        @test is_zero(C())
        @test coeff(C()) == Kzer
        @test even_coeff(C()) == Kzer && even_coeff(C()) == odd_coeff(C())
        @test C(1) == one(C)
        @test is_one(C(1))
        @test coeff(C(1)) == [K(1)]
        @test even_coeff(C(1)) == [K(1)] && odd_coeff(C(1)) == Kzer
      end
      @testset "functions on elements" begin
        x = C([17])
        @test parent(x) == C
        @test parent_type(x) == typeof(C)
        @test typeof(x) == typeof(C())

        @test even_part(x) == C([17])
        @test odd_part(x) == C([0])
        xeven, xodd = even_coeff(x), odd_coeff(x)
        x.even_coeffs, x.odd_coeffs = Kzer, [K(1)]
        @test xeven != even_coeff(x)
        @test xodd != odd_coeff(x)
        _set_even_odd_coeff!(x)
        @test xeven == even_coeff(x) && xodd == odd_coeff(x)

        @test +x == x
        @test -x == C([-17])
        @test divexact(2 * x, 2) == x && divexact(2 * x, QQ(2)) == x
        @test divexact(C(18), OK(2)) == C(9)
        @test_throws ArgumentError divexact(x, 2)
        @test_throws ArgumentError divexact(x, ZZ(2))
        @test_throws ArgumentError divexact(x, QQ(2))
        @test_throws ArgumentError divexact(x, OK(2))
      end
      @testset "defining relations" begin
        @test length(coefficient_ideals(pseudo_basis(C))) == 1
        @test pseudo_basis(C, 1) == pseudo_matrix(identity_matrix(K, 1), [fractional_ideal(OK, one(OK))])
        @test is_empty(gens(C))
        a, b = rand(ZZ, -100:100), rand(ZZ, -100:100)
        @test C(a) * C(b) == C(a * b)
        @test C(a) + C(b) == C(a + b)
        @test (C(a) + C(b)) * (C(a) - C(b)) == C(a)^2 - C(b)^2
      end
      @testset "center and centroid" begin
        @test center(C) == centroid(C)
        @test center(C) == pseudo_matrix(identity_matrix(K, 1), [fractional_ideal(OK, one(OK))])
        @test disq(C) == quadratic_discriminant(C)
        @test disq(C) == (fractional_ideal(OK, one(OK)), K(1))
      end
    end
    @testset "ZZCliffordOrder" begin
      QQzer = [QQ(0)]
      idzero = identity_matrix(QQ, 0)
      qs = quadratic_space(QQ, idzero)
      ls = lattice(qs, idzero)
      C = clifford_order(ls)
      @test is_commutative(C)
      @testset "Construction" begin
        @test typeof(C) == ZZCliffordOrder
        @test elem_type(C) == ZZCliffordOrderElem
        @test elem_type(C) == typeof(C())
        @test base_ring_type(C) == ZZRing
        @test base_ring_type(C) == typeof(base_ring(C))

        @test (base_ring(C), rank(C), lattice(C), gram_matrix(C)) == (ZZ, 1, ls, idzero)
        @test C() == C(0) && C() == zero(C)
        @test is_zero(C())
        @test coeff(C()) == QQzer
        @test even_coeff(C()) == QQzer && even_coeff(C()) == odd_coeff(C())
        @test C(1) == one(C)
        @test is_one(C(1))
        @test coeff(C(1)) == [QQ(1)]
        @test even_coeff(C(1)) == [QQ(1)] && odd_coeff(C(1)) == QQzer
      end
      @testset "functions on elements" begin
        x = C([17])
        @test parent(x) == C
        @test parent_type(x) == typeof(C)
        @test typeof(x) == typeof(C())

        @test even_part(x) == C([17])
        @test odd_part(x) == C([0])
        xeven, xodd = even_coeff(x), odd_coeff(x)
        x.even_coeffs, x.odd_coeffs = QQzer, [QQ(1)]
        @test xeven != even_coeff(x)
        @test xodd != odd_coeff(x)
        _set_even_odd_coeff!(x)
        @test xeven == even_coeff(x) && xodd == odd_coeff(x)

        @test +x == x
        @test -x == C([-17])
        @test divexact(2 * x, 2) == x && divexact(2 * x, QQ(2)) == x
        @test divexact(C(18), ZZ(2)) == C(9)
        @test_throws ArgumentError divexact(x, 2)
        @test_throws ArgumentError divexact(x, ZZ(2))
        @test_throws ArgumentError divexact(x, QQ(2))
        @test_throws ArgumentError divexact(x, 2//1)
      end
      @testset "defining relations" begin
        @test length(basis(C)) == 1
        @test basis(C, 1) == one(C)
        @test is_empty(gens(C))
        a, b = rand(ZZ, -100:100), rand(ZZ, -100:100)
        @test C(a) * C(b) == C(a * b)
        @test C(a) + C(b) == C(a + b)
        @test (C(a) + C(b)) * (C(a) - C(b)) == C(a)^2 - C(b)^2
      end
      @testset "center and centroid" begin
        @test center(C) == centroid(C)
        @test center(C) == [one(C)]
        @test disq(C) == quadratic_discriminant(C)
        @test disq(C) == ZZ(1)
      end
    end
  end
end
