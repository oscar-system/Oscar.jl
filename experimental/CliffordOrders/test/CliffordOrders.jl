
@testset "alltest" begin
  _set_even_odd_coeff!, mul_with_gen = Oscar._set_even_odd_coeff!, Oscar._mul_with_gen
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

  @testset "quadratic_field_sqrt(5)" begin
    K, b = quadratic_field(5)
    O = maximal_order(K)
    ide = fractional_ideal(O, one(O))
    a = 1//2 * (1 + b)
    G = K[2*a 1; 1 2*(1 - a)]
    qsK = quadratic_space(K, G)
    lsK = lattice(qsK)
    C = clifford_order(lsK)
    Kzer = K.([0, 0, 0, 0])
    e(i::Int) = C(matrix(gen(C, i))[1,:])
    
    @test !is_commutative(C)
    
    @testset "construction" begin
      @test typeof(C) == CliffordOrder{elem_type(base_ring_type(C)), typeof(ambient_algebra(C))}
      @test elem_type(C) ==
      CliffordOrderElem{elem_type(base_ring_type(C)), typeof(ambient_algebra(C)), typeof(C)}
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == typeof(base_ring(C))

      @test (base_ring(C), gram_matrix(C), lattice(C), rank(C), coefficient_ideals(C)) == (O, G, lsK, 4, [ide, ide, ide, ide])
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coeff(C()) == Kzer
      @test even_coeff(C()) == Kzer && even_coeff(C()) == odd_coeff(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coeff(C(1)) == K.([1, 0, 0, 0])
      @test even_coeff(C(1)) == K.([1, 0, 0, 0]) && odd_coeff(C(1)) == Kzer
      verr1, verr2 = K.([1, 2, 3]), K.([1, 2, 3, 4, 5])
      @test_throws ArgumentError C(verr1)
      @test_throws ArgumentError C(verr2)
      Ca = C(a)
      @test parent(Ca) == C
      @test Ca == C([a, 0, 0, 0])
      @test Ca == C(K.([a, 0, 0, 0]))
      @test 2 * Ca == C([2 * a, 0, 0, 0])
      @test Ca * 2 == 2 * Ca
      @test 2 * Ca == C(2 * a)
      @test C(2) == C(K(2))
    end
    @testset "functions on elements" begin
      x = C([-1, 1, a, a + 1])
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == typeof(C())

      @test even_part(x) == C([-1, 0, 0, a + 1])
      @test odd_part(x) == C([0, 1, a, 0])
      xeven, xodd = even_coeff(x), odd_coeff(x)
      x.even_coeffs, x.odd_coeffs = Kzer, Kzer
      @test xeven != even_coeff(x) && xodd != odd_coeff(x)
      _set_even_odd_coeff!(x)
      @test xeven == even_coeff(x) && xodd == odd_coeff(x)

      @test +x == x
      @test -x == C([1, -1, -a, -(a + 1)])
      @test_throws ArgumentError divexact(x, 2)
      @test divexact(2 * x, O(2)) == x
      @test divexact(2 * x, K(2)) == x
    end
    @testset "conversion to ambient algebra and vice versa" begin
      CA = ambient_algebra(C)
      x = C([1, a, 1+a, 0])
      @test x == C(CA(x))
      @test CA(x) in C
      @test !(1//2 * CA(x) in C) 
    end
    @testset "defining relations" begin
      @test e(1) == C(matrix(gens(C))[1,:])
      @test e(2) == C(matrix(gens(C))[2,:])
      @test_throws BoundsError gen(C, 0)
      @test_throws BoundsError gen(C, 3)

      @test e(1)^2 == (e(1) * e(1)) && e(1)^2 == C(a)
      @test e(2)^2 == (e(2) * e(2)) && e(2)^2 == C(1 - a)
      @test e(1) * e(2) + e(2) * e(1) == C(1)

      @test e(1) * C(1) == e(1) && C(1) * e(1) == e(1)
      @test e(2) * C(1) == e(2) && C(1) * e(2) == e(2)
      @test is_zero(C(1)^2 - C(1))
      @test e(1) * e(2) == C(matrix(pseudo_basis(C, 4))[1,:])
      @test e(1) * e(2) == C([0, 0, 0, 1])
      @test e(2) * e(1) == C(gram_matrix(C)[1, 2]) - e(1) * e(2)
    end
    @testset "some computations" begin
      x, y = C([-2, 1, -1, 2]), C([2, 4, 6, 8])
      @test is_zero(x - x) && is_zero(x - +x)
      @test even_part(x) == C([-2, 0, 0, 2]) && odd_part(x) == C([0, 1, -1, 0])
      @test x + y == C([0, 5, 5, 10])
      @test x - y == C([-4, -3, -7, -6])
      @test y - x == C([4, 3, 7, 6])
      @test divexact(y, 2) == C([1, 2, 3, 4])
      @test divexact(y, ZZ(2)) == C([1, 2, 3, 4])
      @test divexact(y, 2//1) == C([1, 2, 3, 4])
      @test divexact(y, QQ(2//1)) == C([1, 2, 3, 4])
      @test_throws ArgumentError divexact(y, 4)

      #Cross-checking with MAGMAs results
      @test x^2 == C([8, -2, 2, -4])
      @test y^2 == C([-20 * a + 128, 48, 72, 96])
      @test x * y == C([10 * a + 2, -20 * a + 22, -22, 14])
      @test y * x == C([10 * a + 12, 20 * a - 18, -2, -6])
      @test x * (x + y) * y == x^2 * y + x * y^2
      @test (x + y)^2 == x^2 + x * y + y * x + y^2
      @test (x + y) * (x - y) == x^2 - x * y + y * x - y^2
    end
    @testset "mul_with_gen" begin
      x = C([-a^2, 2, a + 4, -1])
      @test mul_with_gen(coeff(x), 1, gram_matrix(C)) ==
        K.([3 * a + 4, -(a^2 + 1), a, -(a + 4)])
      @test mul_with_gen(coeff(x), 2, gram_matrix(C)) ==
        K.([-(a^2 + 3 * a - 4), -(1 - a), -a^2, 2])
    end
    @testset "center and centroid" begin
      orth = C([-1, 0, 0, 2])
      @test center(C) == pseudo_matrix(K[1 0 0 0], [ide])
      @test centroid(C) == [pseudo_matrix(K[1 0 0 0; 1 0 0 -1], [ide, ide])]
      @test disq(C) == quadratic_discriminant(C)
      @test disq(C) == (fractional_ideal(O, O(5)), K(5))
    end
  end

end
