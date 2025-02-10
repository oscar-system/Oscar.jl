
@testset "all test - Clifford orders" begin
  _set_even_odd_coefficients! = Oscar._set_even_odd_coefficients!
  mul_with_gen = Oscar._mul_with_gen
  @testset "failing constructions" begin
    @testset "CliffordOrder" begin
      K, a = quadratic_field(-5)
      O = maximal_order(K)
      qs = quadratic_space(K, K[1 2; 2 1])
      @test_throws ArgumentError clifford_order(lattice(qs))
      qs1 = quadratic_space(K, K[2 4; 4 2])
      pseu = pseudo_matrix(identity_matrix(K, 2), [fractional_ideal(O, one(K)), fractional_ideal(O, K(1//2))])
      @test_throws ArgumentError clifford_order(lattice(qs1, pseu))
    end
    @testset "ZZCliffordOrder" begin
      @test_throws ArgumentError clifford_order(root_lattice(:I, 1))
      qs = quadratic_space(QQ, QQ[1 2; 2 1])
      @test_throws ArgumentError clifford_order(lattice(qs))
      qs10 = quadratic_space(QQ, identity_matrix(QQ, 10))
      @test_throws ArgumentError clifford_order(lattice(qs10))
    end
  end
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
        @test typeof(C) == CliffordOrder{elem_type(base_ring_type(C)), typeof(algebra(C))}
        @test elem_type(C) == CliffordOrderElem{elem_type(base_ring_type(C)), typeof(algebra(C))}
        @test elem_type(C) == typeof(C())
        @test base_ring_type(C) == typeof(OK)
        @test base_ring_type(C) == typeof(base_ring(C))

        @test (base_ring(C), rank(C), lattice(C), gram_matrix(C), coefficient_ideals(C)) == (OK, 1, ls, idzero, [fractional_ideal(OK, OK(1))])
        @test C() == C(0) && C() == zero(C)
        @test is_zero(C())
        @test C(C()) == C()
        CC = clifford_order(ls)
        @test_throws ArgumentError CC(C())
        @test coefficients(C()) == Kzer
        @test even_coefficients(C()) == Kzer && even_coefficients(C()) == odd_coefficients(C())
        @test C(1) == one(C)
        @test is_one(C(1))
        @test coefficients(C(1)) == [K(1)]
        @test even_coefficients(C(1)) == [K(1)] && odd_coefficients(C(1)) == Kzer
      end
      @testset "equality and wrong parents" begin
        CC = clifford_order(ls)
        x, y = C(), CC()
        @test_throws ErrorException x + y
        @test_throws ErrorException x - y
        @test_throws ErrorException x * y
        @test x != y
        y = C(1)
        x[1] = K(1)
        @test x == y
      end
      @testset "functions on elements" begin
        x = C([17])
        @test parent(x) == C
        @test parent_type(x) == typeof(C)
        @test typeof(x) == typeof(C())

        @test even_part(x) == C([17])
        @test odd_part(x) == C([0])
        xeven, xodd = even_coefficients(x), odd_coefficients(x)
        x.even_coeffs, x.odd_coeffs = Kzer, [K(1)]
        @test xeven != x.even_coeffs
        @test xodd != x.odd_coeffs
        _set_even_odd_coefficients!(x)
        @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

        @test +x == x
        @test -x == C([-17])
        @test divexact(2 * x, 2) == x && divexact(2 * x, QQ(2)) == x
        @test divexact(C(18), OK(2)) == C(9)
        @test_throws ArgumentError divexact(x, 2)
        @test_throws ArgumentError divexact(x, ZZ(2))
        @test_throws ArgumentError divexact(x, QQ(2))
        @test_throws ArgumentError divexact(x, OK(2))
      end
      @testset "conversion to ambient algebra and vice versa" begin
        CA = algebra(C)
        x = C(a)
        @test x == C(CA(x))
        @test CA(x) == CA(C(CA(x)))
        @test CA(x) in C
        @test !(1//2 * CA(x) in C) 
      end
      @testset "defining relations" begin
        @test length(pseudo_basis(C)) == 1
        @test pseudo_basis(C, 1) == (one(C), fractional_ideal(OK, one(OK)))
        @test is_empty(pseudo_gens(C))
        @test_throws BoundsError pseudo_gen(C, 1)
        a, b = rand(ZZ, -100:100), rand(ZZ, -100:100)
        @test C(a) * C(b) == C(a * b)
        @test C(a) + C(b) == C(a + b)
        @test (C(a) + C(b)) * (C(a) - C(b)) == C(a)^2 - C(b)^2
      end
       @testset "getindex" begin
        @test C(17)[1] == K(17)
        @test coeff(C(17), 1) == K(17)
        @test C()[1] == K()
        @test coeff(C(), 1) == K()
        @test_throws BoundsError C()[0]
        @test_throws BoundsError coeff(C(), 0)
        @test_throws BoundsError C()[2]
        @test_throws BoundsError coeff(C(), 2)
      end
      @testset "center and centroid" begin
        @test pseudo_basis_of_center(C) == pseudo_basis_of_centroid(C)
        @test pseudo_basis_of_center(C) == [pseudo_basis(C, 1)]
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
        @test C(C()) == C()
        CC = clifford_order(ls)
        @test_throws ArgumentError CC(C())
        @test coefficients(C()) == QQzer
        @test even_coefficients(C()) == QQzer && even_coefficients(C()) == odd_coefficients(C())
        @test C(1) == one(C)
        @test is_one(C(1))
        @test coefficients(C(1)) == [QQ(1)]
        @test even_coefficients(C(1)) == [QQ(1)] && odd_coefficients(C(1)) == QQzer
      end
      @testset "equality and wrong parents" begin
        CC = clifford_order(ls)
        x, y = C(), CC()
        @test_throws ErrorException x + y
        @test_throws ErrorException x - y
        @test_throws ErrorException x * y
        @test x != y
        y = C(1)
        x[1] = QQ(1)
        @test x == y
      end
      @testset "functions on elements" begin
        x = C([17])
        @test parent(x) == C
        @test parent_type(x) == typeof(C)
        @test typeof(x) == typeof(C())

        @test even_part(x) == C([17])
        @test odd_part(x) == C([0])
        xeven, xodd = even_coefficients(x), odd_coefficients(x)
        x.even_coeffs, x.odd_coeffs = QQzer, [QQ(1)]
        @test xeven != even_coefficients(x)
        @test xodd != odd_coefficients(x)
        _set_even_odd_coefficients!(x)
        @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

        @test +x == x
        @test -x == C([-17])
        @test divexact(2 * x, 2) == x && divexact(2 * x, QQ(2)) == x
        @test divexact(C(18), ZZ(2)) == C(9)
        @test_throws ArgumentError divexact(x, 2)
        @test_throws ArgumentError divexact(x, ZZ(2))
        @test_throws ArgumentError divexact(x, QQ(2))
        @test_throws ArgumentError divexact(x, 2//1)
      end
      @testset "conversion to ambient algebra and vice versa" begin
        CA = algebra(C)
        x = C(17)
        @test x == C(CA(x))
        @test CA(x) == CA(C(CA(x)))
        @test CA(x) in C
        @test !(1//2 * CA(x) in C) 
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
      @testset "getindex" begin
        @test C(17)[1] == QQ(17)
        @test coeff(C(17), 1) == QQ(17)
        @test C()[1] == QQ()
        @test coeff(C(), 1) == QQ()
        @test_throws BoundsError C()[0]
        @test_throws BoundsError coeff(C(), 0)
        @test_throws BoundsError C()[2]
        @test_throws BoundsError coeff(C(), 2)
      end
      @testset "setindex!" begin
        x = C()
        x[1] = QQ(1)
        @test x == C(1)
        x[1] = QQ()
        @test x == C()
        @test_throws BoundsError x[0] = QQ()
        @test_throws BoundsError x[2] = QQ()
      end
      @testset "center and centroid" begin
        @test basis_of_center(C) == basis_of_centroid(C)
        @test basis_of_center(C) == [one(C)]
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
    e(i::Int) = pseudo_gen(C, i)[1]
    
    @test !is_commutative(C)
    
    @testset "construction" begin
      @test typeof(C) == CliffordOrder{elem_type(base_ring_type(C)), typeof(algebra(C))}
      @test elem_type(C) ==
      CliffordOrderElem{elem_type(base_ring_type(C)), typeof(algebra(C))}
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == typeof(base_ring(C))

      @test (base_ring(C), gram_matrix(C), lattice(C), rank(C), coefficient_ideals(C)) == (O, G, lsK, 4, [ide, ide, ide, ide])
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test C(C()) == C()
      CC = clifford_order(lsK)
      @test_throws ArgumentError CC(C())
      @test coefficients(C()) == Kzer
      @test even_coefficients(C()) == Kzer && even_coefficients(C()) == odd_coefficients(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coefficients(C(1)) == K.([1, 0, 0, 0])
      @test even_coefficients(C(1)) == K.([1, 0, 0, 0]) && odd_coefficients(C(1)) == Kzer
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
    @testset "equality and wrong parents" begin
      CC = clifford_order(lsK)
      x, y = C(), CC()
      @test_throws ErrorException x + y
      @test_throws ErrorException x - y
      @test_throws ErrorException x * y
      @test x != y
      y = C(1)
      x[1] = K(1)
      @test x == y
    end
    @testset "functions on elements" begin
      x = C([-1, 1, a, a + 1])
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == typeof(C())

      @test even_part(x) == C([-1, 0, 0, a + 1])
      @test odd_part(x) == C([0, 1, a, 0])
      xeven, xodd = even_coefficients(x), odd_coefficients(x)
      x.even_coeffs, x.odd_coeffs = Kzer, Kzer
      @test xeven != x.even_coeffs && xodd != x.odd_coeffs
      _set_even_odd_coefficients!(x)
      @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

      @test +x == x
      @test -x == C([1, -1, -a, -(a + 1)])
      @test_throws ArgumentError divexact(x, 2)
      @test divexact(2 * x, O(2)) == x
      @test divexact(2 * x, K(2)) == x
    end
    @testset "conversion to ambient algebra and vice versa" begin
      CA = algebra(C)
      x = C([1, a, 1+a, 0])
      @test x == C(CA(x))
      @test CA(x) == CA(C(CA(x)))
      @test CA(x) in C
      @test !(1//2 * CA(x) in C) 
    end
    @testset "defining relations" begin
      @test length(pseudo_gens(C)) == 2
      @test e(1) == pseudo_gens(C)[1][1]
      @test e(2) == pseudo_gens(C)[2][1]
      @test_throws BoundsError pseudo_gen(C, 0)
      @test_throws BoundsError pseudo_gen(C, 3)

      @test e(1)^2 == (e(1) * e(1)) && e(1)^2 == C(a)
      @test e(2)^2 == (e(2) * e(2)) && e(2)^2 == C(1 - a)
      @test e(1) * e(2) + e(2) * e(1) == C(1)

      @test e(1) * C(1) == e(1) && C(1) * e(1) == e(1)
      @test e(2) * C(1) == e(2) && C(1) * e(2) == e(2)
      @test is_zero(C(1)^2 - C(1))
      @test e(1) * e(2) == pseudo_basis(C, 4)[1]
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
    @testset "getindex" begin
      x = C([a, 2 * a, 3 * a, 4 * a])
      @test x[1] == a
      @test coeff(x, 1) == a
      @test x[2] == 2 * a
      @test coeff(x, 2) == 2 * a
      @test x[3] == 3 * a
      @test coeff(x, 3) == 3 * a
      @test x[4] == 4 * a
      @test coeff(x, 4) == 4 * a
      @test_throws BoundsError x[0]
      @test_throws BoundsError coeff(x, 0)
      @test_throws BoundsError x[5]
      @test_throws BoundsError coeff(x, 5)
    end
    @testset "setindex!" begin
      x = C([a, 2 * a, 3 * a, 4 * a])
      for i in 1:4
        x[i] = K()
      end
      x == C()
    end
    @testset "mul_with_gen" begin
      x = C([-a^2, 2, a + 4, -1])
      @test mul_with_gen(coefficients(x), 1, gram_matrix(C)) ==
        K.([3 * a + 4, -(a^2 + 1), a, -(a + 4)])
      @test mul_with_gen(coefficients(x), 2, gram_matrix(C)) ==
        K.([-(a^2 + 3 * a - 4), -(1 - a), -a^2, 2])
    end
    @testset "center and centroid" begin
      @test pseudo_basis_of_center(C) == [pseudo_basis(C, 1)]
      @test pseudo_basis_of_centroid(C) == [pseudo_basis(C, 1), (C([1,0,0,-1]), ide)]
      @test disq(C) == quadratic_discriminant(C)
      @test disq(C) == (fractional_ideal(O, O(5)), K(5))
    end
  end
  
  @testset "hyperbolic plane" begin
    @testset "CliffordOrder" begin
      K, a = quadratic_field(-5)
      O = maximal_order(K)
      ide = fractional_ideal(O, one(O))
      G = K[0 1; 1 0]
      qsK = quadratic_space(K, G)
      lsK = lattice(qsK)
      C = clifford_order(lsK)
      Kzer = K.([0, 0, 0, 0])
      e(i::Int) = pseudo_gen(C, i)[1]
    
      @test !is_commutative(C)
    
      @testset "construction" begin
        @test typeof(C) == CliffordOrder{elem_type(base_ring_type(C)), typeof(algebra(C))}
        @test elem_type(C) ==
        CliffordOrderElem{elem_type(base_ring_type(C)), typeof(algebra(C))}
        @test elem_type(C) == typeof(C())
        @test base_ring_type(C) == typeof(base_ring(C))

        @test (base_ring(C), gram_matrix(C), lattice(C), rank(C), coefficient_ideals(C)) == (O, G, lsK, 4, [ide, ide, ide, ide])
        @test C() == C(0) && C() == zero(C)
        @test is_zero(C())
        @test coefficients(C()) == Kzer
        @test even_coefficients(C()) == Kzer && even_coefficients(C()) == odd_coefficients(C())
        @test C(1) == one(C)
        @test is_one(C(1))
        @test coefficients(C(1)) == K.([1, 0, 0, 0])
        @test even_coefficients(C(1)) == K.([1, 0, 0, 0]) && odd_coefficients(C(1)) == Kzer
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
      @testset "equality and wrong parents" begin
        CC = clifford_order(lsK)
        x, y = C(), CC()
        @test_throws ErrorException x + y
        @test_throws ErrorException x - y
        @test_throws ErrorException x * y
        @test x != y
        y = C(1)
        x[1] = K(1)
        @test x == y
      end
      @testset "functions on elements" begin
        x = C([-1, 1, a, a + 1])
        @test parent(x) == C
        @test parent_type(x) == typeof(C)
        @test typeof(x) == typeof(C())

        @test even_part(x) == C([-1, 0, 0, a + 1])
        @test odd_part(x) == C([0, 1, a, 0])
        xeven, xodd = even_coefficients(x), odd_coefficients(x)
        x.even_coeffs, x.odd_coeffs = Kzer, Kzer
        @test xeven != x.even_coeffs && xodd != x.odd_coeffs
        _set_even_odd_coefficients!(x)
        @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

        @test +x == x
        @test -x == C([1, -1, -a, -(a + 1)])
        @test_throws ArgumentError divexact(x, 2)
        @test divexact(2 * x, O(2)) == x
        @test divexact(2 * x, K(2)) == x
      end
      @testset "conversion to ambient algebra and vice versa" begin
        CA = algebra(C)
        x = C([1, a, 1+a, 0])
        @test x == C(CA(x))
        @test CA(x) == CA(C(CA(x)))
        @test CA(x) in C
        @test !(1//2 * CA(x) in C) 
      end
      @testset "defining relations" begin
        @test length(pseudo_gens(C)) == 2
        @test e(1) == pseudo_gens(C)[1][1]
        @test e(2) == pseudo_gens(C)[2][1]
        @test_throws BoundsError pseudo_gen(C, 0)
        @test_throws BoundsError pseudo_gen(C, 3)

        @test e(1)^2 == (e(1) * e(1)) && e(1)^2 == C()
        @test e(2)^2 == (e(2) * e(2)) && e(2)^2 == C()
        @test e(1) * e(2) + e(2) * e(1) == C(1)

        @test e(1) * C(1) == e(1) && C(1) * e(1) == e(1)
        @test e(2) * C(1) == e(2) && C(1) * e(2) == e(2)
        @test is_zero(C(1)^2 - C(1))
        @test e(1) * e(2) == pseudo_basis(C, 4)[1]
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

        @test x * (x + y) * y == x^2 * y + x * y^2
        @test (x + y)^2 == x^2 + x * y + y * x + y^2
        @test (x + y) * (x - y) == x^2 - x * y + y * x - y^2
      end
      @testset "getindex" begin
        x = C([a, 2 * a, 3 * a, 4 * a])
        @test x[1] == a
        @test coeff(x, 1) == a
        @test x[2] == 2 * a
        @test coeff(x, 2) == 2 * a
        @test x[3] == 3 * a
        @test coeff(x, 3) == 3 * a
        @test x[4] == 4 * a
        @test coeff(x, 4) == 4 * a
        @test_throws BoundsError x[0]
        @test_throws BoundsError coeff(x, 0)
        @test_throws BoundsError x[5]
        @test_throws BoundsError coeff(x, 5)
      end
      @testset "setindex!" begin
        x = C([a, 2 * a, 3 * a, 4 * a])
        for i in 1:4
          x[i] = K()
        end
        x == C()
      end
      @testset "center and centroid" begin
        @test pseudo_basis_of_center(C) == [pseudo_basis(C, 1)]
        @test pseudo_basis_of_centroid(C) == [pseudo_basis(C, 1), (C([1,0,0,-1]), ide)]
        @test disq(C) == quadratic_discriminant(C)
        @test disq(C) == (ide, K(1))
      end
    end
    @testset "ZZCliffordOrder" begin
      G = QQ[0 1; 1 0]
      qs = quadratic_space(QQ, G)
      ls = lattice(qs)
      C = clifford_order(ls)
      QQzer = QQ.([0, 0, 0, 0])
      e(i::Int) = gen(C, i)
    
      @test !is_commutative(C)
    
      @testset "construction" begin
        @test typeof(C) == ZZCliffordOrder
        @test elem_type(C) == ZZCliffordOrderElem
        @test elem_type(C) == typeof(C())
        @test base_ring_type(C) == typeof(base_ring(C))

        @test (base_ring(C), gram_matrix(C), lattice(C), rank(C)) == (ZZ, G, ls, 4)
        @test C() == C(0) && C() == zero(C)
        @test is_zero(C())
        @test C(C()) == C()
        CC = clifford_order(ls)
        @test_throws ArgumentError CC(C())
        @test coefficients(C()) == QQzer
        @test even_coefficients(C()) == QQzer && even_coefficients(C()) == odd_coefficients(C())
        @test C(1) == one(C)
        @test is_one(C(1))
        @test coefficients(C(1)) == QQ.([1, 0, 0, 0])
        @test even_coefficients(C(1)) == QQ.([1, 0, 0, 0]) && odd_coefficients(C(1)) == QQzer
        verr1, verr2 = QQ.([1, 2, 3]), QQ.([1, 2, 3, 4, 5])
        @test_throws ArgumentError C(verr1)
        @test_throws ArgumentError C(verr2)
      end
      @testset "equality and wrong parents" begin
        CC = clifford_order(ls)
        x, y = C(), CC()
        @test_throws ErrorException x + y
        @test_throws ErrorException x - y
        @test_throws ErrorException x * y
        @test x != y
        y = C(1)
        x[1] = QQ(1)
        @test x == y
      end
      @testset "functions on elements" begin
        x = C([-1, 1, 0, 8])
        @test parent(x) == C
        @test parent_type(x) == typeof(C)
        @test typeof(x) == typeof(C())

        @test even_part(x) == C([-1, 0, 0, 8])
        @test odd_part(x) == C([0, 1, 0, 0])
        xeven, xodd = even_coefficients(x), odd_coefficients(x)
        x.even_coeffs, x.odd_coeffs = QQzer, QQzer
        @test xeven != x.even_coeffs && xodd != x.odd_coeffs
        _set_even_odd_coefficients!(x)
        @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

        @test +x == x
        @test -x == C([1, -1, 0, -8])
        @test_throws ArgumentError divexact(x, 2)
        @test divexact(2 * x, 2) == x
        @test divexact(2 * x, ZZ(2)) == x
      end
      @testset "conversion to ambient algebra and vice versa" begin
        CA = algebra(C)
        x = C([1, 8, 9, 0])
        @test x == C(CA(x))
        @test CA(x) == CA(C(CA(x)))
        @test CA(x) in C
        @test !(1//2 * CA(x) in C) 
      end
      @testset "defining relations" begin
        @test length(gens(C)) == 2
        @test e(1) == gens(C)[1]
        @test e(2) == gens(C)[2]
        @test_throws BoundsError gen(C, 0)
        @test_throws BoundsError gen(C, 3)

        @test e(1)^2 == (e(1) * e(1)) && e(1)^2 == C()
        @test e(2)^2 == (e(2) * e(2)) && e(2)^2 == C()
        @test e(1) * e(2) + e(2) * e(1) == C(1)

        @test e(1) * C(1) == e(1) && C(1) * e(1) == e(1)
        @test e(2) * C(1) == e(2) && C(1) * e(2) == e(2)
        @test is_zero(C(1)^2 - C(1))
        @test e(1) * e(2) == basis(C, 4)
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

        @test x * (x + y) * y == x^2 * y + x * y^2
        @test (x + y)^2 == x^2 + x * y + y * x + y^2
        @test (x + y) * (x - y) == x^2 - x * y + y * x - y^2
      end
      @testset "getindex" begin
        x = C([1, 2, 3, 4])
        @test x[1] == QQ(1)
        @test coeff(x, 1) == QQ(1)
        @test x[2] == QQ(2)
        @test coeff(x, 2) == QQ(2)
        @test x[3] == QQ(3)
        @test coeff(x, 3) == QQ(3)
        @test x[4] == QQ(4)
        @test coeff(x, 4) == QQ(4)
        @test_throws BoundsError x[0]
        @test_throws BoundsError coeff(x, 0)
        @test_throws BoundsError x[5]
        @test_throws BoundsError coeff(x, 5)
      end
      @testset "setindex!" begin
        x = C([1, 2, 3, 4])
        for i in 1:4
          x[i] = QQ()
        end
        x == C()
      end
      @testset "center and centroid" begin
        @test basis_of_center(C) == [basis(C, 1)]
        @test basis_of_centroid(C) == [basis(C, 1), C([0,0,0,1])]
        @test disq(C) == quadratic_discriminant(C)
        @test disq(C) == ZZ(1)
      end
    end
  end

  @testset "ZZCliffordE6" begin
    e6 = root_lattice(:E, 6)
    C = clifford_order(e6)
    QQzer = fill(zero(QQ), 2^6)
    QQone = copy(QQzer)
    QQone[1] = one(QQ)
    @test !is_commutative(C)
    @testset "Construction" begin
      @test typeof(C) == ZZCliffordOrder
      @test elem_type(C) == ZZCliffordOrderElem
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == ZZRing
      @test base_ring_type(C) == typeof(base_ring(C))

      @test (base_ring(C), rank(C), lattice(C), gram_matrix(C)) == (ZZ, 2^6, e6, gram_matrix(e6))
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coefficients(C()) == QQzer
      @test even_coefficients(C()) == QQzer && even_coefficients(C()) == odd_coefficients(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coefficients(C(1)) == QQone
      @test even_coefficients(C(1)) == QQone && odd_coefficients(C(1)) == QQzer
    end
    @testset "equality and wrong parents" begin
      CC = clifford_order(e6)
      x, y = C(), CC()
      @test_throws ErrorException x + y
      @test_throws ErrorException x - y
      @test_throws ErrorException x * y
      @test x != y
      y = C(1)
      x[1] = QQ(1)
      @test x == y
    end
    @testset "functions on elements" begin
      x = gen(C, 1) * gen(C, 2)
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == ZZCliffordOrderElem
      @test typeof(x) == typeof(C())

      @test even_part(x) == x
      @test odd_part(x) == zero(C)
      xeven, xodd = even_coefficients(x), odd_coefficients(x)
      x.even_coeffs, x.odd_coeffs = QQzer, QQone
      @test xeven != even_coefficients(x)
      @test xodd != odd_coefficients(x)
      _set_even_odd_coefficients!(x)
      @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

      ngtv = copy(QQzer)
      ngtv[4] = -1
      @test +x == x
      @test -x == C(ngtv)
      @test divexact(2 * x, 2) == x && divexact(2 * x, QQ(2)) == x
      @test divexact(2 * x, ZZ(-2)) == C(ngtv)
      @test_throws ArgumentError divexact(x, 2)
      @test_throws ArgumentError divexact(x, ZZ(2))
      @test_throws ArgumentError divexact(x, QQ(2))
      @test_throws ArgumentError divexact(x, 2//1)
    end
    @testset "conversion to ambient algebra and vice versa" begin
      CA = algebra(C)
      x = C(QQ.(1:64))
      @test x == C(CA(x))
      @test CA(x) == CA(C(CA(x)))
      @test CA(x) in C
      @test !(1//2 * CA(x) in C) 
    end
    @testset "defining relations" begin
      @test length(basis(C)) == 2^6
      @test basis(C, 1) == one(C)
      @test length(gens(C)) == 6
      a, b = rand(ZZ, -100:100), rand(ZZ, -100:100)
      @test C(a) * C(b) == C(a * b)
      @test C(a) + C(b) == C(a + b)
      @test (C(a) + C(b)) * (C(a) - C(b)) == C(a)^2 - C(b)^2
      
      e(i) = gen(C, i)
      for i in 1:6
        @test e(i)^2 == C(divexact(gram_matrix(C)[i, i], 2))
      end
      
      for i in 2:6
        for j in 1:i-1
          @test e(i) * e(j) == C(gram_matrix(C)[i, j]) - e(j) * e(i)
        end
      end

      for i in 1:6
        @test e(i) == basis(C, 2^(i - 1) + 1)
      end
    end
    @testset "getindex" begin
      x = C(QQ.(1:64))
      @test typeof(x) == typeof(C())
      for i in 1:64
        @test coeff(x, i) == i
        @test x[i] == i
      end
      @test_throws BoundsError x[0]
      @test_throws BoundsError coeff(x, 0)
      @test_throws BoundsError x[65]
      @test_throws BoundsError coeff(x, 65)
      @test_throws BoundsError x[rand(65:500)]
    end
    @testset "center and centroid" begin
      @test basis_of_center(C) == [one(C)]
      z_elt = copy(QQzer)
      for i in [1, 4, 25, 34, 37, 49]
        z_elt[i] = QQ(1)
      end
      for i in [28, 40, 52, 58, 61]
        z_elt[i] = QQ(2)
      end
      z_elt[64] = QQ(4)
      z_elt = C(z_elt)
      @test even_part(z_elt) == z_elt
      @test basis_of_centroid(C) == [one(C), z_elt]
      @test is_zero(z_elt^2 - z_elt + one(C))
      orth = 2 * z_elt - one(C)
      @test disq(C) == quadratic_discriminant(C)
      @test disq(C) == ZZ(coefficients(orth^2)[1])
      @test disq(C) == ZZ(-3)
    end
  end
  @testset "ZZCliffordA5" begin
    a5 = root_lattice(:A, 5)
    C = clifford_order(a5)
    QQzer = fill(zero(QQ), 2^5)
    QQone = copy(QQzer)
    QQone[1] = one(QQ)
    @test !is_commutative(C)
    @testset "Construction" begin
      @test typeof(C) == ZZCliffordOrder
      @test elem_type(C) == ZZCliffordOrderElem
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == ZZRing
      @test base_ring_type(C) == typeof(base_ring(C))

      @test (base_ring(C), rank(C), lattice(C), gram_matrix(C)) == (ZZ, 2^5, a5, gram_matrix(a5))
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coefficients(C()) == QQzer
      @test even_coefficients(C()) == QQzer && even_coefficients(C()) == odd_coefficients(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coefficients(C(1)) == QQone
      @test even_coefficients(C(1)) == QQone && odd_coefficients(C(1)) == QQzer
    end
    @testset "equality and wrong parents" begin
      CC = clifford_order(a5)
      x, y = C(), CC()
      @test_throws ErrorException x + y
      @test_throws ErrorException x - y
      @test_throws ErrorException x * y
      @test x != y
      y = C(1)
      x[1] = QQ(1)
      @test x == y
    end
    @testset "functions on elements" begin
      x = gen(C, 1) * gen(C, 2)
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == ZZCliffordOrderElem
      @test typeof(x) == typeof(C())

      @test even_part(x) == x
      @test odd_part(x) == zero(C)
      xeven, xodd = even_coefficients(x), odd_coefficients(x)
      x.even_coeffs, x.odd_coeffs = QQzer, QQone
      @test xeven != even_coefficients(x)
      @test xodd != odd_coefficients(x)
      _set_even_odd_coefficients!(x)
      @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

      ngtv = copy(QQzer)
      ngtv[4] = -1
      @test +x == x
      @test -x == C(ngtv)
      @test divexact(2 * x, 2) == x && divexact(2 * x, QQ(2)) == x
      @test divexact(2 * x, ZZ(-2)) == C(ngtv)
      @test_throws ArgumentError divexact(x, 2)
      @test_throws ArgumentError divexact(x, ZZ(2))
      @test_throws ArgumentError divexact(x, QQ(2))
      @test_throws ArgumentError divexact(x, 2//1)
    end
    @testset "conversion to ambient algebra and vice versa" begin
      CA = algebra(C)
      x = C(QQ.(1:32))
      @test x == C(CA(x))
      @test CA(x) == CA(C(CA(x)))
      @test CA(x) in C
      @test !(1//2 * CA(x) in C) 
    end
    @testset "defining relations" begin
      @test length(basis(C)) == 2^5
      @test basis(C, 1) == one(C)
      @test length(gens(C)) == 5
      a, b = rand(ZZ, -100:100), rand(ZZ, -100:100)
      @test C(a) * C(b) == C(a * b)
      @test C(a) + C(b) == C(a + b)
      @test (C(a) + C(b)) * (C(a) - C(b)) == C(a)^2 - C(b)^2
      
      e(i) = gen(C, i)
      for i in 1:5
        @test e(i)^2 == C(divexact(gram_matrix(C)[i, i], 2))
      end
      
      for i in 2:5
        for j in 1:i-1
          @test e(i) * e(j) == C(gram_matrix(C)[i, j]) - e(j) * e(i)
        end
      end

      for i in 1:5
        @test e(i) == basis(C, 2^(i - 1) + 1)
      end
    end
    @testset "getindex" begin
      x = C(QQ.(1:32))
      @test typeof(x) == typeof(C())
      for i in 1:32
        @test x[i] == i
      end
      @test_throws BoundsError x[0]
      @test_throws BoundsError x[33]
      @test_throws BoundsError x[rand(33:500)]
    end     
    @testset "center and centroid" begin
      @test basis_of_center(C) == basis_of_centroid(C)
      z_elt = copy(QQzer)
      for i in [2, 5, 17]
        z_elt[i] = QQ(1)
      end
      for i in [8, 20, 26, 29]
        z_elt[i] = QQ(2)
      end
      z_elt[32] = QQ(4)
      z_elt = C(z_elt)
      @test odd_part(z_elt) == z_elt
      @test basis_of_centroid(C) == [one(C), z_elt]
      @test disq(C) == quadratic_discriminant(C)
      @test disq(C) == ZZ(coefficients(z_elt^2)[1])
      @test disq(C) == ZZ(3)
    end
  end

end
