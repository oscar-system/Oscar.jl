
@testset "all tests - Clifford algebras" begin
  _dim_qf = Oscar._dim_qf
  _set_even_odd_coefficients! = Oscar._set_even_odd_coefficients! 
  mul_with_gen = Oscar._mul_with_gen 

  @testset "zero-dimensional corner case" begin
    empty_qs = quadratic_space(QQ, identity_matrix(QQ, 0)) #zero-dim quad space over QQ
    C = clifford_algebra(empty_qs)
    QQzer = [QQ(0)]
    @test is_commutative(C)
    @testset "construction" begin
      @test typeof(C) ==
        CliffordAlgebra{typeof(base_ring(C)()),typeof(gram_matrix(empty_qs))}
      @test elem_type(C) ==
        CliffordAlgebraElem{typeof(base_ring(C)()),typeof(gram_matrix(C))}
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == QQField

      @test (base_ring(C), space(C), gram_matrix(C), _dim_qf(C), dim(C)) ==
        (QQ, empty_qs, identity_matrix(QQ, 0), 0, 1)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test C(C()) == C()
      CC = clifford_algebra(empty_qs)
      @test_throws ArgumentError C(CC())
      @test coefficients(C()) == QQzer
      @test even_coefficients(C()) == QQzer && even_coefficients(C()) == odd_coefficients(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coefficients(C(1)) == [QQ(1)]
      @test even_coefficients(C(1)) == [QQ(1)] && odd_coefficients(C(1)) == QQzer
    end
    @testset "equality and wrong parents" begin
      CC = clifford_algebra(empty_qs)
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
      @test divexact(C(18), 2) == C(9)
    end
    @testset "defining relations" begin
      @test length(basis(C)) == 1
      @test basis(C, 1) == C(1)
      @test is_empty(gens(C))
      a, b = rand(QQ, -100:100), rand(QQ, -100:100)
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
    @testset "center and centroid" begin
      @test basis_of_center(C) == basis_of_centroid(C)
      @test basis_of_center(C) == [one(C)]
      @test disq(C) == quadratic_discriminant(C)
      @test disq(C) == 1
    end
  end

  @testset "quadratic_field_sqrt(5)" begin
    K, b = quadratic_field(5)
    a = 1//2 * (1 + b)
    G = K[2 * a 1; 1 2 * (1 - a)]
    qsK = quadratic_space(K, G)
    C = clifford_algebra(qsK)
    @test !is_commutative(C)
    e(i::Int) = gen(C, i) #e(i) returns the i-th generator of C
    Kzer = K.([0, 0, 0, 0])
    @testset "construction" begin
      @test typeof(C) == CliffordAlgebra{typeof(base_ring(C)()),typeof(gram_matrix(qsK))}
      @test elem_type(C) ==
        CliffordAlgebraElem{typeof(base_ring(C)()),typeof(gram_matrix(qsK))}
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == typeof(K)

      @test (base_ring(C), gram_matrix(C), _dim_qf(C), dim(C)) == (K, G, 2, 4)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test C(C()) == C()
      CC = clifford_algebra(qsK)
      @test_throws ArgumentError C(CC())
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
      CC = clifford_algebra(qsK)
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
      @test xeven != even_coefficients(x) && xodd != odd_coefficients(x)
      _set_even_odd_coefficients!(x)
      @test xeven == even_coefficients(x) && xodd == odd_coefficients(x)

      @test +x == x
      @test -x == C([1, -1, -a, -(a + 1)])
      @test divexact(x, 2) == C([-1//2, 1//2, a//2, (a + 1)//2])
      @test divexact(2 * x, 2) == x && divexact(2 * x, K(2)) == x
    end
    @testset "defining relations" begin
      @test e(1) == gens(C)[1]
      @test e(2) == gens(C)[2]
      @test_throws BoundsError gen(C, 0)
      @test_throws BoundsError gen(C, 3)

      @test e(1)^2 == (e(1) * e(1)) && e(1)^2 == C(a)
      @test e(2)^2 == (e(2) * e(2)) && e(2)^2 == C(1 - a)
      @test e(1) * e(2) + e(2) * e(1) == C(1)

      @test e(1) * C(1) == e(1) && C(1) * e(1) == e(1)
      @test e(2) * C(1) == e(2) && C(1) * e(2) == e(2)
      @test is_zero(C(1)^2 - C(1))
      @test e(1) * e(2) == basis(C, 4) && basis(C, 4) == C([0, 0, 0, 1])
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
      @test divexact(y, 4) == C([2//4, 4//4, 6//4, 8//4])

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
      @test x[2] == 2 * a
      @test x[3] == 3 * a
      @test x[4] == 4 * a
      @test_throws BoundsError x[0]
      @test_throws BoundsError x[5]
    end
    @testset "mul_with_gen" begin
      x = C([-a^2, 2, a + 4, -1])
      @test mul_with_gen(coefficients(x), 1, gram_matrix(C)) ==
        K.([3 * a + 4, -(a^2 + 1), a, -(a + 4)])
      @test mul_with_gen(coefficients(x), 2, gram_matrix(C)) ==
        K.([-(a^2 + 3 * a - 4), -(1 - a), -a^2, 2])
    end
    @testset "center and centroid" begin
      orth = C([-1, 0, 0, 2])
      @test basis_of_center(C) == [one(C)]
      @test basis_of_centroid(C) == [one(C), orth]
      @test disq(C) == quadratic_discriminant(C) && disq(C) == coefficients(orth^2)[1]
    end
  end

  @testset "random quaternion alg over rationals" begin
    a, b = QQ(0), QQ(0)
    while is_zero(a * b)
      a, b = rand(QQ, -5:5), rand(QQ, -5:5)
    end
    G = diagonal_matrix(2 * a, 2 * b)
    qs = quadratic_space(QQ, G)
    C = clifford_algebra(qs)
    @test !is_commutative(C)
    e(i::Int) = gen(C, i) #e(i) returns the i-th generator of C
    zer = QQ.([0, 0, 0, 0])
    @testset "construction" begin
      @test typeof(C) == CliffordAlgebra{typeof(base_ring(C)()),typeof(gram_matrix(qs))}
      @test elem_type(C) ==
        CliffordAlgebraElem{typeof(base_ring(C)()),typeof(gram_matrix(qs))}
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == QQField

      @test (base_ring(C), gram_matrix(C), _dim_qf(C), dim(C)) == (QQ, G, 2, 4)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coefficients(C()) == zer
      @test even_coefficients(C()) == zer && even_coefficients(C()) == odd_coefficients(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coefficients(C(1)) == QQ.([1, 0, 0, 0])
      @test even_coefficients(C(1)) == QQ.([1, 0, 0, 0]) && odd_coefficients(C(1)) == zer
    end
    @testset "equality and wrong parents" begin
      CC = clifford_algebra(qs)
      x, y = C(), CC()
      @test x != y
      @test_throws ErrorException x + y
      @test_throws ErrorException x - y
      @test_throws ErrorException x * y
      y = C(1)
      x[1] = QQ(1)
      @test x == y 
    end
    @testset "defining relations" begin
      @test e(1) == gens(C)[1]
      @test e(2) == gens(C)[2]

      @test e(1)^2 == (e(1) * e(1)) && e(1)^2 == C(a)
      @test e(2)^2 == (e(2) * e(2)) && e(2)^2 == C(b)
      @test e(1) * e(2) + e(2) * e(1) == C(zer)

      @test e(1) * C(1) == e(1) && C(1) * e(1) == e(1)
      @test e(2) * C(1) == e(2) && C(1) * e(2) == e(2)
      @test is_zero(C(1)^2 - C(1))
      @test e(1) * e(2) == basis(C, 4) && basis(C, 4) == C([0, 0, 0, 1])
      @test e(2) * e(1) == C(gram_matrix(C)[1, 2]) - e(1) * e(2)
    end
    @testset "some computations" begin
      @test e(2) * (e(1) + e(2)) == C([b, 0, 0, -1])
      @test (1//2 * (1 + e(1)))^2 == 1//4 * C([1 + a, 2, 0, 0])
      @test (e(1) * e(2))^2 == -C(a * b)
      @test C(a * b) == C(a) * C(b)
      @test (e(1) + e(2))^2 == C(a) + C(b)
      @test e(1)^3 == a * e(1)
      @test e(2)^3 == b * e(2)

      x = C([a + b, a, b, a - b])
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == typeof(C())

      @test x^2 == C([
        (a + b)^2 + a^3 + b^3 - (a - b)^2 * a * b,
        2 * a * (a + b),
        2 * b * (a + b),
        2 * (a^2 - b^2),
      ])
      @test x == +x
      @test -x == C(-coefficients(x))
      @test is_zero(x + -x)
    end
    @testset "getindex" begin
      x = C([a, b, a + b, a - b])
      @test x[1] == a
      @test x[2] == b
      @test x[3] == a + b
      @test x[4] == a - b
      @test_throws BoundsError x[0]
      @test_throws BoundsError x[5]
    end
    @testset "mul_with_gen" begin
      x = C([1, 2, 3, 4])
      @test mul_with_gen(coefficients(x), 1, gram_matrix(C)) == QQ.([2 * a, 1, -4 * a, -3])
      @test mul_with_gen(coefficients(x), 2, gram_matrix(C)) == QQ.([3 * b, 4 * b, 1, 2])
    end
    @testset "center and centroid" begin
      @test basis_of_center(C) == [one(C)]
      @test basis_of_centroid(C) == [one(C), C(QQ.([0, 0, 0, 1]))]
      @test disq(C) == quadratic_discriminant(C) && disq(C) == -a * b
      @test disq(C) == coefficients(basis_of_centroid(C)[2]^2)[1]
    end
  end

  @testset "cyclic number field, bigger example" begin
    K, z = cyclotomic_field(5)
    G = 2 * identity_matrix(K, 5)
    for i in 1:5, j in 1:5
      if i != j
        G[i, j] = z^(abs(i - j))
      end
    end
    qs = quadratic_space(K, G)
    C = clifford_algebra(qs)
    @test !is_commutative(C)
    e(i::Int) = gen(C, i) #e(i) returns the i-th generator of C
    zer = fill(K(), 2^5)
    ein = fill(K(), 2^5)
    ein[1] = K(1)
    @testset "construction" begin
      @test typeof(C) == CliffordAlgebra{typeof(base_ring(C)()),typeof(gram_matrix(qs))}
      @test elem_type(C) ==
        CliffordAlgebraElem{typeof(base_ring(C)()),typeof(gram_matrix(qs))}
      @test elem_type(C) == typeof(C())
      @test base_ring_type(C) == typeof(K)

      @test (base_ring(C), gram_matrix(C), _dim_qf(C), dim(C)) == (K, G, 5, 32)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coefficients(C()) == zer
      @test even_coefficients(C()) == zer && even_coefficients(C()) == odd_coefficients(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coefficients(C(1)) == ein
      @test even_coefficients(C(1)) == ein && odd_coefficients(C(1)) == zer
    end
    @testset "equality and wrong parents" begin
      CC = clifford_algebra(qs)
      x, y = C(), CC()
      @test_throws ErrorException x + y
      @test_throws ErrorException x - y
      @test_throws ErrorException x * y
      @test x != y
      y = C(1)
      x[1] = K(1)
      @test x == y 
    end
    @testset "defining relations" begin
      for i in 1:5, j in 1:5
        if i == j
          @test e(i) * e(j) == C(1)
        else
          @test e(i) * e(j) + e(j) * e(i) == C(G[i, j]) && G[i, j] == C(z^abs(i - j))
        end
      end
      for i in 1:5, j in (i + 1):5
        if i != j
          @test e(i) * e(j) == basis(C, 2^(i - 1) + 2^(j - 1) + 1)
        end
      end
    end
    @testset "some computations" begin
      x = sum(map(i -> e(i), 1:5))
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == typeof(C())

      @test x^2 == C(sum(map(i -> i * z^(4 - i), 1:4))) #Easy by-hand computation
      y = C(QQ.(1:32))
      @test y^2 == C([
        -811 * z^3 - 4277 * z^2 - 4364 * z - 1654,
        -838 * z^3 - 4362 * z^2 - 4520 * z - 1772,
        277 * z^3 - 1975 * z^2 - 4912 * z - 3670,
        288 * z^3 - 2124 * z^2 - 5144 * z - 3808,
        522 * z^3 + 3127 * z^2 + 5138 * z + 2420,
        560 * z^3 + 3170 * z^2 + 5108 * z + 2304,
        -804 * z^3 - 3507 * z^2 - 5158 * z - 3512,
        -832 * z^3 - 3496 * z^2 - 5040 * z - 3392,
        856 * z^3 + 2315 * z^2 + 2102 * z + 968,
        904 * z^3 + 2454 * z^2 + 2332 * z + 1088,
        1468 * z^3 + 4657 * z^2 + 5914 * z + 2384,
        1520 * z^3 + 4804 * z^2 + 6024 * z + 2432,
        492 * z^3 - 2241 * z^2 - 4474 * z - 2864,
        552 * z^3 - 2102 * z^2 - 4316 * z - 2752,
        3508 * z^3 + 6941 * z^2 + 6050 * z + 1760,
        3616 * z^3 + 7024 * z^2 + 6000 * z + 1632,
        -1195 * z^3 - 4821 * z^2 - 5180 * z - 1862,
        -1206 * z^3 - 4858 * z^2 - 5224 * z - 1868,
        949 * z^3 + 153 * z^2 - 2000 * z - 2134,
        976 * z^3 + 36 * z^2 - 2264 * z - 2304,
        1706 * z^3 + 4295 * z^2 + 4562 * z + 852,
        1808 * z^3 + 4466 * z^2 + 4628 * z + 832,
        876 * z^3 - 131 * z^2 - 1958 * z - 2008,
        928 * z^3 + 8 * z^2 - 1744 * z - 1920,
        1320 * z^3 - 645 * z^2 - 2186 * z - 1816,
        1384 * z^3 - 474 * z^2 - 1988 * z - 1728,
        3452 * z^3 + 6721 * z^2 + 5082 * z + 176,
        3600 * z^3 + 6836 * z^2 + 5032 * z + 64,
        284 * z^3 - 1201 * z^2 - 1722 * z - 1424,
        424 * z^3 - 934 * z^2 - 1468 * z - 1344,
        4852 * z^3 + 7405 * z^2 + 4290 * z - 1248,
        4960 * z^3 + 7520 * z^2 + 4304 * z - 1312,
      ]) #According to Magma
      @test divexact(y, z) == inv(z) * y
      @test divexact(y, K(2)) == K(1//2) * y
      @test divexact(y, 2) == QQ(1//2) * y
      @test divexact(y, QQ(2)) == QQ(1//2) * y
    end
    @testset "getindex" begin
      x = C(QQ.(1:32))
      @test typeof(x) == typeof(C())
      for i in 1:32
        @test x[i] == i
        @test coeff(x, i) == i
      end
      @test_throws BoundsError x[0]
      @test_throws BoundsError x[33]
      @test_throws BoundsError x[rand(33:500)]
    end
    @testset "center and centroid" begin
      @test basis_of_center(C) == basis_of_centroid(C)
      orth = C([0, z^2, -z^3, 0, z^2, 0, 0, -2 * z, -z^3,
        0, 0, 2 * z^2, 0, -2 * z^3, -2 * z^3 - 2 * z^2 - 2 * z - 2,
        0, z^2, 0, 0, -2 * z, 0, 2 * z^2, -2 * z^3, 0,
        0, -2 * z, 2 * z^2, 0, -2 * z, 0, 0, 4])
      @test basis_of_centroid(C) == [one(C), orth]
      @test disq(C) == coefficients(orth^2)[1]
      @test disq(C) == quadratic_discriminant(C)
      @test disq(C) == -3 * z^3 - 19 * z^2 - 3 * z + 13
    end
  end

  @testset "multiplication helpers on bitvectors" begin
    
    @testset "test shift_entries!" begin
      shift_entries! = Oscar._shift_entries!
      A3, s3 = QQ.([1, 2, 3, 4, 5, 6]), 0
      A4, s4 = QQ.([1, 2, 3, 4, 5, 0]), 1
      A5, s5 = QQ.([1//4, -5//6, 13//2, 0, 0, 0, 0]), 3
      @test shift_entries!(A3, s3) == A3
      @test shift_entries!(A4, s4) == QQ.([0, 1, 2, 3, 4, 5])
      @test shift_entries!(A5, s5) == QQ.([0, 0, 0, 1//4, -5//6, 13//2, 0])
    end

    @testset "test baseelt_with_gen" begin
      mul_baseelt_with_gen = Oscar._mul_baseelt_with_gen
      @testset "mul_baseelt_with_gen orthogonal" begin
        gram1, gram2 = QQ[4 0; 0 -6], QQ[1 0; 0 1]
        @test mul_baseelt_with_gen(1, 1, gram1) == QQ.([0, 1, 0, 0])
        @test mul_baseelt_with_gen(1, 2, gram1) == QQ.([0, 0, 1, 0])
        @test mul_baseelt_with_gen(2, 1, gram1) == QQ.([2, 0, 0, 0])
        @test mul_baseelt_with_gen(2, 2, gram1) == QQ.([0, 0, 0, 1])
        @test mul_baseelt_with_gen(3, 1, gram1) == QQ.([0, 0, 0, -1])
        @test mul_baseelt_with_gen(3, 2, gram1) == QQ.([-3, 0, 0, 0])
        @test mul_baseelt_with_gen(4, 1, gram1) == QQ.([0, 0, -2, 0])
        @test mul_baseelt_with_gen(4, 2, gram1) == QQ.([0, -3, 0, 0])

        @test mul_baseelt_with_gen(4, 1, gram2) == QQ.([0, 0, -1//2, 0])
        @test mul_baseelt_with_gen(4, 2, gram2) == QQ.([0, 1//2, 0, 0])
      end
      @testset "mul_baseelt_with_gen non-orthogonal" begin
        H2 = QQ[0 1; 1 0]
        @test mul_baseelt_with_gen(2, 1, H2) == QQ.([0, 0, 0, 0])
        @test mul_baseelt_with_gen(3, 2, H2) == QQ.([0, 0, 0, 0])

        A3 = QQ[2 1 0; 1 2 1; 0 1 2]
        @test mul_baseelt_with_gen(6, 1, A3) == QQ.([0, 0, 0, 0, -1, 0, 0, 0])
        @test mul_baseelt_with_gen(6, 2, A3) == QQ.([0, 1, 0, 0, 0, 0, 0, -1])
        @test mul_baseelt_with_gen(6, 3, A3) == QQ.([0, 1, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(7, 1, A3) == QQ.([0, 0, 0, 0, -1, 0, 0, 1])
        @test mul_baseelt_with_gen(7, 2, A3) == QQ.([0, 0, 1, 0, -1, 0, 0, 0])
        @test mul_baseelt_with_gen(7, 3, A3) == QQ.([0, 0, 1, 0, 0, 0, 0, 0])

        X = QQ[2 -4 -1 -1; -4 4 -3 1; -1 -3 -2 3; -1 1 3 2]
        @test mul_baseelt_with_gen(14, 1, X) ==
          QQ.([0, 0, 0, 0, 0, X[1, 4], 0, 0, 0, -X[1, 3], 0, 0, X[1, 1] / 2, 0, 0, 0])
        @test mul_baseelt_with_gen(14, 2, X) ==
          QQ.([0, 0, 0, 0, 0, X[2, 4], 0, 0, 0, -X[2, 3], 0, 0, 0, 0, 0, 1])
        @test mul_baseelt_with_gen(14, 3, X) ==
          QQ.([0, 0, 0, 0, 0, X[3, 4], 0, 0, 0, -X[3, 3] / 2, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(14, 4, X) ==
          QQ.([0, 0, 0, 0, 0, X[4, 4] / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        @test mul_baseelt_with_gen(4, 1, X) ==
          QQ.([0, X[1, 2], -X[1, 1] / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(4, 2, X) ==
          QQ.([0, X[2, 2] / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(4, 3, X) ==
          QQ.([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(4, 4, X) ==
          QQ.([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0])

        K, a = quadratic_field(3)
        G = K[1 a; a 1]
        @test mul_baseelt_with_gen(4, 1, G) == [0, a, (-1//2), 0]
        @test mul_baseelt_with_gen(4, 2, G) == [0, (1//2), 0, 0]
      end
    end
  end
end
