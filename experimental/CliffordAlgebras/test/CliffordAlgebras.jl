
@testset "all tests" begin
  coeffs, dim_qf, set_even_odd_coeffs!, mul_with_gen = Oscar.coeffs,
  Oscar.dim_qf, Oscar.set_even_odd_coeffs!,
  Oscar._mul_with_gen

  @testset "some corner cases" begin
    ma = matrix_algebra(QQ, 2)
    Gfail = matrix([basis(ma)[1] basis(ma)[2]; basis(ma)[2] basis(ma)[4]])
    @test_throws ArgumentError clifford_algebra(Gfail)
    G = identity_matrix(ZZ, 0) #Empty matrix
    C = clifford_algebra(G)
    zer = [ZZ(0)]
    @test is_commutative(C)
    @testset "construction" begin
      @test typeof(C) == AlgClf{typeof(base_ring(C)())}
      @test elem_type(C) == AlgClfElem{typeof(base_ring(C)()),typeof(C)}
      @test elem_type(C) == typeof(C())

      @test (base_ring(C), gram(C), dim_qf(C), dim(C)) == (ZZ, G, 0, 1)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coeffs(C()) == zer
      @test even_coeffs(C()) == zer && even_coeffs(C()) == odd_coeffs(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coeffs(C(1)) == [ZZ(1)]
      @test even_coeffs(C(1)) == [ZZ(1)] && odd_coeffs(C(1)) == zer
    end
    @testset "functions on elements" begin
      x = C([17])
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == typeof(C())

      @test even_part(x) == C([17])
      @test odd_part(x) == C([0])
      xeven, xodd = even_coeffs(x), odd_coeffs(x)
      x.even_coeffs, x.odd_coeffs = zer, [ZZ(1)]
      @test xeven != even_coeffs(x)
      @test xodd != odd_coeffs(x)
      set_even_odd_coeffs!(x)
      @test xeven == even_coeffs(x) && xodd == odd_coeffs(x)

      @test +x == x
      @test -x == C([-17])
      @test_throws ArgumentError divexact(x, 2)
      @test divexact(2 * x, 2) == x && divexact(2 * x, ZZ(2)) == x
      @test divexact(C(18), 2) == C(9)
    end
    @testset "defining relations" begin
      @test length(basis(C)) == 1
      @test basis(C, 1) == C(1)
      @test is_empty(gens(C))
      a, b = rand(ZZ, -100:100), rand(ZZ, -100:100)
      @test C(a) * C(b) == C(a * b)
      @test C(a) + C(b) == C(a + b)
      @test (C(a) + C(b)) * (C(a) - C(b)) == C(a)^2 - C(b)^2
    end
  end

  @testset "ring_of_integers_sqrt(5)" begin
    K, _ = quadratic_field(5)
    OK = ring_of_integers(K)
    a = Oscar.basis(OK)[2] #1//2 * (1 + sqrt(5))
    G = OK[2*a 1; 1 2*(1 - a)]
    C = clifford_algebra(G)
    @test !is_commutative(C)
    e(i::Int) = gen(C, i) #e(i) returns the i-th generator of C
    zerOK = OK.([0, 0, 0, 0])
    @testset "construction" begin
      @test typeof(C) == AlgClf{typeof(base_ring(C)())}
      @test elem_type(C) == AlgClfElem{typeof(base_ring(C)()),typeof(C)}
      @test elem_type(C) == typeof(C())

      Gfail1, Gfail2 = OK[2*a 5; 1 2*(1 - a)], OK[a 1; 1 (1-a)]
      @test_throws ArgumentError clifford_algebra(Gfail1)
      @test_throws ArgumentError clifford_algebra(Gfail2)

      @test (base_ring(C), gram(C), dim_qf(C), dim(C)) == (OK, G, 2, 4)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coeffs(C()) == zerOK
      @test even_coeffs(C()) == zerOK && even_coeffs(C()) == odd_coeffs(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coeffs(C(1)) == OK.([1, 0, 0, 0])
      @test even_coeffs(C(1)) == OK.([1, 0, 0, 0]) && odd_coeffs(C(1)) == zerOK
      verr1, verr2 = OK.([1, 2, 3]), OK.([1, 2, 3, 4, 5])
      @test_throws ArgumentError C(verr1)
      @test_throws ArgumentError C(verr2)
      Ca = C(a)
      @test parent(Ca) == C
      @test Ca == C([a, 0, 0, 0])
      @test Ca == C(OK.([a, 0, 0, 0]))
      @test 2 * Ca == C([2 * a, 0, 0, 0])
      @test Ca * 2 == 2 * Ca
      @test 2 * Ca == C(2 * a)
      @test C(2) == C(OK(2))
    end
    @testset "functions on elements" begin
      x = C([-1, 1, a, a + 1])
      @test parent(x) == C
      @test parent_type(x) == typeof(C)
      @test typeof(x) == typeof(C())

      @test even_part(x) == C([-1, 0, 0, a + 1])
      @test odd_part(x) == C([0, 1, a, 0])
      xeven, xodd = even_coeffs(x), odd_coeffs(x)
      x.even_coeffs, x.odd_coeffs = zerOK, zerOK
      @test xeven != even_coeffs(x) && xodd != odd_coeffs(x)
      set_even_odd_coeffs!(x)
      @test xeven == even_coeffs(x) && xodd == odd_coeffs(x)

      @test +x == x
      @test -x == C([1, -1, -a, -(a + 1)])
      @test_throws ArgumentError divexact(x, 2)
      @test divexact(2 * x, 2) == x && divexact(2 * x, OK(2)) == x
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
      @test e(2) * e(1) == C(gram(C)[1, 2]) - e(1) * e(2)
    end
    @testset "some computations" begin
      x, y = C([-2, 1, -1, 2]), C([2, 4, 6, 8])
      @test is_zero(x - x) && is_zero(x - +x)
      @test even_part(x) == C([-2, 0, 0, 2]) && odd_part(x) == C([0, 1, -1, 0])
      @test x + y == C([0, 5, 5, 10])
      @test x - y == C([-4, -3, -7, -6])
      @test y - x == C([4, 3, 7, 6])

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
      @test mul_with_gen(coeffs(x), 1, gram(C)) == OK.([3 * a + 4, -(a^2 + 1), a, -(a + 4)])
      @test mul_with_gen(coeffs(x), 2, gram(C)) ==
        OK.([-(a^2 + 3 * a - 4), -(1 - a), -a^2, 2])
    end
  end

  @testset "random quaternion alg over rationals" begin
    a, b = rand(QQ, -5:5), rand(QQ, -5:5)
    G = diagonal_matrix(2 * a, 2 * b)
    C = clifford_algebra(G)
    @test !is_commutative(C)
    e(i::Int) = gen(C, i) #e(i) returns the i-th generator of C
    zer = QQ.([0, 0, 0, 0])
    @testset "construction" begin
      @test typeof(C) == AlgClf{typeof(base_ring(C)())}
      @test elem_type(C) == AlgClfElem{typeof(base_ring(C)()),typeof(C)}
      @test elem_type(C) == typeof(C())

      @test (base_ring(C), gram(C), dim_qf(C), dim(C)) == (QQ, G, 2, 4)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coeffs(C()) == zer
      @test even_coeffs(C()) == zer && even_coeffs(C()) == odd_coeffs(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coeffs(C(1)) == QQ.([1, 0, 0, 0])
      @test even_coeffs(C(1)) == QQ.([1, 0, 0, 0]) && odd_coeffs(C(1)) == zer
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
      @test e(2) * e(1) == C(gram(C)[1, 2]) - e(1) * e(2)
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
      @test -x == C(-coeffs(x))
      @test is_zero(x + -x)
    end
    @testset "mul_with_gen" begin
      x = C([1, 2, 3, 4])
      @test mul_with_gen(coeffs(x), 1, gram(C)) == QQ.([2 * a, 1, -4 * a, -3])
      @test mul_with_gen(coeffs(x), 2, gram(C)) == QQ.([3 * b, 4 * b, 1, 2])
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
    C = clifford_algebra(G)
    @test !is_commutative(C)
    e(i::Int) = gen(C, i) #e(i) returns the i-th generator of C
    zer = fill(K(), 2^5)
    ein = fill(K(), 2^5)
    ein[1] = K(1)
    @testset "construction" begin
      @test typeof(C) == AlgClf{typeof(base_ring(C)())}
      @test elem_type(C) == AlgClfElem{typeof(base_ring(C)()),typeof(C)}
      @test elem_type(C) == typeof(C())

      @test (base_ring(C), gram(C), dim_qf(C), dim(C)) == (K, G, 5, 32)
      @test C() == C(0) && C() == zero(C)
      @test is_zero(C())
      @test coeffs(C()) == zer
      @test even_coeffs(C()) == zer && even_coeffs(C()) == odd_coeffs(C())
      @test C(1) == one(C)
      @test is_one(C(1))
      @test coeffs(C(1)) == ein
      @test even_coeffs(C(1)) == ein && odd_coeffs(C(1)) == zer
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
      y = C(ZZ.(1:32))
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
      ]) #Magma
    end
  end

  @testset "multiplication helpers on bitvectors" begin
    bitvec_to_int = Oscar._bitvec_to_int
    @testset "test bitvec_to_int" begin
      @test bitvec_to_int(BitVector([0])) == 0
      @test bitvec_to_int(BitVector([0, 0, 0, 0])) == 0
      @test bitvec_to_int(BitVector([1, 0, 0, 0])) == 1
      @test bitvec_to_int(BitVector([0, 1, 0, 0])) == 2
      @test bitvec_to_int(BitVector([0, 0, 1, 0])) == 4
      @test bitvec_to_int(BitVector([0, 0, 0, 1])) == 8
      @test bitvec_to_int(BitVector([1, 1, 1, 1])) == 15
      @test bitvec_to_int(BitVector([0, 1, 0, 1])) == 10
      @test bitvec_to_int(BitVector([1, 1, 1, 0, 1, 0, 0, 1])) == 1 + 2 + 4 + 16 + 128
    end
    @testset "test getbasis_index" begin
      getbasis_index = Oscar._get_basis_index
      b = BitVector([1, 0, 0, 1, 1, 0, 1])
      @test getbasis_index(b) == bitvec_to_int(b) + 1
    end
    @testset "test getfirst and getlast" begin
      #getfirst = Oscar._get_first
      getlast = Oscar._get_last
      #@test getfirst(BitVector([0, 0, 0, 0, 0])) == 6
      @test getlast(BitVector([0, 0, 0, 0, 0])) == 0

      #@test getfirst(BitVector([0, 1, 1, 0, 1, 0])) == 2
      @test getlast(BitVector([0, 1, 1, 0, 1, 0])) == 5

      #@test getfirst(BitVector([1])) == 1
      @test getlast(BitVector([1])) == 1
    end
    @testset "test shift_entries!" begin
      shift_entries! = Oscar._shift_entries!
      A3, s3 = ZZ.([1, 2, 3, 4, 5, 6]), 0
      A4, s4 = ZZ.([1, 2, 3, 4, 5, 0]), 1
      A5, s5 = QQ.([1//4, -5//6, 13//2, 0, 0, 0, 0]), 3
      @test shift_entries!(A3, s3) == A3
      @test shift_entries!(A4, s4) == ZZ.([0, 1, 2, 3, 4, 5])
      @test shift_entries!(A5, s5) == QQ.([0, 0, 0, 1//4, -5//6, 13//2, 0])
    end

    @testset "test baseelt_with_gen" begin
      mul_baseelt_with_gen = Oscar._mul_baseelt_with_gen
      @testset "mul_baseelt_with_gen orthogonal" begin
        gram1, gram2 = ZZ[4 0; 0 -6], QQ[1 0; 0 1]
        @test mul_baseelt_with_gen(BitVector([0, 0]), 1, gram1) == ZZ.([0, 1, 0, 0])
        @test mul_baseelt_with_gen(BitVector([0, 0]), 2, gram1) == ZZ.([0, 0, 1, 0])
        @test mul_baseelt_with_gen(BitVector([1, 0]), 1, gram1) == ZZ.([2, 0, 0, 0])
        @test mul_baseelt_with_gen(BitVector([1, 0]), 2, gram1) == ZZ.([0, 0, 0, 1])
        @test mul_baseelt_with_gen(BitVector([0, 1]), 1, gram1) == ZZ.([0, 0, 0, -1])
        @test mul_baseelt_with_gen(BitVector([0, 1]), 2, gram1) == ZZ.([-3, 0, 0, 0])
        @test mul_baseelt_with_gen(BitVector([1, 1]), 1, gram1) == ZZ.([0, 0, -2, 0])
        @test mul_baseelt_with_gen(BitVector([1, 1]), 2, gram1) == ZZ.([0, -3, 0, 0])

        @test mul_baseelt_with_gen(BitVector([1, 1]), 1, gram2) == QQ.([0, 0, -1//2, 0])
        @test mul_baseelt_with_gen(BitVector([1, 1]), 2, gram2) == QQ.([0, 1//2, 0, 0])
      end
      @testset "mul_baseelt_with_gen non-orthogonal" begin
        H2 = ZZ[0 1; 1 0]
        @test mul_baseelt_with_gen(BitVector([1, 0]), 1, H2) == ZZ.([0, 0, 0, 0])
        @test mul_baseelt_with_gen(BitVector([0, 1]), 2, H2) == ZZ.([0, 0, 0, 0])

        A3 = QQ[2 1 0; 1 2 1; 0 1 2]
        b101 = BitVector([1, 0, 1])
        @test mul_baseelt_with_gen(b101, 1, A3) == QQ.([0, 0, 0, 0, -1, 0, 0, 0])
        @test mul_baseelt_with_gen(b101, 2, A3) == QQ.([0, 1, 0, 0, 0, 0, 0, -1])
        @test mul_baseelt_with_gen(b101, 3, A3) == QQ.([0, 1, 0, 0, 0, 0, 0, 0])
        b011 = BitVector([0, 1, 1])
        @test mul_baseelt_with_gen(b011, 1, A3) == QQ.([0, 0, 0, 0, -1, 0, 0, 1])
        @test mul_baseelt_with_gen(b011, 2, A3) == QQ.([0, 0, 1, 0, -1, 0, 0, 0])
        @test mul_baseelt_with_gen(b011, 3, A3) == QQ.([0, 0, 1, 0, 0, 0, 0, 0])

        X = ZZ[2 -4 -1 -1; -4 4 -3 1; -1 -3 -2 3; -1 1 3 2]
        b1011 = BitVector([1, 0, 1, 1])
        @test mul_baseelt_with_gen(b1011, 1, X) ==
          ZZ.([0, 0, 0, 0, 0, X[1, 4], 0, 0, 0, -X[1, 3], 0, 0, X[1, 1] / 2, 0, 0, 0])
        @test mul_baseelt_with_gen(b1011, 2, X) ==
          ZZ.([0, 0, 0, 0, 0, X[2, 4], 0, 0, 0, -X[2, 3], 0, 0, 0, 0, 0, 1])
        @test mul_baseelt_with_gen(b1011, 3, X) ==
          ZZ.([0, 0, 0, 0, 0, X[3, 4], 0, 0, 0, -X[3, 3] / 2, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(b1011, 4, X) ==
          ZZ.([0, 0, 0, 0, 0, X[4, 4] / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        b1100 = BitVector([1, 1, 0, 0])
        @test mul_baseelt_with_gen(b1100, 1, X) ==
          ZZ.([0, X[1, 2], -X[1, 1] / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(b1100, 2, X) ==
          ZZ.([0, X[2, 2] / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(b1100, 3, X) ==
          ZZ.([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
        @test mul_baseelt_with_gen(b1100, 4, X) ==
          ZZ.([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0])

        K, a = quadratic_field(3)
        G = K[1 a; a 1]
        @test mul_baseelt_with_gen(BitVector([1, 1]), 1, G) == [0, a, (-1//2), 0]
        @test mul_baseelt_with_gen(BitVector([1, 1]), 2, G) == [0, (1//2), 0, 0]
      end
    end
  end
end
