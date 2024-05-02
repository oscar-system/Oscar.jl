@testset "AbelianClosure" begin
  @testset "Interface" begin
    K, z = abelian_closure(QQ)
    ConformanceTests.test_Field_interface(K)
  end

  @testset "Creation" begin
    K, z = abelian_closure(QQ)
    @inferred abelian_closure(QQ)
    @test K === abelian_closure(QQ)[1]
    @test K isa QQAbField
    @test elem_type(K) === QQAbFieldElem{AbsSimpleNumFieldElem}
    @test elem_type(typeof(K)) === QQAbFieldElem{AbsSimpleNumFieldElem}
    @test parent_type(QQAbFieldElem{AbsSimpleNumFieldElem}) === QQAbField{AbsSimpleNumField}
    @test parent_type(one(K)) === QQAbField{AbsSimpleNumField}

    a = @inferred K()
    @test a isa QQAbFieldElem
    @test parent(a) === K

    a = @inferred K(1)
    @test parent(a) === K
    @test a isa QQAbFieldElem
    @test isone(a)
    @test isone(one(a))
    @test !iszero(a)

    a = @inferred K(0)
    @test parent(a) === K
    @test a isa QQAbFieldElem
    @test iszero(a)
    @test iszero(zero(a))
    @test !isone(a)

    c = one(K) + one(K)
    for T in [Int, BigInt, ZZRingElem, QQFieldElem, Rational{Int}, Rational{BigInt}]
      b = @inferred K(T(2))
      @test b == c
    end

    a = @inferred z(4)
    @test isone(a^4) && !isone(a) && !isone(a^2) && !isone(a^3)
    a = @inferred root_of_unity(K, 4)
    @test isone(a^4) && !isone(a) && !isone(a^2) && !isone(a^3)

    @test a === K(a)

    @test copy(a) == a
    b = deepcopy(a)
    @test a !== b
    @test a == b
  end

  @testset "Printing" begin
    K, z = abelian_closure(QQ)

    @test AbstractAlgebra.PrettyPrinting.detailed(K) == "Abelian closure of rational field"
    @test AbstractAlgebra.PrettyPrinting.oneline(K) == "Abelian closure of rational field"
    @test AbstractAlgebra.PrettyPrinting.supercompact(K) == "Abelian closure of QQ"

    @test AbstractAlgebra.PrettyPrinting.detailed(z) == "Generator of abelian closure of rational field"
    @test AbstractAlgebra.PrettyPrinting.oneline(z) == "Generator of abelian closure of rational field"
    @test AbstractAlgebra.PrettyPrinting.supercompact(z) == "Generator of abelian closure of QQ"

    orig = get_variable(K)
    @test orig == "zeta"
    s = sprint(show, "text/plain", z)

    a = z(1)
    sprint(show, "text/plain", a) == "1"
    a = z(4)
    sprint(show, "text/plain", a) == "z(4)"
    Oscar.with_unicode() do
      # The following holds only if `set_variable!`
      # has not been called before for `K`.
      @test get_variable(K) == "ζ"
      s = sprint(show, "text/plain", z)

      a = z(1)
      sprint(show, "text/plain", a) == "1"
      a = z(4)
      sprint(show, "text/plain", a) == "ζ(4)"
    end

    t = set_variable!(K, "zz")
    @test t == "zeta"
    @test get_variable(K) == "zz"
    sprint(show, "text/plain", a) == "zz(4)"

    zz = gen(K, "ω")
    @test get_variable(K) == "ω"
    sprint(show, "text/plain", zz(4)) == "ω(4)"

    a = @inferred root_of_unity(K, 4)
    @test isone(a^4) && !isone(a) && !isone(a^2)

    # reset variable for any subsequent (doc-)tests
    @test set_variable!(K, orig) == "ω"
  end

  @testset "Coercion" begin
    K, z = abelian_closure(QQ)
    @test Oscar.AbelianClosure.is_conductor(3)
    @test !Oscar.AbelianClosure.is_conductor(2)
    a, b = z(4), z(3)
    c, d = @inferred Oscar.AbelianClosure.make_compatible(a, b)
    fa = minpoly(Oscar.AbelianClosure.data(a))
    fb = minpoly(Oscar.AbelianClosure.data(b))
    fc = minpoly(Oscar.AbelianClosure.data(c))
    fd = minpoly(Oscar.AbelianClosure.data(d))
    @test fa == minpoly(a)
    @test iszero(fa(c)) && iszero(fc(a))
    @test iszero(fb(d)) && iszero(fd(b))
    @test_throws Hecke.NotImplemented Oscar.AbelianClosure.coerce_down(Hecke.rationals_as_number_field()[1], 1, z(2))
  end

  @testset "Conversion to ZZRingElem and QQFieldElem" begin
    K, z = abelian_closure(QQ)
    x = z(5)
    y = ZZ(x^5)
    @test y isa ZZRingElem
    @test y == 1
    y = QQ(x^5)
    @test y isa QQFieldElem
    @test y == 1
    @test_throws ErrorException ZZ(x)
    @test_throws ErrorException QQ(x)
  end

  @testset "Conversion to QQBarFieldElem" begin
    K, z = abelian_closure(QQ)
    F = QQBarField()
    x = z(5)
    y = F(x)
    @test F(x) == QQBarFieldElem(x)
    @test F(x^2) == y^2
    @test y == root_of_unity(F, 5)
    @test is_one(y^5)
    a = x + x^4
    b = x^2 + x^3
    n = zero(K)
    @test a > n && F(a) > F(n)
    @test a > 0
    @test a > BigInt(0)
    @test a > 0 // 1
    @test a > ZZ(0)
    @test a > QQ(0)
    @test b < n && F(b) < F(n)
    @test b < 0
    @test b < BigInt(0)
    @test b < 0 // 1
    @test b < ZZ(0)
    @test b < QQ(0)
    @test b < a && F(b) < F(a)
    @test is_positive(a)
    @test is_negative(b)
    @test_throws DomainError x < n
    @test_throws DomainError y < F(n)

    # compatibility of roots of unity in `K` and `F`
    i = F(z(4))
    for n in 1:20
      @test z(n) == root_of_unity(K, n)
      @test all(d -> z(n)^d == z(div(n, d)), divisors(n))
      @test F(z(n)) == cospi(F(2)/n) + i*sinpi(F(2)/n)
    end
  end

  @testset "Conversion to Float64 and ComplexF64" begin
    K, z = abelian_closure(QQ)
    x = z(5)
    a = x + x^4
    @test Float64(a) > 0
    @test_throws InexactError Float64(x)
    @test is_one(ComplexF64(one(K)))
  end

  @testset "Promote rule" begin
    @test Oscar.AbstractAlgebra.promote_rule(QQAbFieldElem, Int) == QQAbFieldElem
    @test Oscar.AbstractAlgebra.promote_rule(QQAbFieldElem, ZZRingElem) == QQAbFieldElem
    @test Oscar.AbstractAlgebra.promote_rule(QQAbFieldElem, QQFieldElem) == QQAbFieldElem
  end

  @testset "Test for integrality and rationality" begin
    for sparse in (true, false)
      K, z = abelian_closure(QQ; sparse)
      z_3 = z(3)
      @test isinteger(z_3^3)
      @test !isinteger(z_3)
      @test is_rational(1//2*z_3^3)
      @test !is_rational(1//2*z_3)
      @test is_algebraic_integer(z_3)
      @test !is_algebraic_integer(1//2*z_3)
    end
  end

  @testset "Arithmetic" begin
    K, z = abelian_closure(QQ)
    @test is_unit(z(1))
    @test !is_unit(zero(K))
    rand_elem() = begin n = rand([3, 4, 5]); sum(rand(-1:1) * z(n) for i in 1:3) end
    a = one(K)
    @test isone(inv(a))
    for i in 1:10
      a = rand_elem()
      b = rand_elem()
      c = rand_elem()
      @test iszero(a - a)
      @test iszero(a + (-a))
      @test a * (b + c) == a*b + a*c
      @test (b + c) * a == b*a + c*a
      @test a * (b - c) == a*b - a*c
      @test (b - c) * a == b*a - c*a
      if !iszero(a)
        @test isone(a * inv(a))
        @test divexact(a * b, a) == b
        @test //(a * b, a) == b
        @test div(a * b, a) == b
      end
      @test a^10 == reduce(*, fill(a, 10))
      @test a^ZZRingElem(10) == reduce(*, fill(a, 10))
    end
  end

  @testset "Reduction to finite fields" for gf_fun in [GF, Nemo.Native.GF]
    K, z = abelian_closure(QQ)
    F = gf_fun(2, 1)
    @test reduce(z(2), F) == one(F)
    @test reduce(one(K), F) == one(F)
    @test reduce(zero(K), F) == zero(F)
    @test_throws ErrorException reduce(z(3), F)
    @test_throws ErrorException reduce(z(4), F)
    F = gf_fun(2, 2)
    red = reduce(z(3), F)
    @test red != one(F) && red^3 == one(F)
    @test reduce(one(K), F) == one(F)
    @test reduce(zero(K), F) == zero(F)
    @test_throws ErrorException reduce(z(7), F)
    F = gf_fun(2, 3)
    red = reduce(z(7), F)
    @test red != one(F) && red^7 == one(F)
    @test reduce(one(K), F) == one(F)
    @test reduce(zero(K), F) == zero(F)
    F = gf_fun(3, 2)
    @test reduce(one(K), F) == one(F)
    @test reduce(zero(K), F) == zero(F)
    @test_throws ErrorException reduce(z(3), F)
  end

  @testset "Unsafe operations" begin
    K, z = abelian_closure(QQ)
    rand_elem() = begin n = rand([3, 4, 5]); sum(rand(-1:1) * z(n) for i in 1:3) end
    for i in 1:10
      a = rand_elem()
      b = rand_elem()
      c = rand_elem()
      aa = deepcopy(a)
      bb = deepcopy(b)
      @test add!(a, b) == aa + bb
      @test b == bb
      aa = deepcopy(a)
      @test Oscar.AbelianClosure.neg!(a) == -aa

      a = rand_elem()
      b = rand_elem()
      c = rand_elem()
      aa = deepcopy(a)
      bb = deepcopy(b)
      @test mul!(c, a, b) == a * b
    end
  end

  @testset "Ad hoc operations" begin
    K, z = abelian_closure(QQ)
    rand_elem() = begin n = rand([3, 4, 5]); sum(rand(-1:1) * z(n) for i in 1:3) end
    for i in 1:10
      a = rand_elem()
      _b = rand(-10:10)
      for T in [Int, BigInt, ZZRingElem, QQFieldElem, Rational{Int}, Rational{BigInt}]
        b = (!iszero(_b) && (T === QQFieldElem || T <: Rational)) ? inv(T(_b)) : T(_b)
        @test a * b == a * K(b)
        @test b * a == a * K(b)
        @test a + b == a + K(b)
        @test b + a == a + K(b)
        @test a - b == a - K(b)
        @test b - a == K(b) - a
        if !iszero(b)
          @test (a * b) // b == a
          @test (a // b) * b == a

          @test (a * b) / b == a
          @test (a / b) * b == a
        end

        if !iszero(a)
          @test (b * a) // a == b
          @test (b // a) * a == b

          @test (b * a) / a == b
          @test (b / a) * a == b
        end
      end
    end
  end

  @testset "Comparison" begin
    K, z = abelian_closure(QQ)
    a = z(2)
    b = z(4)^2
    @test @inferred a == b

    for T in [Int, BigInt, ZZRingElem, QQFieldElem, Rational{Int}, Rational{BigInt}]
      @test @inferred b == T(-1)
      @test @inferred T(-1) == b
    end

    a = z(5)
    b = z(15)^3
    @test @inferred a == b
    @test hash(a, zero(UInt)) == hash(b, zero(UInt))
  end

  @testset "Roots" begin
    K, z = abelian_closure(QQ)
    b = z(5)
    @test root(b, 2)^2 == b
    for n in [3,4,5]
      bs = @inferred roots(b, n)
      @test length(bs) == n
      @test all(c -> c^n == b, bs)

      @test is_root_of_unity(z(5))
      @test !is_root_of_unity(z(5) + 1)
    end

    @test length(roots(8*b, 3)) == 3
    Kx, x = polynomial_ring(K)
    @test length(roots(x^15-2^15)) == 15

    @test order(z(5)) == 5

    @test cbrt(K(8))^3 == K(8)
    @test_throws ErrorException("Element 4 does not have a cube root") cbrt(K(4))

  end

  @testset "Automorphism" begin
    K, z = abelian_closure(QQ)
    a = z(5)
    f = hom(K, K, 3)
    @test f(a) == a^3
    @test f(z(9)) == z(9)

    @test isone(conj(a)*a)
    @test !isreal(a)
    @test isreal(a+conj(a))
  end

  @testset "Square roots" begin
    K, z = abelian_closure(QQ)

    sqrts = 
    [ 4*z(17)^14 + 4*z(17)^12 + 4*z(17)^11 + 4*z(17)^10 + 4*z(17)^7 + 4*z(17)^6 + 4*z(17)^5 + 4*z(17)^3 + 3, 
      20*z(13)^11 + 20*z(13)^8 + 20*z(13)^7 + 20*z(13)^6 + 20*z(13)^5 + 20*z(13)^2 + 23,
      -2*z(8)^3 - 2*z(8),
      -2*z(11)^9 - 2*z(11)^5 - 2*z(11)^4 - 2*z(11)^3 - 2*z(11) - 1,
      -4*z(12)^2 + 2,
      2*z(12)^3 - 4*z(12),
      -2*z(15)^7 + 2*z(15)^5 - 4*z(15)^4 + 2*z(15)^3 - 2*z(15)^2 - 4*z(15) + 3,
     ]

    for a in sqrts
      x, y, n = Oscar.AbelianClosure.quadratic_irrationality_info(a)
      @test ((a - x)//y)^2 == n
    end

    @test sqrt(K(2))^2 == K(2)
    @test_throws ErrorException("Element zeta(4) + 1 does not have a square root") sqrt(z(4)+1)

    @test Oscar.AbelianClosure.quadratic_irrationality_info(z(5)) === nothing

    @testset for d in -50:50
      if d == 0
        N = 1
      elseif mod(d, 4) == 1
        N = abs(d)
      else
        N = 4*abs(d)
      end
      x = Oscar.AbelianClosure.square_root_in_cyclotomic_field(K, d, N)
      @test x^2 == d
    end

    @test Oscar.AbelianClosure.square_root_in_cyclotomic_field(K, 7, 7) === nothing
    @test Oscar.AbelianClosure.square_root_in_cyclotomic_field(K, 7, 5) === nothing
  end

  @testset "Natural embeddings" begin
    K, _ = abelian_closure(QQ)
    @testset "Natural embeddings of cyclotomic fields" for N in 1:50
      F, z = cyclotomic_field(N)
      x = K(z)
      @test order(x) == N
    end

    @testset "Natural embeddings of quadratic fields" for d in -50:50
      if ! is_square(d)
        F, z = quadratic_field(d)
        x = K(z)
        @test x^2 == d
      end
    end

    R, x = polynomial_ring(QQ)
    F, z = number_field(x^3-2)
    @test_throws ArgumentError K(z)
  end

  @testset "Polynomial" begin
    K, z = abelian_closure(QQ)
    Kx, x = K[:x]
    @test (x^2 + 1)(z(4)) == z(4)^2 + 1
  end

  @testset "Matrix" begin
    K, z = abelian_closure(QQ)
    M = matrix(K, 2, 2, [z(3), z(4), z(5), z(3)])
    @test M^2 == M * M
  end

  @testset "Singular ring" begin
    K, z = abelian_closure(QQ)
    iso = Oscar.iso_oscar_singular_coeff_ring(K)
    L = codomain(iso)
    a = z(4)
    @test preimage(iso, iso(a)) == a
  end

  K, z = abelian_closure(QQ)
  S = [z(3)]
  @test degree(number_field(QQ, S)[1]) == 2

  @testset "Syntax to create number fields from QabElem" begin
    K, z = abelian_closure(QQ)
    F, a = QQ[z(5)]
    @test F isa AbsSimpleNumField
    @test dim(F) == 4
    F, a = QQ[[z(5), z(3)]]
    @test F isa AbsSimpleNumField
    @test dim(F) == 8
    F, a = QQ[z(5), z(3)]
    @test F isa AbsSimpleNumField
    @test dim(F) == 8
  end

  @testset "Atlas irrationalities" begin
    K, z = abelian_closure(QQ)
    vals = [z(7) + z(7)^2 + z(7)^4,
            z(8) + z(8)^3,
            z(9) + z(9)^-1,
            z(3) - z(3)^2]
    for val in vals
      @test atlas_irrationality(atlas_description(val)) == val
    end
  end
end
