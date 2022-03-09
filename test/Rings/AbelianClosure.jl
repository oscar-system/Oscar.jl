@testset "AbelianClousre" begin

  @testset "Creation" begin
    K, z = abelian_closure(QQ)
    @inferred abelian_closure(QQ)
    @test K === abelian_closure(QQ)[1]
    @test K isa QabField
    @test elem_type(K) === QabElem{nf_elem}
    @test elem_type(typeof(K)) === QabElem{nf_elem}
    @test parent_type(QabElem) === QabField{AnticNumberField}
    @test parent_type(one(K)) === QabField{AnticNumberField}

    a = @inferred K()
    @test a isa QabElem
    @test parent(a) === K

    a = @inferred K(1)
    @test parent(a) === K
    @test a isa QabElem
    @test isone(a)
    @test isone(one(a))
    @test !iszero(a)

    a = @inferred K(0)
    @test parent(a) === K
    @test a isa QabElem
    @test iszero(a)
    @test iszero(zero(a))
    @test !isone(a)

    c = one(K) + one(K)
    for T in [Int, BigInt, fmpz, fmpq, Rational{Int}, Rational{BigInt}]
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
    @test sprint(show, "text/plain", K) == "Abelian closure of Q"
    @test get_variable(K) == "ζ"
    s = sprint(show, "text/plain", z)

    a = z(1)
    sprint(show, "text/plain", a) == "1"
    a = z(4)
    sprint(show, "text/plain", a) == "ζ(4)"

    t = set_variable!(K, "z")
    @test t == "ζ"
    @test get_variable(K) == "z"
    sprint(show, "text/plain", a) == "z(4)"

    zz = gen(K, "ω")
    @test get_variable(K) == "ω"
    sprint(show, "text/plain", zz(4)) == "ω(4)"

    a = @inferred root_of_unity(K, 4)
    @test isone(a^4) && !isone(a) && !isone(a^2)
  end

  @testset "Coercion" begin
    K, z = abelian_closure(QQ)
    @test Oscar.AbelianClosure.isconductor(3)
    @test !Oscar.AbelianClosure.isconductor(2)
    a, b = z(4), z(3)
    c, d = @inferred Oscar.AbelianClosure.make_compatible(a, b)
    fa = minpoly(Oscar.AbelianClosure.data(a))
    fb = minpoly(Oscar.AbelianClosure.data(b))
    fc = minpoly(Oscar.AbelianClosure.data(c))
    fd = minpoly(Oscar.AbelianClosure.data(d))
    @test iszero(fa(c)) && iszero(fc(a))
    @test iszero(fb(d)) && iszero(fd(b))
    @test_throws Hecke.NotImplemented Oscar.AbelianClosure.coerce_down(Hecke.rationals_as_number_field()[1], 1, z(2))
  end

  @testset "Promote rule" begin
    @test Oscar.AbstractAlgebra.promote_rule(QabElem, Int) == QabElem
    @test Oscar.AbstractAlgebra.promote_rule(QabElem, fmpz) == QabElem
    @test Oscar.AbstractAlgebra.promote_rule(QabElem, fmpq) == QabElem
  end

  @testset "Arithmetic" begin
    K, z = abelian_closure(QQ)
    @test isunit(z(1))
    @test !isunit(zero(K))
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
      @test a^fmpz(10) == reduce(*, fill(a, 10))
    end
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
      @test addeq!(a, b) == aa + bb
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
      for T in [Int, BigInt, fmpz, fmpq, Rational{Int}, Rational{BigInt}]
        b = (!iszero(_b) && (T === fmpq || T <: Rational)) ? inv(T(_b)) : T(_b)
        @test a * b == a * K(b)
        @test b * a == a * K(b)
        @test a + b == a + K(b)
        @test b + a == a + K(b)
        @test a - b == a - K(b)
        @test b - a == K(b) - a
        if !iszero(b)
          @test //(a * b, b) == a
        end
      end
    end
  end

  @testset "Comparison" begin
    K, z = abelian_closure(QQ)
    a = z(2)
    b = z(4)^2
    @test @inferred a == b

    for T in [Int, BigInt, fmpz, fmpq, Rational{Int}, Rational{BigInt}]
      @test @inferred b == T(-1)
      @test @inferred T(-1) == b
    end
  end

  @testset "Roots" begin
    K, z = abelian_closure(QQ)
    b = z(5)
    @test root(b, 2)^2 == b
    for n in [3,4,5]
      bs = @inferred roots(b, n)
      @test length(bs) == n
      @test all(c -> c^n == b, bs)

      @test isroot_of_unity(z(5))
      @test !isroot_of_unity(z(5) + 1)
    end

    @test length(roots(8*b, 3)) == 3
    Kx, x = PolynomialRing(K)
    @test length(roots(x^15-2^15)) == 15

    @test order(z(5)) == 5
  end

  @testset "Automorphism" begin
    K, z = abelian_closure(QQ)
    a = z(5)
    f = hom(K, K, 3)
    @test f(a) == a^3
    @test f(z(9)) == z(9)
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

    @test Oscar.AbelianClosure.quadratic_irrationality_info(z(5)) === nothing
  end

  @testset "Polynomial" begin
    K, z = abelian_closure(QQ)
    Kx, x = K["x"]
    @test (x^2 + 1)(z(4)) == z(4)^2 + 1
  end

  @testset "Matrix" begin
    K, z = abelian_closure(QQ)
    M = matrix(K, 2, 2, [z(3), z(4), z(5), z(3)])
    @test M^2 == M * M
  end

  @testset "Singular ring" begin
    K, z = abelian_closure(QQ)
    L = Oscar.singular_ring(K)
    a = z(4)
    @test K(L(a)) == a
  end
end
