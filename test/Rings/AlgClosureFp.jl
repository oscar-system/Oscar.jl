function test_elem(K::AlgClosure{T}) where T <: FinField
  d = rand(1:8)
  F = ext_of_degree(K, d)
  return K(rand(F))
end

@testset "AlgClosureFp" begin
  @testset "Interface" begin
    K = algebraic_closure(GF(3,1))
    test_Field_interface(K)
  end

  @testset "Creation" begin
    F = GF(3, 1)
    K = algebraic_closure(F)
    @test base_field(K) === F
    @test characteristic(K) == characteristic(F)

    @inferred algebraic_closure(F)
    @test K === algebraic_closure(F)
    @test K isa AlgClosure
    @test elem_type(K) === AlgClosureElem{typeof(F)}
    @test elem_type(typeof(K)) === AlgClosureElem{typeof(F)}
    @test parent_type(AlgClosureElem{typeof(F)}) === AlgClosure{typeof(F)}
    @test parent_type(one(K)) === AlgClosure{typeof(F)}

    a = @inferred K()
    @test a isa AlgClosureElem
    @test parent(a) === K

    a = @inferred K(1)
    @test parent(a) === K
    @test a isa AlgClosureElem
    @test isone(a)
    @test isone(one(a))
    @test !iszero(a)

    a = @inferred K(0)
    @test parent(a) === K
    @test a isa AlgClosureElem
    @test iszero(a)
    @test iszero(zero(a))
    @test !isone(a)

    c = one(K) + one(K)
    for T in [Int, BigInt, ZZRingElem]
      b = @inferred K(T(2))
      @test b == c
    end

    b = deepcopy(a)
    @test a !== b
    @test a == b
  end

  @testset "Printing" begin
    F = GF(3, 1)
    K = algebraic_closure(F)
    @test sprint(show, "text/plain", K) == "Algebraic Closure of Finite field of degree 1 over GF(3)"

    for x in F
      @test sprint(show, "text/plain", K(x)) == sprint(show, "text/plain", x)
    end
  end

  @testset "Coercion and conversion" begin
    F1 = GF(3, 1)
    K = algebraic_closure(F1)
    F2 = ext_of_degree(K, 2)
    a = gen(F2)
    c = K(a)
    @test F2(c) == a
    F3 = ext_of_degree(K, 3)
    b = gen(F3)
    d = K(b)
    @test F3(d) == b
    @test_throws ErrorException F2(d)
  end

  @testset "Arithmetic" begin
    F1 = GF(3, 1)
    K = algebraic_closure(F1)
    @test is_unit(K(gen(F1)))
    @test !is_unit(zero(K))
    a = one(K)
    @test isone(inv(a))
    for i in 1:10
      a = test_elem(K)
      b = test_elem(K)
      c = test_elem(K)
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

  @testset "Ad hoc operations" begin
    F1 = GF(3, 1)
    K = algebraic_closure(F1)
    for i in 1:10
      a = test_elem(K)
      for T in [Int, BigInt, ZZRingElem]
        b = rand(-10:10)
        @test a * b == a * K(b)
        @test b * a == a * K(b)
        @test a + b == a + K(b)
        @test b + a == a + K(b)
        @test a - b == a - K(b)
        @test b - a == K(b) - a
        if !iszero(b % characteristic(K))
          @test //(a * b, b) == a
        end
      end
    end
  end

  @testset "Comparison" begin
    p = 3
    F1 = GF(p, 1)
    K = algebraic_closure(F1)
    F2 = ext_of_degree(K, 2)
    F3 = ext_of_degree(K, 3)
    F4 = ext_of_degree(K, 4)
    a = K(gen(F1))
    b = K(gen(F2))
    c = K(gen(F4))

#TODO: make the setup type stable
#   @test @inferred a == b^(p+1)
    @test a == b^(p+1)
    @test b == c^(p^2+1)

    @test hash(a, zero(UInt)) == hash(b^(p+1), zero(UInt))
  end

  @testset "Polynomial" begin
    p = 3
    F1 = GF(p, 1)
    K = algebraic_closure(F1)
    Kx, x = K["x"]
    @test (x^2 + 1)(K(1)) == 2*K(1)
  end

  @testset "Matrix" begin
    p = 3
    F1 = GF(p, 1)
    K = algebraic_closure(F1)
    z = zero(K)
    o = one(K)
    M = matrix(K, 2, 2, [o, o, z, o])
    @test M^2 == M * M
  end

  @testset "Singular ring" begin
    p = 3
    F1 = GF(p, 1)
    K = algebraic_closure(F1)
    L = Oscar.singular_coeff_ring(K)
    a = K(gen(F1))
    @test K(L(a)) == a
  end
end
