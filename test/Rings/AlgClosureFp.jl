@testset "AlgClosureFp" begin
  @testset "Interface for $F" for F in [GF(3, 1), Nemo.Native.GF(3, 1)]
    K = algebraic_closure(GF(3,1))
    ConformanceTests.test_Field_interface(K)
  end

  @testset "Creation for $F" for F in [GF(3, 1), Nemo.Native.GF(3, 1)]
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

  @testset "Printing for $F" for F in [GF(3, 1), Nemo.Native.GF(3, 1)]
    K = algebraic_closure(F)
    if F isa FqField
      @test AbstractAlgebra.PrettyPrinting.detailed(K) == "Algebraic closure of prime field of characteristic 3"
      @test AbstractAlgebra.PrettyPrinting.oneline(K) == "Algebraic closure of prime field of characteristic 3"
      @test AbstractAlgebra.PrettyPrinting.supercompact(K) == "Algebraic closure of GF(3)"
    elseif F isa fqPolyRepField
      @test AbstractAlgebra.PrettyPrinting.detailed(K) == "Algebraic closure of finite field of degree 1 over GF(3)"
      @test AbstractAlgebra.PrettyPrinting.oneline(K) == "Algebraic closure of finite field of degree 1 over GF(3)"
      @test AbstractAlgebra.PrettyPrinting.supercompact(K) == "Algebraic closure of GF(3)"
    else
      error("unreachable")
    end

    for x in F
      @test sprint(show, "text/plain", K(x)) == sprint(show, "text/plain", x)
    end
  end

  @testset "Coercion and conversion for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
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

  @testset "Arithmetic for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
    K = algebraic_closure(F1)
    @test is_unit(K(one(F1)))
    @test !is_unit(zero(K))
    a = one(K)
    @test isone(inv(a))
    for i in 1:10
      a = ConformanceTests.generate_element(K)
      b = ConformanceTests.generate_element(K)
      c = ConformanceTests.generate_element(K)
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
      x = 5
      @test a // x == a // ZZRingElem(x)
    end
  end

  @testset "Ad hoc operations for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
    K = algebraic_closure(F1)
    for i in 1:10
      a = ConformanceTests.generate_element(K)
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

  @testset "Comparison for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
    K = algebraic_closure(F1)
    F2 = ext_of_degree(K, 2)
    F3 = ext_of_degree(K, 3)

#TODO: make the setup type stable
#   @test @inferred K(one(F1)) == K(one(F2))
    @test K(one(F1)) == K(one(F2))
    @test K(one(F2)) == K(one(F3))

    @test hash(K(one(F1)), zero(UInt)) == hash(K(one(F2)), zero(UInt))
    @test hash(K(one(F2)), zero(UInt)) == hash(K(one(F3)), zero(UInt))
  end

  @testset "Embeddings for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
    K = algebraic_closure(F1)
    F2 = ext_of_degree(K, 2)
    F3 = ext_of_degree(K, 3)
    emb = embedding(F2, K)
    for x in [one(F2), gen(F2)]
      img = emb(x)
      @test K(x) == img
      @test preimage(emb, img) == x
    end
    @test has_preimage_with_preimage(emb, K(one(F3)))[1]
    @test preimage(emb, K(one(F3))) == one(F2)
    @test !has_preimage_with_preimage(emb, K(gen(F3)))[1]
  end

  @testset "Polynomial for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
    p = characteristic(F1)
    K = algebraic_closure(F1)
    Kx, x = K[:x]
    @test (x^2 + 1)(K(1)) == 2*K(1)

    r = roots(x^4 -1)
    @test sort!(map(degree, r)) == [1, 1, 2, 2]
  end

  @testset "Matrix for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
    p = characteristic(F1)
    K = algebraic_closure(F1)
    z = zero(K)
    o = one(K)
    M = matrix(K, 2, 2, [o, o, z, o])
    @test M^2 == M * M
  end

  @testset "Singular ring for $F1" for F1 in [GF(3, 1), Nemo.Native.GF(3, 1)]
    p = characteristic(F1)
    K = algebraic_closure(F1)
    L = Oscar.singular_coeff_ring(K)
    a = K(gen(F1))
    @test K(L(a)) == a
  end

  R = algebraic_closure(GF(3,1))
  Kt, t = rational_function_field(R, "t")
  @test sprint(show, t) isa String
  @test is_perfect(R)
end
