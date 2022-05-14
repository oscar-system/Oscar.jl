@testset "Number field" begin

  Qx, x = FlintQQ["x"]
  k, _ = number_field(x^2 + 1)
  ku, u = k["u1", "u2"]
  Ik = ideal(ku, [u[1]^3 + u[2]^3 - 3, u[1]^5 + u[2]^5 - 5])

  Qy, y = FlintQQ["y1", "y2"]
  IQ = ideal(Qy, [y[1]^3 + y[2]^3 - 3, y[1]^5 + y[2]^5 - 5])

  for (Bk, Pk, I) in [(k, ku, Ik), (FlintQQ, Qy, IQ)]
    gg = gens(Pk)
    @test_throws ErrorException number_field(ideal([gg[1]]))
    K,  = @inferred number_field(I, [:a1, :a2])
    @assert symbols(K) == [:a1, :a2]
    @test_throws ArgumentError number_field(I, [:a1])
    @test_throws ArgumentError number_field(I, [:a1, :a2, :a3])

    K, = @inferred number_field(I, ["a1", "a2"])
    @test symbols(K) == [:a1, :a2]
    @test_throws ArgumentError number_field(I, ["a1"])
    @test_throws ArgumentError number_field(I, ["a1", "a2", "a3"])

    K, = @inferred number_field(I, "a")
    @test symbols(K) == [:a1, :a2]

    K, = @inferred number_field(I)

    K, a = @inferred number_field(I, :a)
    @test symbols(K) == [:a1, :a2]

    # parent and element type
    @inferred parent_type(a[1])
    @inferred elem_type(K)
    @test parent_type(elem_type(typeof(K))) == typeof(K)
    @test elem_type(parent_type(a[1])) == typeof(a[1])

    # field access
    @test Bk == @inferred base_field(K)
    @test Pk == @inferred Oscar.polynomial_ring(K)
    @test I == @inferred Oscar.defining_ideal(K)
    @test 12 == @inferred degree(K)
    @test K.(gens(Pk)) == @inferred gens(K)
    @test K(gen(Pk, 1)) == @inferred gen(K, 1)
    @test K(gen(Pk, 2)) == @inferred gen(K, 2)
    @test gen(Pk, 1) == @inferred Hecke.data(a[1])
    @test K === @inferred parent(a[1])
    @test 2 == @inferred ngens(K)
    @test [:a1, :a2] == @inferred symbols(K)

    # string i/o
    s = sprint(show, "text/plain", K)
    @test s isa String
    s = sprint(show, "text/plain", a[1])
    @test s isa String

    @test (@inferred check_parent(a[1], a[1])) === nothing
    KK, aa = number_field(I, [:a1, :a2])
    @test_throws ArgumentError check_parent(a[1], aa[1])

    # Ring interface

    b = @inferred one(K)
    @test @inferred isone(b)
    b = @inferred zero(K)
    @test @inferred iszero(b)
    @test a[1] == @inferred canonical_unit(a[1])

    # Arithmetic

    for i in 1:10
      b = rand(K, -2:2)
      c = rand(K, -2:2)
      d = rand(K, -2:2)
      e = @inferred b + c
      @test e == c + b
      f = @inferred e * d
      @test f == d * e
      @test d * b + d * c == f
      e = @inferred b - c
      f = @inferred e * d
      @test d * b - d * c == f

      @test one(K) * b == b
      @test zero(K) * b == zero(K)
      @test iszero(b + @inferred(-b))

      @test isone(@inferred b^0)
      @test b == b^1
      @test b * b == b^2
      @test b^10 == reduce(*, [b for j in 1:10])
      @test b^11 == reduce(*, [b for j in 1:11])

      while iszero(b)
        b = rand(K, -2:2)
      end
      bi = @inferred inv(b)
      @test isone(bi * b)

      e = @inferred divexact(c, b)
      @test e * b == c

      e = @inferred c//b
      @test e * b == c

    end

    b = zero(K)
    @test iszero(b^1)
    @test_throws ArgumentError inv(b)

    # Ad hoc arithmetic

    for i in 1:10
      b = rand(K, -1:1)
      for R in Any[Bk, Int32, Int, BigInt, fmpz,
                   Base.Rational{Int}, Base.Rational{BigInt}, fmpq]
        @test b == @inferred (b + R(0))
        @test b == @inferred (R(0) + b)
        @test b == @inferred (b * R(1))
        @test b == @inferred (R(1) * b)
        @test b == @inferred b//R(1)
        @test b == @inferred divexact(b, R(1))
        if !iszero(b)
          @test inv(b) == @inferred R(1)//b
          @test inv(b) == @inferred divexact(R(1), b)
        end
      end
    end

    b = rand(K, -1:1)
    bb = deepcopy(b)
    @test parent(bb) === K
    @test b !== bb

    # In place operations
    for i in 1:10
      b = rand(K, -2:2)
      c = rand(K, -2:2)
      d = rand(k, -2:2)
      e = zero(K)
      @inferred mul!(e, b, c)
      @test e == b * c
      @inferred add!(e, b, c)
      @test e == b + c

      e = deepcopy(b)
      @inferred addeq!(b, c)
      @test b == e + c
    end

    # comparison
    for i in 1:10
      b = rand(K, -2:2)
      @test (@inferred b == b)
      @test (@inferred b != (b + 1))
      for R in Any[Bk, Int32, Int, BigInt, fmpz,
                   Base.Rational{Int}, Base.Rational{BigInt}, fmpq]
        o = one(K)
        @test (@inferred o == R(1))
        @test (@inferred o != R(0))
        @test (@inferred R(1) == o)
        @test (@inferred R(0) != o)

      end
    end

    # parent call overloading
    b = @inferred K()
    @test parent(b) === K
    b = @inferred K(gen(Pk, 1))
    @test parent(b) === K

    PP, _x = Bk["x1", "x2", "x3"]
    @test_throws ErrorException K(_x[1])

    for R in Any[Bk, Int, BigInt, fmpz,
                Base.Rational{Int}, Base.Rational{BigInt}, fmpq]
      b = @inferred K(R(1))
    end

    # denominator
    if Bk === FlintQQ
      b = rand(K, -2:2)
      d = @inferred denominator(b)
      @test d isa fmpz
    end

    # basis
    B = @inferred basis(K)
    @test length(B) == degree(K)
    BA = @inferred absolute_basis(K)
    @test length(BA) == absolute_degree(K)

    for i in 1:10
      b = rand(K, -2:2)
      v = @inferred coordinates(b)
      @test b == sum(v[i] * B[i] for i in 1:degree(K))

      v = @inferred absolute_coordinates(b)
      @test length(v) == absolute_degree(K)
      @test b == sum(v[i] * BA[i] for i in 1:length(BA))
    end

    # basis matrix
    if Bk == FlintQQ
      for i in 1:10
        BB = [rand(K, -2:2) for j in 1:rand(1:10)]
        M = @inferred basis_matrix(BB, FakeFmpqMat)
        @test nrows(M) == length(BB)
        @test ncols(M) == degree(K)
        for n in 1:nrows(M)
          @test BB[n] == sum(M[n, m] * B[m] for m in 1:length(B))
        end
       end
    end

    # representation matrix
    for i in 1:10
      b = rand(K, -2:2)
      M = @inferred representation_matrix(b)
      @test nrows(M) == degree(K)
      @test ncols(M) == degree(K)
      @test (M isa dense_matrix_type(Bk))
      for n in 1:length(B)
        @test b * B[n] == sum(M[n, m] * B[m] for m in 1:length(B))
      end

      if Bk == FlintQQ
        M, d = @inferred representation_matrix_q(b)
        @test nrows(M) == degree(K)
        @test ncols(M) == degree(K)
        @test (M isa fmpz_mat)
        for n in 1:length(B)
          @test b * B[n] == sum(M[n, m]//d * B[m] for m in 1:length(B))
          c = @inferred Oscar.Hecke.elem_from_mat_row(K, M, n, d)
          @test b * B[n] == c
        end

        # elem_to_mat_row!
        c = rand(-10:10)
        while iszero(c)
          c = rand(-10:10)
        end
        b = b//c
        MM = zero_matrix(FlintZZ, nrows(M), ncols(M))
        dd = fmpz()
        j = rand(1:nrows(MM))
        Oscar.Hecke.elem_to_mat_row!(MM, j, dd, b)
        @test b == sum(MM[j, m]//dd * B[m] for m in 1:length(B))
      end
    end

    # minpoly and charpoly
    for i in 1:10
      b = rand(K, -2:2)
      f = @inferred minpoly(b)
      @test iszero(f(b))
      @test is_irreducible(f)

      f = @inferred charpoly(b)
      @test iszero(f(b))
      @test degree(f) == degree(K)
    end

    # trace and norm
    for i in 1:10
      b = rand(K, -2:2)

      f = @inferred charpoly(b)
      @test iszero(f(b))
      @test degree(f) == degree(K)

      t = @inferred tr(b)
      @test parent(t) === Bk
      @test t == -coeff(f, degree(K) - 1)

      t = @inferred norm(b)
      @test parent(t) === Bk
      @test t == (isodd(degree(K)) ? -1 : 1) * coeff(f, 0)

      b = rand(Bk, -10:10)
      t = @inferred tr(K(b))
      @test t == degree(K) * b

      t = @inferred norm(K(b))
      @test t == b^degree(K)
    end

    # (Maximal) order needs some adjustments on the Hecke side
    #@test_broken maximal_order(K)

    # Random
    b = @inferred rand(K, -1:1)
    @test parent(b) === K

    z = @inferred Oscar.primitive_element(K)
    @test degree(minpoly(z)) == degree(K)

    # simple extension
    Ks, KstoK = simple_extension(K)
    for i in 1:10
      b = rand(Ks, -10:10)
      c = rand(Ks, -10:10)
      @inferred KstoK(b)
      @test parent(KstoK(b)) === K
      @test KstoK(b + c) == KstoK(b) + KstoK(c)
      @test KstoK(b * c) == KstoK(b) * KstoK(c)

      @test KstoK\(KstoK(b)) == b
      b = rand(K, -10:10)
      c = rand(K, -10:10)
      @inferred KstoK\(b)
      @test parent(KstoK\(b)) === Ks
      @test KstoK\(b + c) == (KstoK\b) + (KstoK\c)
      @test KstoK\(b * c) == (KstoK\b) * (KstoK\c)
      @test KstoK(KstoK\(b)) == b
    end

    # simple extension
    if Bk == FlintQQ
      Ks, KstoK = simple_extension(K, simplify = true)
    else
      Ks, KstoK = simple_extension(K)
    end

    for i in 1:10
      b = rand(Ks, -10:10)
      c = rand(Ks, -10:10)
      @inferred KstoK(b)
      @test parent(KstoK(b)) === K
      @test KstoK(b + c) == KstoK(b) + KstoK(c)
      @test KstoK(b * c) == KstoK(b) * KstoK(c)

      if Bk == FlintQQ
        @test KstoK\(KstoK(b)) == b
        b = rand(K, -10:10)
        c = rand(K, -10:10)
        @inferred KstoK\(b)
        @test parent(KstoK\(b)) === Ks
        @test KstoK\(b + c) == (KstoK\b) + (KstoK\c)
        @test KstoK\(b * c) == (KstoK\b) * (KstoK\c)
        @test KstoK(KstoK\(b)) == b
      end
    end

    # map
    f = hom(K, K, gens(K))
    @test K === @inferred domain(f)
    @test K === @inferred codomain(f)
    for i in 1:10
      b = rand(K, -10:10)
      @test b == @inferred f(b)
    end

    f = id_hom(K)
    @test f * f == f

    g = hom(K, K, [gen(K, 2), gen(K, 1)])
    @test K === @inferred domain(g)
    @test K === @inferred codomain(g)
    for i in 1:10
      b = rand(K, -10:10)
      @test b == @inferred g(g(b))
    end

    for i in 1:10
      b = rand(K, -10:10)
      @test b == @inferred g\(g(b))
    end

    @test @inferred f == f
    @test @inferred g == g
    @test @inferred f != g
    @test f == @inferred g * g

    f = id_hom(K)
    @test f * f == f
  end
end
