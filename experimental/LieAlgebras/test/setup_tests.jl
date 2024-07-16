if !isdefined(Main, :GAPWrap)
  import Oscar: GAPWrap
end

if !isdefined(Main, :lie_algebra_conformance_test)
  function lie_algebra_conformance_test(
    L::LieAlgebra{C}, parentT::DataType, elemT::DataType; num_random_tests::Int=10
  ) where {C<:FieldElem}
    @testset "basic manipulation" begin
      x = L(rand(-10:10, dim(L)))

      @test parentT <: LieAlgebra{C}
      @test elemT <: LieAlgebraElem{C}
      @test L isa parentT
      @test x isa elemT

      @test parent_type(elemT) == parentT
      @test elem_type(parentT) == elemT

      @test parent(x) === L

      @test coefficient_ring(x) === coefficient_ring(L)
      @test elem_type(coefficient_ring(L)) == C

      @test characteristic(L) == characteristic(coefficient_ring(L))

      # this block stays only as long as `ngens` and `gens` are not specialized for Lie algebras
      @test dim(L) == ngens(L)
      @test basis(L) == gens(L)
      @test all(i -> basis(L, i) == gen(L, i), 1:dim(L))

      @test dim(L) == length(basis(L))
      @test all(i -> basis(L, i) == basis(L)[i], 1:dim(L))

      @test dim(L) == length(symbols(L))

      @test iszero(zero(L))

      @test coefficients(x) == [coeff(x, i) for i in 1:dim(L)]
      @test all(i -> coeff(x, i) == x[i], 1:dim(L))
      @test sum(x[i] * basis(L, i) for i in 1:dim(L); init=zero(L)) == x

      @test x == x
      @test deepcopy(x) == x
      @test hash(deepcopy(x)) == hash(x)
    end

    @testset "parent object call overload" begin
      @test L() == zero(L) == L(zeros(coefficient_ring(L), dim(L)))

      for _ in 1:num_random_tests
        coeffs = rand(-10:10, dim(L))
        x1 = L(coeffs)
        x2 = L(coefficient_ring(L).(coeffs))
        x3 = L(matrix(coefficient_ring(L), 1, dim(L), coeffs))
        x4 = L(sparse_row(matrix(coefficient_ring(L), 1, dim(L), coeffs)))
        x5 = L(x1)
        @test x1 == x2
        @test x1 == x3
        @test x1 == x4
        @test x1 == x5
      end
    end

    @testset "vector space axioms" begin
      for _ in 1:num_random_tests
        x = L(rand(-10:10, dim(L)))
        y = L(rand(-10:10, dim(L)))
        z = L(rand(-10:10, dim(L)))

        @test x + y == y + x
        @test x + (y + z) == (x + y) + z

        @test x + zero(L) == x
        @test zero(L) + x == x

        @test -x + x == zero(L)
        @test x + (-x) == zero(L)

        @test x - y == x + (-y)

        @test x * 0 == zero(L)
        @test 0 * x == zero(L)

        @test 2 * x == x + x
        @test x * 2 == x + x
        @test coefficient_ring(L)(2) * x == x + x
        @test x * coefficient_ring(L)(2) == x + x
      end
    end

    @testset "Lie algebra axioms" begin
      for _ in 1:num_random_tests
        x = L(rand(-10:10, dim(L)))
        y = L(rand(-10:10, dim(L)))
        z = L(rand(-10:10, dim(L)))

        @test x * y == bracket(x, y)

        @test (x + y) * z == x * z + y * z
        @test x * (y + z) == x * y + x * z

        @test x * x == zero(L)
        @test x * y == -(y * x)

        @test (x * (y * z)) + (y * (z * x)) + (z * (x * y)) == zero(L)
      end
    end
  end
end

if !isdefined(Main, :lie_algebra_module_conformance_test)
  function lie_algebra_module_conformance_test(
    L::LieAlgebra{C},
    V::LieAlgebraModule{C},
    parentT::DataType=LieAlgebraModule{C},
    elemT::DataType=LieAlgebraModuleElem{C};
    num_random_tests::Int=10,
  ) where {C<:FieldElem}
    @testset "basic manipulation" begin
      v = V(rand(-10:10, dim(V)))

      # @test parentT <: LieAlgebraModule{C}
      # @test elemT <: LieAlgebraModuleElem{C}
      @test V isa parentT
      @test v isa elemT

      @test parent_type(elemT) == parentT
      @test elem_type(parentT) == elemT

      @test parent(v) === V

      @test coefficient_ring(v) === coefficient_ring(V)
      @test elem_type(coefficient_ring(V)) == C

      @test base_lie_algebra(V) === L

      # this block stays only as long as `ngens` and `gens` are not specialized for Lie algebra modules
      @test dim(V) == ngens(V)
      @test basis(V) == gens(V)
      @test all(i -> basis(V, i) == gen(V, i), 1:dim(V))

      @test dim(V) == length(basis(V))
      @test all(i -> basis(V, i) == basis(V)[i], 1:dim(V))

      @test iszero(zero(V))

      @test coefficients(v) == [coeff(v, i) for i in 1:dim(V)]
      @test all(i -> coeff(v, i) == v[i], 1:dim(V))
      @test sum(v[i] * basis(V, i) for i in 1:dim(V); init=zero(V)) == v

      @test v == v
      @test deepcopy(v) == v
      @test hash(deepcopy(v)) == hash(v)
    end

    @testset "parent object call overload" begin
      @test V() == zero(V) == V(zeros(coefficient_ring(V), dim(V)))

      for _ in 1:num_random_tests
        coeffs = rand(-10:10, dim(V))
        v1 = V(coeffs)
        v2 = V(coefficient_ring(V).(coeffs))
        v3 = V(matrix(coefficient_ring(V), 1, dim(V), coeffs))
        v4 = V(sparse_row(matrix(coefficient_ring(V), 1, dim(V), coeffs)))
        v5 = V(v1)
        @test v1 == v2
        @test v1 == v3
        @test v1 == v4
        @test v1 == v5
      end
    end

    @testset "vector space axioms" begin
      for _ in 1:num_random_tests
        v = V(rand(-10:10, dim(V)))
        w = V(rand(-10:10, dim(V)))
        w2 = V(rand(-10:10, dim(V)))

        @test v + w == w + v
        @test v + (w + w2) == (v + w) + w2

        @test v + zero(V) == v
        @test zero(V) + v == v

        @test -v + v == zero(V)
        @test v + (-v) == zero(V)

        @test v - w == v + (-w)

        @test v * 0 == zero(V)
        @test 0 * v == zero(V)

        @test 2 * v == v + v
        @test v * 2 == v + v
        @test coefficient_ring(V)(2) * v == v + v
        @test v * coefficient_ring(V)(2) == v + v
      end
    end

    @testset "Lie algebra action axioms" begin
      for _ in 1:num_random_tests
        x = L(rand(-10:10, dim(L)))
        y = L(rand(-10:10, dim(L)))
        v = V(rand(-10:10, dim(V)))
        w = V(rand(-10:10, dim(V)))

        @test (x * v) isa elemT
        @test parent(x * v) == parent(v)

        @test (x + y) * v == x * v + y * v
        @test x * (v + w) == x * v + x * w

        @test (x * y) * v == x * (y * v) - y * (x * v)
      end
    end
  end
end
