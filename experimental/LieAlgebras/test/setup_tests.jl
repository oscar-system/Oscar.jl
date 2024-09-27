include(joinpath(Oscar.oscardir, "test", "Serialization", "setup_tests.jl"))

if !isdefined(Main, :GAPWrap)
  import Oscar: GAPWrap
end

if !isdefined(Main, :lie_algebra_conformance_test) || isinteractive()
  function lie_algebra_conformance_test(
    L::LieAlgebra{C}, parentT::DataType, elemT::DataType; num_random_tests::Int=10
  ) where {C<:FieldElem}
    @testset "basic manipulation" begin
      x = L(rand(-10:10, dim(L)))

      @test parentT <: LieAlgebra{C}
      @test elemT <: LieAlgebraElem{C}
      @test L isa parentT
      @test x isa elemT

      @test (@inferred parent_type(elemT)) == parentT
      @test (@inferred elem_type(parentT)) == elemT

      @test (@inferred parent(x)) === L

      @test (@inferred coefficient_ring(x)) === (@inferred coefficient_ring(L))
      @test (@inferred elem_type(coefficient_ring(L))) == C

      @test (@inferred characteristic(L)) == characteristic(coefficient_ring(L))

      # this block stays only as long as `ngens` and `gens` are not specialized for Lie algebras
      @test (@inferred dim(L)) == (@inferred ngens(L))
      @test (@inferred basis(L)) == (@inferred gens(L))
      dim(L) >= 1 && @test (@inferred basis(L, 1)) == (@inferred gen(L, 1))
      @test all(i -> basis(L, i) == gen(L, i), 1:dim(L))

      @test dim(L) == length(basis(L))
      @test all(i -> basis(L, i) == basis(L)[i], 1:dim(L))

      @test dim(L) == length(symbols(L))

      z = @inferred zero(L)
      @test @inferred iszero(z)

      @test (@inferred coefficients(x)) == [@inferred coeff(x, i) for i in 1:dim(L)]
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
        x1 = @inferred L(coeffs)
        x2 = @inferred L(coefficient_ring(L).(coeffs))
        x3 = @inferred L(matrix(coefficient_ring(L), 1, dim(L), coeffs))
        x4 = @inferred L(sparse_row(matrix(coefficient_ring(L), 1, dim(L), coeffs)))
        x5 = @inferred L(x1)
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

    @testset "Root systems" begin
      if !(L isa DirectSumLieAlgebra) # TODO: make root_system work for DirectSumLieAlgebra
        try
          root_system(L)
        catch e
          e isa ArgumentError || rethrow()
          @test any(
            msg -> contains(e.msg, msg),
            ["Killing form is degenerate", "Cartan subalgebra is not split"],
          )
        end
        if has_root_system(L)
          rs = root_system(L)
          @test rs isa RootSystem
          @test dim(L) == n_roots(rs) + n_simple_roots(rs)
          chev = @inferred chevalley_basis(L)
          @test length(chev) == 3
          @test length(chev[1]) == length(chev[2])
          @test dim(L) == sum(length, chev; init=0)
          H = cartan_subalgebra(L)
          @test all(in(H), chev[3])
          @test all(xs -> bracket(xs...) in H, zip(chev[1], chev[2]))
        end
      end
    end

    @testset "Serialization" begin
      mktempdir() do path
        is_abelian(L) # call something that puts an attribute on L

        test_save_load_roundtrip(
          path,
          L;
          with_attrs=false,
          check_func=loaded -> !has_attribute(loaded, :is_abelian),
        ) do loaded
          # nothing, cause `L === loaded` anyway
        end

        test_save_load_roundtrip(
          path,
          L;
          with_attrs=true,
          check_func=loaded -> all((
            has_attribute(loaded, :is_abelian),
            get_attribute(loaded, :is_abelian) == get_attribute(L, :is_abelian),
            sprint(show, "text/plain", loaded) == sprint(show, "text/plain", L) ||
              occursin(
                "cyclotomic field",
                lowercase(sprint(print, "text/plain", coefficient_ring(L))),
              ), # cyclotomic fields are printed as number fields after (de)serialization          
          )),
        ) do loaded
          # nothing, cause `L === loaded` anyway
        end

        if dim(L) >= 1
          x = basis(L, 1)
          test_save_load_roundtrip(path, x) do loaded
            @test parent(loaded) === L
            @test coefficients(loaded) == coefficients(x)
          end
        end

        if dim(L) >= 1 # TODO: remove this condition once deserializing empty vectors keeps the type (https://github.com/oscar-system/Oscar.jl/issues/3983)
          test_save_load_roundtrip(path, basis(L)) do loaded
            @test length(loaded) == dim(L)
            @test all(
              coefficients(loaded[i]) == coefficients(basis(L, i)) for i in 1:dim(L)
            )
          end
        end
      end
    end
  end
end

if !isdefined(Main, :lie_algebra_module_conformance_test) || isinteractive()
  function lie_algebra_module_conformance_test(
    L::LieAlgebra{C},
    V::LieAlgebraModule{C},
    parentT::DataType=LieAlgebraModule{C,elem_type(L)},
    elemT::DataType=LieAlgebraModuleElem{C,elem_type(L)};
    num_random_tests::Int=10,
  ) where {C<:FieldElem}
    @testset "basic manipulation" begin
      v = V(rand(-10:10, dim(V)))

      # @test parentT <: LieAlgebraModule{C}
      # @test elemT <: LieAlgebraModuleElem{C}
      @test V isa parentT
      @test v isa elemT

      @test (@inferred parent_type(elemT)) == parentT
      @test (@inferred elem_type(parentT)) == elemT

      @test (@inferred parent(v)) === V

      @test (@inferred coefficient_ring(v)) === (@inferred coefficient_ring(V))
      @test (@inferred elem_type(coefficient_ring(V))) == C

      @test (@inferred base_lie_algebra(V)) === L

      # this block stays only as long as `ngens` and `gens` are not specialized for Lie algebra modules
      @test (@inferred dim(V)) == (@inferred ngens(V))
      @test (@inferred basis(V)) == (@inferred gens(V))
      dim(V) >= 1 && @test (@inferred basis(V, 1)) == (@inferred gen(V, 1))
      @test all(i -> basis(V, i) == gen(V, i), 1:dim(V))

      @test dim(V) == length(basis(V))
      @test all(i -> basis(V, i) == basis(V)[i], 1:dim(V))

      z = @inferred zero(V)
      @test @inferred iszero(z)

      @test (@inferred coefficients(v)) == [@inferred coeff(v, i) for i in 1:dim(V)]
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
        v1 = @inferred V(coeffs)
        v2 = @inferred V(coefficient_ring(V).(coeffs))
        v3 = @inferred V(matrix(coefficient_ring(V), 1, dim(V), coeffs))
        v4 = @inferred V(sparse_row(matrix(coefficient_ring(V), 1, dim(V), coeffs)))
        v5 = @inferred V(v1)
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

    if dim(V) <= 50 # for better test runtimes
      @testset "Serialization" begin
        mktempdir() do path
          test_save_load_roundtrip(
            path,
            V;
            with_attrs=false,
          ) do loaded
            # nothing, cause `V === loaded` anyway
          end

          test_save_load_roundtrip(
            path,
            V;
            with_attrs=true,
            check_func=loaded -> all((
              sprint(show, "text/plain", loaded) == sprint(show, "text/plain", V) ||
                occursin(
                  "cyclotomic field",
                  lowercase(sprint(print, "text/plain", coefficient_ring(V))),
                ), # cyclotomic fields are printed as number fields after (de)serialization
            )),
          ) do loaded
            # nothing, cause `V === loaded` anyway
          end

          if dim(V) >= 1
            v = basis(V, 1)
            test_save_load_roundtrip(path, v) do loaded
              @test parent(loaded) === V
              @test coefficients(loaded) == coefficients(v)
            end
          end

          if dim(V) >= 1 # TODO: remove this condition once deserializing empty vectors keeps the type (https://github.com/oscar-system/Oscar.jl/issues/3983)
            test_save_load_roundtrip(path, basis(V)) do loaded
              @test length(loaded) == dim(V)
              @test all(
                coefficients(loaded[i]) == coefficients(basis(V, i)) for i in 1:dim(V)
              )
            end
          end
        end
      end
    end
  end
end
