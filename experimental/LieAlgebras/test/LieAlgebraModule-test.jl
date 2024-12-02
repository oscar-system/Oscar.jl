
@testset "LieAlgebras.LieAlgebraModule" begin
  @testset "conformance tests" begin
    @testset "0-dim module of sl_2(QQ)" begin
      L = special_linear_lie_algebra(QQ, 2)
      V = trivial_module(L, 0)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "3-dim trivial module of so_3(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 2)
      V = trivial_module(L, 3)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V of sl_2(QQ) using structure constants" begin
      sc = Matrix{sparse_row_type(elem_type(QQ))}(undef, 3, 2)
      sc[1, 1] = sparse_row(QQ, [1, 2], [0, 0])
      sc[1, 2] = sparse_row(QQ, [1, 2], [1, 0])
      sc[2, 1] = sparse_row(QQ, [1, 2], [0, 1])
      sc[2, 2] = sparse_row(QQ, [1, 2], [0, 0])
      sc[3, 1] = sparse_row(QQ, [1, 2], [1, 0])
      sc[3, 2] = sparse_row(QQ, [1, 2], [0, -1])

      L = special_linear_lie_algebra(QQ, 2)
      V = abstract_module(L, 2, sc)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "λ = [1,1,0] of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = simple_module(L, [1, 1, 0])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "λ = [1,1,0,0] of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = simple_module(L, [1, 1, 0, 0])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "λ = [0,1] of A_2(QQ)" begin
      L = lie_algebra(QQ, :A, 2)
      V = simple_module(L, [0, 1])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "λ = [0,1,1] of B_3(QQ)" begin
      L = lie_algebra(QQ, :B, 3)
      V = simple_module(L, [0, 1, 1])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = standard_module(L)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = standard_module(L)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V of sl_4(GF(3))" begin
      L = special_linear_lie_algebra(GF(3), 4)
      V = standard_module(L)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V of gl_4(GF(2^3))" begin
      L = general_linear_lie_algebra(GF(2, 3), 4)
      V = standard_module(L)
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V^* of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = dual(standard_module(L))
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "(T^2 S^2 V)^* of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = dual(tensor_power(symmetric_power(standard_module(L), 2)[1], 2)[1])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⊕(T^2 V) of sl_4(GF(5))" begin
      L = special_linear_lie_algebra(GF(5), 4)
      V = direct_sum(tensor_power(standard_module(L), 2)[1])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 V ⊕ V^* of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = direct_sum(symmetric_power(standard_module(L), 2)[1], dual(standard_module(L)))
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V ⊕ V ⊕ V of sl_4(CF(4))" begin
      L = special_linear_lie_algebra(cyclotomic_field(4)[1], 4)
      V = direct_sum(standard_module(L), standard_module(L), standard_module(L))
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⊗(T^2 V) of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = tensor_product(tensor_power(standard_module(L), 2)[1])
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 V ⊗ V^* of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = tensor_product(
        symmetric_power(standard_module(L), 2)[1], dual(standard_module(L))
      )
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "V ⊗ V ⊗ V of sl_4(GF(3^2))" begin
      L = special_linear_lie_algebra(GF(3, 2), 4)
      V = tensor_product(standard_module(L), standard_module(L), standard_module(L))
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⋀^2 S^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = exterior_power(symmetric_power(standard_module(L), 2)[1], 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⋀^2 T^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = exterior_power(tensor_power(standard_module(L), 2)[1], 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "⋀^2 V of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = exterior_power(standard_module(L), 3)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 ⋀^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = symmetric_power(exterior_power(standard_module(L), 2)[1], 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 T^2 V of so_4(QQ)" begin
      L = special_orthogonal_lie_algebra(QQ, 4)
      V = symmetric_power(tensor_power(standard_module(L), 2)[1], 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "S^2 V of so_4(CF(4))" begin
      L = special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)
      V = symmetric_power(standard_module(L), 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "T^2 ⋀^2 V of sl_4(QQ)" begin
      L = special_linear_lie_algebra(QQ, 4)
      V = tensor_power(exterior_power(standard_module(L), 2)[1], 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "T^2 S^2 V of sl_4(CF(4))" begin
      L = special_linear_lie_algebra(cyclotomic_field(4)[1], 4)
      V = tensor_power(symmetric_power(standard_module(L), 2)[1], 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end

    @testset "T^2 V of sl_4(GF(2^3))" begin
      L = special_linear_lie_algebra(GF(2, 3), 4)
      V = tensor_power(standard_module(L), 2)[1]
      lie_algebra_module_conformance_test(L, V)
    end
  end

  @testset "module constructions" begin
    module_type_bools(V) = (
      Oscar._is_standard_module(V),
      Oscar._is_dual(V)[1],
      Oscar._is_direct_sum(V)[1],
      Oscar._is_tensor_product(V)[1],
      Oscar._is_exterior_power(V)[1],
      Oscar._is_symmetric_power(V)[1],
      Oscar._is_tensor_power(V)[1],
    )

    @testset for R in [QQ, cyclotomic_field(4)[1], GF(3), GF(2, 3)]
      @testset for L in [special_linear_lie_algebra(R, 3)]
        @testset "standard_module" begin
          V = standard_module(L)
          @test V === standard_module(L)
          @test V !== standard_module(L; cached=false)
          @test dim(V) == L.n
          @test length(repr(V)) < 10^4 # outputs tend to be excessively long due to recursion

          @test module_type_bools(V) == (true, false, false, false, false, false, false) # standard_module

          x = L(rand(-10:10, dim(L)))
          v = V(rand(-10:10, dim(V)))
          @test x * v == V(matrix_repr(x) * coefficients(v))
        end

        @testset "dual" begin
          V = direct_sum(standard_module(L), tensor_power(standard_module(L), 2)[1])
          type_V = module_type_bools(V)

          dual_V = dual(V)
          @test dual_V === dual(V)
          @test dual_V !== dual(V; cached=false)
          @test type_V == module_type_bools(V) # construction of dual_V should not change type of V
          @test Oscar._is_dual(dual_V) === (true, V)
          @test dim(dual_V) == dim(V)
          @test length(repr(dual_V)) < 10^4 # outputs tend to be excessively long due to recursion

          @test module_type_bools(dual_V) ==
            (false, true, false, false, false, false, false) # dual

          dual_dual_V = dual(dual_V)
          @test dim(dual_dual_V) == dim(V)
          @test all(
            i ->
              Oscar.LieAlgebras.transformation_matrix(dual_dual_V, i) ==
              Oscar.LieAlgebras.transformation_matrix(V, i),
            1:dim(L),
          )
        end

        @testset "direct_sum" begin
          V = tensor_power(standard_module(L), 2)[1]
          type_V = module_type_bools(V)

          for k in 1:3
            ds_V = direct_sum([V for _ in 1:k]...)
            @test_broken ds_V === direct_sum([V for _ in 1:k]...)
            @test_broken ds_V !== direct_sum([V for _ in 1:k]...; cached=false)
            @test type_V == module_type_bools(V) # construction of ds_V should not change type of V
            @test Oscar._is_direct_sum(ds_V) == (true, [V for _ in 1:k])
            @test all(
              x -> x[1] === x[2], zip(Oscar._is_direct_sum(ds_V)[2], [V for _ in 1:k])
            )
            @test dim(ds_V) == k * dim(V)
            @test length(repr(ds_V)) < 10^4 # outputs tend to be excessively long due to recursion

            @test module_type_bools(ds_V) ==
              (false, false, true, false, false, false, false) # direct_sum

            x = L(rand(-10:10, dim(L)))
            a = [V(rand(-10:10, dim(V))) for _ in 1:k]
            @test ds_V([x * v for v in a]) == x * ds_V(a)
          end

          V1 = dual(standard_module(L))
          V2 = tensor_power(standard_module(L), 2)[1]
          type_V1 = module_type_bools(V1)
          type_V2 = module_type_bools(V2)

          ds_V = direct_sum(V1, V2)
          @test type_V1 == module_type_bools(V1) # construction of ds_V should not change type of V1
          @test type_V2 == module_type_bools(V2) # construction of ds_V should not change type of V2
          @test Oscar._is_direct_sum(ds_V) == (true, [V1, V2])
          @test dim(ds_V) == dim(V1) + dim(V2)
          @test length(repr(ds_V)) < 10^4 # outputs tend to be excessively long due to recursion

          @test module_type_bools(ds_V) == (false, false, true, false, false, false, false) # direct_sum

          x = L(rand(-10:10, dim(L)))
          a = [Vi(rand(-10:10, dim(Vi))) for Vi in [V1, V2]]
          @test ds_V([x * v for v in a]) == x * ds_V(a)
        end

        @testset "tensor_product" begin
          V = direct_sum(standard_module(L), dual(standard_module(L)))
          type_V = module_type_bools(V)

          for k in 1:3
            tp_V = tensor_product([V for _ in 1:k]...)
            @test_broken tp_V === tensor_product([V for _ in 1:k]...)
            @test_broken tp_V !== tensor_product([V for _ in 1:k]...; cached=false)
            @test type_V == module_type_bools(V) # construction of tp_V should not change type of V
            @test Oscar._is_tensor_product(tp_V) == (true, [V for _ in 1:k])
            @test all(
              x -> x[1] === x[2], zip(Oscar._is_tensor_product(tp_V)[2], [V for _ in 1:k])
            )
            @test dim(tp_V) == dim(V)^k
            @test length(repr(tp_V)) < 10^4 # outputs tend to be excessively long due to recursion

            @test module_type_bools(tp_V) ==
              (false, false, false, true, false, false, false) # tensor_product

            x = L(rand(-10:10, dim(L)))
            a = [V(rand(-10:10, dim(V))) for _ in 1:k]
            @test sum(tp_V([i == j ? x * v : v for (j, v) in enumerate(a)]) for i in 1:k) ==
              x * tp_V(a)

            @test tp_V == tensor_power(V, k)[1]
          end

          V1 = dual(standard_module(L))
          V2 = tensor_power(standard_module(L), 2)[1]
          type_V1 = module_type_bools(V1)
          type_V2 = module_type_bools(V2)

          tp_V = tensor_product(V1, V2)
          @test type_V1 == module_type_bools(V1) # construction of tp_V should not change type of V1
          @test type_V2 == module_type_bools(V2) # construction of tp_V should not change type of V2
          @test Oscar._is_tensor_product(tp_V) == (true, [V1, V2])
          @test dim(tp_V) == dim(V1) * dim(V2)
          @test length(repr(tp_V)) < 10^4 # outputs tend to be excessively long due to recursion

          @test module_type_bools(tp_V) == (false, false, false, true, false, false, false) # tensor_product

          x = L(rand(-10:10, dim(L)))
          a = [Vi(rand(-10:10, dim(Vi))) for Vi in [V1, V2]]
          @test sum(tp_V([i == j ? x * v : v for (j, v) in enumerate(a)]) for i in 1:2) ==
            x * tp_V(a)
        end

        @testset "exterior_power" begin
          if characteristic(R) != 0
            @test_throws ArgumentError symmetric_power(standard_module(L), 2)[1]
          else
            V = symmetric_power(standard_module(L), 2)[1]
            type_V = module_type_bools(V)

            for k in 1:3
              E, map = exterior_power(V, k)
              @test E === exterior_power(V, k)[1]
              @test E !== exterior_power(V, k; cached=false)[1]
              @test type_V == module_type_bools(V) # construction of E should not change type of V
              @test Oscar._is_exterior_power(E) === (true, V, k)
              @test dim(E) == binomial(dim(V), k)
              @test length(repr(E)) < 10^4 # outputs tend to be excessively long due to recursion

              @test module_type_bools(E) == (false, false, false, false, true, false, false) # exterior_power

              if k == 1
                x = L(rand(-10:10, dim(L)))
                a = V(rand(-10:10, dim(V)))
                # @test E(x * a) == x * E(a)  # TODO: fix this
                @test E([x * a]) == x * E([a])
              elseif k == 2
                a = V()
                b = V()
                while iszero(a) ||
                        iszero(b) ||
                        can_solve(
                          Oscar.LieAlgebras._matrix(a),
                          Oscar.LieAlgebras._matrix(b);
                          side=:left,
                        )
                  a = V(rand(-10:10, dim(V)))
                  b = V(rand(-10:10, dim(V)))
                end
                @test !iszero(E(a, b))
                @test iszero(E(a, b) + E(b, a))
                @test !iszero(E(a, b) - E(b, a))
                @test iszero(E(a, a))
              end

              as = NTuple{k,elem_type(V)}(V(rand(-10:10, dim(V))) for _ in 1:k)
              @test E(as) == map(as)

              T = get_attribute(E, :embedding_tensor_power)
              E_to_T = get_attribute(E, :embedding_tensor_power_embedding)
              T_to_E = get_attribute(E, :embedding_tensor_power_projection)
              @test T == tensor_power(V, k)[1]
              @test domain(E_to_T) === E
              @test codomain(E_to_T) === T
              @test domain(T_to_E) === T
              @test codomain(T_to_E) === E
              @test compose(E_to_T, T_to_E) == identity_map(E)
            end
          end
        end

        @testset "symmetric_power" begin
          if characteristic(R) != 0
            @test_throws ArgumentError exterior_power(standard_module(L), 2)[1]
          else
            V = exterior_power(standard_module(L), 2)[1]
            type_V = module_type_bools(V)

            for k in 1:3
              S, map = symmetric_power(V, k)
              @test S === symmetric_power(V, k)[1]
              @test S !== symmetric_power(V, k; cached=false)[1]
              @test type_V == module_type_bools(V) # construction of S should not change type of V
              @test Oscar._is_symmetric_power(S) === (true, V, k)
              @test dim(S) == binomial(dim(V) + k - 1, k)
              @test length(repr(S)) < 10^4 # outputs tend to be excessively long due to recursion

              @test module_type_bools(S) == (false, false, false, false, false, true, false) # symmetric_power

              if k == 1
                x = L(rand(-10:10, dim(L)))
                a = V(rand(-10:10, dim(V)))
                # @test S(x * a) == x * S(a)  # TODO: fix this
                @test S([x * a]) == x * S([a])
              elseif k == 2
                a = V()
                b = V()
                while iszero(a) ||
                        iszero(b) ||
                        Oscar.can_solve(
                          Oscar.LieAlgebras._matrix(a),
                          Oscar.LieAlgebras._matrix(b);
                          side=:left,
                        )
                  a = V(rand(-10:10, dim(V)))
                  b = V(rand(-10:10, dim(V)))
                end
                @test !iszero(S(a, b))
                @test !iszero(S(a, b) + S(b, a))
                @test iszero(S(a, b) - S(b, a))
                @test !iszero(S(a, a))
              end

              as = NTuple{k,elem_type(V)}(V(rand(-10:10, dim(V))) for _ in 1:k)
              @test S(as) == map(as)

              T = get_attribute(S, :embedding_tensor_power)
              S_to_T = get_attribute(S, :embedding_tensor_power_embedding)
              T_to_S = get_attribute(S, :embedding_tensor_power_projection)
              @test T == tensor_power(V, k)[1]
              @test domain(S_to_T) === S
              @test codomain(S_to_T) === T
              @test domain(T_to_S) === T
              @test codomain(T_to_S) === S
              @test compose(S_to_T, T_to_S) == identity_map(S)
            end
          end
        end

        @testset "tensor_power" begin
          V = standard_module(L)
          type_V = module_type_bools(V)

          for k in 1:3
            T, map = tensor_power(V, k)
            @test T === tensor_power(V, k)[1]
            @test T !== tensor_power(V, k; cached=false)[1]
            @test type_V == module_type_bools(V) # construction of T should not change type of V
            @test Oscar._is_tensor_power(T) === (true, V, k)
            @test dim(T) == dim(V)^k
            @test length(repr(T)) < 10^4 # outputs tend to be excessively long due to recursion

            @test module_type_bools(T) == (false, false, false, false, false, false, true) # tensor_power

            if k == 1
              x = L(rand(-10:10, dim(L)))
              a = V(rand(-10:10, dim(V)))
              # @test T(x * a) == x * T(a)  # TODO: fix this
              @test T([x * a]) == x * T([a])
            elseif k == 2
              a = V()
              b = V()
              while iszero(a) ||
                      iszero(b) ||
                      can_solve(
                        Oscar.LieAlgebras._matrix(a),
                        Oscar.LieAlgebras._matrix(b);
                        side=:left,
                      )
                a = V(rand(-10:10, dim(V)))
                b = V(rand(-10:10, dim(V)))
              end
              @test !iszero(T(a, b))
              @test !iszero(T(a, b) + T(b, a))
              @test !iszero(T(a, b) - T(b, a))
              @test !iszero(T(a, a))
            end

            as = NTuple{k,elem_type(V)}(V(rand(-10:10, dim(V))) for _ in 1:k)
            @test T(as) == map(as)
          end
        end
      end
    end
  end

  @testset "so_n correctness regression" begin
    function lie_algebra_module_struct_const(
      L::LieAlgebra{C}, V::LieAlgebraModule{C}
    ) where {C<:FieldElem}
      R = coefficient_ring(L)
      dimL = dim(L)
      dimV = dim(V)
      struct_const_V = Matrix{Vector{Tuple{elem_type(R),Int}}}(undef, dimL, dimV)
      for (i, xi) in enumerate(basis(L)), (j, vj) in enumerate(basis(V))
        struct_const_V[i, j] = [
          (c, k) for (k, c) in enumerate(coefficients(xi * vj)) if !iszero(c)
        ]
      end
      return struct_const_V
    end

    L = special_orthogonal_lie_algebra(QQ, 3, identity_matrix(QQ, 3))

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 3)], [], [(QQ(1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    @test lie_algebra_module_struct_const(L, standard_module(L)) == struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 3)], [], [(QQ(1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 1)[1]) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 6)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-2), 2)],
      [(QQ(1), 1), (QQ(-1), 4)],
      [(QQ(-1), 5)],
      [(QQ(2), 2)],
      [(QQ(1), 3)],
      [],
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-2), 3)],
      [(QQ(-1), 5)],
      [(QQ(1), 1), (QQ(-1), 6)],
      [],
      [(QQ(1), 2)],
      [(QQ(2), 3)],
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [],
      [(QQ(-1), 3)],
      [(QQ(1), 2)],
      [(QQ(-2), 5)],
      [(QQ(1), 4), (QQ(-1), 6)],
      [(QQ(2), 5)],
    ]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 2)[1]) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 10)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-3), 2)],
      [(QQ(1), 1), (QQ(-2), 4)],
      [(QQ(-2), 5)],
      [(QQ(2), 2), (QQ(-1), 7)],
      [(QQ(1), 3), (QQ(-1), 8)],
      [(QQ(-1), 9)],
      [(QQ(3), 4)],
      [(QQ(2), 5)],
      [(QQ(1), 6)],
      [],
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-3), 3)],
      [(QQ(-2), 5)],
      [(QQ(1), 1), (QQ(-2), 6)],
      [(QQ(-1), 8)],
      [(QQ(1), 2), (QQ(-1), 9)],
      [(QQ(2), 3), (QQ(-1), 10)],
      [],
      [(QQ(1), 4)],
      [(QQ(2), 5)],
      [(QQ(3), 6)],
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [],
      [(QQ(-1), 3)],
      [(QQ(1), 2)],
      [(QQ(-2), 5)],
      [(QQ(1), 4), (QQ(-1), 6)],
      [(QQ(2), 5)],
      [(QQ(-3), 8)],
      [(QQ(1), 7), (QQ(-2), 9)],
      [(QQ(2), 8), (QQ(-1), 10)],
      [(QQ(3), 9)],
    ]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 4)[1]) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 3)], [], [(QQ(1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 2)[1]) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 3, 3)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[[], [(QQ(-1), 3)], [(QQ(1), 2)]]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(1), 3)], [], [(QQ(-1), 1)]]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[[(QQ(-1), 2)], [(QQ(1), 1)], []]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 3)[1]) ==
      struct_const_V

    L = special_orthogonal_lie_algebra(QQ, 4, identity_matrix(QQ, 4))

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], []
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [], [(QQ(1), 2)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    @test lie_algebra_module_struct_const(L, standard_module(L)) == struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], []
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [], [(QQ(1), 2)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    @test lie_algebra_module_struct_const(L, symmetric_power(standard_module(L), 1)[1]) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], []
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [], [(QQ(1), 2)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 1)[1]) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 6)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 4)], [(QQ(-1), 5)], [(QQ(1), 2)], [(QQ(1), 3)], []
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 4)], [], [(QQ(-1), 6)], [(QQ(-1), 1)], [], [(QQ(1), 3)]
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 5)], [(QQ(1), 6)], [], [], [(QQ(-1), 1)], [(QQ(-1), 2)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], [], [(QQ(-1), 6)], [(QQ(1), 5)]
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 3)], [], [(QQ(1), 1)], [(QQ(1), 6)], [], [(QQ(-1), 4)]
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], [(QQ(-1), 5)], [(QQ(1), 4)], []
    ]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 2)[1]) ==
      struct_const_V

    struct_const_V = Matrix{Vector{Tuple{QQFieldElem,Int64}}}(undef, 6, 4)
    struct_const_V[1, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [], [(QQ(-1), 4)], [(QQ(1), 3)]
    ]
    struct_const_V[2, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(1), 4)], [], [(QQ(-1), 2)]
    ]
    struct_const_V[3, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 4)], [], [], [(QQ(1), 1)]
    ]
    struct_const_V[4, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [], [(QQ(-1), 3)], [(QQ(1), 2)], []
    ]
    struct_const_V[5, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(1), 3)], [], [(QQ(-1), 1)], []
    ]
    struct_const_V[6, :] = Vector{Tuple{QQFieldElem,Int64}}[
      [(QQ(-1), 2)], [(QQ(1), 1)], [], []
    ]
    @test lie_algebra_module_struct_const(L, exterior_power(standard_module(L), 3)[1]) ==
      struct_const_V
  end

  @testset "dim_of_simple_module" begin
    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let L = lie_algebra(QQ, :A, 6)
      dim = @inferred dim_of_simple_module(
        Int, L, [1, 3, 5, 0, 1, 0]
      )
      @test dim isa Int
      @test dim == 393513120
    end

    let R = root_system(:B, 7)
      dim = @inferred dim_of_simple_module(
        ZZRingElem, R, [7, 2, 5, 1, 0, 2, 6]
      )
      @test dim isa ZZRingElem
      @test dim == 307689492858882008424585750
    end

    let L = lie_algebra(QQ, :C, 3)
      dim = @inferred dim_of_simple_module(
        L, [3, 3, 3]
      )
      @test dim isa Int
      @test dim == 262144
    end

    let L = lie_algebra(QQ, :D, 5)
      dim = @inferred dim_of_simple_module(
        Int128, L, [1, 2, 3, 4, 5]
      )
      @test dim isa Int128
      @test dim == 591080490000
    end

    let L = lie_algebra(QQ, :E, 6)
      dim = @inferred dim_of_simple_module(
        ZZRingElem, L, [6, 5, 4, 3, 2, 1]
      )
      @test dim isa ZZRingElem
      @test dim == 53947263633682628459250
    end

    let L = lie_algebra(QQ, :F, 4)
      dim = @inferred dim_of_simple_module(
        BigInt, L, [2, 4, 1, 2]
      )
      @test dim isa BigInt
      @test dim == 5989283015625
    end

    let L = lie_algebra(QQ, :G, 2)
      dim = @inferred dim_of_simple_module(L, ZZ.([2, 2]))
      @test dim isa Int
      @test dim == 729
    end

    let L = special_linear_lie_algebra(QQ, 2) # type A_1 but without known root system
      dim = @inferred dim_of_simple_module(Int8, L, [15])
      @test dim isa Int8
      @test dim == 16
    end

    let L = special_orthogonal_lie_algebra(QQ, 7) # type B_3 but without known root system
      dim = @inferred dim_of_simple_module(L, ZZ.([1, 2, 1]))
      @test dim isa Int
      @test dim == 2800
    end
  end

  @testset "dominant_character" begin
    function check_dominant_character(
      LR::Union{LieAlgebra,RootSystem}, hw::Vector{<:Oscar.IntegerUnion}
    )
      domchar = @inferred dominant_character(LR, hw)
      R = LR isa RootSystem ? LR : root_system(LR)
      @test domchar isa Dict{WeightLatticeElem,Int}
      @test domchar[WeightLatticeElem(R, hw)] == 1
      @test issetequal(keys(domchar), dominant_weights(LR, hw))
      @test all(is_dominant, keys(domchar))
      @test all(>=(1), values(domchar))

      domchar_ZZ = @inferred dominant_character(ZZRingElem, LR, hw)
      @test domchar_ZZ isa Dict{WeightLatticeElem,ZZRingElem}
      @test domchar_ZZ == Dict{WeightLatticeElem,ZZRingElem}(domchar)

      return domchar
    end

    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let R = root_system(Tuple{Symbol,Int}[]), hw = Int[]
      domchar = check_dominant_character(R, hw)
      @test domchar == Dict(WeightLatticeElem(R, Int[]) => 1)
    end

    let L = lie_algebra(QQ, :A, 2), hw = [0, 0]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(WeightLatticeElem(R, [0, 0]) => 1)
    end

    let L = lie_algebra(QQ, :A, 3), hw = ZZ.([1, 1, 1])
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 1, 1]) => 1,
        WeightLatticeElem(R, [2, 0, 0]) => 2,
        WeightLatticeElem(R, [0, 0, 2]) => 2,
        WeightLatticeElem(R, [0, 1, 0]) => 4,
      )
    end

    let L = lie_algebra(QQ, :C, 3), hw = [2, 0, 1]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [2, 0, 1]) => 1,
        WeightLatticeElem(R, [0, 1, 1]) => 1,
        WeightLatticeElem(R, [3, 0, 0]) => 1,
        WeightLatticeElem(R, [1, 1, 0]) => 3,
        WeightLatticeElem(R, [0, 0, 1]) => 6,
        WeightLatticeElem(R, [1, 0, 0]) => 7,
      )
    end

    let L = lie_algebra(QQ, :D, 4), hw = [0, 3, 1, 0]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [0, 3, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 1, 2, 1]) => 1,
        WeightLatticeElem(R, [1, 2, 0, 1]) => 2,
        WeightLatticeElem(R, [2, 0, 3, 0]) => 1,
        WeightLatticeElem(R, [2, 0, 1, 2]) => 2,
        WeightLatticeElem(R, [0, 0, 3, 2]) => 1,
        WeightLatticeElem(R, [2, 1, 1, 0]) => 3,
        WeightLatticeElem(R, [0, 1, 3, 0]) => 2,
        WeightLatticeElem(R, [0, 1, 1, 2]) => 3,
        WeightLatticeElem(R, [0, 2, 1, 0]) => 7,
        WeightLatticeElem(R, [3, 0, 0, 1]) => 4,
        WeightLatticeElem(R, [1, 0, 2, 1]) => 8,
        WeightLatticeElem(R, [1, 0, 0, 3]) => 4,
        WeightLatticeElem(R, [1, 1, 0, 1]) => 12,
        WeightLatticeElem(R, [2, 0, 1, 0]) => 16,
        WeightLatticeElem(R, [0, 0, 3, 0]) => 12,
        WeightLatticeElem(R, [0, 0, 1, 2]) => 16,
        WeightLatticeElem(R, [0, 1, 1, 0]) => 26,
        WeightLatticeElem(R, [1, 0, 0, 1]) => 36,
        WeightLatticeElem(R, [0, 0, 1, 0]) => 56,
      )
    end

    let R = root_system(:E, 6), hw = [1, 0, 1, 0, 1, 0]
      domchar = check_dominant_character(R, hw)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 0, 1, 0, 1, 0]) => 1,
        WeightLatticeElem(R, [0, 0, 0, 1, 1, 0]) => 2,
        WeightLatticeElem(R, [2, 1, 0, 0, 0, 1]) => 3,
        WeightLatticeElem(R, [0, 1, 1, 0, 0, 1]) => 6,
        WeightLatticeElem(R, [2, 0, 1, 0, 0, 0]) => 10,
        WeightLatticeElem(R, [1, 0, 0, 0, 1, 1]) => 16,
        WeightLatticeElem(R, [1, 2, 0, 0, 0, 0]) => 15,
        WeightLatticeElem(R, [0, 0, 2, 0, 0, 0]) => 20,
        WeightLatticeElem(R, [1, 0, 0, 1, 0, 0]) => 44,
        WeightLatticeElem(R, [0, 1, 0, 0, 0, 2]) => 36,
        WeightLatticeElem(R, [0, 1, 0, 0, 1, 0]) => 92,
        WeightLatticeElem(R, [2, 0, 0, 0, 0, 1]) => 104,
        WeightLatticeElem(R, [0, 0, 1, 0, 0, 1]) => 204,
        WeightLatticeElem(R, [1, 1, 0, 0, 0, 0]) => 425,
        WeightLatticeElem(R, [0, 0, 0, 0, 0, 2]) => 416,
        WeightLatticeElem(R, [0, 0, 0, 0, 1, 0]) => 836,
        WeightLatticeElem(R, [1, 0, 0, 0, 0, 0]) => 1600,
      )
    end

    let L = lie_algebra(QQ, :G, 2), hw = [1, 2]
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 2]) => 1,
        WeightLatticeElem(R, [4, 0]) => 1,
        WeightLatticeElem(R, [2, 1]) => 2,
        WeightLatticeElem(R, [0, 2]) => 2,
        WeightLatticeElem(R, [3, 0]) => 3,
        WeightLatticeElem(R, [1, 1]) => 5,
        WeightLatticeElem(R, [2, 0]) => 7,
        WeightLatticeElem(R, [0, 1]) => 7,
        WeightLatticeElem(R, [1, 0]) => 10,
        WeightLatticeElem(R, [0, 0]) => 10,
      )
    end

    let L = special_linear_lie_algebra(QQ, 2), hw = [7] # type A_1 but without known root system
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [7]) => 1,
        WeightLatticeElem(R, [5]) => 1,
        WeightLatticeElem(R, [3]) => 1,
        WeightLatticeElem(R, [1]) => 1,
      )
    end

    let L = special_orthogonal_lie_algebra(QQ, 7), hw = ZZ.([1, 2, 0]) # type B_3 but without known root system
      domchar = check_dominant_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(
        WeightLatticeElem(R, [1, 2, 0]) => 1,
        WeightLatticeElem(R, [2, 0, 2]) => 1,
        WeightLatticeElem(R, [2, 1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, 2]) => 2,
        WeightLatticeElem(R, [0, 2, 0]) => 2,
        WeightLatticeElem(R, [3, 0, 0]) => 2,
        WeightLatticeElem(R, [1, 0, 2]) => 3,
        WeightLatticeElem(R, [1, 1, 0]) => 6,
        WeightLatticeElem(R, [2, 0, 0]) => 6,
        WeightLatticeElem(R, [0, 0, 2]) => 9,
        WeightLatticeElem(R, [0, 1, 0]) => 9,
        WeightLatticeElem(R, [1, 0, 0]) => 15,
        WeightLatticeElem(R, [0, 0, 0]) => 15,
      )
    end
  end

  @testset "character" begin
    function check_character(
      LR::Union{LieAlgebra,RootSystem}, hw::Vector{<:Oscar.IntegerUnion}
    )
      char = @inferred character(LR, hw)
      R = LR isa RootSystem ? LR : root_system(LR)
      @test char isa Dict{WeightLatticeElem,Int}
      @test char[WeightLatticeElem(R, hw)] == 1
      @test all(>=(1), values(char))
      @test sum(values(char)) == dim_of_simple_module(LR, hw)
      domchar = @inferred dominant_character(LR, hw)
      @test all(w -> domchar[w] == char[w], keys(domchar))

      char_ZZ = @inferred character(ZZRingElem, LR, hw)
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}
      @test char_ZZ == Dict{WeightLatticeElem,ZZRingElem}(char)

      return char
    end

    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let R = root_system(Tuple{Symbol,Int}[]), hw = Int[]
      domchar = check_character(R, hw)
      @test domchar == Dict(WeightLatticeElem(R, Int[]) => 1)
    end

    let L = lie_algebra(QQ, :A, 2), hw = [0, 0]
      domchar = check_character(L, hw)
      R = root_system(L)
      @test domchar == Dict(WeightLatticeElem(R, [0, 0]) => 1)
    end

    let L = lie_algebra(QQ, :A, 3), hw = ZZ.([1, 1, 0])
      char = check_character(L, hw)
      R = root_system(L)
      @test char == Dict(
        WeightLatticeElem(R, [1, 1, 0]) => 1,
        WeightLatticeElem(R, [2, -1, 1]) => 1,
        WeightLatticeElem(R, [-1, 2, 0]) => 1,
        WeightLatticeElem(R, [2, 0, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 2,
        WeightLatticeElem(R, [1, -2, 2]) => 1,
        WeightLatticeElem(R, [0, 1, -1]) => 2,
        WeightLatticeElem(R, [-2, 1, 1]) => 1,
        WeightLatticeElem(R, [1, -1, 0]) => 2,
        WeightLatticeElem(R, [-1, -1, 2]) => 1,
        WeightLatticeElem(R, [-2, 2, -1]) => 1,
        WeightLatticeElem(R, [1, 0, -2]) => 1,
        WeightLatticeElem(R, [-1, 0, 0]) => 2,
        WeightLatticeElem(R, [0, -2, 1]) => 1,
        WeightLatticeElem(R, [-1, 1, -2]) => 1,
        WeightLatticeElem(R, [0, -1, -1]) => 1,
      )
    end

    let L = lie_algebra(QQ, :C, 3), hw = [2, 0, 1]
      char = check_character(L, hw)
    end

    let R = root_system(:D, 4), hw = [0, 3, 1, 0]
      char = check_character(R, hw)
    end

    let L = lie_algebra(QQ, :E, 6), hw = [1, 0, 1, 0, 1, 0]
      char = check_character(L, hw)
    end

    let L = lie_algebra(QQ, :G, 2), hw = [3, 2]
      char = check_character(L, hw)
    end

    let L = special_linear_lie_algebra(QQ, 2), hw = [7] # type A_1 but without known root system
      char = check_character(L, hw)
      R = root_system(L)
      @test char == Dict(
        WeightLatticeElem(R, [7]) => 1,
        WeightLatticeElem(R, [5]) => 1,
        WeightLatticeElem(R, [3]) => 1,
        WeightLatticeElem(R, [1]) => 1,
        WeightLatticeElem(R, [-1]) => 1,
        WeightLatticeElem(R, [-3]) => 1,
        WeightLatticeElem(R, [-5]) => 1,
        WeightLatticeElem(R, [-7]) => 1,
      )
    end

    let L = special_orthogonal_lie_algebra(QQ, 7), hw = ZZ.([1, 2, 1]) # type B_3 but without known root system
      char = check_character(L, hw)
    end
  end

  @testset "tensor_product_decomposition" begin
    function test_tensor_product_decomposition(
      LR::Union{LieAlgebra,RootSystem},
      hw1::Vector{<:Oscar.IntegerUnion},
      hw2::Vector{<:Oscar.IntegerUnion},
    )
      dec = @inferred tensor_product_decomposition(LR, hw1, hw2)
      @test dec == @inferred tensor_product_decomposition(LR, hw2, hw1)
      @test multiplicity(dec, Int.(hw1 + hw2)) == 1
      dim_prod = dim_of_simple_module(LR, hw1) * dim_of_simple_module(LR, hw2)
      dim_dec = sum(multiplicity(dec, w) * dim_of_simple_module(LR, w) for w in unique(dec))
      @test dim_prod == dim_dec
      return dec
    end

    # All concrete test results have been computed using the LiE CAS (http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) v2.2.2
    let R = root_system(Tuple{Symbol,Int}[]), hw = Int[]
      dec = test_tensor_product_decomposition(R, hw, hw)
      @test dec == multiset(Dict(Int[] => 1))
    end

    let L = lie_algebra(QQ, :A, 3), hw1 = [2, 1, 2], hw2 = [1, 0, 1]
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [3, 1, 3] => 1,
          [3, 2, 1] => 1,
          [1, 2, 3] => 1,
          [4, 0, 2] => 1,
          [2, 0, 4] => 1,
          [1, 3, 1] => 1,
          [2, 1, 2] => 3,
          [2, 2, 0] => 1,
          [0, 2, 2] => 1,
          [3, 0, 1] => 1,
          [1, 0, 3] => 1,
          [1, 1, 1] => 1,
        ),
      )
    end

    let L = lie_algebra(QQ, :B, 4), hw1 = [1, 1, 0, 0], hw2 = [0, 1, 0, 1]
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [1, 2, 0, 1] => 1,
          [2, 0, 1, 1] => 1,
          [0, 1, 1, 1] => 1,
          [2, 1, 0, 1] => 1,
          [1, 0, 0, 3] => 1,
          [0, 2, 0, 1] => 1,
          [1, 0, 1, 1] => 2,
          [3, 0, 0, 1] => 1,
          [1, 1, 0, 1] => 3,
          [0, 0, 0, 3] => 1,
          [0, 0, 1, 1] => 2,
          [2, 0, 0, 1] => 2,
          [0, 1, 0, 1] => 2,
          [1, 0, 0, 1] => 2,
          [0, 0, 0, 1] => 1,
        ),
      )
    end

    let L = lie_algebra(QQ, :C, 2), hw1 = [2, 2], hw2 = ZZ.([2, 0])
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [4, 2] => 1,
          [2, 3] => 1,
          [4, 1] => 1,
          [0, 4] => 1,
          [2, 2] => 2,
          [4, 0] => 1,
          [0, 3] => 1,
          [2, 1] => 1,
          [0, 2] => 1,
        ),
      )
    end

    let L = lie_algebra(QQ, :D, 5), hw1 = ZZ.([1, 1, 3, 0, 2]), hw2 = [2, 1, 0, 2, 0]
      dec = test_tensor_product_decomposition(L, hw1, hw2)
    end

    let R = root_system(:E, 6), hw1 = [1, 1, 0, 0, 1, 2], hw2 = [2, 0, 1, 1, 0, 0]
      dec = test_tensor_product_decomposition(R, hw1, hw2)
    end

    let R = root_system(:G, 2), hw1 = [1, 3], hw2 = [5, 2]
      dec = test_tensor_product_decomposition(R, hw1, hw2)
    end

    let L = special_linear_lie_algebra(QQ, 2), hw1 = [7], hw2 = [2] # type A_1 but without known root system
      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(Dict([9] => 1,
        [7] => 1,
        [5] => 1,
      ))
    end

    let L = special_orthogonal_lie_algebra(QQ, 7) # type B_3 but without known root system
      hw1 = ZZ.([1, 1, 0])
      hw2 = ZZ.([0, 1, 1])

      dec = test_tensor_product_decomposition(L, hw1, hw2)
      @test dec == multiset(
        Dict(
          [1, 2, 1] => 1,
          [2, 0, 3] => 1,
          [2, 1, 1] => 1,
          [0, 1, 3] => 1,
          [0, 2, 1] => 1,
          [3, 0, 1] => 1,
          [1, 0, 3] => 2,
          [1, 1, 1] => 3,
          [2, 0, 1] => 2,
          [0, 0, 3] => 2,
          [0, 1, 1] => 2,
          [1, 0, 1] => 2,
          [0, 0, 1] => 1,
        ),
      )
    end
  end

  @testset "demazure_character" begin
    function demazure_character_trivial_tests(R::RootSystem, w::WeightLatticeElem)
      W = weyl_group(R)

      char = demazure_character(R, w, longest_element(W))
      @test char == character(R, w)
      @test char isa Dict{WeightLatticeElem,Int}

      char = demazure_character(R, w, one(W))
      @test char == Dict(w => 1)
      @test char isa Dict{WeightLatticeElem,Int}
    end
    # The "correct" solution in the below tests has been computed using Magma V2.28-7.
    # Note that Magma considers right actions of Weyl groups on weights,
    # so to reproduce one needs to reverse all Weyl words when copying a test case to Magma.

    @testset "Type A" begin
      for n in 2:4
        R = root_system(:A, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end

      R = root_system(:A, 2)
      W = weyl_group(R)

      #x = epsilon
      char = demazure_character(R, fundamental_weight(R, 1), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 1]) => 1
      )

      #x = s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[1]*s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [0, -1]) => 1,
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [1, -2]) => 1,
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [0, 0]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[2]*s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 0]) => 1,
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [-2, 1]) => 1,
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [0, 0]) => 1,
        WeightLatticeElem(R, [1, 1]) => 1,
      )

      #x = s[1]*s[2]*s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1]) => 1,
        WeightLatticeElem(R, [0, -1]) => 1,
        WeightLatticeElem(R, [1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 0]) => 1,
        WeightLatticeElem(R, [1, -1]) => 1,
        WeightLatticeElem(R, [0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 2, 1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2]) => 1,
        WeightLatticeElem(R, [-2, 1]) => 1,
        WeightLatticeElem(R, [2, -1]) => 1,
        WeightLatticeElem(R, [1, -2]) => 1,
        WeightLatticeElem(R, [-1, -1]) => 1,
        WeightLatticeElem(R, [0, 0]) => 2,
        WeightLatticeElem(R, [1, 1]) => 1,
      )
    end

    @testset "Type B" begin
      for n in 2:4
        R = root_system(:B, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end

      R = root_system(:B, 3)
      W = weyl_group(R)

      #x = epsilon
      char = demazure_character(R, fundamental_weight(R, 1), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), one(W))
      @test char == Dict(
        WeightLatticeElem(R, [1, 1, 1]) => 1
      )

      #x = s[1]
      char = demazure_character(R, fundamental_weight(R, 1), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), W([1]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2, 1]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
      )

      #x = s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [0, 1, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 1]) => 1
      )

      char = demazure_character(R, weyl_vector(R), W([2]))
      @test char == Dict(
        WeightLatticeElem(R, [2, -1, 3]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
      )

      #x = s[3]
      char = demazure_character(R, fundamental_weight(R, 1), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([3]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 2, -1]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
      )

      #x = s[1]*s[3]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 1, 1]) => 1,
        WeightLatticeElem(R, [-1, 2, 1]) => 1,
        WeightLatticeElem(R, [-1, 3, -1]) => 1,
        WeightLatticeElem(R, [1, 2, -1]) => 1,
      )

      #x=s[1]*s[2]*s[3]
      char = demazure_character(R, fundamental_weight(R, 1), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 0, 0]) => 1,
        WeightLatticeElem(R, [0, 1, -2]) => 1,
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [0, -1, 2]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1,
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [1, 1, -2]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([1, 2, 3]))
      @test char == Dict(
        WeightLatticeElem(R, [2, 0, 1]) => 1,
        WeightLatticeElem(R, [1, 2, -3]) => 1,
        WeightLatticeElem(R, [1, 1, -1]) => 1,
        WeightLatticeElem(R, [1, 1, 1]) => 1,
        WeightLatticeElem(R, [0, 3, -3]) => 1,
        WeightLatticeElem(R, [0, 1, 1]) => 1,
        WeightLatticeElem(R, [-1, 3, -1]) => 1,
        WeightLatticeElem(R, [2, -1, 3]) => 1,
        WeightLatticeElem(R, [0, 0, 3]) => 1,
        WeightLatticeElem(R, [1, -2, 5]) => 1,
        WeightLatticeElem(R, [0, 2, -1]) => 1,
        WeightLatticeElem(R, [-1, 2, 1]) => 1,
        WeightLatticeElem(R, [1, -1, 3]) => 1,
        WeightLatticeElem(R, [2, 2, -3]) => 1,
        WeightLatticeElem(R, [1, 2, -1]) => 1,
        WeightLatticeElem(R, [1, 0, 1]) => 1,
        WeightLatticeElem(R, [1, 3, -5]) => 1,
        WeightLatticeElem(R, [2, 1, -1]) => 1,
      )

      #x=s[3]*s[2]*s[3]*s[1]*s[2]*s[3]*s[1]*s[2]
      char = demazure_character(R, fundamental_weight(R, 1), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [0, -1, 2]) => 1,
        WeightLatticeElem(R, [0, 0, 0]) => 1,
        WeightLatticeElem(R, [1, -1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, -2]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 2), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [-1, 2, -2]) => 1,
        WeightLatticeElem(R, [1, 1, -2]) => 1,
        WeightLatticeElem(R, [-2, 1, 0]) => 1,
        WeightLatticeElem(R, [2, -1, 0]) => 1,
        WeightLatticeElem(R, [0, -1, 0]) => 1,
        WeightLatticeElem(R, [-1, -1, 2]) => 1,
        WeightLatticeElem(R, [-1, 1, 0]) => 1,
        WeightLatticeElem(R, [1, -2, 2]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [-1, 0, 2]) => 1,
        WeightLatticeElem(R, [1, -1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, 0]) => 1,
        WeightLatticeElem(R, [0, 1, -2]) => 1,
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [-1, 1, -2]) => 1,
        WeightLatticeElem(R, [0, -1, 2]) => 1,
        WeightLatticeElem(R, [0, 0, 0]) => 3,
        WeightLatticeElem(R, [1, 0, -2]) => 1,
        WeightLatticeElem(R, [-1, 0, 0]) => 1,
      )

      char = demazure_character(R, fundamental_weight(R, 3), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test char == Dict(
        WeightLatticeElem(R, [0, -1, 1]) => 1,
        WeightLatticeElem(R, [0, 0, 1]) => 1,
        WeightLatticeElem(R, [1, 0, -1]) => 1,
        WeightLatticeElem(R, [-1, 0, 1]) => 1,
        WeightLatticeElem(R, [0, 0, -1]) => 1,
        WeightLatticeElem(R, [-1, 1, -1]) => 1,
        WeightLatticeElem(R, [0, 1, -1]) => 1,
        WeightLatticeElem(R, [1, -1, 1]) => 1,
      )

      char = demazure_character(R, weyl_vector(R), W([3, 2, 3, 1, 2, 3, 1, 2]))
      @test length(char) == 124
      @test char[WeightLatticeElem(R, [-3, 3, -3])] == 1
      @test char[WeightLatticeElem(R, [0, 2, -3])] == 7
      @test char[WeightLatticeElem(R, [0, 1, -1])] == 12
      @test char[WeightLatticeElem(R, [3, 0, -3])] == 2
      @test char[WeightLatticeElem(R, [0, 1, -3])] == 6
    end

    @testset "Type C" begin
      for n in 2:4
        R = root_system(:C, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end
    end

    @testset "Type D" begin
      for n in 4:4
        R = root_system(:D, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end
    end

    @testset "Type E" begin
      for n in 6:6
        R = root_system(:E, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        #demazure_character_trivial_tests(R, weyl_vector(R)) #takes too long for CI
        #demazure_character_trivial_tests(R, 2*weyl_vector(R)) #takes too long for CI
      end
    end

    @testset "Type F" begin
      for n in 4:4
        R = root_system(:F, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        #demazure_character_trivial_tests(R, 2*weyl_vector(R)) #takes too long for CI
      end
    end

    @testset "Type G" begin
      for n in 2:2
        R = root_system(:G, n)
        for i in 1:n
          w = fundamental_weight(R, i)
          demazure_character_trivial_tests(R, w)
        end
        demazure_character_trivial_tests(R, weyl_vector(R))
        demazure_character_trivial_tests(R, 2 * weyl_vector(R))
      end
    end

    @testset "Input Conversion" begin
      R = root_system(:B, 3)
      W = weyl_group(R)
      L = lie_algebra(QQ, R)
      w_int = [0, 1, 0]
      w_weight = WeightLatticeElem(R, w_int)
      x = W([1, 2, 3])
      reduced_expression = word(x)
      result = Dict(
        WeightLatticeElem(R, [0, 1, 0]) => 1,
        WeightLatticeElem(R, [1, -1, 2]) => 1,
        WeightLatticeElem(R, [1, 0, 0]) => 1,
        WeightLatticeElem(R, [1, 1, -2]) => 1,
      )

      char = @inferred demazure_character(R, w_weight, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_weight, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(R, w_weight, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_weight, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(R, w_int, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_int, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(R, w_int, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, R, w_int, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_weight, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_weight, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_weight, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_weight, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_int, x)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_int, x)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}

      char = @inferred demazure_character(L, w_int, reduced_expression)
      @test char == result
      @test char isa Dict{WeightLatticeElem,Int}

      char_ZZ = @inferred demazure_character(ZZRingElem, L, w_int, reduced_expression)
      @test char_ZZ == result
      @test char_ZZ isa Dict{WeightLatticeElem,ZZRingElem}
    end
  end
end
