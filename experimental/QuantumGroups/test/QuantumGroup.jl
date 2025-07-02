@testset "QuantumGroups.QuantumGroup" begin
  @testset "Conformance Tests" begin
    #ConformanceTests.test_NCRing_interface(quantum_group(:A, 5))
    #ConformanceTests.test_NCRing_interface(quantum_group(:B, 3))
    #ConformanceTests.test_NCRing_interface(quantum_group(:C, 3))
    #ConformanceTests.test_NCRing_interface(quantum_group(:D, 4))
    #ConformanceTests.test_NCRing_interface(quantum_group(:G, 2))
  end

  @testset "negative_chevalley_gens" begin
    for (fam, rk) in [(:A, 5), (:B, 3), (:C, 3), (:D, 4), (:G, 2)]
      U = quantum_group(fam, rk)
      R = root_system(U)
      F = negative_chevalley_gens(U)
      @test length(F) == rank(R)
      for i in 1:rank(R)
        @test is_gen_with_index(F[i]) == (true, U.cvx[i])
      end
    end
  end

  @testset "bar_automorphism" begin
    for (fam, rk) in [(:A, 5), (:B, 3), (:C, 3), (:D, 4), (:G, 2)]
      U = quantum_group(fam, rk)
      bar = bar_automorphism(U)
      for g in gens(U)
        @test bar(bar(g)) == g
      end
      for g in negative_chevalley_gens(U)
        @test bar(g) == g
      end
    end
  end

  @testset "canonical_basis_elem" begin
    for (fam, rk) in [(:A, 5), (:B, 3), (:C, 3), (:D, 4), (:G, 2)]
      U = quantum_group(fam, rk)
      bar = bar_automorphism(U)

      for _ in 1:5
        empty!(U.canonical_basis)
        b = zeros(Int, ngens(U))
        for _ in 1:2
          b[rand(1:ngens(U))] = rand(0:2)
        end
        elem = canonical_basis_elem(U, b)
        @test bar(elem) == elem
        @test exponent_vector(elem, length(elem)) == b
        @test coeff(elem, length(elem)) == one(coefficient_ring(U))

        for i in 1:(length(elem) - 1)
          cf = coeff(elem, i).d
          exp = exponent_vector(elem, i)

          k = 0
          while iszero(coeff(numerator(cf), k))
            k += 1
          end
          @test k > degree(denominator(cf))
          @test is_monomial(denominator(cf))
        end
      end
    end

    # TODO: once bilinear form is implemented, test (b, b) in 1 + qA
  end

  @testset "monomial" begin end

  @testset "canonical_basis_expansion" begin
    # property tests
    for (fam, rk) in [(:A, 5), (:B, 3), (:C, 3), (:D, 4), (:G, 2)]
      U = quantum_group(fam, rk)
      # - expansion of a PBW monomial
      for _ in 1:5
        empty!(U.canonical_basis)
        b = zeros(Int, ngens(U))
        for _ in 1:2
          b[rand(1:ngens(U))] = rand(0:2)
        end
        F = pbw_monomial(U, b)
        exp = canonical_basis_expansion(F)
        @test first(exp) == (one(coefficient_ring(U)), b)
        @test F == sum(c * canonical_basis_elem(U, b) for (c, b) in exp)
      end

      # - expansion of a canonical basis element
      for _ in 1:5
        empty!(U.canonical_basis)
        b = zeros(Int, ngens(U))
        for _ in 1:2
          b[rand(1:ngens(U))] = rand(0:2)
        end
        elem = canonical_basis_elem(U, b)
        exp = canonical_basis_expansion(elem)
        @test only(exp) == (one(coefficient_ring(U)), b)
      end
    end

    # example tests
    U = quantum_group(:A, 4)
    F = negative_chevalley_gens(U)
    q = gen(coefficient_ring(U))

    f = monomial(U, [1, 2, 3, 4, 1, 2, 3, 4], [4, 3, 2, 1, 1, 2, 3, 4])
    exp = canonical_basis_expansion(f)
    @test length(exp) == 12
    @test exp[1] == (one(q), [4, 1, 2, 0, 2, 2, 0, 0, 1, 4])
    @test exp[2] == ((q^2 + 1)//q, [4, 1, 2, 0, 2, 3, 0, 0, 0, 5])
    @test exp[3] == ((q^6 + 2 * q^4 + 2 * q^2 + 1)//q^3, [4, 1, 3, 0, 1, 3, 0, 0, 1, 4])
    @test exp[10] == (
      (q^14 + 4 * q^12 + 8 * q^10 + 11 * q^8 + 11 * q^6 + 8 * q^4 + 4 * q^2 + 1)//q^7,
      [5, 0, 4, 0, 1, 4, 0, 0, 0, 5],
    )
    @test exp[11] == (
      (q^14 + 3 * q^12 + 6 * q^10 + 8 * q^8 + 8 * q^6 + 6 * q^4 + 3 * q^2 + 1)//q^7,
      [5, 0, 5, 0, 0, 4, 0, 0, 1, 4],
    )
    @test exp[12] == (
      (
        q^20 + 4 * q^18 + 9 * q^16 + 15 * q^14 + 20 * q^12 + 22 * q^10 + 20 * q^8 +
        15 * q^6 + 9 * q^4 + 4 * q^2 + 1
      )//q^10,
      [5, 0, 5, 0, 0, 5, 0, 0, 0, 5],
    )

    U = quantum_group(:G, 2)
    F = negative_chevalley_gens(U)
    q = gen(coefficient_ring(U))

    f = monomial(U, [1, 2, 1, 2, 1], [1, 1, 2, 1, 1])
    exp = canonical_basis_expansion(f)
    @test length(exp) == 2
    @test exp[1] == (one(q), [1, 0, 1, 0, 1, 0])
    @test exp[2] == (one(q), [2, 0, 0, 0, 2, 0])
  end
end
