@testset "QuantumGroups.QuantumGroup" begin
  @testset "negative_chevalley_gens" begin
    R = root_system(U)
    F = negative_chevalley_gens(U)
    @test length(F) == rank(R)
    for i in 1:rank(R)
      @test is_gen_with_index(F[i]) == true, U.cvx[i]
    end
  end
  
  @testset "bar_automorphism" begin
    bar = bar_automorphism(U)
    for g in gens(U)
      @test bar(bar(g)) == g
    end
    for g in negative_chevalley_gens(U)
      @test bar(g) == g
    end
  end

  @testset "canonical_basis_elem" begin
    # test bar invariance
    # test coefficients
    # test PBW monomial
    
    # TODO: once bilinear form is implemented, test (b, b) = 1 + qA
  end
  
  @testset "canonical_basis_expansion" begin
    # property tests
    # - expansion of a PBW monomial

    # - expansion of a canonical basis element
    b = rand(1:3, ngens(U))
    elem = canonical_basis_elem(U, b)
    exp = canonical_basis_expansion(elem)
    @test only(exp) == (one(coefficient_ring(U)), b)
    
    # example tests
    U = quantum_group(:A, 4)
    F = negative_chevalley_gens(U)
    q = gen(coefficient_ring(U))
    
    f = F[1]^4*F[2]^3*F[3]^2*F[4]*F[1]*F[2]^2*F[3]^3*F[4]^4
    exp = canonical_basis_expansion(f)
    @test length(exp) == 12
    @test exp[1] == (one(q), [4, 1, 2, 0, 2, 2, 0, 0, 1, 4])
    @test exp[2] == ((q^2+1)//q, [5, 0, 5, 0, 0, 5, 0, 0, 0, 5])
    @test exp[3] == ((q^6 + 2*q^4 + 2*q^2 + 1)//q^3, [4, 1, 3, 0, 1, 3, 0, 0, 1, 4])
    @test exp[10] == ((q^14 + 4*q^12 + 8*q^10 + 11*q^8 + 11*q^6 + 8*q^4 + 4*q^2 + 1)//q^7, [5, 0, 4, 0, 1, 4, 0, 0, 0, 5])
    @test exp[11] == ((q^14 + 3*q^12 + 6*q^10 + 8*q^8 + 8*q^6 + 6*q^4 + 3*q^2 + 1)//q^7, [5, 0, 5, 0, 0, 4, 0, 0, 1, 4])
    @test exp[12] == ((q^20 + 4*q^18 + 9*q^16 + 15*q^14 + 20*q^12 + 22*q^10 + 20*q^8 + 15*q^6 + 9*q^4 + 4*q^2 + 1)//q^10, [5, 0, 5, 0, 0, 5, 0, 0, 0, 5])
    
    U = quantum_group(:G, 2)
    F = negative_chevalley_gens(U)
    q = gen(coefficient_ring(U))
    
    f = F[1]*F[2]^2*F[1]^2*F[2]^2*F[1]
    exp = canonical_basis_expansion(f)
    @test length(exp) == 2
    @test exp[1] == (one(q), [1, 0, 1, 0, 1, 0])
    @test exp[2] == (one(q), [2, 0, 0, 0, 2, 0])
  end
end
