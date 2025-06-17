@testset "QuantumGroups.QuantumField" begin
  @testset "Conformance Tests" begin
    QF, _ = quantum_field()
    ConformanceTests.test_Field_interface(QF)
  end

  @testset "q_integer" begin
    QF, q = quantum_field()
    for _ in 1:10
      n = rand(0:10)
      @test evaluate(q_integer(n, q).d, 1) == n
    end

    # also covered by doctests
  end

  @testset "q_factorial" begin
    QF, q = quantum_field()
    for _ in 1:10
      n = rand(0:10)
      @test evaluate(q_factorial(n, q).d, 1) == factorial(n)
    end

    # also covered by doctests
  end

  @testset "q_binomial" begin
    for _ in 1:10
      n = rand(0:10)
      k = rand(0:n)
      @test evaluate(q_binomial(n, k, q).d, 1) == binomial(n, k)
    end

    # also covered by doctests
  end
end
