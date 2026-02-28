@testset "LAG Group conformace test for $(LAGname)" for (LAGname, LAG) in [
  ("A4", linear_algebraic_group(root_system(:A, 4), GF(7))),
  ("A5", linear_algebraic_group(:A, 5, GF(5))),
  #("A2", linear_algebraic_group(root_system(:A, 2) ,QQ)),
]
  ConformanceTests.test_Group_interface(LAG)
  ConformanceTests.test_GroupElem_interface(rand(LAG, 2)...)
end

@testset "LAG subgroups test $(LAGname)" for (LAGname, LAG) in [
  ("A4", linear_algebraic_group(root_system(:A, 4), GF(7))),
  ("A5", linear_algebraic_group(:A, 5, GF(5))),
  ("A2", linear_algebraic_group(root_system(:A, 2), GF(3))),
]
  for alpha in positive_roots(root_system(LAG))
    U = root_subgroup(LAG, alpha)
    b, _ = is_subgroup(U, LAG)
    @test b
  end
  T = maximal_torus(LAG)
  B = borel_subgroup(LAG)
  b, _ = is_subgroup(B, LAG)
  @test b
  b, _ = is_subgroup(T, LAG)
  @test b
end

@testset "Bruhat decomposition" begin
  LAG = linear_algebraic_group(root_system(:A, 2), GF(3))
  bruh = bruhat_decomp(LAG)
  W = weyl_group(root_system(LAG))
  @test length(bruh) == order(W)
  check = 0
  for cell in bruh
    check += order(cell)
  end
  @test check == order(LAG)
end
