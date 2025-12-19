@testset "LAG Group conformace test for $(LAGname)" for (LAGname, LAG) in [  
  ("A4", linear_algebraic_group(root_system(:A, 4) ,GF(7))),
  ("A5", linear_algebraic_group(:A, 5,GF(5))),
  #("A2", linear_algebraic_group(root_system(:A, 2) ,QQ)),
]
  ConformanceTests.test_Group_interface(LAG)
  ConformanceTests.test_GroupElem_interface(rand(LAG, 2)...)
end
