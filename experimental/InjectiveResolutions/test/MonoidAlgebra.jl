@testset "Ring interface for MonoidAlgebra" begin
  a = 2
  function ConformanceTests.generate_element(A::MonoidAlgebra{<:FieldElem, T}) where {T<:MPolyQuoRing}
    return A(rand(base_ring(A.algebra), -a:a, 0:a, 0:a))
  end

  function ConformanceTests.generate_element(A::MonoidAlgebra{<:FieldElem, T}) where {T<:MPolyRing}
    return A(rand(A.algebra, -a:a, 0:a, 0:a))
  end

  S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
  J = ideal(S, [x*z-y^2])
  R_Q, phi = quo(S, J)

  # get MonoidAlgebra
  kQ = Oscar.MonoidAlgebra(R_Q)
  Oscar.ConformanceTests.test_Ring_interface(kQ)
  # Oscar.ConformanceTests.test_Ring_interface_recursive(kQ) # always times out somehow...
end

