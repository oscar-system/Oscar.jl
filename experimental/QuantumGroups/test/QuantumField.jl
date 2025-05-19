@testset "QuantumGroups.QuantumField" begin
  function ConformanceTests.generate_element(QF::QuantumField)
    return QuantumFieldElem(ConformanceTests.generate_element(QF.d))
  end

  QF, _ = quantum_field()
  #AbstractAlgebra.ConformanceTests.test_Field_interface(QF)
end
