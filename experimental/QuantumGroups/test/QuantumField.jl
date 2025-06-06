@testset "QuantumGroups.QuantumField" begin
  QF, _ = quantum_field()
  AbstractAlgebra.ConformanceTests.test_Field_interface(QF)
end
