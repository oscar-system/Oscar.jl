if !isdefined(Main, :test_Field_interface) || !isdefined(Main, :test_NCRing_interface)
  import Oscar.AbstractAlgebra
  include(joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"))
end
