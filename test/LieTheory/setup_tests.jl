if !isdefined(Main, :test_mutating_op_like_zero)
  import Oscar.AbstractAlgebra
  include(
    joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl")
  )
end

if !isdefined(Main, :test_Group_interface)
  import Oscar.AbstractAlgebra
  import Oscar.AbstractAlgebra: Group
  include(
    joinpath(pathof(AbstractAlgebra), "..", "..", "test", "Groups-conformance-tests.jl")
  )
end

if !isdefined(Main, :GAPWrap)
  import Oscar: GAPWrap
end
