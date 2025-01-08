if !isdefined(Main, :test_mutating_op_like_zero)
  include(
    joinpath(
      pathof(Oscar.Nemo.AbstractAlgebra), "..", "..", "test", "Rings-conformance-tests.jl"
    ),
  )
end

if !isdefined(Main, :GAPWrap)
  import Oscar: GAPWrap
end
