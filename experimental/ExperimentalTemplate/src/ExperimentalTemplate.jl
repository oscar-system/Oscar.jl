# Add your new types, functions, and methods here.

mutable struct ExampleStruct
  i::Int
end

@doc raw"""
    my_access_func(S::ExampleStruct)

This is a dummy sample function to teach the use of docstrings.
"""
function my_access_func(S::ExampleStruct)
  return S.i
end
