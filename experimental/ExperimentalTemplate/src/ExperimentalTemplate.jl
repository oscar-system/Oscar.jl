# Add your new types, functions, and methods here.

mutable struct ExampleStruct
  i::Int
end

@doc raw"""
    content(S::ExampleStruct)

This is a dummy sample function to teach the use of docstrings.
"""
function content(S::ExampleStruct)
  return S.i
end

