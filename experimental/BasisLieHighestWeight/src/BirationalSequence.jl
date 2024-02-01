struct BirationalSequence
  operators::Vector{GAP.Obj}
  operators_vectors::Vector{Vector{Any}}
  weights_w::Vector{Vector{ZZRingElem}}
  weights_alpha::Vector{Vector{QQFieldElem}}
end

function Base.show(io::IO, birational_sequence::BirationalSequence)
  println(io, "BirationalSequence")
  println(io, "Operators: ", birational_sequence.operators)
  print(io, "Weights in alpha_i:", birational_sequence.weights_alpha)
end
