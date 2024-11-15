struct BirationalSequence
  operator_roots::Vector{RootSpaceElem}
  operator_weights::Vector{WeightLatticeElem}
end

function Base.show(io::IO, birational_sequence::BirationalSequence)
  println(io, "BirationalSequence")
  println(io, "Operators: ", birational_sequence.operator_roots)
  print(io, "Operators as weights:", birational_sequence.operator_weights)
end
