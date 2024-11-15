struct BirationalSequence
  operator_roots::Vector{RootSpaceElem}
  operator_weights::Vector{WeightLatticeElem}

  function BirationalSequence(
    operator_roots::Vector{RootSpaceElem}, operator_weights::Vector{WeightLatticeElem}
  )
    @req length(operator_roots) == length(operator_weights) "Different lengths"
    return new(operator_roots, operator_weights)
  end
end

function birational_sequence(
  operator_roots::Vector{RootSpaceElem}, operator_weights::Vector{WeightLatticeElem}
)
  return BirationalSequence(operator_roots, operator_weights)
end

function birational_sequence(operator_roots::Vector{RootSpaceElem})
  return birational_sequence(operator_roots, WeightLatticeElem.(operator_roots))
end

function birational_sequence(operator_weights::Vector{WeightLatticeElem})
  return birational_sequence(RootSpaceElem.(operator_weights), operator_weights)
end

function Base.show(io::IO, birational_seq::BirationalSequence)
  println(io, "BirationalSequence")
  println(io, "Operators: ", operators_as_roots(birational_seq))
  print(io, "Operators as weights:", operators_as_weights(birational_seq))
end

function length(birational_seq::BirationalSequence)
  return length(birational_seq.operator_roots)
end

function operators_as_roots(birational_seq::BirationalSequence)
  return birational_seq.operator_roots
end

function operator_as_root(birational_seq::BirationalSequence, i::Int)
  return birational_seq.operator_roots[i]
end

function operators_as_weights(birational_seq::BirationalSequence)
  return birational_seq.operator_weights
end

function operator_as_weight(birational_seq::BirationalSequence, i::Int)
  return birational_seq.operator_weights[i]
end
