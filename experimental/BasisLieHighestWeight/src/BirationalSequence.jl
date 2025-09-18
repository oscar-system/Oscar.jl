struct BirationalSequence
  operator_roots::Vector{RootSpaceElem}
  operator_weights::Vector{WeightLatticeElem}
  root_system::RootSystem

  function BirationalSequence(
    operator_roots::Vector{RootSpaceElem}, operator_weights::Vector{WeightLatticeElem},
    root_sys::RootSystem,
  )
    @req length(operator_roots) == length(operator_weights) "Different lengths"
    @req all(rootspace_elem -> root_system(rootspace_elem) == root_sys, operator_roots) "Different root systems"
    @req all(weight_elem -> root_system(weight_elem) == root_sys, operator_weights) "Different root systems"
    return new(operator_roots, operator_weights, root_sys)
  end
end

function birational_sequence(
  operator_roots::Vector{RootSpaceElem}, operator_weights::Vector{WeightLatticeElem},
  root_system::RootSystem,
)
  return BirationalSequence(operator_roots, operator_weights, root_system)
end

function birational_sequence(operator_roots::Vector{RootSpaceElem}, root_system::RootSystem)
  return birational_sequence(
    operator_roots, WeightLatticeElem.(operator_roots), root_system
  )
end

function birational_sequence(
  operator_weights::Vector{WeightLatticeElem}, root_system::RootSystem
)
  return birational_sequence(
    RootSpaceElem.(operator_weights), operator_weights, root_system
  )
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

function root_system(birational_seq::BirationalSequence)
  return birational_seq.root_system
end
