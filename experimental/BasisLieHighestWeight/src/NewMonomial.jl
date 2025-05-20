function weight(mon::ZZMPolyRingElem, birational_seq::BirationalSequence)
  @assert length(birational_seq) == nvars(parent(mon))
  return sum(
    exp * weight for
    (exp, weight) in zip(degrees(mon), operators_as_weights(birational_seq));
    init=zero(weight_lattice(root_system(birational_seq))),
  )
end

function calc_vec(
  v0::SRow{ZZRingElem},
  mon::ZZMPolyRingElem,
  matrices_of_operators::Vector{<:SMat{ZZRingElem}},
)
  v = v0
  degree_mon = degrees(mon)
  for i in length(degree_mon):-1:1
    for _ in 1:degree_mon[i]
      v = v * matrices_of_operators[i]
    end
  end
  return v
end
