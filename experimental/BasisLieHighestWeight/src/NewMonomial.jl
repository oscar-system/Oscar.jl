function calc_weight(
  mon::ZZMPolyRingElem, weights_w::Vector{Vector{ZZRingElem}}
)::Vector{ZZRingElem}
  """
  calculates weight associated with monomial mon
  """
  degree_mon = degrees(mon)
  weight_w = [ZZ(0) for i in 1:length(weights_w[1])]
  for i in 1:length(degree_mon)
    weight_w .+= degree_mon[i] * weights_w[i]
  end
  return weight_w
end

function calc_vec(
  v0::SRow{ZZRingElem},
  mon::ZZMPolyRingElem,
  matrices_of_operators::Union{
    Vector{SMat{ZZRingElem,Hecke.ZZRingElem_Array_Mod.ZZRingElem_Array}},
    Vector{SMat{ZZRingElem}},
  },
)::SRow{ZZRingElem}
  """
  calculates vector associated with monomial mon
  """
  vec = v0
  degree_mon = degrees(mon)
  for i in length(degree_mon):-1:1
    for j in 1:degree_mon[i]
      # currently there is no sparse matrix * vector mult
      # this is also the line that takes up almost all the computation time for big examples
      vec = mul(vec, transpose(matrices_of_operators[i]))
    end
  end
  return vec
end
