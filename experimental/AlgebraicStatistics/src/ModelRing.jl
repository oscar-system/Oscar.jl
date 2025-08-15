struct ModelRing
  R::Ring
  index_to_gen::Dict{T, RingElem}
  gen_to_index::Dict{RingElem, T}

  function ModelRing(S::Ring, varnames::Vector{VarName})
    
  end
end

