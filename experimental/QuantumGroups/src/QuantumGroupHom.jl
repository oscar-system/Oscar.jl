function bar_involution(U::QuantumGroup)
  if isdefined(U, :bar_involution)
    return U.bar_involution
  end
end
