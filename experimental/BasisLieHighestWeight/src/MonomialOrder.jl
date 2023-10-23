function get_monomial_ordering_lt(
  ordering_input::Union{Symbol,Function}, ZZx::ZZMPolyRing
)::Function
  """
  Returns the desired monomial_ordering function less than, i.e. return true <=> mon1 < mon2
  """
  if isa(ordering_input, Function)
    choosen_monomial_order = ordering_input
  else
    if ordering_input == :oplex
      return oplex_lt
    else
      choosen_monomial_order = monomial_ordering(ZZx, ordering_input)
    end
  end
  return (mon1::ZZMPolyRingElem, mon2::ZZMPolyRingElem) ->
    (cmp(choosen_monomial_order, mon1, mon2) < 0)
end

function oplex_lt(mon1::ZZMPolyRingElem, mon2::ZZMPolyRingElem)
  """
  Less-than function for monomials in oplex order
  (mon1, mon2) -> (mon1 < mon2)
  """
  deg1 = degrees(mon1)
  deg2 = degrees(mon2)

  # Start comparing, starting with the first degree
  for i in 1:length(deg1)
    diff = deg1[i] - deg2[i]

    if diff != 0
      return diff > 0 # return mon1 < mon2 if first non-zero of difference is positive
    end
  end

  return false # mon1 == mon2 and therefore not <
end

#function oplex_lt(ZZx::ZZMPolyRing, mon1::ZZMPolyRingElem, mon2::ZZMPolyRingElem)
#    # opposite of lex, return true if first non-zero if - is positive.
#    if degrees(mon1) == degrees(mon2)
#        return false
#    else
#        x = gens(ZZx)
#        lex_order = eval(Symbol("lex"))(x)
#        return (cmp(lex_order, mon1, mon2) == 1)
#    end
#    return false
#end
