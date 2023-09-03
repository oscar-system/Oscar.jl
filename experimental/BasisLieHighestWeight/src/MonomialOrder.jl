function get_monomial_order_lt(monomial_order::Union{String, Function}, ZZx::ZZMPolyRing)::Function
    """
    Returns the desired monomial_order function less than, i.e. return true <=> mon1 < mon2
    """
    if isa(monomial_order, Function)
       choosen_monomial_order = monomial_order
    else
        # Ensure that `monomial_order` is a valid function before trying to call it
        if isdefined(Main, Symbol(monomial_order))
            x = gens(ZZx)
            choosen_monomial_order = eval(Symbol(monomial_order))(x)
        elseif monomial_order == "oplex"
            return oplex_lt
        else
            error("No monomial_order: $monomial_order")
        end
    end
    return (mon1, mon2) -> (cmp(choosen_monomial_order, mon1, mon2) == -1)
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