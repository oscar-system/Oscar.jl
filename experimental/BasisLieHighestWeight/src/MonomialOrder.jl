function get_monomial_order_lt(monomial_order::Union{String, Function}, ZZx::ZZMPolyRing)::Function
    """
    Returns the desired monomial_order function less than
    """
    if isa(monomial_order, Function)
       choosen_monomial_order = monomial_order
    else
        # Ensure that `monomial_order` is a valid function before trying to call it
        if isdefined(Main, Symbol(monomial_order))
            x = gens(ZZx)
            choosen_monomial_order = eval(Symbol(monomial_order))(x)
        else
            error("No monomial_order: $monomial_order")
        end
    end
    return (mon1, mon2) -> (cmp(choosen_monomial_order, mon1, mon2) == -1)
end
