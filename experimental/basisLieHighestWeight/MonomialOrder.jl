# manages methods for different orders of monomials

function highest_calc_sub_monomial(mon::Vector{Int}, calc_monomials)
    """
    returns the key in calc_monomials that can be extended by the least amount of left-operations to mon
    """
    sub_mon = copy(mon)

    m = length(mon)
    for i in 1:m
        while sub_mon[i] > 0
            #println(sub_mon)
            #println(calc_monomials)
            if haskey(calc_monomials, sub_mon)
                return sub_mon
            else
                sub_mon[i] -= 1
            end
        end
    end
    #println(sub_mon)
    return sub_mon # [0 for i in 1:m]
end

function lt_monomial_order(monomial_order)
    """
    Returns the desired monomial_order function
    """
    if isa(monomial_order, Function)
        return monomial_order
    elseif monomial_order == "GRevLex"
        return lt_grevlex
    elseif monomial_order == "RevLex"
        return lt_rlex
    elseif monomial_order == "Lex"
        return lt_lex
    elseif monomial_order == "GLex"
        return lt_glex
    else
        println("no suitable order picked")
    end
end

function lt_grevlex(mon1, mon2)
    """
    graded reverse lexicographic order
    returns true if mon1 is less than mon2
    """
    if sum(mon1) == sum(mon2)
        for i=1:length(mon1)
            if mon1[i] < mon2[i]
                return false
            elseif mon1[i] > mon2[i]
                return true
            end
        end
        return false
    else
        return sum(mon1) < sum(mon2)
    end
end

function lt_glex(mon1, mon2)
    """
    graded lexicographic order
    returns true if mon1 is less than mon2
    """
    if sum(mon1) == sum(mon2) 
        for i=length(mon1):-1:1
            if mon1[i] < mon2[i]
                return true
            elseif mon1[i] > mon2[i]
                return false
            end
        end
        return false
    else
        return sum(mon1) < sum(mon2)
    end
end

function lt_lex(mon1, mon2)
    """
    lexicographic order
    returns true if mon1 is less than mon2
    """
    for i=length(mon1):-1:1
        if mon1[i] < mon2[i]
            return true
        elseif mon1[i] > mon2[i]
            return false
        end
    end
    return false
end

function lt_rlex(mon1, mon2)
    """
    reverse lexicographic order
    returns true if mon1 is less than mon2
    """
    for i=1:length(mon1)
        if mon1[i] < mon2[i]
            return true
        elseif mon1[i] > mon2[i]
            return false
        end
    end
    return false
end