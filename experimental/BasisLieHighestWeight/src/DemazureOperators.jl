# lambda, mu in w_i
# beta in alpha_i

# ZZ_x, x = PolynomialRing(ZZ, 3)

function monomial_from_degrees(ZZ_x::ZZMPolyRing, degs::Vector{Int})
    vars = gens(ZZ_x)


    if length(degs) != length(vars)
        throw(ArgumentError("Length of degree vector must match the number of variables in the polynomial ring!"))
    end

    monomial = prod(v^d for (v, d) in zip(vars, degs))

    return monomial
end

function is_monomial(p::Oscar.ZZMPolyRingElem)
    return length(terms(p)) == 1
end

function demazure_scalar_prod(beta::Int, lambda::Vector{Int})
    return lambda[beta]
end 

function demazure_s(beta::Int, lambda::Vector{Int})
    new_lambda = copy(lambda)
    new_lambda[beta] = 0
    return new_lambda
end

function demazure_operator_monom(ZZ_x::ZZMPolyRing, beta::Int, e_lambda::Oscar.ZZMPolyRingElem)    
    lambda = degrees(e_lambda)
    scalar_prod = demazure_scalar_prod(beta, lambda)
    if scalar_prod >= 0
        result = ZZ_x(1)
        for i in 0:lambda[beta]
            result += monomial_from_degrees(ZZ_x, lambda)
            lambda[beta] -= 1
        end
    elseif scalar_prod == -1
        result = ZZ_x(0)
    else
        result = ZZ_x(1)
        for i in 0:lambda[beta]
            result += monomial_from_degrees(ZZ_x, lambda)
            lambda[beta] -= 1
        end
    end

    return result
end

function demazure_operator(ZZ_x::ZZMPolyRing, beta::Int, e_lambda::ZZMPolyRingElem)
    monoms = terms(e_lambda)
    result_poly = zero(e_lambda)
    
    for monom in monoms
        result_poly += demazure_operator_monom(ZZ_x, beta, monom)
    end
    
    return result_poly
end


function demazure_operator_monom_sum(ZZ_x::ZZMPolyRing, beta::Int, e_lambda::ZZMPolyRingElem, e_mu::ZZMPolyRingElem)
    if !is_monomial(e_lambda) || !is_monomial(e_mu)
        throw(ArgumentError("The input must be a monomial!"))
    end
    e_s_beta_mu = monomial_from_degrees(ZZ_x, demazure_scalar_prod(beta, degrees(e_mu)))
    return e_lambda*demazure_operator_monom(ZZ_x, beta, e_mu) 
            + e_s_beta_mu*demazure_operator(ZZ_x, beta, e_mu)
            - e_s_beta_mu*e_lambda
end