# lambda, mu in w_i
# beta in alpha_i

# ZZ_x, x = PolynomialRing(ZZ, 3)

function monomial_from_degrees(
    ZZ_x::AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}, 
    degs::Vector{Int}
    )
    vars = gens(ZZ_x)


    if length(degs) != length(vars)
        throw(ArgumentError("Length of degree vector must match the number of variables in the polynomial ring!"))
    end
    
    monomial = prod(v^d for (v, d) in zip(vars, degs))
    return monomial
end

function is_monomial(p::AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing})
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

function demazure_operator_monom(
    ZZ_x::AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}, 
    beta::Int,
    beta_wi::Vector{Int},
    e_lambda::AbstractAlgebra.Generic.LaurentMPolyWrap{ZZRingElem, ZZMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}}
    )
    lambda = leading_exponent_vector(e_lambda)
    scalar_prod = demazure_scalar_prod(beta, lambda)
    if scalar_prod >= 0
        result = ZZ_x(0)
        for i in 0:lambda[beta]
            result += monomial_from_degrees(ZZ_x, lambda)
            lambda -= beta_wi
        end
    elseif scalar_prod == -1
        result = ZZ_x(0)
    else
        result = ZZ_x(0)
        for i in 0:lambda[beta]
            result += monomial_from_degrees(ZZ_x, lambda)
            lambda -= beta_wi
        end
    end
    return result
end

function demazure_operator(
    ZZ_x::AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}, 
    beta::Int,
    beta_wi::Vector{Int},
    e_lambda::AbstractAlgebra.Generic.LaurentMPolyWrap{ZZRingElem, ZZMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}}
    )
    monoms = terms(e_lambda)
    result_poly = zero(e_lambda)
    
    for monom in monoms
        result_poly += demazure_operator_monom(ZZ_x, beta, beta_wi, monom)
    end
    
    return result_poly
end

function demazure_operators_summary(
    type::String,
    rank::Int,
    lambda::Vector{Int},
    weyl_word::Vector{Int}
    )
    alpha_list = [ [i == j ? 1 : 0 for i in 1:rank] for j in 1:rank] # [..., [0, .., 0, 1, 0, ..., 0], ...]
    alpha_wi_list = [alpha_to_w(type, rank, alpha_i) for alpha_i in alpha_list]
    ZZ_x, x = LaurentPolynomialRing(ZZ, length(lambda))
    sub_word = []
    p = monomial_from_degrees(ZZ_x, lambda)
    for alpha_i in reverse(weyl_word)
        append!(sub_word, alpha_i)
        p = demazure_operator(ZZ_x, alpha_i, alpha_wi_list[alpha_i], p)
        println("")
        println("sub_word: ", sub_word)
        println("p: ", p)
        for term in terms(p)
            println(leading_exponent_vector(term), ": ", leading_coefficient(term))
        end
    end

    println("")
    println("Dimension calculated through basis_lie_highest_weight for full word:")
    lie_algebra = LieAlgebraStructure(type, rank)
    println(get_dim_weightspace(lie_algebra, lambda))
end

function demazure_operator_monom_sum(
    ZZ_x::AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}, 
    beta::Int, 
    e_lambda::ZZMPolyRingElem, 
    e_mu::ZZMPolyRingElem
    )
    if !is_monomial(e_lambda) || !is_monomial(e_mu)
        throw(ArgumentError("The input must be a monomial!"))
    end
    e_s_beta_mu = monomial_from_degrees(ZZ_x, demazure_scalar_prod(beta, degrees(e_mu)))
    return e_lambda*demazure_operator_monom(ZZ_x, beta, e_mu) 
            + e_s_beta_mu*demazure_operator(ZZ_x, beta, e_mu)
            - e_s_beta_mu*e_lambda
end