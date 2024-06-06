################################################################################
#
#  Tropicalization of binomial ideals
#
#  WARNING: assumes without test that `gens(I)` is binomial
#
################################################################################

function tropical_variety_binomial(I::MPolyIdeal,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)
    ###
    # Construct matrix of exponent vector differences
    # and vector of coefficient valuation differences
    ###
    G = gens(I)
    A = matrix(ZZ,[first(expv)-last(expv) for expv in collect.(exponents.(G))])
    b = [ QQ(nu(last(coeff)/first(coeff))) for coeff in collect.(coefficients.(G))]

    ###
    # Compute tropical variety multiplicity
    ###
    snfAdiag = elementary_divisors(A)
    weight = abs(prod([m for m in snfAdiag if !iszero(m)]))

    ###
    # Constructing tropical variety set-theoretically
    ###
    A = matrix(QQ, A)
    can_solve, V, L = can_solve_with_solution_and_kernel(transpose(A), matrix(QQ,[b]); side=:left)
    @req can_solve "tropical variety cannot be empty"
    SigmaV = polyhedral_complex(IncidenceMatrix([[1]]), V, nothing, L)

    ###
    # Assemble tropical variety
    ###
    TropV = tropical_variety(SigmaV,[weight],convention(nu))
    if !weighted_polyhedral_complex_only
        set_attribute!(TropV,:algebraic_ideal,I)
        set_attribute!(TropV,:tropical_semiring_map,nu)
    end
    return TropV
end
