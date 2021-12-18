###
# Valuations over exact fields for tropical geometry
# ==================================================
###

struct ValuationMap{typeofValuedField,typeofUniformizer}
  valued_field::typeofValuedField
  uniformizer::typeofUniformizer
  residue_field
  residue_map
end
export ValuationMap


###
# p-adic valuation on QQ
###

# Constructor:
function ValuationMap(Q::FlintRationalField,p::fmpz)
    residue_map(c) = FiniteField(p)(c)
    return ValuationMap{typeof(Q),typeof(p)}(Q,p,FiniteField(p),residue_map)
end
ValuationMap(Q::FlintRationalField,p) = ValuationMap(Q,ZZ(p)) # for other types of `p` such as `Int`

# Evaluation:
(val::ValuationMap{FlintRationalField,fmpz})(c) = valuation(c, val.uniformizer)



###
# Laurent valuation on K(t)
###

# Constructor:
function t_adic_valuation(t::AbstractAlgebra.Generic.Rat,c)
    num = numerator(c)
    nom = denominator(c)
    valnum = first(i for i in 0:degree(num) if !iszero(coeff(num, i)))
    valnom = first(i for i in 0:degree(nom) if !iszero(coeff(nom, i)))
    return valnum-valnom
end

function ValuationMap(Kt::AbstractAlgebra.Generic.RationalFunctionField,t::AbstractAlgebra.Generic.Rat)
    function residue_map(c)
        valc = t_adic_valuation(t,c)
        if (valc<0)
            error("residue_map: input has negative valuation, not in valuation ring")
        end
        return base_ring(Kt)(evaluate(c,0))
    end
    return ValuationMap{typeof(Kt),typeof(t)}(Kt,t,base_ring(Kt),residue_map)
end

# Evaluation:
(val::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K},AbstractAlgebra.Generic.Rat{K}} where {K})(c) = t_adic_valuation(val.uniformizer,c)
