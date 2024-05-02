function check_char(A::Ring, B::Ring)
  @req characteristic(A) == characteristic(B) "wrong characteristic"
end

################################################################################
#
#  Some ad hoc promote rules
#
################################################################################

AbstractAlgebra.promote_rule(::Type{AbsSimpleNumFieldElem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{QQFieldElem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{fpFieldElem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{FpFieldElem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{FqFieldElem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{FqPolyRepFieldElem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{fqPolyRepFieldElem}, ::Type{Singular.n_transExt}) = Singular.n_transExt
