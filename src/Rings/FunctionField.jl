function check_char(A::Ring, B::Ring)
  @req characteristic(A) == characteristic(B) "wrong characteristic"
end

# conversion between fraction_field and Singular function field
# univariate
function Oscar.singular_coeff_ring(F::AbstractAlgebra.Generic.FracField{T}) where T <: Union{QQPolyRingElem, fpPolyRingElem, FqPolyRingElem}
  R = base_ring(F)
  return Singular.FunctionField(singular_coeff_ring(base_ring(R)), [string(R.S)])[1]
end

function (F::Singular.N_FField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{QQPolyRingElem, fpPolyRingElem, FqPolyRingElem}
  check_char(F, parent(x))
  @req Singular.transcendence_degree(F) == 1 "wrong number of generators"
  a = Singular.transcendence_basis(F)[1]
  numerator(x)(a) // denominator(x)(a)
end

function (F::AbstractAlgebra.Generic.FracField{T})(x::Singular.n_transExt) where T <: Union{QQPolyRingElem, fpPolyRingElem, FqPolyRingElem}
  check_char(F, parent(x))
  @req Singular.transcendence_degree(parent(x)) == 1 "wrong number of generators"
  n, d = Singular.n_transExt_to_spoly.([numerator(x), denominator(x)])
  F(n) // F(d)
end

# multivariate
function Oscar.singular_coeff_ring(F::AbstractAlgebra.Generic.FracField{T}) where T <: Union{QQMPolyRingElem, fpMPolyRingElem, FqMPolyRingElem}
  R = base_ring(F)
  return Singular.FunctionField(singular_coeff_ring(base_ring(R)), _variables_for_singular(symbols(R)))[1]
end

function (F::Singular.N_FField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{QQMPolyRingElem, fpMPolyRingElem, FqMPolyRingElem}
  check_char(F, parent(x))
  @req Singular.transcendence_degree(F) == ngens(base_ring(parent(x))) "wrong number of generators"
  a = Singular.transcendence_basis(F)
  numerator(x)(a...) // denominator(x)(a...)
end

function (F::AbstractAlgebra.Generic.FracField{T})(x::Singular.n_transExt) where T <: Union{QQMPolyRingElem, fpMPolyRingElem, FqMPolyRingElem}
  check_char(F, parent(x))
  @req ngens(base_ring(F)) == Singular.transcendence_degree(parent(x)) "wrong number of generators"
  n, d = Singular.n_transExt_to_spoly.([numerator(x), denominator(x)])
  F(n) // F(d)
end

# conversion between Singular and AbstractAlgebra function fields
function Oscar.singular_coeff_ring(F::AbstractAlgebra.Generic.RationalFunctionField{T}) where T <: FieldElem
  return singular_coeff_ring(F.fraction_field) # TODO: this is not correct
end

function (F::Singular.N_FField)(x::AbstractAlgebra.Generic.RationalFunctionFieldElem{T}) where T <: FieldElem
  return F(x.d)
end

function (F::AbstractAlgebra.Generic.RationalFunctionField{T})(x::Singular.n_transExt) where T <: FieldElem
  return F(F.fraction_field(x))
end

function (Ox::PolyRing)(f::Singular.spoly)
  O = base_ring(Ox)
  @assert ngens(parent(f)) == 1
  Ox(O.(coefficients_of_univariate(f)))
end


################################################################################
#
#  Some ad hoc promote rules
#
################################################################################

AbstractAlgebra.promote_rule(::Type{AbsSimpleNumFieldElem}, ::Type{C}) where {C <: Singular.n_transExt} = C

AbstractAlgebra.promote_rule(::Type{fpFieldElem}, ::Type{C}) where {C <: Singular.n_transExt} = C

AbstractAlgebra.promote_rule(::Type{FpFieldElem}, ::Type{C}) where {C <: Singular.n_transExt} = C

AbstractAlgebra.promote_rule(::Type{FqFieldElem}, ::Type{C}) where {C <: Singular.n_transExt} = C

AbstractAlgebra.promote_rule(::Type{FqPolyRepFieldElem}, ::Type{C}) where {C <: Singular.n_transExt} = C

AbstractAlgebra.promote_rule(::Type{fqPolyRepFieldElem}, ::Type{C}) where {C <: Singular.n_transExt} = C
