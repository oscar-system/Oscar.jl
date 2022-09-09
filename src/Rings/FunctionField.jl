function check_char(A::Ring, B::Ring)
  characteristic(A) == characteristic(B) || throw(ArgumentError("wrong characteristic"))
end

# conversion between FractionField and Singular function field
# univariate
function Oscar.singular_coeff_ring(F::AbstractAlgebra.Generic.FracField{T}) where T <: Union{fmpq_poly, gfp_poly}
  R = base_ring(F)
  return Singular.FunctionField(singular_coeff_ring(base_ring(R)), [string(R.S)])[1]
end

function (F::Singular.N_FField)(x::AbstractAlgebra.Generic.Frac{T}) where T <: Union{fmpq_poly, gfp_poly}
  check_char(F, parent(x))
  Singular.transcendence_degree(F) == 1 || throw(ArgumentError("wrong number of generators"))
  a = Singular.transcendence_basis(F)[1]
  numerator(x)(a) // denominator(x)(a)
end

function (F::AbstractAlgebra.Generic.FracField{T})(x::Singular.n_transExt) where T <: Union{fmpq_poly, gfp_poly}
  check_char(F, parent(x))
  Singular.transcendence_degree(parent(x)) == 1 || throw(ArgumentError("wrong number of generators"))
  n, d = Singular.n_transExt_to_spoly.([numerator(x), denominator(x)])
  F(n) // F(d)
end

# multivariate
function Oscar.singular_coeff_ring(F::AbstractAlgebra.Generic.FracField{T}) where T <: Union{fmpq_mpoly, gfp_mpoly}
  R = base_ring(F)
  return Singular.FunctionField(singular_coeff_ring(base_ring(R)), string.(R.S))[1]
end

function (F::Singular.N_FField)(x::AbstractAlgebra.Generic.Frac{T}) where T <: Union{fmpq_mpoly, gfp_mpoly}
  check_char(F, parent(x))
  Singular.transcendence_degree(F) == ngens(base_ring(parent(x))) || throw(ArgumentError("wrong number of generators"))
  a = Singular.transcendence_basis(F)
  numerator(x)(a...) // denominator(x)(a...)
end

function (F::AbstractAlgebra.Generic.FracField{T})(x::Singular.n_transExt) where T <: Union{fmpq_mpoly, gfp_mpoly}
  check_char(F, parent(x))
  ngens(base_ring(F)) == Singular.transcendence_degree(parent(x)) || throw(ArgumentError("wrong number of generators"))
  n, d = Singular.n_transExt_to_spoly.([numerator(x), denominator(x)])
  F(n) // F(d)
end

# conversion between Singular and AbstractAlgebra function fields
function Oscar.singular_coeff_ring(F::AbstractAlgebra.Generic.RationalFunctionField{T}) where T <: FieldElem
  return singular_coeff_ring(F.fraction_field) # TODO: this is not correct
end

function (F::Singular.N_FField)(x::AbstractAlgebra.Generic.Rat{T}) where T <: FieldElem
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

# coercion
function (F::FlintRationalField)(x::AbstractAlgebra.Generic.Frac{T}) where T <: Union{fmpq_poly, fmpq_mpoly}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:fmpq, fmpq, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::FlintRationalField)(x::AbstractAlgebra.Generic.Rat{fmpq})
  return F(x.d)
end

function (F::Nemo.GaloisField)(x::AbstractAlgebra.Generic.Frac{T}) where T <: Union{gfp_poly, gfp_mpoly}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:gfp_elem, gfp_elem, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::Nemo.GaloisField)(x::AbstractAlgebra.Generic.Rat{gfp_elem})
  return F(x.d)
end

################################################################################
#
#  Some ad hoc promote rules
#
################################################################################

AbstractAlgebra.promote_rule(::Type{nf_elem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{fmpq}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{gfp_elem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{gfp_fmpz_elem}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{fq}, ::Type{Singular.n_transExt}) = Singular.n_transExt

AbstractAlgebra.promote_rule(::Type{fq_nmod}, ::Type{Singular.n_transExt}) = Singular.n_transExt
