# conversion between Singular and AbstractAlgebra function fields
function Oscar.singular_ring(F::AbstractAlgebra.Generic.RationalFunctionField{T}) where T <: FieldElem
  return Singular.FunctionField(singular_ring(base_ring(F)), [string(F.S)])[1]
end

function (F::Singular.N_FField)(x::AbstractAlgebra.Generic.Rat{T}) where T <: FieldElem
  @assert Singular.transcendence_degree(F) == 1
  a = Singular.transcendence_basis(F)[1]
  numerator(x)(a) // denominator(x)(a)
end

function (F::AbstractAlgebra.Generic.RationalFunctionField{T})(x::Singular.n_transExt) where T <: FieldElem
  @assert Singular.transcendence_degree(parent(x)) == 1
  n, d = Singular.n_transExt_to_spoly.([numerator(x), denominator(x)])
  F(n) // F(d)
end

function (F::FlintRationalField)(x::AbstractAlgebra.Generic.Rat{fmpq})
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError) end
  F(cst_num) // F(cst_denom)
end

function (F::Singular.N_FField)(x::gfp_elem)
  @assert characteristic(F) == characteristic(parent(x))
  F(lift(x))
end

# These should probably be put in Singular.jl
function (F::Singular.N_FField)(x::fmpq)
  @assert characteristic(F) == 0
  F(numerator(x)) // F(denominator(x))
end

function (Ox::PolyRing)(f::Singular.spoly)
  O = base_ring(Ox)
  @assert ngens(parent(f)) == 1
  Ox(O.(coefficients_of_univariate(f)))
end
