###############################################################################
#
#   Constructors
#
###############################################################################

# We need to make this constructor internal for now,
# since QuaGroup does not allow arbitrary Cartan matrices
# and bilinear forms.

function _quantum_group(
  R::RootSystem;
  bilinear_form::ZZMatrix=Oscar.bilinear_form(R),
  w0::Vector{<:Integer}=word(longest_element(weyl_group(R))),
)
  QF, _ = quantum_field()
  return _quantum_group(QF, R; bilinear_form=bilinear_form, w0=w0)
end

raw"""
    quantum_group(R::RootSystem; bilinear_form::ZZMatrix=Oscar.bilinear_form(R),
      w0::Vector{<:Integer}=word(longest_element(weyl_group(R))),
    ) -> QuantumGroup
"""
function _quantum_group(
  QF::QuantumField,
  R::RootSystem;
  bilinear_form::ZZMatrix=Oscar.bilinear_form(R),
  w0::Vector{<:Integer}=word(longest_element(weyl_group(R))),
)
  q = gen(QF)
  vars = Symbol[]
  for i in 1:length(w0)
    push!(vars, Symbol("F$i"))
  end
  #=
  for i in 1:rank(R)
    push!(vars, Symbol("[K$i; 1]"))
    push!(vars, Symbol("K$i"))
  end
  for i in 1:length(w0)
    push!(vars, Symbol("E$i"))
  end
  =#

  P, theta = polynomial_ring(QF, vars)

  # for now we rely on the QuaGroup package to compute the PBW relations
  GAP.Packages.load("QuaGroup")
  gapR = GAP.Globals.Objectify(
    GAP.Globals.NewType(
      GAP.Globals.NewFamily(GAP.GapObj("RootSystemFam"), GAP.Globals.IsObject),
      GAP.evalstr("IsAttributeStoringRep and IsRootSystem"),
    ),
    GAP.evalstr("rec()"),
  )

  gapSim = GAP.GapObj(map(r -> GAP.GapObj(coefficients(r))[1], (simple_roots(R))))
  gapPos = GAP.GapObj(map(r -> GAP.GapObj(coefficients(r))[1], (positive_roots(R))))
  gapBil = GAP.GapObj(bilinear_form)

  GAP.Globals.SetPositiveRoots(gapR, gapPos)
  GAP.Globals.SetNegativeRoots(gapR, -gapPos)
  GAP.Globals.SetSimpleSystem(gapR, gapSim)
  GAP.Globals.SetCartanMatrix(gapR, GAP.Obj(transpose(cartan_matrix(R))))
  GAP.Globals.SetBilinearFormMat(gapR, gapBil)
  GAP.Globals.SetPositiveRootsNF(gapR, gapPos)
  GAP.Globals.SetSimpleSystemNF(gapR, gapSim)
  GAP.Globals.SetBilinearFormMatNF(gapR, gapBil)
  GAP.Globals.SetTypeOfRootSystem(
    gapR, collect(Iterators.flatmap(t -> GAP.GapObj.(t), root_system_type(R)))
  )
  GAP.Globals.SetLongestWeylWord(gapR, GAP.GapObj(GAP.GapObj.(w0)))

  gapU = GAP.Globals.QuantizedUEA(gapR)
  gapF = GAP.Globals.GeneratorsOfAlgebra(gapU)

  # set rels
  npos = number_of_positive_roots(R)
  rels = Matrix{elem_type(P)}(undef, length(vars), length(vars))

  term = one(P)
  for i in 1:length(vars), j in (i + 1):length(vars)
    rep = Oscar.GAPWrap.ExtRepOfObj(gapF[j] * gapF[i])
    rels[i, j] = zero(P)
    for n in 1:2:length(rep)
      for m in 1:2:length(rep[n])
        if !GAP.Globals.IsList(rep[n][m])
          if rep[n][m] <= npos
            for _ in 1:rep[n][m + 1]
              term = mul!(term, theta[rep[n][m]])
            end
          else
            for _ in 1:rep[n][m + 1]
              term = mul!(term, theta[rank(R) + rep[n][m]])
            end
          end
        else
          nk = 2 * (rep[n][m][1] - npos - 1) + npos + 1
          for _ in 1:rep[n][m][2]
            term = mul!(term, theta[nk + 1])
          end
          for _ in 1:rep[n][m + 1]
            term = mul!(term, theta[nk])
          end
        end
        term = mul!(term, inv(q_factorial(rep[n][m + 1], q)))
      end
      coeffRep = GAP.Globals.CoefficientsOfLaurentPolynomial(rep[n + 1])
      coeff = zero(QF)
      for n in 1:length(coeffRep[1])
        coeff = addmul!(coeff, q^(n + coeffRep[2] - 1), coeffRep[1][n])
      end
      rels[i, j] = addmul!(rels[i, j], term, coeff)
      term = one!(term)
    end
  end

  cvx = zeros(Int, npos)
  for i in 1:npos
    beta = w0[i]
    for j in (i - 1):-1:1
      beta = weyl_group(R).refl[w0[j], beta]
    end
    cvx[beta] = i
  end

  alg = pbw_algebra(P, rels)
  return QuantumGroup(
    alg,
    R,
    bilinear_form,
    w0,
    cvx,
    [q^div(r * gapBil * r, 2) for r in GAP.Globals.PositiveRootsInConvexOrder(gapR)],
    Dict{Vector{Int},PBWAlgebraElem}(),
    _bar_automorphism(alg, R, cvx),
  )
end

@doc raw"""
    quantum_group(fam::Symbol, rk::Int; w0::Union{Vector{<:Integer},Nothing}=nothing) -> QuantumGroup
    
Return a quantum group for a finite root system with family `fam` with rank `rk`.
Optionally, a reduced decomposition `w0` can be provided to choose a convex order for the positive roots.
This choice determines the PBW basis of the quantum group.
"""
function quantum_group(fam::Symbol, rk::Int; w0::Union{Vector{<:Integer},Nothing}=nothing)
  QF, _ = quantum_field()
  return quantum_group(QF, fam, rk; w0=w0)
end

@doc raw"""
    quantum_group(QF::QuantumField, fam::Symbol, rk::Int; w0::Union{Vector{<:Integer},Nothing}=nothing) -> QuantumGroup
"""
function quantum_group(
  QF::QuantumField, fam::Symbol, rk::Int; w0::Union{Vector{<:Integer},Nothing}=nothing
)
  R = root_system(fam, rk)
  if !is_finite(weyl_group(R))
    error("Quantum groups are currently only supported for finite root systems.")
  end

  if isnothing(w0)
    w0 = word(longest_element(weyl_group(R)))
  end
  return _quantum_group(QF, R; w0=w0)
end

###############################################################################
#
#   Accessors
#
###############################################################################

function coefficient_ring(U::QuantumGroup)
  return coefficient_ring(U.algebra)
end

@doc raw"""
    root_system(U::QuantumGroup) -> RootSystem
    
Return the root system of the quantum group `U`.
"""
function root_system(U::QuantumGroup)
  return U.root_system
end

function parent(x::QuantumGroupElem)
  return x.parent
end

###############################################################################
#
#   Type system
#
###############################################################################

function elem_type(::Type{QuantumGroup})
  return QuantumGroupElem
end

function parent_type(::Type{QuantumGroupElem})
  return QuantumGroup
end

function is_exact_type(::Type{QuantumGroupElem})
  return is_exact_type(QuantumFieldElem)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function ngens(U::QuantumGroup)
  return ngens(U.algebra)
end

function gen(U::QuantumGroup, i::Int)
  return QuantumGroupElem(U, gen(U.algebra, i))
end

function gens(U::QuantumGroup)
  return [gen(U, i) for i in 1:ngens(U)]
end

function is_gen(x::QuantumGroupElem)
  return is_gen(x.elem)
end

function is_gen_with_index(x::QuantumGroupElem)
  return is_gen_with_index(x.elem)
end

@doc raw"""
    negative_chevalley_gens(U::QuantumGroup) -> Vector{QuantumGroupElem}
"""
function negative_chevalley_gens(U::QuantumGroup)
  return [gen(U, U.cvx[i]) for i in 1:rank(root_system(U))]
end

function one(U::QuantumGroup)
  return QuantumGroupElem(U, one(U.algebra))
end

function zero(U::QuantumGroup)
  return QuantumGroupElem(U, zero(U.algebra))
end

function coeff(x::QuantumGroupElem, i::Int)
  exp = exponent_vector(x.elem, i)
  cf = deepcopy(coeff(x.elem, i))
  for j in 1:length(exp)
    cf = mul!(cf, q_factorial(exp[j], parent(x).qi[j]))
  end
  return cf
end

function length(x::QuantumGroupElem)
  return length(x.elem)
end

function exponent_vector(x::QuantumGroupElem, i::Int)
  return exponent_vector(x.elem, i)
end

function set_exponent_vector!(x::QuantumGroupElem, i::Int, exp::Vector{Int})
  set_exponent_vector!(x.elem, i, exp)
  return x
end

###############################################################################
#
#   Copying / Hashing
#
###############################################################################

function deepcopy_internal(x::QuantumGroupElem, dict::IdDict)
  return get!(dict, x) do
    QuantumGroupElem(
      x.parent,
      deepcopy_internal(x.elem, dict),
    )
  end
end

function Base.hash(x::QuantumGroupElem, h::UInt)
  b = 0xe2cb30215ca391a1 % UInt
  h = hash(parent(x), h)
  h = hash(x.elem, h)

  return xor(h, b)
end

###############################################################################
#
#   Printing
#
###############################################################################

function Base.show(io::IO, U::QuantumGroup)
  print(io, "Quantum Group over ", root_system(U))
end

#function Base.show(io::IO, x::QuantumGroupElem)
#  show(io, x.elem)
#end

function expressify(x::QuantumGroupElem; context=nothing)
  expr = Expr(:call, :+)

  U = parent(x)
  exp = Vector{Int}(undef, ngens(U))
  for i in 1:length(x)
    exponent_vector!(exp, x.elem, i)
    expr2 = Expr(:call, :*, expressify(coeff(x, i); context=context))
    for j in 1:length(exp)
      if exp[j] != 0
        push!(expr2.args, string("F[$j]", exp[j] > 1 ? "^($(exp[j]))" : ""))
      end
    end
    push!(expr.args, expr2)
  end
  return expr
end

@enable_all_show_via_expressify QuantumGroupElem

###############################################################################
#
#   Arithmetic
#
###############################################################################

function add!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = add!(z.elem, x.elem, y.elem)
  return z
end

function div!(z::QuantumGroupElem, x::QuantumGroupElem, a::QuantumFieldElem)
  z.elem = div!(z.elem, x.elem, a)
  return z
end

function mul!(z::QuantumGroupElem, x::QuantumGroupElem, a::QuantumFieldElem)
  z.elem = mul!(z.elem, x.elem, a)
  return z
end

function mul!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = mul!(z.elem, x.elem, y.elem)
  return z
end

function neg!(z::QuantumGroupElem, x::QuantumGroupElem)
  z.elem = neg!(z.elem, x.elem)
  return z
end

function sub!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = sub!(z.elem, x.elem, y.elem)
  return z
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function Base.:+(x::QuantumGroupElem, y::QuantumGroupElem)
  check_parent(x, y)
  return add!(zero(x), x, y)
end

function Base.:-(x::QuantumGroupElem, y::QuantumGroupElem)
  check_parent(x, y)
  return sub!(zero(x), x, y)
end

function Base.:-(x::QuantumGroupElem)
  return neg!(zero(x), x)
end

function Base.:^(x::QuantumGroupElem, n::Int)
  return QuantumGroupElem(parent(x), x.elem^n)
end

function Base.:*(x::QuantumGroupElem, a::QuantumFieldElem)
  return mul!(zero(x), x, a)
end

function Base.:*(a::QuantumFieldElem, x::QuantumGroupElem)
  return mul!(zero(x), x, a)
end

function Base.:*(x::QuantumGroupElem, y::QuantumGroupElem)
  check_parent(x, y)
  return mul!(zero(x), x, y)
end

function Base.:/(x::QuantumGroupElem, a::QuantumFieldElem)
  return divexact_right(x, a)
end

function divexact_right(x::QuantumGroupElem, a::QuantumFieldElem)
  return QuantumGroupElem(parent(x), divexact_right(x.elem, a))
end

function divexact_right(x::QuantumGroupElem, y::QuantumGroupElem)
  check_parent(x, y)
  return QuantumGroupElem(parent(x), divexact(x.elem, y.elem))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(x::QuantumGroupElem, y::QuantumGroupElem)
  check_parent(x, y)
  return x.elem == y.elem
end

function isone(x::QuantumGroupElem)
  return isone(x.elem)
end

function iszero(x::QuantumGroupElem)
  return iszero(x.elem)
end

###############################################################################
#
#   Monomials
#
###############################################################################

@doc raw"""
    monomial(F::Vector{QuantumGroupElem}, m::Vector{Tuple{Int,Int}}) -> QuantumGroupElem
    
For a vector `F` of generators of a quantum group `U` return the monomial
for the word `i` with divided powers `n`.
"""
function monomial(U::QuantumGroup, i::Vector{Int}, n::Vector{Int})
  f = one(U)
  for (i, n) in Iterators.zip(i, n)
    f = mul!(f, gen(U, U.cvx[i])^n)
    f = div!(f, q_factorial(n, U.qi[U.cvx[i]]))
  end

  return f
end

function pbw_monomial(U, n::Vector{Int})
  @req ngens(U) == length(n) "length must match number of generators"

  f = one(U)
  set_exponent_vector!(f.elem, 1, n)
  for i in 1:length(n)
    f = div!(f, q_factorial(n[i], U.qi[i]))
  end

  return f
end

###############################################################################
#
#   Conformance 
#
###############################################################################

function ConformanceTests.generate_element(U::QuantumGroup)
  return QuantumGroupElem(U, ConformanceTests.generate_element(U.algebra))
end
