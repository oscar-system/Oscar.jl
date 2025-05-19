###############################################################################
#
#   Constructors
#
###############################################################################

@doc raw"""
    quantum_group(R::RootSystem; bilinear_form::ZZMatrix=Oscar.bilinear_form(R),
      w0::Vector{<:Integer}=word(longest_element(weyl_group(R))),
    ) -> QuantumGroup
"""
function quantum_group(
  R::RootSystem;
  bilinear_form::ZZMatrix=Oscar.bilinear_form(R),
  w0::Vector{<:Integer}=word(longest_element(weyl_group(R))),
)
  QF, q = quantum_field()
  P, theta = polynomial_ring(QF, :F => 1:length(w0))

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
  rels = Matrix{elem_type(P)}(undef, npos, npos) # zero_matrix(P, npos, npos)

  term = one(P)
  for i in 1:npos, j in (i + 1):npos
    rep = Oscar.GAPWrap.ExtRepOfObj(gapF[j] * gapF[i])
    rels[i, j] = zero(P)
    for n in 1:2:length(rep)
      for m in 1:2:length(rep[n])
        for _ in 1:rep[n][m + 1]
          term = mul!(term, theta[rep[n][m]])
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

  alg = pbw_algebra(P, rels) # lex(theta)
  return QuantumGroup(
    alg,
    R,
    bilinear_form,
    w0,
    cvx,
    # [q^div(r * gapBil * r, 2) for r in GAP.Globals.PositiveRootsInConvexOrder(gapR)],
    Dict{Vector{Int},PBWAlgebraElem}(),
    #_bar_involution(),
  )
end

@doc raw"""
    quantum_group(fam::Symbol, rk::Int) -> QuantumGroup
"""
function quantum_group(fam::Symbol, rk::Int)
  return quantum_group(root_system(fam, rk))
end

###############################################################################
#
#   Accessors
#
###############################################################################

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

function ngens(U::QuantumGroup)
  return ngens(U.algebra)
end

function gen(U::QuantumGroup, i::Int)
  return QuantumGroupElem(U, gen(U.algebra, i))
end

function gens(U::QuantumGroup)
  return [gen(U, i) for i in 1:ngens(U)]
end

function chevalley_gens(U::QuantumGroup)
  return [gen(U, U.cvx[i]) for i in 1:rank(root_system(U))]
end

function one(U::QuantumGroup)
  return QuantumGroupElem(U, one(U.algebra))
end

function zero(U::QuantumGroup)
  return QuantumGroupElem(U, zero(U.algebra))
end

###############################################################################
#
#   
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

function Base.show(io::IO, x::QuantumGroupElem)
  show(io, x.elem)
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function add!(z::QuantumGroupElem, x::QuantumGroupElem, y::QuantumGroupElem)
  z.elem = add!(z.elem, x.elem, y.elem)
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
