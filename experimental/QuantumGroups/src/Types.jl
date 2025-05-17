const pbwAlg_multGrow = 5

struct MPolyRing{T} <: AbstractAlgebra.MPolyRing{T}
  N::Int
  coefficient_ring::Ring
  graded::Bool
  vars::Vector{Symbol}
end

mutable struct MPolyRingElem{T} <: AbstractAlgebra.MPolyRingElem{T}
  parent::MPolyRing{T}
  coeffs::Memory{T}
  exps::Memory{Int}
  len::Int
end

function MPolyRingElem(R::MPolyRing{T}) where {T}
  return MPolyRingElem(
    R,
    #UInt32(8),
    Memory{T}(),
    Memory{Int}(),
    0,
  )
end

function MPolyRingElem(R::MPolyRing{T}, coeff::T, exp::Vector{Int}) where {T}
  x = MPolyRingElem(R)
  fit!(x, 1)
  x.coeffs[1] = coeff
  copyto!(x.exps, exp)
  x.len = 1

  return x
end

struct PBWAlgebra{T} <: NCRing
  R::MPolyRing{T}
  mult::Vector{Matrix{MPolyRingElem{T}}}
end

function PBWAlgebra(
  R::MPolyRing{T}, rels::Vector{MPolyRingElem{T}}
) where {T<:FieldElem}
  mult = Vector{Matrix{MPolyRingElem{T}}}(undef, length(rels))
  for i in 1:length(rels)
    if length(rels[i]) == 1 # quasi-commuative case
      mult[i] = Matrix{MPolyRingElem{T}}(undef, 1, 1)
    else
      mult[i] = Matrix{MPolyRingElem{T}}(undef, pbwAlg_multGrow, pbwAlg_multGrow)
    end
    mult[i][1, 1] = rels[i]
  end
  return PBWAlgebra(R, mult)
end

mutable struct PBWAlgebraElem{T} # <: NCRingElem
  parent::PBWAlgebra{T}
  poly::MPolyRingElem{T}
end

function PBWAlgebraElem(A::PBWAlgebra{T}) where {T}
  return PBWAlgebraElem(A, zero(A.R))
end

struct PBWAlgebraHom{T}
  domain::PBWAlgebra{T}
  codomain::PBWAlgebra{T}
  img::Vector{PBWAlgebraElem{T}}
end

const QuantumFieldElem = FracFieldElem{
  LaurentPolyWrap{
    ZZRingElem,
    ZZPolyRingElem,
    AbstractAlgebra.Generic.LaurentPolyWrapRing{ZZRingElem,ZZPolyRing},
  },
}

mutable struct QuantumGroup
  algebra::PBWAlgebra{QuantumFieldElem}

  root_system::RootSystem
  bilinear_form::ZZMatrix
  w0::Vector{Int}

  cvx::Vector{Int} # permutation of roots (by index) to convex order

  # cache
  canonical_basis::Dict{
    Vector{Int},PBWAlgebraElem{QuantumFieldElem}
  }

  bar_involution::PBWAlgebraHom{QuantumFieldElem}
end

mutable struct QuantumGroupElem <: NCRingElem
  parent::QuantumGroup
  elem::PBWAlgebraElem{FracFieldElem{QuantumFieldElem}}
end

struct QuantumGroupHom
  domain::QuantumGroup
  codomain::QuantumGroup
  img::Vector{QuantumGroupElem}
end