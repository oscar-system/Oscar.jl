const pbwAlg_multGrow = 5

struct MPolyRing{T} <: AbstractAlgebra.MPolyRing{T}
  N::Int
  coefficient_ring::Ring
  graded::Bool
  vars::Vector{Symbol}
end

mutable struct MPolyRingElem{T} <: AbstractAlgebra.MPolyRingElem{T}
  parent::MPolyRing{T}
  coeffs::Vector{Union{T,Nothing}}
  exps::Vector{Int}
  len::Int
end

function MPolyRingElem(R::MPolyRing{T}) where {T}
  return MPolyRingElem(
    R,
    Vector{Union{T,Nothing}}(),
    Vector{Int}(),
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

mutable struct PBWAlgebraElem{T} <: NCRingElem
  parent::PBWAlgebra{T}
  poly::MPolyRingElem{T}
end

function PBWAlgebraElem(A::PBWAlgebra{T}) where {T}
  return PBWAlgebraElem(A, zero(A.R))
end

struct PBWAlgebraHom{T,S,H<:Map} <:
       Map{PBWAlgebra{T},PBWAlgebra{S},Hecke.HeckeMap,PBWAlgebraHom{T,S,H}}
  domain::PBWAlgebra{T}
  codomain::PBWAlgebra{S}
  field_automorphism::H
  img::Vector{PBWAlgebraElem{T}}
end

struct QuantumField <: Field
  d::AbstractAlgebra.Generic.RationalFunctionField{QQFieldElem,QQPolyRingElem}
  q_factorial::Dict
end

mutable struct QuantumFieldElem <: FieldElem
  parent::QuantumField
  d::AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem,QQPolyRingElem}
end

struct _BarAutomorphism <: Map{QuantumField,QuantumField,Hecke.HeckeMap,_BarAutomorphism} end

mutable struct QuantumGroup <: NCRing
  algebra::PBWAlgebra{QuantumFieldElem}

  root_system::RootSystem
  bilinear_form::ZZMatrix
  w0::Vector{Int}

  cvx::Vector{Int} # permutation of roots (by index) to convex order
  qi::Vector{QuantumFieldElem}

  # cache
  canonical_basis::Dict{
    Vector{Int},PBWAlgebraElem{QuantumFieldElem}
  }

  bar_automorphism::PBWAlgebraHom{QuantumFieldElem,QuantumFieldElem,_BarAutomorphism}
end

mutable struct QuantumGroupElem <: NCRingElem
  parent::QuantumGroup
  elem::PBWAlgebraElem{QuantumFieldElem}
end

struct QuantumGroupHom <: Map{QuantumGroup,QuantumGroup,Hecke.HeckeMap,QuantumGroupHom}
  domain::QuantumGroup
  codomain::QuantumGroup
  hom::PBWAlgebraHom{QuantumFieldElem,QuantumFieldElem}
end
