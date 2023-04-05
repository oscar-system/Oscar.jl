################################################################################
# Lower case constructors
################################################################################

projective_scheme(S::MPolyDecRing) = ProjectiveScheme(S)

projective_scheme(S::MPolyDecRing, I::MPolyIdeal{T}) where {T<:MPolyDecRingElem} = ProjectiveScheme(S, I)

projective_scheme(I::MPolyIdeal{<:MPolyDecRingElem}) = ProjectiveScheme(base_ring(I), I)

projective_scheme(Q::MPolyQuoRing{MPolyDecRingElem{T, PT}}) where {T, PT<:MPolyRingElem{T}} = ProjectiveScheme(Q)

################################################################################
# Subschemes
################################################################################

function subscheme(P::AbsProjectiveScheme, f::RingElem)
  S = homogeneous_coordinate_ring(P)
  parent(f) === S || return subscheme(P, S(f))
  Q, _ = quo(S, ideal(S, [f]))
  result = projective_scheme(Q)
  if isdefined(P, :Y) 
    set_base_scheme!(result, base_scheme(P))
  end
  return result
end

function subscheme(
    P::AbsProjectiveScheme, 
    f::Vector{T}
  ) where {T<:RingElem}
  length(f) == 0 && return P #TODO: Replace P by an honest copy!
  S = homogeneous_coordinate_ring(P)
  for i in 1:length(f)
    parent(f[i]) === S || return subscheme(P, S.(f))
  end
  Q, _ = quo(S, ideal(S, f))
  result = projective_scheme(Q)
  if isdefined(P, :Y) 
    set_base_scheme!(result, base_scheme(P))
  end
  return result
end

function subscheme(P::AbsProjectiveScheme,
    I::Ideal{T}
  ) where {T<:RingElem}
  S = homogeneous_coordinate_ring(P)
  base_ring(I) === S || error("ideal does not belong to the correct ring")
  Q, _ = quo(S, I)
  result = projective_scheme(Q)
  if isdefined(P, :Y) 
    set_base_scheme!(result, base_scheme(P))
  end
  return result
end

################################################################################
# Projective space
################################################################################

function projective_space(
    A::CoeffRingType, 
    var_symb::Vector{Symbol}
  ) where {CoeffRingType<:Ring}
  n = length(var_symb)
  R, _ = polynomial_ring(A, var_symb)
  S, _ = grade(R, [1 for i in 1:n ])
  return projective_scheme(S)
end

projective_space(
                 A::CoeffRingType, 
                 var_names::Vector{String}
                ) where {CoeffRingType<:Ring} = projective_space(A, Symbol.(var_names))


function projective_space(
    A::CoeffRingType, 
    r::Int; 
    var_name::String="s"
  ) where {CoeffRingType<:Ring}
  R, _ = polynomial_ring(A, [var_name*"$i" for i in 0:r])
  S, _ = grade(R, [1 for i in 0:r ])
  return projective_scheme(S)
end

function projective_space(
    W::Union{<:SpecOpen, <:AbsSpec}, 
    r::Int; 
    var_name::String="s"
  ) 
  P = projective_space(OO(W), r, var_name=var_name)
  set_base_scheme!(P, W)
  return P
end

function projective_space(
    W::Union{<:SpecOpen, <:AbsSpec}, 
    var_names::Vector{Symbol}
  ) 
  P = projective_space(OO(W), var_names)
  set_base_scheme!(P, W)
  return P
end

function projective_space(
    W::Union{<:SpecOpen, <:AbsSpec}, 
    var_names::Vector{String}
  ) 
  P = projective_space(OO(W), var_names)
  set_base_scheme!(P, W)
  return P
end

################################################################################
# reduced scheme
################################################################################

function reduced_scheme(X::ProjectiveScheme)
  I = defining_ideal(X)
  Irad = radical(I)
  Xred = subscheme(ambient_space(X), I)
  set_attribute!(Xred, :is_reduced=>true)
  return Xred
end
