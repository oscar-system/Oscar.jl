
function subscheme(P::AbsProjectiveScheme, f::RingElemType) where {RingElemType<:MPolyDecRingElem}
  S = ambient_coordinate_ring(P)
  parent(f) == S || error("ring element does not belong to the correct ring")
  Q = ProjectiveScheme(S, ideal(S, vcat(gens(defining_ideal(P)), [f])))
  if isdefined(P, :Y) 
    set_base_scheme!(Q, base_scheme(P))
  end
  return Q
end

function subscheme(
    P::AbsProjectiveScheme, 
    f::Vector{RingElemType}
  ) where {RingElemType<:MPolyDecRingElem}
  length(f) == 0 && return P #TODO: Replace P by an honest copy!
  S = ambient_coordinate_ring(P)
  for i in 1:length(f)
    parent(f[i]) == S || error("ring element does not belong to the correct ring")
  end
  Q = ProjectiveScheme(S, ideal(S, vcat(gens(defining_ideal(P)),f)))
  if isdefined(P, :Y) 
    set_base_scheme!(Q, base_scheme(P))
  end
  return Q
end

function subscheme(P::AbsProjectiveScheme, 
    I::MPolyIdeal{T}
  ) where {T<:RingElem}
  S = ambient_coordinate_ring(P)
  base_ring(I) == S || error("ideal does not belong to the correct ring")
  Q = ProjectiveScheme(S, ideal(S, vcat(gens(I), gens(defining_ideal(P)))))
  if isdefined(P, :Y) 
    set_base_scheme!(Q, base_scheme(P))
  end
  return Q
end

function projective_space(
    A::CoeffRingType, 
    var_symb::Vector{Symbol}
  ) where {CoeffRingType<:Ring}
  n = length(var_symb)
  R, _ = polynomial_ring(A, var_symb)
  S, _ = grade(R, [1 for i in 1:n ])
  I = ideal(S, [zero(S)])
  return ProjectiveScheme(S, I)
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
  I = ideal(S, [zero(S)])
  return ProjectiveScheme(S, I)
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

