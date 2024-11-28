# pushforward of algebraic cycles for proper maps

@doc raw"""
    pushforward(p::BlowupMorphism, W::AbsAlgebraicCycle)

Given a `BlowupMorphism` ``p : X → Y`` and an `AbsAlgebraicCycle` ``W ∈ Z_k(X)``,
compute ``p_*(W) ∈ Z_k(Y)``.
"""
function pushforward(p::BlowupMorphism, W::AbsAlgebraicCycle)
  error("not implemented")
end

include("fixes.jl")
include("erroxefunctions.jl")

