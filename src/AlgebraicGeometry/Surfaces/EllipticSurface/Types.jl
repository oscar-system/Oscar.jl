@doc raw"""
    EllipticSurface{BaseField<:Field, BaseCurveFieldType} <: AbsCoveredScheme{BaseField}

A relatively minimal elliptic surface defined as follows.
    
A genus ``1``-fibration is a proper map
```math
\pi \colon X \to C
```
from a smooth projective surface ``X`` to a smooth projective curve ``C`` whose generic fiber is a curve of (arithmetic) genus ``1``.

The fibration is relatively minimal if its fibers do not contain any ``(-1)``-curves.
We call the fibration elliptic if it is relatively minimal and comes equipped with a section ``\sigma_0\colon \mathbb{P}^1 \to X``.
This turns the generic fiber of ``\pi`` into an elliptic curve ``E/k(C)`` where
``k(C)`` is the function field of the curve ``C``.
  
We further require that ``\pi`` has at least one singular fiber and that the base field ``k`` is perfect. 

For now functionality is restricted to ``C = \mathbb{P}^1``.

This datatype stores a subgroup of the Mordell-Weil group ``E(k(C))``.
It is referred to as the  Mordell-Weil subgroup of `X`.
"""
@attributes mutable struct EllipticSurface{BaseField<:Field, BaseCurveFieldType} <: AbsCoveredSurface{BaseField}
  Y::CoveredScheme{BaseField}  # the underlying_scheme
  E::EllipticCurve{BaseCurveFieldType}
  MWL::Vector{EllipticCurvePoint{BaseCurveFieldType}} # basis for the mordell weil group
  MWLtors::Vector{EllipticCurvePoint{BaseCurveFieldType}} # torsion sections
  Weierstrasschart::AbsAffineScheme
  Weierstrassmodel::CoveredScheme
  inc_Weierstrass::CoveredClosedEmbedding # inclusion of the weierstrass chart in its ambient projective bundle
  inc_Y::CoveredClosedEmbedding # inclusion of Y in its ambient blown up projective bundle
  euler_characteristic::Int
  resolution_strategy::Symbol
  # the following are temporary until we have a dedicated type for
  # iterated blow ups
  blowup::AbsCoveredSchemeMorphism
  blowups::Vector{<:AbsCoveredSchemeMorphism}
  # exceptionals not used for now
  ambient_blowups::Vector{<:BlowupMorphism}
  ambient_exceptionals::Vector{<:EffectiveCartierDivisor}
  fibration::AbsCoveredSchemeMorphism # the projection to IP^1
  fibration_weierstrass_model::AbsCoveredSchemeMorphism # the projection from the Weierstrass model

  function EllipticSurface(
      generic_fiber::EllipticCurve{F}, 
      euler_characteristic::Int, 
      mwl_basis::Vector{<:EllipticCurvePoint};
      resolution_strategy::Symbol=:iterative
    ) where F
    B = typeof(coefficient_ring(base_ring(base_field(generic_fiber))))
    S = new{B,F}()
    S.E = generic_fiber
    S.MWL = mwl_basis
    S.euler_characteristic = euler_characteristic
    set_attribute!(S, :is_irreducible=>true)
    set_attribute!(S, :is_reduced=>true)
    set_attribute!(S, :is_integral=>true)
    set_attribute!(S, :is_equidimensional=>true)
    S.resolution_strategy = resolution_strategy
    return S
  end
end

@doc raw"""
    EllipticSurfaceSection <: AbsWeilDivisor
    
A section of an elliptic fibration represented as a Weil divisor.
"""
@attributes mutable struct EllipticSurfaceSection{
    CoveredSchemeType<:AbsCoveredScheme, 
    CoefficientRingType<:AbstractAlgebra.Ring, 
    CoefficientRingElemType<:AbstractAlgebra.RingElem
   } <: AbsWeilDivisor{CoveredSchemeType, CoefficientRingType}
  D::WeilDivisor{CoveredSchemeType, CoefficientRingType, CoefficientRingElemType}
  P::EllipticCurvePoint

  function EllipticSurfaceSection(X::EllipticSurface, P::EllipticCurvePoint; coefficient_ring::Ring=ZZ)
    @vprint :EllipticSurface 3 "Computing a section from a point on the generic fiber\n"
    weierstrass_contraction(X) # trigger required computations
    PX = _section_on_weierstrass_ambient_space(X, P)
    for f in X.ambient_blowups
      PX = strict_transform(f , PX)
    end
    PY = pullback(X.inc_Y, PX)
    set_attribute!(PY, :name, string("section: (",P[1]," : ",P[2]," : ",P[3],")"))
    set_attribute!(PY, :_self_intersection, -euler_characteristic(X))
    W =  WeilDivisor(PY, check=false)
    set_attribute!(W, :is_prime=>true)
    I = first(components(W))
    set_attribute!(I, :is_prime=>true)
    return new{typeof(X), typeof(coefficient_ring), elem_type(coefficient_ring)}(W, P)
  end
end

underlying_divisor(D::EllipticSurfaceSection) = D.D
rational_point(D::EllipticSurfaceSection) = D.P
