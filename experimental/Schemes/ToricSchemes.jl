export ToricSpec

export affine_normal_toric_variety

export ToricCoveredScheme

export normal_toric_variety

@attributes mutable struct ToricSpec{BRT, RT, SpecType<:Spec} <: AbsSpec{BRT, RT}
  X::SpecType
  antv::AffineNormalToricVariety

  function ToricSpec(antv::AffineNormalToricVariety; 
      R::MPolyRing=base_ring(toric_ideal(antv))
    )
    I = toric_ideal(R, antv)
    X = Spec(R, I)

    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X, antv)
  end
end

underlying_scheme(X::ToricSpec) = X.X
affine_normal_toric_variety(X::ToricSpec) = X.antv
antv = affine_normal_toric_variety
cone(X::ToricSpec) = cone(antv(X))
hilbert_basis(X::ToricSpec) = matrix(ZZ, hilbert_basis(cone(X)))

@Markdown.doc """
    torus_inclusions(X::ToricSpec)

For an affine toric variety ``X`` this returns a list `l` 
containing the inclusions ``Tʳ⁽ⁱ⁾ ↪ X`` of the different 
tori. 
"""
function torus_inclusions(X::ToricSpec)::Vector{<:AbsSpecMor}
end

@Markdown.doc """
    torus_action(X::ToricSpec)

For an affine toric variety ``X`` with a dense open torus ``T``
this returns a quintuple of morphisms `(pT, pX, incX, mult)` 
consisting of 

 * the projection ``T × X → T`` of the product with the torus ``T`` to ``T``
 * the projection ``T × X → X``
 * the inclusion ``X ↪ T × X`` taking ``x`` to ``(1, x)``
 * the group action ``T × X → X``.
"""
function torus_action(X::ToricSpec)::AbsSpecMor
end

function Base.show(io::IO, X::ToricSpec) 
  print(io, "Spec of a toric variety with Hilbert basis $(hilbert_basis(X))")
end

@attributes mutable struct ToricCoveredScheme{BRT, 
                                              CoveredSchemeType<:CoveredScheme
                                             } <: AbsCoveredScheme{BRT}
  X::CoveredSchemeType
  ntv::NormalToricVariety

  function ToricCoveredScheme(ntv::NormalToricVariety)
    F = fan(ntv)
    C = maximal_cones(F)

    antv = affine_open_covering(ntv)
    rings = [PolynomialRing(coefficient_ring(antv), ["x_$(i)_$(k)" for i in 1:ngens(base_ring(toric_ideal(antv[k])))], cached=false)[1] for k in 1:length(antv)]
    U = [ToricSpec(U, R=R) for (U, R) in zip(antv, rings)]
    cov = Covering(U)
    ### insert the glueings
    X = CoveredScheme(cov)
    return new{typeof(base_ring(X)), typeof(X)}(X, ntv)
  end
end

underlying_scheme(X::ToricCoveredScheme) = X.X
normal_toric_variety(X::ToricCoveredScheme) = X.ntv
ntv = normal_toric_variety
