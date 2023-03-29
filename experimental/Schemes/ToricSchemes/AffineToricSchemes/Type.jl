@attributes mutable struct ToricSpec{BRT, RT} <: AbsSpec{BRT, RT}
  antv::AffineNormalToricVariety
  var_names::Vector{Symbol}
  X::Spec
  dual_cone::Cone
  hb::ZZMatrix
  function ToricSpec(antv::AffineNormalToricVariety; 
      var_names::Vector{Symbol} = Vector{Symbol}()
    )
    # R = polynomial_ring(QQ, var_names)
    # Cdual = polarize(cone(antv))
    # hb = matrix(ZZ, hilbert_basis(Cdual))
    # I = toric_ideal(R, hb)
    # X = Spec(R, I)
    return new{typeof(QQ), MPolyQuoRing}(antv, var_names)
  end
end

function Base.show(io::IO, X::ToricSpec) 
  print(io, "Spec of an affine toric variety with cone spanned by $(rays(affine_normal_toric_variety(X)))")
end

