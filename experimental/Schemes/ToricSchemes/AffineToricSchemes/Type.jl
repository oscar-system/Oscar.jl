export ToricSpec

@attributes mutable struct ToricSpec{BRT, RT, SpecType<:Spec} <: AbsSpec{BRT, RT}
  X::SpecType
  antv::AffineNormalToricVariety
  dual_cone::Cone
  hb::fmpz_mat
  function ToricSpec(antv::AffineNormalToricVariety; R::MPolyRing=base_ring(toric_ideal(antv)))
    Cdual = polarize(cone(antv))
    hb = matrix(ZZ, hilbert_basis(Cdual))
    I = toric_ideal(R, hb)
    X = Spec(R, I)
    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X, antv, Cdual, hb)
  end
end

function Base.show(io::IO, X::ToricSpec) 
  print(io, "Spec of an affine toric variety with cone spanned by $(rays(affine_normal_toric_variety(X)))")
end
