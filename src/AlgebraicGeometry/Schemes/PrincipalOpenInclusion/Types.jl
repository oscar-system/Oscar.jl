# This type is needed to handle the morphisms in refinements of coverings efficiently. 
# In principal, it does not much more than `PrincipalOpenSubset` with it's implicit 
# embedding into its ambient scheme. But this allows for an actual morphism of rings 
# in the background, so that the identification of the domain with a principal open 
# subset of its codomain can also be realized for schemes with different `ambient_coordinate_ring`s. 
# This is important for the refinements, since we also allow `SimplifiedAffineScheme`.
@doc raw"""
    PrincipalOpenEmbedding{DomainType, CodomainType, PullbackType} <: AbsAffineSchemeMor

An open inclusion ``ι : U ↪ X`` of one affine scheme ``U`` into another 
one ``X`` such that the image is a principal open subset. 
"""
@attributes mutable struct PrincipalOpenEmbedding{DomainType, CodomainType, PullbackType} <: AbsAffineSchemeMor{DomainType, CodomainType, PullbackType, OpenInclusion, Nothing}
  inc::AffineSchemeMor{DomainType, CodomainType, PullbackType}
  complement_equations::Vector{<:RingElem}
  image::PrincipalOpenSubset
  inverse_on_image::AbsAffineSchemeMor

  function PrincipalOpenEmbedding(f::AbsAffineSchemeMor, complement_equations::Vector{T}; check::Bool=true) where {T<:RingElem}
    X = domain(f)
    Y = codomain(f)
    @assert all(x->parent(x) === OO(Y), complement_equations)
    @check all(h->is_unit(pullback(f)(h)), complement_equations) "image of the map is not contained in the complent of the zero locus of the given equations"
    @check begin
      U = PrincipalOpenSubset(Y, complement_equations)
      f_res = morphism(X, U, pullback(f).(gens(OO(Y))), check=check)
      is_isomorphism(f_res) || error("restriction is not an isomorphism")
    end

    return new{typeof(X), typeof(Y), pullback_type(f)}(f, complement_equations)
  end
end

