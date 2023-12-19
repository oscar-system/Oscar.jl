@doc raw"""
    PrincipalOpenEmbedding{DomainType, CodomainType, PullbackType} <: AbsSpecMor

An open inclusion ``ι : U ↪ X`` of one affine scheme ``U`` into another 
one ``X`` such that the image is a principal open subset. 
"""
@attributes mutable struct PrincipalOpenEmbedding{DomainType, CodomainType, PullbackType} <: AbsSpecMor{DomainType, CodomainType, PullbackType, OpenInclusion, Nothing}
  inc::SpecMor{DomainType, CodomainType, PullbackType}
  complement_equations::Vector{<:RingElem}
  image::PrincipalOpenSubset
  inverse_on_image::AbsSpecMor

  function PrincipalOpenEmbedding(f::AbsSpecMor, complement_equations::Vector{T}; check::Bool=true) where {T<:RingElem}
    X = domain(f)
    Y = codomain(f)
    @assert all(x->parent(x) === OO(Y), complement_equations)
    @check all(h->is_unit(pullback(f)(h)), complement_equations) "image of the map is not contained in the complent of the zero locus of the given equations"
    @check begin
      U = PrincipalOpenSubset(Y, complement_equations)
      f_res = SpecMor(X, U, pullback(f).(gens(OO(Y))), check=check)
      # is_isomorphism does not work for all required constellations.
      #is_isomorphism(f) || error("restriction is not an isomorphism")
    end

    return new{typeof(X), typeof(Y), pullback_type(f)}(f, complement_equations)
  end
end

