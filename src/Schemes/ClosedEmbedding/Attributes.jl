export image_ideal, complement

########################################################################
# Attributes for ClosedEmbedding                                       #
########################################################################
underlying_morphism(f::ClosedEmbedding) = f.inc

@Markdown.doc """
    image_ideal(f::ClosedEmbedding)

For a closed embedding ``f : X → Y`` of affine schemes ``X = Spec(S)`` 
into ``Y = Spec(R)`` such that ``S ≅ R/I`` via ``f`` for some ideal 
``I ⊂ R`` this returns ``I``.
"""
image_ideal(f::ClosedEmbedding) = f.I::ideal_type(OO(codomain(f)))

function complement(f::ClosedEmbedding)
  if !isdefined(f, :U)
    U = SpecOpen(codomain(f), image_ideal(f))
    f.U = U
  end
  return f.U
end

ideal_type(::Type{RT}) where {RT<:MPolyRing} = MPolyIdeal{elem_type(RT)}
ideal_type(::Type{RT}) where {PolyType, RT<:MPolyQuo{PolyType}} = MPolyQuoIdeal{PolyType}
ideal_type(R::Ring) = ideal_type(typeof(R))


