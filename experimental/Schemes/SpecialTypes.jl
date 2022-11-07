export underlying_morphism, complement_ideal, complement_scheme

export image_ideal

export ideal_type

function is_constant(a::MPolyLocalizedRingElem) 
  reduce_fraction(a)
  return is_constant(numerator(a)) && is_constant(denominator(a))
end

### Already implemented in AA -- but probably buggy?
#function is_zero_divisor(f::MPolyElem)
#  iszero(f) && return true
#  if is_constant(f)
#    c = coefficients(f)[1]
#    return is_zero_divisor(c)
#  end
#  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
#end

function is_zero_divisor(f::MPolyLocalizedRingElem)
  iszero(f) && return true
  if is_constant(f)
    c = first(coefficients(numerator(f)))
    return is_zero_divisor(c)
  end
  return is_zero_divisor(numerator(f))
end

### The following method is only implemented when the coefficient ring is a field.
# The code should be valid generically, but the Singular backend needed for the 
# ideal quotient is probably buggy for non-fields.
function is_zero_divisor(f::MPolyQuoElem{<:MPolyElem{<:FieldElem}})
  iszero(f) && return true
  b = simplify(f)
  # The next block is basically useless when the coefficient ring is 
  # a field, because it is merely another `is_zero`-check. However, 
  # once more functionality is working, it will actually do stuff and 
  # the above signature can be widened.
  if is_constant(lift(b))
    c = first(coefficients(lift(b)))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end

function is_zero_divisor(f::MPolyQuoLocalizedRingElem{<:Field})
  iszero(f) && return true
  # The next block is basically useless when the coefficient ring is 
  # a field, because it is merely another `is_zero`-check. However, 
  # once more functionality is working, it will actually do stuff and 
  # the above signature can be widened.
  if is_constant(lifted_numerator(f)) && is_constant(lifted_denominator(f))
    c = first(coefficients(lift(numerator(f))))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end



########################################################################
# Methods for ClosedEmbedding                                          #
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


