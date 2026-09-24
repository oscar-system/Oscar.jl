function domain(hom::PBWAlgebraHom)
  return hom.domain
end

function codomain(hom::PBWAlgebraHom)
  return hom.codomain
end

function image!(
  z::PBWAlgebraElem{S}, hom::PBWAlgebraHom{T,S}, x::PBWAlgebraElem{T}
) where {T,S}
  z = zero!(z)
  t1 = one(hom.codomain)
  t2 = zero(hom.codomain)
  cf = zero(coefficient_ring(parent(x)))

  for i in 1:length(x)
    for j in 1:ngens(hom.domain)
      for _ in 1:exponent(x, i, j)
        t2 = mul!(t2, t1, hom.img[j])
        t1 = swap!(t1, t2)
      end
    end

    cf = image!(cf, hom.field_automorphism, coeff(x, i))
    z = addmul!(z, t1, cf)
    t1 = one!(t1)
  end

  return z
end

function image(hom::PBWAlgebraHom{T,S}, x::PBWAlgebraElem{T}) where {T,S}
  @req parent(x) === hom.domain "parent mismatch"
  return image!(zero(hom.codomain), hom, x)
end

function (hom::PBWAlgebraHom{T})(x::PBWAlgebraElem{T}) where {T}
  return image(hom, x)
end
