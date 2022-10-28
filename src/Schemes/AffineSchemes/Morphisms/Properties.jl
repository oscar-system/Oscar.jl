export is_isomorphism, is_inverse_of, is_identity_map


########################################
# (1) Isomorphism, inverse and identity
########################################

@attr Bool function is_isomorphism(f::AbsSpecMor)
  has_attribute(f, :inverse) && return true
  is_isomorphism(pullback(f)) || return false
  set_attribute!(f, :inverse, SpecMor(codomain(f), domain(f), inverse(pullback(f))))
  return true
end

function is_inverse_of(f::S, g::T) where {S<:AbsSpecMor, T<:AbsSpecMor}
  return is_isomorphism(f) && (inverse(f) == g)
end

is_identity_map(f::AbsSpecMor) = (domain(f) == codomain(f)) && all(x->(pullback(f)(x) == x), gens(OO(domain(f))))
