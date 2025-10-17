mutable struct RationalSection{SheafType<:AbsCoherentSheaf}
  sheaf::SheafType
  cache::IdDict{AbsAffineScheme, Tuple{<:RingElem, <:RingElem}}

  function RationalSection(F::AbsCoherentSheaf)
    return new{typeof(F)}(F, IdDict{AbsAffineScheme, Tuple{<:RingElem, <:RingElem}}())
  end
end

sheaf(s::RationalSection) = s.sheaf
scheme(s::RationalSection) = scheme(sheaf(s))

function (s::RationalSection)(U::AbsAffineScheme)
  return get!(s.cache, U) do
    produce_rational_section(s, U)
  end
end

function produce_rational_section(s::RationalSection, U::AbsAffineScheme)
  @assert has_ancestor(x->any(y===x for y in patches(trivializing_covering(sheaf(s)))), U)
  if isempty(s.cache)
    num = den = one(OO(U))
    return num, den
  end

  F = sheaf(s)
  @assert is_one(ngens(F(U))) "local representatives of rational sections can only be given on charts for which the sheaf is free of rank one"
  X = scheme(s)

  # if an ancestor has a value already stored, restrict from there.
  for (V, (num, den)) in s.cache
    if has_ancestor(x->(x===V), U)
      @assert is_one(ngens(F(V)))
      e = first(gens(F(V)))
      res = F(V, U)
      new_num = res(num*e)[1]
      new_den = res(den*e)[1]
      return new_num, new_den
    end
  end

  triv_cov = trivializing_covering(F)
  if any(x->(x===U), patches(triv_cov))
    for (V, (num, den)) in s.cache
      if any(x->x===V, patches(triv_cov))
        glue = triv_cov[V, U]
        VU, UV = gluing_domains(glue)
        trans = F(V, UV)
        e = first(gens(F(V)))
        new_num = trans(num*e)[1]
        new_den = trans(den*e)[1]
        frac = (lifted_numerator(new_num)*lifted_denominator(new_den))//(lifted_denominator(new_num)*lifted_numerator(new_den))
        res = F(U, UV)
        comp = matrix(res)[1, 1]
        frac = frac*(lifted_denominator(comp)//lifted_numerator(comp))
        return OO(U)(numerator(frac)), OO(U)(denominator(frac))
      else
        VV = __find_chart(V, triv_cov)
        res = F(VV, V)
        e = first(gens(F(V)))
        new_num = preimage(res, num*e)[1]
        new_den = preimage(res, den*e)[1]
        s.cache[VV] = (new_num, new_den)
        return produce_rational_section(s, U)
      end
    end
  end

  V = __find_chart(U, triv_cov)
  num, den = s(V)
  # cache has been filled by the previous call so that the next 
  # one is caught above
  return produce_rational_section(s, U)
end
      
