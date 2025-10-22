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
      new_den = OO(X)(V, U)(den)
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
        new_den = OO(X)(V, UV)(den)
        frac = (lifted_numerator(new_num)*lifted_denominator(new_den))//(lifted_denominator(new_num)*lifted_numerator(new_den))
        res = F(U, UV)
        comp = matrix(res)[1, 1]
        frac = frac*(lifted_denominator(comp)//lifted_numerator(comp))
        return OO(U)(numerator(frac)), OO(U)(denominator(frac))
      else
        VV = __find_chart(V, triv_cov)
        res = F(VV, V)
        mat_res = matrix(res)[1, 1]
        num_res = OO(X)(VV, V)(num)
        den_res = OO(X)(VV, V)(den)
        frac = (lifted_numerator(num_res)*lifted_denominator(den_res))//(lifted_denominator(num_res)*lifted_numerator(den_res))
        frac = (numerator(frac)*lifted_denominator(mat_res))//(denominator(frac)*lifted_numerator(mat_res))
        new_num = numerator(frac)
        new_den = denominator(frac)
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

function weil_divisor(s::RationalSection; ring::Ring=ZZ)
  M = sheaf(s)
  X = scheme(M)
  cov = trivializing_covering(M)
  #decomp_info = decomposition_info(cov)
  result = weil_divisor(X, ring)

  for U in patches(cov)
    num, den = s(U)
    num_ideal = ideal(num)
    den_ideal = ideal(den)
    num_primes = [P for P in minimal_primes(num_ideal) if dim(P) == dim(X) - 1]
    den_primes = [P for P in minimal_primes(den_ideal) if dim(P) == dim(X) - 1]
    num_mults = Dict{Ideal, Int}(P=>_colength_in_localization(num_ideal, P) for P in num_primes)
    den_mults = Dict{Ideal, Int}(P=>_colength_in_localization(den_ideal, P) for P in den_primes)
    for (P, k) in den_mults
      if haskey(num_mults, P)
        e = num_mults[P]
        e = e - k
        if is_zero(e)
          delete!(num_mults, P)
        else
          num_mults[P] = e
        end
      else
        num_mults[P] = -k
      end
    end

    for (P, k) in num_mults
      found = false
      for Q in components(result) 
        if Q(U) == P
          @assert k == result[Q]
          found = true
          break
        end
      end
      found && continue
      PP = PrimeIdealSheafFromChart(X, U, P)
      inc = k*weil_divisor(PP, ring)
      result = result + inc
    end
  end
  return result
end

function check_codim(U,I)
    n = dim(U)
    @assert base_ring(I) == OO(U) "The provided ideal is not an ideal of the provided ring"
    n - dim(I) == 1
end
    

function div(s::Oscar.RationalSection)
    X = Oscar.scheme(s)
    F = Oscar.sheaf(s)
    triv_cov = trivializing_covering(F)
    simp_num = algebraic_cycle(X,ZZ)
    simp_denom = algebraic_cycle(X,ZZ)
    for U in triv_cov
        num = s(U)[1]
        denom = s(U)[2]
        if !is_unit(num)
            N = minimal_primes(ideal(OO(X)(U),num))
            for p in N
                if check_codim(U,p)
                    if !any(map(t -> t(U) == p,components(simp_num)))
                        stalk = localization(OO(U),complement_of_prime_ideal(p))
                        fp = ideal(stalk[1],stalk[2](num))
                        multp = length(quotient_ring_as_module(quo(stalk[1],fp)[1]))
                        simp_num = simp_num + multp*algebraic_cycle(Oscar.PrimeIdealSheafFromChart(X,U,p))
                    end
                end
            end
        end
        if !is_unit(denom)
            D = minimal_primes(ideal(OO(X)(U),denom))
            for p in D
                if check_codim(U,p)
                    if !any(map(t -> t(U) == p,components(simp_denom)))
                        stalk = localization(OO(U),complement_of_prime_ideal(p))
                        fp = ideal(stalk[1],stalk[2](denom))
                        multp = length(quotient_ring_as_module(quo(stalk[1],fp)[1]))
                        simp_denom = simp_denom + multp*algebraic_cycle(Oscar.PrimeIdealSheafFromChart(X,U,p))
                    end
                end
            end
        end
    end
    simp_num - simp_denom
end
