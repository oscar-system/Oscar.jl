@doc raw"""
    order_omega_mod_N(d::IntegerUnion, q::IntegerUnion, N::IntegerUnion) -> Pair{Bool, Bool}

Return `(flag_plus, flag_minus)` where `flag_plus` and `flag_minus`
are `true` or `false`, depending on whether `N` divides the order
of the orthogonal groups $\Omega^+(d, q)$ and $\Omega^-(d, q)$.

# Examples
```jldoctest
julia> Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 60)
(false, true)

julia> Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 5, 60)
(true, true)
```
"""
function order_omega_mod_N(d::IntegerUnion, q::IntegerUnion, N::IntegerUnion)
  @req is_even(d) "d must be even"
  m = div(d, 2)
  exp, N = remove(N, q)
  facts = collect(factor(q))
  p = facts[1][1]
  if mod(N, p) == 0
    exp = exp + 1
    _, N = remove(N, p)
  end
  if m*(m-1) < exp
    # A group of order `N` does not embed in any candidate.
    return (false, false)
  end

  q2 = ZZ(q)^2
  q2i = ZZ(1)
  for i in  1:(m-1)
    q2i = q2 * q2i
    if i == 1 && is_odd(q)
      g = gcd(N, div(q2i-1, 2))
    else
      g = gcd(N, q2i-1)
    end
    N = div(N, g)
    if N == 1
      # A group of order N may embed in both candidates.
      return (true, true)
    end
  end

  # embeds in + type?, embeds in - type?
  return (mod(q^m-1, N) == 0, mod(q^m+1, N) == 0)
end


@doc raw"""
    reduce_mod_squares(val::nf_elem)

Return an element of `F = parent(val)` that is equal to `val`
modulo squares in `F`.

If `val` describes an integer then the result corresponds to the
squarefree part of this integer.
Otherwise the coefficients of the result have a squarefree g.c.d.

# Examples
```jldoctest
julia> F, z = cyclotomic_field(4);

julia> Oscar.OrthogonalDiscriminants.reduce_mod_squares(4*z^0)
1

julia> Oscar.OrthogonalDiscriminants.reduce_mod_squares(-8*z^0)
-2
```
"""
function reduce_mod_squares(val::nf_elem)
  is_zero(val) && return val
  d = denominator(val)
  if ! isone(d)
    val = val * d^2
  end
  if is_integer(val)
    intval = ZZ(val)
    sgn = sign(intval)
    good = [x[1] for x in collect(factor(intval)) if is_odd(x[2])]
    F = parent(val)
    return F(prod(good, init = sgn))
  end
  # Just get rid of the square part of the gcd of the coefficients.
  c = map(numerator, coefficients(val))
  s = 1
  for (p, e) in collect(factor(gcd(c)))
    if iseven(e)
      s = s * p^e
    elseif e > 1
      s = s * p^(e-1)
    end
  end
  return val//s
end
