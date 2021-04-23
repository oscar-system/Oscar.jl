##############################################################################
#
# implementations for Q(t)[x]
#
##############################################################################

# from Q(t1, t2)[x1, x2] to S = Q[x1, x2, t1, t2]
# the product of the two outputs should be equal to the input
function _remove_denominators(S::MPolyRing, f::MPolyElem)
  R = parent(f)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  den = one(R2)
  cont = zero(R2)
  for (c1, e1) in zip(coefficients(f), exponent_vectors(f))
    cont = gcd(cont, numerator(c1))
    den = lcm(den, denominator(c1))
  end
  num = MPolyBuildCtx(S)
  for (c1, e1) in zip(coefficients(f), exponent_vectors(f))
    cf = divexact(den, denominator(c1))
    if !iszero(cont)
      cf *= divexact(numerator(c1), cont)
    end
    for (c2, e2) in zip(coefficients(cf), exponent_vectors(cf))
      push_term!(num, c2, vcat(e1, e2))
    end
  end
  return (finish(num), cont//den)
end

# from Q[x1, x2, t1, t2] to R = Q(t1, t2)[x1, x2]
# the output is equal to the input
function _restore_numerators(R::MPolyRing, g::MPolyElem)
  S = parent(g)                   # Q[x1, x2, t1, t2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  n = nvars(R)
  d = Dict{Vector{Int}, typeof(MPolyBuildCtx(R2))}()
  for (c, e) in zip(coefficients(g), exponent_vectors(g))
    e1 = e[1:n]
    e2 = e[n+1:end]
    if !haskey(d, e1)
      d[e1] = MPolyBuildCtx(R2)
    end
    push_term!(d[e1], c, e2)
  end
  z = MPolyBuildCtx(R)
  for (e, c) in d
    push_term!(z, base_ring(R)(finish(c)), e)
  end
  return finish(z)
end

# if R2 is a singular ring, keep it that way since limited types can be factored
function _add_variables(R2::MPolyRing, n::Int)
   n += nvars(R2)
   x = map(i->"x"*string(i), 1:n)
   if isa(R2, Singular.PolyRing)
      return Singular.PolynomialRing(base_ring(R2), x)[1]
   else
      return PolynomialRing(base_ring(R2), x)[1]
   end
end


function Oscar.gcd(
  a::AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{T}},
  b::AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{T}}
) where T <: MPolyElem

  R = parent(a)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  R == parent(b) || error("parents do not match")
  S = _add_variables(R2, nvars(R))
  A = _remove_denominators(S, a)[1]
  B = _remove_denominators(S, b)[1]
  return _restore_numerators(R, gcd(A, B))    # not nec monic
end


function _convert_fac(R, u, fac)
  Rfac = Fac{elem_type(R)}()
  Rfac.unit = R(u)
  for (f, e) in fac
    t = _restore_numerators(R, f)
    if isconstant(t)
      Rfac.unit *= t^e
    else
      Rfac[t] = e
    end
  end
  return Rfac
end

function Oscar.factor(
  a::AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{T}}
) where T <: MPolyElem

  R = parent(a)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  S = _add_variables(R2, nvars(R))
  A, u = _remove_denominators(S, a)
  return _convert_fac(R, u, factor(A))
end

function Oscar.factor_squarefree(
  a::AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{T}}
) where T <: MPolyElem

  R = parent(a)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  S = _add_variables(R2, nvars(R))
  A, u = _remove_denominators(S, a)
  return _convert_fac(R, u, factor_squarefree(A))
end


