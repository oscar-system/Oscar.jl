##############################################################################
#
# Various special cases for nested rings. Currently:
#  1. frac(R[t1, t2, ...])[x1, x2, ...]
#  2. R[z1, ...]...[y1, ...][x1, ...] and univariate possibilities
#
##############################################################################

export denest, renest

##############################################################################
#
# 1. implementations for Q(t)[x]
#
##############################################################################

# from Q(t1, t2)[x1, x2] to S = Q[x1, x2, t1, t2]
# the product of the two outputs should be equal to the input
function _remove_denominators(
  S::Union{MPolyRing, PolyRing},
  f::Union{MPolyElem, PolyElem})

  R = parent(f)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  den = one(R2)
  cont = zero(R2)
  for (c1, e1) in zip(coefficients(f), exponent_vectors(f))
    cont = gcd(cont, numerator(c1))
    den = lcm(den, denominator(c1))
  end
  num = MPolyBuildCtx(S)  # PolyElem should satisfy the exponent vector interface
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
function _restore_numerators(
  R::Union{MPolyRing, PolyRing},
  g::Union{MPolyElem, PolyElem})

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
function _add_variables(R2::Singular.PolyRing, n::Int)
  n += nvars(R2)
  x = map(i->"x"*string(i), 1:n)
  return Singular.PolynomialRing(base_ring(R2), x, cached=false)[1]
end

function _add_variables(R2::Union{MPolyRing, PolyRing}, n::Int)
  n += nvars(R2)
  x = map(i->"x"*string(i), 1:n)
  return PolynomialRing(base_ring(R2), x, cached=false)[1]
end


function Oscar.gcd(
  a::AbstractAlgebra.MPolyElem{AbstractAlgebra.Generic.Frac{T}},
  b::AbstractAlgebra.MPolyElem{AbstractAlgebra.Generic.Frac{T}}
) where T <: Union{PolyElem, MPolyElem}

  R = parent(a)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  R == parent(b) || error("parents do not match")
  S = _add_variables(R2, nvars(R))
  A = _remove_denominators(S, a)[1]
  B = _remove_denominators(S, b)[1]
  return _restore_numerators(R, gcd(A, B))    # not nec monic
end

function Oscar.gcd(
  a::AbstractAlgebra.PolyElem{AbstractAlgebra.Generic.Frac{T}},
  b::AbstractAlgebra.PolyElem{AbstractAlgebra.Generic.Frac{T}}
) where T <: Union{PolyElem, MPolyElem}

  R = parent(a)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  R == parent(b) || error("parents do not match")
  S = _add_variables(R2, nvars(R))
  A = _remove_denominators(S, a)[1]
  B = _remove_denominators(S, b)[1]
  return _restore_numerators(R, gcd(A, B))    # not nec monic
end

function _convert_frac_fac(R, u, fac)
  Rfac = Fac{elem_type(R)}()
  Rfac.unit = R(u)*_restore_numerators(R, fac.unit)
  for (f, e) in fac
    t = _restore_numerators(R, f)
    if is_constant(t)
      Rfac.unit *= t^e
    else
      Rfac[t] = e
    end
  end
  return Rfac
end

function Oscar.factor(
  a::Union{AbstractAlgebra.MPolyElem{AbstractAlgebra.Generic.Frac{T}},
           AbstractAlgebra.PolyElem{AbstractAlgebra.Generic.Frac{T}}}
) where T <: Union{PolyElem, MPolyElem}

  R = parent(a)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  S = _add_variables(R2, nvars(R))
  A, u = _remove_denominators(S, a)
  return _convert_frac_fac(R, u, factor(A))
end

function Oscar.factor_squarefree(
  a::Union{AbstractAlgebra.MPolyElem{AbstractAlgebra.Generic.Frac{T}},
           AbstractAlgebra.PolyElem{AbstractAlgebra.Generic.Frac{T}}}
) where T <: Union{PolyElem, MPolyElem}

  R = parent(a)                   # Q(t1, t2)[x1, x2]
  R2 = base_ring(base_ring(R))    # Q[t1, t2]
  S = _add_variables(R2, nvars(R))
  A, u = _remove_denominators(S, a)
  return _convert_frac_fac(R, u, factor_squarefree(A))
end

##############################################################################
#
# 2. implementations for Q[z][y][x] ect
#
##############################################################################

# denest for poly ring

function _denest_recursive(BR::Union{MPolyRing, PolyRing},
                           R::Union{MPolyRing, PolyRing}, n::Int)
  return _denest_recursive(base_ring(BR), BR, n + nvars(R))
end

function _denest_recursive(BR, R::Union{MPolyRing, PolyRing}, n::Int)
  return R, n
end

@doc Markdown.doc"""
    denest(R::Union{PolyRing, MPolyRing})

Return a multivariate polynomial ring resulting from denesting an iterated
polynomial ring `R`.
"""
function denest(R::Union{PolyRing, MPolyRing})
  R2, n = _denest_recursive(base_ring(R), R, 0)
  return _add_variables(R2, n)
end

# denest for poly elem

function _denest_recursive(r::MPolyBuildCtx, f::PolyElem{T},
                           e0::Vector{Int}) where T <: Union{PolyElem, MPolyElem}
  for i in degree(f):-1:0
    _denest_recursive(r, coeff(f, i), vcat(e0, [i]))
  end
end

function _denest_recursive(r::MPolyBuildCtx, f::PolyElem, e0::Vector{Int})
  for i in degree(f):-1:0
    push_term!(r, coeff(f, i), vcat(e0, [i]))
  end
end

function _denest_recursive(r::MPolyBuildCtx, f::MPolyElem{T},
                           e0::Vector{Int}) where T <: Union{PolyElem, MPolyElem}
  for (c1, e1) in zip(coefficients(f), exponent_vectors(f))
    _denest_recursive(r, c1, vcat(e0, e1))
  end
end

function _denest_recursive(r::MPolyBuildCtx, f::MPolyElem, e0::Vector{Int})
  for (c1, e1) in zip(coefficients(f), exponent_vectors(f))
    push_term!(r, c1, vcat(e0, e1))
  end
end

@doc Markdown.doc"""
    denest(S::MPolyRing, f::Union{PolyElem, MPolyElem})

Return an element of `S` resulting from denesting a element `f` of an
iterated polynomial ring. The ring `S` should have the same base ring and
number of variables as `denest(parent(f))`.
"""
function denest(S::MPolyRing, f::Union{PolyElem, MPolyElem})
  r = MPolyBuildCtx(S)
  _denest_recursive(r, f, Int[])
  return finish(r)
end

# renest for poly

function _renest_recursive_coeff(R::Union{MPolyRing{T}, PolyRing{T}},
                                 off::Int, idxs::Vector{Int},
                                 gcoeffs, gexps) where T <: Union{PolyElem, MPolyElem}
  return _renest_recursive(base_ring(R), off, idxs, gcoeffs, gexps)
end

function _renest_recursive_coeff(R::Union{MPolyRing, PolyRing}, off::Int,
                                 idxs::Vector{Int}, gcoeffs, gexps)
  @assert length(idxs) == 1
  return gcoeffs[idxs[1]]
end

function _renest_recursive(R::PolyRing, off::Int, idxs::Vector{Int}, gcoeffs, gexps)
  d = Dict{Int, Vector{Int}}()   # exp => indices
  for i in idxs
    e1 = gexps[i][off]
    if !haskey(d, e1)
      d[e1] = [i]
    else
      append!(d[e1], i)
    end
  end
  # TODO once PolyElem has a better required constructor
  z = zero(R)
  for (e1, ii) in d
    z += _renest_recursive_coeff(R, off + 1, ii, gcoeffs, gexps)*gen(R)^e1
  end
  return z
end

function _renest_recursive(R::MPolyRing, off::Int, idxs::Vector{Int}, gcoeffs, gexps)
  n = nvars(R)
  d = Dict{Vector{Int}, Vector{Int}}()  # exp => indices
  for i in idxs
    e1 = gexps[i][off:off+n-1]
    if !haskey(d, e1)
      d[e1] = [i]
    else
      append!(d[e1], i)
    end
  end
  z = MPolyBuildCtx(R)
  for (e1, ii) in d
    push_term!(z, _renest_recursive_coeff(R, off + n, ii, gcoeffs, gexps), e1)
  end
  return finish(z)
end

@doc Markdown.doc"""
    renest(R::Union{PolyRing, MPolyRing}, g::MPolyElem)

Return an element of iterated polynomial ring `R` from its denested
counterpart. The return satisfies `f == renest(R, denest(denest(R), f))`.
"""
function renest(R::Union{PolyRing, MPolyRing}, g::MPolyElem)
  gcoeffs = collect(coefficients(g))
  gexps = collect(exponent_vectors(g))
  return _renest_recursive(R, 1, collect(1:length(gcoeffs)), gcoeffs, gexps)
end

function _nested_gcd(a, b)
  R = parent(a)
  R == parent(b) || error("parents do not match")
  S = denest(R)
  return renest(R, gcd(denest(S, a), denest(S, b)))
end

function Oscar.gcd(
   a::AbstractAlgebra.Generic.Poly{T},
   b::AbstractAlgebra.Generic.Poly{T}
) where T <: Union{PolyElem, MPolyElem}
   return _nested_gcd(a, b)
end

function Oscar.gcd(
   a::AbstractAlgebra.Generic.MPoly{T},
   b::AbstractAlgebra.Generic.MPoly{T}
) where T <: Union{PolyElem, MPolyElem}
   return _nested_gcd(a, b)
end

function _convert_iter_fac(R, fac::Fac)
  Rfac = Fac{elem_type(R)}()
  Rfac.unit = renest(R, fac.unit)
  for (f, e) in fac
    Rfac[renest(R, f)] = e
  end
  return Rfac
end

function Oscar.factor(a::Union{PolyElem{T}, MPolyElem{T}}) where
                                                T <: Union{PolyElem, MPolyElem}
  A = denest(denest(parent(a)), a)
  return _convert_iter_fac(parent(a), factor(A))
end

function Oscar.factor_squarefree(a::Union{PolyElem{T}, MPolyElem{T}}) where
                                                T <: Union{PolyElem, MPolyElem}
  A = denest(denest(parent(a)), a)
  return _convert_iter_fac(parent(a), factor_squarefree(A))
end
