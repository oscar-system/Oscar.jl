
# We are given a flag bundle F = FB(E)_{_1, ..., r_l} --> X, where the
# vector bundle E on X has rank n, and the positive integers r_i, the ranks
# of the subquotient bundles Q_i, sum up to n. The Chow ring of F is the
# quotient ring A[c_ij's]/I, where the coefficient ring A = K[x's]/IA is the
# Chow ring of X and the c_ij are the Chern classes cj( Q_i) with associated
# weight j. See the docstrings of the `flag_bundle` constructor in `main.jl`
# for the relations collected in I.

# Shared presentation and public interface for the two flag algorithms.
# Engines are initialized independently: :gysin never constructs a reducer.
mutable struct _FlagBundleData{N,A,S,C}
  coefficient_ring::A
  ring::S # = A[c_ij's], reduction happens there
  ranks::Vector{Int} # the r_i, the ranks of the subquotien bundles
  starts::Vector{Int} # ends-ranks
  ends::Vector{Int} # the ranks of the subbundles (= dims in flag bundle constructor)
  chern::Vector{C} # the Chern classes of E
  default_algorithm::Symbol
  reducer::Union{Nothing,_FlagBundleReducer}
  gysin::Union{Nothing,_FlagBundleGysin}
  cache_lock::ReentrantLock
end

function _flag_check_algorithm(algorithm::Symbol)
  @req algorithm in (:monic, :gysin) "Expected algorithm=:monic or algorithm=:gysin"
  return nothing
end

function _FlagBundleData(A::Ring, ranks::Vector{Int}, chern::Vector;
                         names::Vector{String}=["c$i" for i in 1:sum(ranks)],
                         algorithm::Symbol=:monic)
  _flag_check_algorithm(algorithm)
  @req !isempty(ranks) && all(>(0), ranks) "Expected positive subquotient ranks"
  n = sum(ranks)
  @req length(chern) == n "Expected c_1(E),...,c_n(E), including trailing zeros"
  @req length(names) == n "Expected one variable name per Chern class"
  P = polynomial_ring(A, names; cached=false)[1]
  C = elem_type(A)
  b = C[_flag_coefficient(A(c)) for c in chern]
  ends = cumsum(ranks)
  return _FlagBundleData{n,typeof(A),typeof(P),C}(
    A, P, copy(ranks), ends - ranks, ends, b, algorithm,
    nothing, nothing, ReentrantLock(),
  )
end

function _flag_engine(D::_FlagBundleData, algorithm::Symbol)
  # sets up engines for the two algorithms
  _flag_check_algorithm(algorithm)
  return lock(D.cache_lock) do
    if algorithm === :monic
      if isnothing(D.reducer)
        D.reducer = _FlagBundleReducer(
          D.coefficient_ring, D.ranks, D.chern; relative_ring=D.ring
        )
      end
      return D.reducer
    else
      if isnothing(D.gysin)
        D.gysin = _FlagBundleGysin(D)
      end
      return D.gysin
    end
  end
end

function _FlagStandardExponents(D::_FlagBundleData{N}) where {N}
  return _FlagStandardExponents{N}(D.starts, D.ends)
end

# The adapter is deliberately separate from the coefficient-ring-independent
# reducer. The ambient Oscar quotient keeps its existing ordering and APIs.
_flag_polynomial(a::MPolyRingElem) = a
_flag_polynomial(a::MPolyDecRingElem) = a.f   # forget decoration
_flag_polynomial(a::MPolyQuoRingElem) = _flag_polynomial(lift(a))

function _flag_relative_polynomial(D::_FlagBundleData{N}, f) where {N}
  # K[x's, c_ij's]--> S without grading
  raw = _flag_polynomial(f)
  A = D.coefficient_ring
  x = gens(A)
  @req ngens(parent(raw)) == N + length(x) "Unexpected ambient polynomial ring"
  ctx = MPolyBuildCtx(D.ring) # D.ring = S without grading
  for (a, e) in zip(coefficients(raw), AbstractAlgebra.exponent_vectors(raw))
    b = A(a)
    for j in eachindex(x)
      e[N + j] == 0 && continue
      b *= x[j]^e[N + j]
    end
    push_term!(ctx, _flag_coefficient(b), e[1:N])
  end
  return finish(ctx)
end

function _flag_flatten(Fl::AbstractVariety, D::_FlagBundleData, f)
  # S without grading --> K[x's, c_ij's]
  R = base_ring(chow_ring(Fl))
  ctx = MPolyBuildCtx(R.R) # R.R: forget grading
  for (a, e) in zip(coefficients(f), AbstractAlgebra.exponent_vectors(f))
    raw = _flag_polynomial(_flag_coefficient(a))
    for (b, v) in zip(coefficients(raw), AbstractAlgebra.exponent_vectors(raw))
      push_term!(ctx, b, vcat(e, v))
    end
  end
  return R(finish(ctx))
end

function _flag_bundle_data(Fl::AbstractVariety)
  # getter function
  D = get_attribute(Fl, :flag_bundle_data)
  @req D isa _FlagBundleData "Expected a variety constructed by flag_bundle"
  @req D.coefficient_ring === chow_ring(codomain(structure_map(Fl))) "The base Chow ring has changed"
  return D
end

# The input may be any representative, not only a quotient normal form.
function _flag_bundle_pushforward(Fl::AbstractVariety, a::MPolyDecRingOrQuoElem,
                                  algorithm::Symbol)
  # gatekeeper w.r.t. algorithm
  _flag_check_algorithm(algorithm)
  D = _flag_bundle_data(Fl)
  Q = chow_ring(Fl)
  @req parent(a) === Q || parent(a) === base_ring(Q) "Class belongs to a different variety"
  f = _flag_relative_polynomial(D, a)
  return _flag_pushforward(_flag_engine(D, algorithm), f)
end

function _flag_pushforward_override(p::AbstractVarietyMap,
                                    a::MPolyDecRingOrQuoElem, algorithm::Symbol)
  _flag_check_algorithm(algorithm)
  Fl = domain(p)
  @req get_attribute(Fl, :flag_bundle_data) isa _FlagBundleData "Algorithm selection is only supported for flag-bundle structure maps"
  @req p === structure_map(Fl) "Algorithm selection is only supported for the flag-bundle structure map"
  return _flag_bundle_pushforward(Fl, a, algorithm)
end

@doc raw"""
    flag_bundle_normal_form(Fl::AbstractVariety, a::MPolyDecRingOrQuoElem)

Return the relative standard-monomial representative of `a` for a variety
constructed by [`flag_bundle`](@ref). The result belongs to
`base_ring(chow_ring(Fl))`. The input may be a class in `chow_ring(Fl)` or a
polynomial in that ambient graded polynomial ring.

For successive subquotient ranks `r[i]`, set `s[i] = sum(r[1:i-1])`. In the
result, every block of Chern variables satisfies
`sum(exponents_in_block_i) <= s[i]`. This is an **unweighted** exponent bound;
reduction uses the product of weighted degree reverse lexicographic orders,
with weights `1:r[i]` in block `i`, and block `1` the most significant.
Coefficients are reduced in the actual Chow ring of the base.

The result need not agree with `lift(simplify(a))`: the ambient quotient may
use a different ordering. No fraction field of the base is introduced.
See [GSS22](@cite), Corollaries 2.13 and 2.16. This function always uses
monic reduction, even when the flag's pushforward algorithm is `:gysin`.

# Examples
```jldoctest
julia> pt = abstract_point();

julia> Fl = flag_bundle(4*OO(pt), 2);

julia> Q = last(tautological_bundles(Fl));

julia> u, v = chern_class(Q, 1), chern_class(Q, 2);

julia> nf = flag_bundle_normal_form(Fl, u^4);

julia> nf == 2*lift(v)^2
true

julia> chow_ring(Fl)(nf) == u^4
true
```
"""
function flag_bundle_normal_form(Fl::AbstractVariety, a::MPolyDecRingOrQuoElem)
  D = _flag_bundle_data(Fl)
  Q = chow_ring(Fl)
  @req parent(a) === Q || parent(a) === base_ring(Q) "Class belongs to a different variety"
  f = _flag_relative_polynomial(D, a)
  return _flag_flatten(Fl, D, _flag_normal_form(_flag_engine(D, :monic), f))
end

@doc raw"""
    flag_bundle_basis(Fl::AbstractVariety)

Return a lazy iterator over a basis of `chow_ring(Fl)` as a free module over
the Chow ring of its base. The variety must have been constructed by
[`flag_bundle`](@ref). The iterator yields classes in `chow_ring(Fl)`; call
`collect` to materialize them. Neither construction nor pushforward needs
to enumerate this basis.

For successive ranks `r[i]`, its monomials are exactly those satisfying
`sum(exponents_in_block_i) <= sum(r[1:i-1])`. There are
`factorial(sum(r)) / prod(factorial.(r))` of them. The first monomial is `1`;
no other ordering of the output is part of the interface.

# Examples
```jldoctest
julia> pt = abstract_point();

julia> Fl = flag_bundle(4*OO(pt), 2);

julia> B = collect(flag_bundle_basis(Fl));

julia> length(B)
6

julia> p = structure_map(Fl);

julia> count(b -> pushforward(p, b) == 1, B)
1
```
"""
function flag_bundle_basis(Fl::AbstractVariety)
  D = _flag_bundle_data(Fl)
  R = base_ring(chow_ring(Fl))
  c = gens(R)
  return (chow_ring(Fl)(prod(c[j]^e[j] for j in eachindex(e)))
          for e in _FlagStandardExponents(D))
end
