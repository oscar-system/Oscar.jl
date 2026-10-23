# Relative monic reduction for ordinary flag bundles.
# Grayson--Seceleanu--Stillman, GSS22, Corollaries 2.13, 2.16 and 5.3.
# Our convention is C_1 >> ... >> C_l, with weights 1,...,r_i in C_i.
# The coefficient ring is the ACTUAL Chow ring of the base, not its fraction
# field. In particular, no coefficient division is allowed in this file.

struct _FlagBundleReducer{N,A,S,P,C}
  coefficient_ring::A
  ring::S
  ranks::Vector{Int}
  starts::Vector{Int}
  ends::Vector{Int}
  stage_coefficients::Vector{Vector{P}}
  top_exponent::NTuple{N,Int}
  relative_dimension::Int
  rules::Dict{NTuple{N,Int},Vector{Pair{NTuple{N,Int},C}}}
  cache_lock::ReentrantLock
end

_flag_coefficient(a) = a
_flag_coefficient(a::MPolyQuoRingElem) = simplify(deepcopy(a))

function _FlagBundleReducer(A::Ring, ranks::Vector{Int}, chern::Vector;
                            names::Vector{String}=["c$i" for i in 1:sum(ranks)],
                            relative_ring=nothing)
  @req !isempty(ranks) && all(>(0), ranks) "Expected positive subquotient ranks"
  n = sum(ranks)
  @req length(chern) == n "Expected c_1(E),...,c_n(E), including trailing zeros"
  @req length(names) == n "Expected one variable name per Chern class"
  P = isnothing(relative_ring) ? polynomial_ring(A, names; cached=false)[1] : relative_ring
  @req base_ring(P) === A && ngens(P) == n "Unexpected relative polynomial ring"
  c = gens(P)
  T, C = elem_type(P), elem_type(A)
  ends = cumsum(ranks)
  starts = ends - ranks
  stages = Vector{Vector{T}}(undef, length(ranks))
  a = vcat(T[one(P)], T[P(_flag_coefficient(A(b))) for b in chern])
  for i in reverse(eachindex(ranks)) # find stage_coefficients
    stages[i] = a
    q, m = starts[i], ranks[i]
    # Monic polynomial division, written as truncated Chern convolution.
    quotient = T[one(P)]
    for k in 1:q
      t = a[k + 1]
      for j in 1:min(k, m)
        t -= c[q + j] * quotient[k - j + 1]
      end
      push!(quotient, t)
    end
    a = quotient
  end
  top = zeros(Int, n)
  for i in eachindex(ranks)
    top[ends[i]] = starts[i]
  end
  rules = Dict{NTuple{n,Int},Vector{Pair{NTuple{n,Int},C}}}()
  return _FlagBundleReducer{n,typeof(A),typeof(P),T,C}(
    A, P, copy(ranks), starts, ends, stages, Tuple(top),
    sum(ranks .* starts), rules, ReentrantLock(),
  )
end

function _flag_monomial(P, e)
  ctx = MPolyBuildCtx(P)
  push_term!(ctx, one(base_ring(P)), collect(e))
  return finish(ctx)
end

# Strict product weighted-degrevlex comparison; base coefficients are ignored.
function _flag_exp_lt(a, b, D::_FlagBundleReducer)
  for i in eachindex(D.ranks)
    s, m = D.starts[i], D.ranks[i]
    wa = sum(j * a[s + j] for j in 1:m)
    wb = sum(j * b[s + j] for j in 1:m)
    wa != wb && return wa < wb
    for j in m:-1:1
      a[s + j] != b[s + j] && return a[s + j] > b[s + j]
    end
  end
  return false
end

function _flag_weight(e, D::_FlagBundleReducer)
  return sum(sum(j * e[D.starts[i] + j] for j in 1:D.ranks[i])
             for i in eachindex(D.ranks))
end

# The determinant is indexed by the desired leading monomial, so no minors
# without a useful leading monomial are generated.
function _flag_relation(D::_FlagBundleReducer, i::Int, beta)
  P, q, m = D.ring, D.starts[i], D.ranks[i]
  c, a = gens(P), D.stage_coefficients[i]
  if q == 0
    j = findfirst(j -> beta[j] != 0, 1:m)
    return c[j] - a[j + 1]
  elseif m == 1
    # (-1)^(q+1) f_i(-c): the projective/complete-flag fast path.
    h = one(P)
    for k in 1:(q + 1)
      h = h * c[q + 1] + (-1)^k * a[k + 1]
    end
    return h
  end
  js = Int[]
  for j in 1:m
    append!(js, fill(j, beta[q + j]))
  end
  @assert length(js) == q + 1
  M = zero_matrix(P, q + 1, q + 1)
  for row in 1:(q + 1)
    h = js[row] + row - 1
    for col in 1:(q + 1)
      j = h - col + 1
      M[row, col] = j == 0 ? one(P) : 1 <= j <= m ? c[q + j] : zero(P)
    end
    M[row, 1] -= a[h + 1]
  end
  # Unlike pivoted or fraction-free elimination, this algorithm requires no
  # inverses or exact divisions, even when A has zero divisors.
  return Oscar.AbstractAlgebra.det_clow(M)
end

function _flag_rule(D::_FlagBundleReducer{N,A,S,P,C}, i::Int,
                    beta::NTuple{N,Int}) where {N,A,S,P,C}
  return lock(D.cache_lock) do
    get!(D.rules, beta) do
      tail = _flag_monomial(D.ring, beta) - _flag_relation(D, i, beta)
      result = Pair{NTuple{N,Int},C}[]
      for (a, e) in zip(coefficients(tail), exponent_vectors(tail))
        a = _flag_coefficient(a)
        iszero(a) && continue
        t = Tuple(e)::NTuple{N,Int}
        @assert _flag_exp_lt(t, beta, D)
        push!(result, t => a)
      end
      result
    end
  end
end

# A small max-heap keeps reduction deterministic and combines equal pending
# monomials before rewriting them. Cancelled entries may remain in the heap;
# the pending dictionary is the authority for whether an entry is still live.
function _flag_heap_push!(heap, e, D)
  push!(heap, e)
  k = length(heap)
  while k > 1
    p = k >> 1
    _flag_exp_lt(heap[p], heap[k], D) || break
    heap[p], heap[k] = heap[k], heap[p]
    k = p
  end
  return nothing
end

function _flag_heap_pop!(heap, D)
  result = heap[1]
  last_entry = pop!(heap)
  isempty(heap) && return result
  heap[1] = last_entry
  k = 1
  while 2k <= length(heap)
    j = 2k
    if j < length(heap) && _flag_exp_lt(heap[j], heap[j + 1], D)
      j += 1
    end
    _flag_exp_lt(heap[k], heap[j], D) || break
    heap[k], heap[j] = heap[j], heap[k]
    k = j
  end
  return result
end

function _flag_enqueue!(pending, heap, D, e, a, top_only)
  top_only && _flag_weight(e, D) < D.relative_dimension && return nothing
  iszero(a) && return nothing
  if haskey(pending, e)
    b = pending[e] + a
    if iszero(b)
      delete!(pending, e)
    else
      pending[e] = b
    end
  else
    pending[e] = a
    _flag_heap_push!(heap, e, D)
  end
  return nothing
end

function _flag_divisor(D::_FlagBundleReducer{N}, e::NTuple{N,Int}) where {N}
  for i in eachindex(D.ranks)
    q, m = D.starts[i], D.ranks[i]
    sum(e[q + j] for j in 1:m) <= q && continue
    beta = zeros(Int, N)
    remaining = q + 1
    for j in 1:m
      beta[q + j] = min(e[q + j], remaining)
      remaining -= beta[q + j]
      remaining == 0 && break
    end
    return i, Tuple(beta)
  end
  return 0, e
end

function _flag_reduce(D::_FlagBundleReducer{N,A,S,P,C}, f::P;
                      top_only::Bool=false) where {N,A,S,P,C}
  @req parent(f) === D.ring "Polynomial belongs to a different relative ring"
  pending = Dict{NTuple{N,Int},C}()
  answer = Dict{NTuple{N,Int},C}()
  heap = NTuple{N,Int}[]
  for (a, e) in zip(coefficients(f), AbstractAlgebra.exponent_vectors(f))
    _flag_enqueue!(pending, heap, D, Tuple(e), _flag_coefficient(a), top_only)
  end
  while !isempty(heap)
    e = _flag_heap_pop!(heap, D)
    haskey(pending, e) || continue
    a = pop!(pending, e)
    i, beta = _flag_divisor(D, e)
    if i == 0
      answer[e] = _flag_coefficient(a)
      continue
    end
    for (v, b) in _flag_rule(D, i, beta)
      t = ntuple(j -> e[j] - beta[j] + v[j], N)
      _flag_enqueue!(pending, heap, D, t, a * b, top_only)
    end
  end
  return answer
end

function _flag_normal_form(D::_FlagBundleReducer, f)
  ctx = MPolyBuildCtx(D.ring)
  for (e, a) in _flag_reduce(D, f)
    push_term!(ctx, a, collect(e))
  end
  return finish(ctx)
end

function _flag_pushforward(D::_FlagBundleReducer, f)
  # Each rule weakly decreases the total Chern-variable weight. A term of
  # weight < dim(F/X) therefore cannot contribute to the top coefficient.
  nf = _flag_reduce(D, f; top_only=true)
  return get(nf, D.top_exponent, zero(D.coefficient_ring))
end

struct _FlagStandardExponents{N}
  starts::Vector{Int}
  ends::Vector{Int}
end

function _FlagStandardExponents(D::_FlagBundleReducer{N}) where {N}
  return _FlagStandardExponents{N}(D.starts, D.ends)
end

Base.IteratorSize(::Type{<:_FlagStandardExponents}) = Base.SizeUnknown()
Base.eltype(::Type{_FlagStandardExponents{N}}) where {N} = NTuple{N,Int}

function Base.iterate(I::_FlagStandardExponents{N}) where {N}
  e = zeros(Int, N)
  return Tuple(e), e
end

function Base.iterate(I::_FlagStandardExponents, state::Vector{Int})
  e = copy(state)
  for i in reverse(eachindex(I.starts))
    s, t = I.starts[i], I.ends[i]
    total = sum(@view e[(s + 1):t])
    for j in t:-1:(s + 1)
      if total < s
        e[j] += 1
        return Tuple(e), e
      end
      total -= e[j]
      e[j] = 0
    end
  end
  return nothing
end
