# other

@doc raw"""
    adapted_string(b::KashiwaraCrystalElem, i::Vector{Int})

Return the adapted string of `b` with respect to the $i$-string `i`.
"""
function adapted_string(p::KashiwaraCrystalElem, rdec::Vector{Int})
  s = zero(rdec)
  b = deepcopy(p)
  for i in 1:length(rdec)
    s[i] = Base.eps(b, rdec[i]) # TODO: this will not work for non normal crystals
    KashiwaraCrystal.e!(b, rdec[i], s[i])
  end
  return s
end

###############################################################################
#
#   LS path model
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", P::LSPathModel)
  #@show_name(io, P)
  #@show_special(io, mime, P)
  #io = pretty(io)
  print(io, "LS path model for ")
  show(io, mime, P.wt)
end

function Base.show(io::IO, P::LSPathModel)
  #@show_name(io, R)
  #@show_special(io, R)
  if is_terse(io)
    print(io, "LS path model")
  else
    print(io, "LS path model for $(P.wt)")
  end
end

@doc raw"""
    ls_path_model(wt::WeightLatticeElem) -> LSPathModel

Return the LS path model for the dominant weight `wt`.
"""
function ls_path_model(wt::WeightLatticeElem)
  return LSPathModel(wt)
end

@doc raw"""
    ls_path_model(R::RootSystem, wt::Vector{<:IntegerUnion}) -> LSPathModel

Return the LS path model for the root system `R` and dominant weight `wt`.
"""
function ls_path_model(R::RootSystem, wt::Vector{<:IntegerUnion})
  return ls_path_model(WeightLatticeElem(R, wt))
end

function _extremal_weight(P::LSPathModel, w::WeylGroupElem)
  wt = get(P.ext, word(w), nothing)
  if isnothing(wt)
    # we need to make a copy, because w may be modified
    return get!(P.ext, deepcopy(word(w)), w * P.wt)
  end

  return wt
end

# TODO: the input here is counterintuitive
function (P::LSPathModel)(wv::Vector{WeylGroupElem}, tv::Vector{<:RationalUnion})
  return LSPathModelElem(P, [LSPathSegment(t, w) for (t, w) in Iterators.zip(tv, wv)])
end

function (P::LSPathModel)(p::YTPathModelElem)
  @req P.wt == parent(p).wt "dominant weights must match"

  pp = dominant_path(P)
  i = Int.(word(longest_element(weyl_group(root_system(P.wt)))))
  KashiwaraCrystal.f!(pp, i, adapted_string(p, i))
  return pp
end

function root_system(P::LSPathModel)
  return root_system(P.wt)
end

function Base.hash(s::LSPathSegment, h::UInt)
  b = 0x73a37b46bd4d49b4 % UInt
  h = hash(s.t, h)
  h = hash(s.w, h)

  return xor(h, b)
end

function Base.:(==)(s1::LSPathSegment, s2::LSPathSegment)
  return s1.t == s2.t && s1.w == s2.w
end

function parent(p::LSPathModelElem)
  return p.parent
end

function Base.:(==)(p::LSPathModelElem, q::LSPathModelElem)
  return parent(p) === parent(q) && p.s == q.s
end

function Base.:(<)(p::LSPathModelElem, q::LSPathModelElem)
  for i in 1:min(length(p.s), length(q.s))
    if p.s[i].w < q.s[i].w
      return true
    elseif p.s[i].w == q.s[i].w
      if p.s[i].t < q.s[i].t
        return true
      elseif p.s[i].t > q.s[i].t
        return false
      end
    else
      return false
    end
  end
  return false
end

function Base.deepcopy_internal(p::LSPathModelElem, dict::IdDict)
  if haskey(dict, p)
    return dict[p]
  end

  p2 = LSPathModelElem(parent(p), deepcopy_internal(p.s, dict))
  dict[p] = p2
  return p2
end

function Base.hash(p::LSPathModelElem, h::UInt)
  b = 0xd379e91be1f36479 % UInt
  h = hash(parent(p), h)
  h = hash(p.s, h)

  return xor(h, b)
end

function expressify(p::LSPathModelElem, s=:e; context=nothing)
  sum = Expr(:call, :+)
  for seg in p.s
    push!(
      sum.args,
      Expr(:call, :*, expressify(seg.t; context=context), "$s($(seg.w))"),
    )
  end
  return sum
end
@enable_all_show_via_expressify LSPathModelElem

@doc raw"""
    max(p::LSPathModelElem) -> WeylGroupElem

Return the maximal element in the support of `p`.
"""
function Base.max(p::LSPathModelElem)
  return p.s[1].w
end

@doc raw"""
    ls_sequence(p::LSPathModelElem) -> Tuple{Vector{WeylGroupElem}, Vector{QQFieldElem}}

Return the LS sequence for `p`.
"""
function ls_sequence(p::LSPathModelElem)
  return map(s -> s.w, p.s), [zero(QQ); accumulate((t, s) -> t + s.t, p.s; init=zero(QQ))]
end

@doc raw"""
    dominant_path(P::LSPathModel) -> LSPathModelElem

Return the dominant path for `P`.
"""
function dominant_path(P::LSPathModel)
  return LSPathModelElem(P, [LSPathSegment(one(QQ), one(weyl_group(root_system(P))))])
end

@doc raw"""
    halpha(p::LSPathModelElem, i::Int) -> Vector{QQFieldElem}

Return the rational turning points of the function $h_\alpha$ for `p` where $\alpha$ is the `i`th simple root.
"""
function halpha(p::LSPathModelElem, i::Int)
  h = sizehint!([zero(QQ)], length(p.s) + 1)
  for k in 1:length(p.s)
    push!(h, h[k] + p.s[k].t * _extremal_weight(parent(p), p.s[k].w)[i])
  end
  return h
end

# implement KashiwaraCrystalElem methods

function KashiwaraCrystal.e!(p::LSPathModelElem, i::Int)
  h = halpha(p, i)

  j = argmin(h) # first time the global min is assumed
  k = findprev(>=(h[j] + 1), h, j - 1)

  if isnothing(k)
    empty!(p.s)
    return p
  end

  of = (h[j] + 1 - h[k])//_extremal_weight(parent(p), p.s[k].w)[i]
  # of > 0, so we need to split this segment
  if !iszero(of)
    insert!(p.s, k, LSPathSegment(of, deepcopy(p.s[k].w)))
    sub!(p.s[k + 1].t, p.s[k + 1].t, of)
    k += 1
    j += 1
  end
  if j <= length(p.s) && length(p.s[j - 1].w) == length(p.s[j].w) + 1
    add!(p.s[j].t, p.s[j].t, p.s[j - 1].t)
    deleteat!(p.s, j - 1)
    j -= 1
  end

  for l in k:(j - 1)
    lmul!(p.s[l].w, i)
  end

  return p
end

function KashiwaraCrystal.f!(p::LSPathModelElem, i::Int)
  h = halpha(p, i)

  # find the last time the global minimum is assumed
  j = 1
  for l in 2:length(h)
    if h[l] <= h[j]
      j = l
    end
  end

  k = findnext(>=(h[j] + 1), h, j + 1)
  if isnothing(k)
    empty!(p.s)
    return p
  end

  of = (h[k] - h[j] - 1)//_extremal_weight(parent(p), p.s[k - 1].w)[i]
  # of > 0, we need to cut the segment
  if !iszero(of)
    insert!(p.s, k, LSPathSegment(of, deepcopy(p.s[k - 1].w)))
    sub!(p.s[k - 1].t, p.s[k - 1].t, of)
  end
  if j > 1 && length(p.s[j - 1].w) == length(p.s[j].w) + 1
    add!(p.s[j - 1].t, p.s[j - 1].t, p.s[j].t)
    deleteat!(p.s, j)
    k -= 1
  end

  for l in j:(k - 1)
    lmul!(p.s[l].w, i)
  end

  return p
end

function KashiwaraCrystal.eps(p::LSPathModelElem, i::Int)
  return -Int(minimum(halpha(p, i)))
end

function KashiwaraCrystal.phi(p::LSPathModelElem, i::Int)
  h = halpha(p, i)
  return Int(h[end]) - Int(minimum(h))
end

function KashiwaraCrystal.weight(p::LSPathModelElem)
  v = sum(s -> mul!(QQ.(coefficients(_extremal_weight(parent(p), s.w))), s.t), p.s)
  return WeightLatticeElem(root_system(parent(p)), ZZ.(v))
end

# LSPathModel specific methods

function iszero(p::LSPathModelElem)
  return isempty(p.s)
end

function (P::LSPathModel)(v::Vector{<:IntegerUnion})
  return P(WeightLatticeElem(root_system(P), v))
end

function (P::LSPathModel)(w::WeightLatticeElem)
  nf = inv(QQ.(cartan_matrix(root_system(P)))) * coefficients(P.wt - w)
  if !all(is_integral, nf)
    return LSPathModelElem[]
  end
  nf = [Int(nf[i]) for i in 1:length(nf)]

  W = weyl_group(root_system(P))
  w0 = longest_element(W)
  path = zeros(Int, length(w0))
  a = dominant_path(P)

  points = LSPathModelElem[]
  num = zeros(Int, length(nf)) # number of times acted with fi

  G = gens(W)

  i = length(w0)
  ret = false
  while true
    s = Int(w0[i])

    mi = min(phi(a, s), nf[s] - num[s])
    KashiwaraCrystal.f!(a, s, mi)
    num[s] += mi
    path[i] += mi

    ok = num == nf
    if i == 1 || ok
      if ok && !(a in points)
        push!(points, deepcopy(a))
      end

      KashiwaraCrystal.e!(a, s, path[i])
      num[s] -= path[i]
      path[i] = 0

      while path[i] == 0
        i += 1
        if i > length(w0)
          @goto done
        end
      end

      s = Int(w0[i])
      KashiwaraCrystal.e!(a, s)
      num[s] -= 1
      path[i] -= 1
    end
    i -= 1
  end
  @label done

  return points
end
