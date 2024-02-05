# TODO: ascending_partitions should be renamed to ascending_compositions
# according to U. Thiel, see the discussion
# https://github.com/oscar-system/Oscar.jl/pull/3159#discussion_r1446198249
@doc raw"""
    ascending_partitions(n::IntegerUnion; algorithm::Symbol=:ks)

Instead of encoding a partition of an integer ``n ≥ 0`` as a *descending*
sequence (which is our convention), one can also encode it as an *ascending*
sequence. In the papers Kelleher & O'Sullivan (2014) and Merca (2012) it is
said that generating the list of all ascending partitions is more efficient
than generating descending ones. To test this, we have implemented the
algorithms given in the papers:
1. `:ks` (*default*) is the algorithm "AccelAsc" (Algorithm 4.1) in [KO14](@cite).
2. `:m` is Algorithm 6 in [Mer12](@cite). This is actually similar to `:ks`.

The ascending partitions are stored here as arrays and are not of type
`Partition` since the latter are descending by our convention. We are using `:ks`
as default since it looks slicker.

# Comparison

We don't see a significant speed difference to the descending encoding:
```julia-repl
julia> @btime partitions(Int8(90));
3.376 s (56634200 allocations: 6.24 GiB)

julia> @btime ascending_partitions(Int8(90), algorithm=:ks);
3.395 s (56634200 allocations: 6.24 GiB)

julia> @btime ascending_partitions(Int8(90), algorithm=:m);
3.451 s (56634200 allocations: 6.24 GiB)
```
"""
function ascending_partitions(n::IntegerUnion; algorithm::Symbol=:ks)
  if algorithm === :ks
    return ascending_partitions_kelleher_osullivan(n)
  elseif algorithm === :m
    return ascending_partitions_merca(n)
  else
    error("algorithm must be either :ks or :m")
  end
end

function ascending_partitions_kelleher_osullivan(n::T) where {T <: IntegerUnion}
  #Argument checking
  @req n >= 0 "n >= 0 required"

  # Some trivial cases
  if n==0
    return Vector{T}[ [] ]
  elseif n==1
    return Vector{T}[ [1] ]
  end

  # Now, the algorithm starts
  P = Vector{T}[]    #this will be the array of all partitions
  a = zeros(T, n)
  k = 2
  y = n-1
  while k != 1
    k -= 1
    x = a[k] + 1
    while 2*x <= y
      a[k] = x
      y -= x
      k += 1
    end
    l = k + 1
    while x <= y
      a[k] = x
      a[l] = y
      push!(P, a[1:l])
      x += 1
      y -= 1
    end
    y += x - 1
    a[k] = y + 1
    push!(P, a[1:k])
  end
  return P
end

function ascending_partitions_merca(n::T) where {T <: IntegerUnion}
  #Argument checking
  @req n >= 0 "n >= 0 required"

  # Some trivial cases
  if n==0
    return Vector{T}[ [] ]
  elseif n==1
    return Vector{T}[ [1] ]
  end

  P = Vector{T}[]    #this will be the array of all partitions
  a = zeros(T, n)
  k = 1
  x = 1
  y = n-1
  while true
    while 3*x <= y
      a[k] = x
      y = y-x
      k += 1
    end
    t = k + 1
    u = k + 2
    while 2*x <= y
      a[k] = x
      a[t] = x
      a[u] = y - x
      push!(P, a[1:u])
      p = x + 1
      q = y - p
      while p <= q
        a[t] = p
        a[u] = q
        push!(P, a[1:u])
        p += 1
        q -= 1
      end
      a[t] = y
      push!(P, a[1:t])
      x += 1
      y -= 1
    end
    while x<=y
      a[k] = x
      a[t] = y
      push!(P, a[1:t])
      x += 1
      y -= 1
    end
    y += x-1
    a[k] = y+1
    push!(P, a[1:k])
    k -= 1

    # I think there's a mistake in the publication
    # because here k could be zero and then we access
    # a[k].
    # That's why I do a while true and check k > 0 here.
    if k == 0
      break
    else
      x = a[k] + 1
    end
  end
  return P
end
