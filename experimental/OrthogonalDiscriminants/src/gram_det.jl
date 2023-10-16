
# The following is a translation to Julia of GAP code written by Frank Lübeck.
# The idea is that for a partition $\lambda$ of $n$,
# the Specht module $L(\lambda)$ of the symmetric group $S_n$
# has a basis that is parameterized by standard tableaux of shape $\lambda$.
# The determinant of the Gram matrix with respect to this basis can be
# computed combinatorically using a formula of Jantzen and Schaper
# (5.33 in the book by Andrew Mathas).

# note that partitions are transposed, compared to Mathas' book, 5.33
function schaper_decomposition(lambda::Partition{T}, p::Int) where T <: IntegerUnion
  n = sum(lambda)
  sgn = (-1)^div(n, 2)
  mu = conjugate(lambda)

  # hook lengths, col-wise
  hooklens = []
  mux = collect(mu) - collect(0:(length(mu)-1))
  lax = collect(lambda) - collect(1:length(lambda))
  for i in 1:length(mu)
    push!(hooklens, mux[i] .+ lax[1:mu[i]])
  end

  # beta numbers for mu, see 5.20
  bmu = collect(lambda)
  for i in (length(bmu)+1):n
    push!(bmu, 0)
  end
  bmu = bmu + collect(n-1:-1:0)

  # collect coeffs and partitions, formula 5.33
  coeffs = ZZRingElem[]
  parts = Partition{T}[]
  for i in 1:length(mu)
    # nu_p(h_x,i)
    nups = []
    for a in hooklens[i]
      b = 0
      while mod(a, p) == 0
        b = b+1;
        a = div(a, p)
      end
      push!(nups, b)
    end
    for j in 1:(mu[i]-1)
      for k in (j+1):mu[i]
        d = nups[j] - nups[k]
        # only term if resulting list gives beta numbers of a partition, 
        # see before 5.32
        if d != 0
          a = bmu[j] + hooklens[i][k]
          if !(a in bmu)
            beta = copy(bmu)
            beta[j] = a
            a = beta[k] - hooklens[i][k]
            if a >= 0 && !(a in beta)
              beta[k] = a
              b = sign(Perm(sortperm(beta))) * sgn
              push!(coeffs, b*d)
              beta = sort!(beta)[n:-1:1] - collect((n-1):-1:0)
              pos = 1
              l = 1
              while l <= length(beta) && beta[l] != 0
                l = l + 1
              end
              beta = beta[1:(l-1)]
              # note again, that partitions are transposed here
              push!(parts, partition(beta))
            end
          end
        end
      end
    end
  end
  return coeffs, parts
end

@doc raw"""
    dimension_specht_module(mu::Partition{T}) where T <: IntegerUnion -> ZZRingElem

Return the dimension of the Specht module for `mu`.

# Examples
```jldoctest
julia> print([dimension_specht_module(p) for p in partitions(4)])
ZZRingElem[1, 3, 2, 3, 1]
```
"""
function dimension_specht_module(mu::Partition{T}) where T <: IntegerUnion
  la = conjugate(mu)
  d = ZZ(1)
  f = ZZ(1)
  m = ZZ(1)
  mux = collect(mu) - collect(0:(length(mu)-1))
  lax = collect(la) - collect(1:length(la))
  for i in 1:length(mu)
    for j in 1:mu[i]
      d = d * (mux[i] + lax[j])
      f = f * m
      m = m + 1
    end
  end
  return div(f, d)
end


@doc raw"""
    gram_determinant_specht_module(mu::Partition{T}) where T <: IntegerUnion

Return the determinant of the Gram matrix for the Specht module for `mu`,
in factorized collected form.

# Examples
```jldoctest
julia> print(gram_determinant_specht_module(partition([4, 3, 2, 1])))
Vector{ZZRingElem}[[3, 1152], [5, 768], [7, 384]]
```
"""
function gram_determinant_specht_module(mu::Partition{T}) where T <: IntegerUnion
  n = sum(mu)
  det = Vector{ZZRingElem}[]
  for p in 2:n
    if is_prime(p)
      cfs, parts = schaper_decomposition(mu, p)
      if length(cfs) > 0
        push!(det, [p, dot(cfs, map(dimension_specht_module, parts))])
      end
    end
  end
  return det
end
