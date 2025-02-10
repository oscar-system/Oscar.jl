################################################################################
# Schur Polynomials
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
# by Tom Schmit and Ulrich Thiel; OSCAR-ified by Claudia He Yun and Matthias Zach.
################################################################################

@doc raw"""
    schur_polynomial([R::ZZMPolyRing], lambda::Partition, n::Int = length(lambda))

Return the Schur polynomial `s` of the partition `lambda` in `n` variables.

The ambient ring of `s` may optionally be supplied as a first argument.

# Examples
```jldoctest
julia> R, _ = ZZ[:a, :b, :c];

julia> schur_polynomial(R, partition([2, 1]))
a^2*b + a*b^2

julia> schur_polynomial(R, partition([2, 1]), 3)
a^2*b + a^2*c + a*b^2 + 2*a*b*c + a*c^2 + b^2*c + b*c^2

julia> schur_polynomial(partition([2]))
x1^2
```
"""
function schur_polynomial(lambda::Partition{T}, n::Int = length(lambda)) where T <: IntegerUnion
  @req n >= 0 "n >= 0 required"
  R, _ = polynomial_ring(ZZ, n, cached = false)
  return schur_polynomial(R, lambda, n)
end

function schur_polynomial(R::ZZMPolyRing, lambda::Partition{T}, n::Int = length(lambda)) where T <: IntegerUnion
  @req n >= 0 "n >= 0 required"
  if n == 0 || n < length(lambda)
    if isempty(lambda)
      return one(R)
    else
      return zero(R)
    end
  end

  @req n <= nvars(R) "n <= nvars(R) required"

  if n >= 10
    return schur_polynomial_combinat(R, lambda, n)
  end
  #decide which Algorithm to use if n < 10
  bo = sum(lambda) <= [140, 50, 17, 11, 10, 10, 11, 13, 14][n]
  if bo
    return schur_polynomial_combinat(R, lambda, n) #Combinatorial formula
  else
    return schur_polynomial_cbf(R, lambda, n) #Cauchy's bialternant formula
  end
end

#returning the schur polynomial in the first n generators of R using Cauchy's bialternant formula.
function schur_polynomial_cbf(R::ZZMPolyRing, lambda::Partition{T}, n::Int = length(lambda))  where T <: IntegerUnion
  # the corner cases are handled in the calling function
  @req n > 0 "number of variables must be > 0"
  @req n >= length(lambda) "number of variables must be at least the length of the partition"
  x = gens(R)[1:n]

  #=
  To calculate the determinant we use the Laplace expansion along the last row.
  Furthermore we use the fact that each column consist of the same variable with decreasing powers.

  This allows us to factorize by the smallest power, thus our last row of the
  minors is always 1, this reduces the amount of polynomial multiplications we
  have to do.

  To avoid calculating the same minors multiple times, we calculate them from
  1x1 up to nxn, always storing them in a Dictionary.
  Keep in mind that each minor consists of k columns and the top k rows.
  =#

  # initializing a few helpful Variables
  exponents = Int[getindex_safe(lambda, i) + n - i for i = 1:n] #the exponents from the Matrix read from top to bottom
  exp_incr = zeros(Int, n) #the increment with which exponents increase
  for i = 1:n - 1
    exp_incr[i] = exponents[i] - exponents[i + 1]
  end
  exp_incr[n] = exponents[n]

  factor = one(R)
  d = R()

  # Initialize Dictionaries (calculating all possible combinations of k = 1...n columns)
  sub_dets = [Dict{BitArray{1}, ZZMPolyRingElem}() for i = 1:n]
  # sub_dets[i] holds all the minors of size i, i.e.
  # sub_dets[2][[false, true, true, false, false]] is the minor of the 2x2 matrix
  # consisting of the first two rows intersected with the 2nd and 3rd columns.
  b = falses(n)
  ntrues = 0
  pointer = 1
  while true
    if b[pointer]
      while pointer <= n && b[pointer]
        b[pointer] = false
        ntrues -= 1
        pointer += 1
      end
      if pointer > n
        break
      end
    end
    b[pointer] = true
    ntrues += 1
    pointer = 1
    sub_dets[ntrues][copy(b)] = R()
  end

  #initialize sub_dets[1]
  for (columnview, ) in sub_dets[1]
    sub_dets[1][columnview] = x[columnview][1]^(exp_incr[1])
  end

  #calculate sub_dets[2:n] using Laplace extension
  exp = zeros(Int, ngens(R))
  for i = 2:n
    for (columnview, ) in sub_dets[i]
      d = R() #the alternating sum of minors
      s = ((i + n)%2 == 1 ? 1 : -1) #the alternating sign
      for j in (1:n)[columnview]
        columnview[j] = false
        d += s*sub_dets[i - 1][columnview]
        columnview[j] = true
        s *= -1 #change sign
      end

      #multiply by the factorized term
      factor = MPolyBuildCtx(R)
      exp = zeros(Int, ngens(R))
      for j in 1:n
        columnview[j] || continue
        exp[j] = exp_incr[i]
      end
      push_term!(factor, one(ZZ), exp)
      sub_dets[i][columnview] = mul!(d, d, finish(factor))
    end
  end
  sp = sub_dets[n][trues(n)]

  # divide by the product
  for i = 1:n - 1
    for j = i + 1:n
      sp = divexact(sp, x[i] - x[j])
    end
  end

  return sp
end

# returning the schur polynomial in the first k generators of R using the
# Combinatorial formula.
function schur_polynomial_combinat(R::ZZMPolyRing, lambda::Partition{T}, k::Int = length(lambda)) where T <: IntegerUnion
  if isempty(lambda)
    return one(R)
  end

  S = base_ring(R)
  sf = MPolyBuildCtx(R)

  #version of the function semistandard_tableaux(shape::Vector{T}, max_val = sum(shape))
  len = length(lambda)
  tab = [(fill(i, lambda[i])) for i = 1:len]
  m = len
  n = lambda[m]

  count = zeros(Int, nvars(R))
  valid = true
  while true
    count .= 0
    for i = 1:len
      for j = 1:lambda[i]
        if tab[i][j] <= k
          count[tab[i][j]] += 1
        else
          valid = false
          break
        end
      end
    end
    if valid
      push_term!(sf, S(1), count)
    end
    #raise one element by 1
    while !(tab[m][n] < k &&
            (n == lambda[m] || tab[m][n] < tab[m][n + 1]) &&
            (m == len || lambda[m + 1] < n || tab[m][n] + 1 < tab[m + 1][n]))
      if n > 1
        n -= 1
      elseif m > 1
        m -= 1
        n = lambda[m]
      else
        return finish(sf)
      end
    end

    tab[m][n] += 1

    #minimize trailing elements
    if n < lambda[m]
      i = m
      j = n + 1
    else
      i = m + 1
      j = 1
    end
    while (i <= len && j <= lambda[i])
      if i == 1
        tab[1][j] = tab[1][j - 1]
      elseif j == 1
        tab[i][1] = tab[i - 1][1] + 1
      else
        tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1)
      end
      if j < lambda[i]
        j += 1
      else
        j = 1
        i += 1
      end
    end
    m = len
    n = lambda[len]
  end #while true
end
