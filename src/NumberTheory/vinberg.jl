###################################################################
# Vinberg's algorithm
###################################################################

@doc raw"""
    _check_v0(Q::ZZMatrix, v0::ZZMatrix) -> ZZMatrix

Check if the inner product of the given `v0` is positive or generate such a `v0` if `v0==0`
"""
function _check_v0(Q::ZZMatrix, v0::ZZMatrix)
  if iszero(v0) # generate `v0` by finding the negative eigenvalue `neg_ev` of `Q`. We set `v0` as a primitive generator of the eigenspace of `neg_ev`
    V = quadratic_space(QQ,Q)
    diag, trafo = diagonal_with_transform(V)
    @req count(is_positive, diag) == 1 "Q must be of signature (1,n)"
    @req count(is_zero, diag) == 0 "Q must be of signature (1,n)"
    @req count(is_negative, diag) > 0 "Q must be of signature (1,n)"
    i = findfirst(is_positive, diag)
    v0 = trafo[i:i,:]
    v0 = ZZ.(denominator(v0)*v0)
    return _rescale_primitive(v0)
  else
    @req (v0*Q*transpose(v0))[1, 1] > 0 "v0 has non positive inner product"
    v0 = _rescale_primitive(v0)
    return v0
  end
end

@doc raw"""
    _all_root_lengths(Q::ZZMatrix) -> Vector{ZZRingElem}

If the user does not want any specific root lengths, we take all of them.

For `Q` the matrix representing the hyperbolic reflection lattice and `l` a length to check,
it is necessary that `l` divides $2 i$ for `l` being a possible root length,  
with `i` being the level of `Q`, i.e. the last invariant (biggest entry) of the Smith normal form of `Q`.
Proof: `r` being a root implies that $\frac{2.r.L}{r^2}$ is a subset of $\Z$
     `r` being primitive implies that $2.i.\Z$ is a subset of $2.r.L$,
     which is a subset of $(r^2).\Z$, which implies that $r^2$ divides $2.i$
"""
function _all_root_lengths(Q::ZZMatrix)
  l = nrows(Q)
  S = snf(Q)
  return (-1) * divisors(2 * S[l, l])
end

@doc raw"""
    _check_root_lengths(Q::ZZMatrix, root_lengths::Vector{ZZRingElem}) -> Vector{ZZRingElem}

Check whether the given lengths are possible root lengths.
Unnecessary roots are sorted out and reported if this is wanted.
Return only the possible root lengths.
"""
function _check_root_lengths(Q::ZZMatrix, root_lengths::Vector{ZZRingElem})
  l = nrows(Q)
  S = snf(Q)
  possible_root_lengths = divisors(2 * S[l, l])
  if isempty(root_lengths)
    return possible_root_lengths
  end
  result = ZZRingElem[]
  for r in root_lengths
    if abs(r) in possible_root_lengths
      push!(result, r)
    else
      @vprintln :Vinberg 1 "$r is not a possible root length"
    end
  end
  @req !isempty(result) "No possible root lengths found"
  return result
end

@doc raw"""
    _distance_indices(upper_bound, root_lengths::Vector{ZZRingElem}) -> Vector{Tuple{Int, ZZRingElem}}

Returns all possible tuples $(n, l)$ with `n` a natural number and `l` contained in `root_lengths`,
sorted by increase of the value $\frac{n^2}{l}$, stopping by `upper_bound`.
"""
function _distance_indices(upper_bound, root_lengths::Vector{ZZRingElem})

  result = Tuple{Int,ZZRingElem}[]

  for l in root_lengths
    # consider for the `upper_bound` $u$ $\frac{n^2}{l} < u \implies n < \sqrt{l*u}$
    x = isqrt(abs(l * upper_bound))
    for n in 1:x
      push!(result, (n, l))
    end
  end

  sort!(result; by=x -> (x[1]^2) // abs(x[2]))
  return result
end

@doc raw"""
    _rescale_primitive(v::ZZMatrix) -> ZZMatrix

Check if `v` is primitive, which is equivalent to checking whether the greatest common divisor `g` of the entries of `v`
is equal to 1. If `v` is not primitive, we divide `v` by `g`.
"""
function _rescale_primitive(v::ZZMatrix)
  g = reduce(gcd, v)
  return isone(g) ? v : v / g
end

@doc raw"""
    _crystallographic_condition(Q::ZZMatrix, v::ZZMatrix)

Check if the reflection by `v` preserves the lattice, i.e. check if for `v` a row vector
$\frac{2.v.Q}{v^2}$ is an integer matrix. 
"""
function _crystallographic_condition(Q::ZZMatrix, v::ZZMatrix)
  A = Q * transpose(v)
  b = (v * A)[1, 1]
  return all(iszero(mod(2*x, b)) for x in A)
end
  
@doc raw"""
    _crystallographic_condition(Q::ZZMatrix, v::ZZMatrix)

Check if the reflection by `v` preserves the lattice, i.e. check if for `v` a row vector
$\frac{2 v Q}{vQv^t}$ is an integer matrix
where `k` is $vQv^t$.
"""
@inline function _crystallographic_condition(Qv::ZZMatrix, k::ZZRingElem)
  return all(iszero(mod(2*x, k)) for x in Qv)
end

@doc raw"""
    _check_coorientation(Q::ZZMatrix, roots::Vector{ZZMatrix}, v::ZZMatrix, v0::ZZMatrix)

First check whether `v` has non-obtuse angles with all roots `r` already found, by checking that $v.r \geq 0$.
"""
function _has_non_obtuse_angles(Q::ZZMatrix, roots::Vector{ZZMatrix}, v::ZZMatrix)
  Qv = Q * transpose(v)
  for r in roots
    if (r * Qv)[1, 1] < 0
      return false
    end
  end
  return true
end

# non-allocating version
function _has_non_obtuse_angles!(tmp::ZZMatrix, Qv::ZZMatrix, roots::Vector{ZZMatrix})
  for r in roots
    #tmp = r*Qv
    mul!(tmp, r, Qv)
    if !is_zero_entry(tmp,1, 1) && !is_positive_entry(tmp, 1, 1)
      return false
    end
  end
  return true
end

@doc raw"""
    _distance_0(Q::ZZMatrix, v0::ZZMatrix, root_lengths::Vector{ZZRingElem}, direction_vector::ZZMatrix) -> Vector{ZZMatrix}

Return the fundamental roots which are orthogonal to `v0`.

After gathering all primitive vectors `v` which are orthogonal to `v0` with length contained in `root_lengths`
satisfying the crystallographic condition, we check for those `v` in order of increase of the value $\frac{(v.v_1)^2}{v^2}$ 
if they are roots, e.g. if they also satisfy
- $(v.v_1) > 0$ for the `direction_vector` $v_1$
- $(v.r) \geq 0$ for all roots `r` already found
"""
function _distance_0(Q::ZZMatrix, v0::ZZMatrix, root_lengths::Vector{ZZRingElem}, direction_vector::ZZMatrix)
  roots = ZZMatrix[]
  possible_vec = ZZMatrix[]
  # _short_vectors_gram is more efficient than short_vectors_affine
  # Thus we have to make some preparations
  bm = kernel(Q * transpose(v0); side=:left)
  QI = bm * Q * transpose(bm)
  QI = QI[1, 1] > 0 ? QI : -QI
  for d in root_lengths # gather all vectors which are orthogonal on `v0` with length contained in `root_lengths`
    d = abs(d)
    V = Hecke._short_vectors_gram(Hecke.LatEnumCtx, map_entries(QQ, QI), d, d, ZZRingElem) # `QI` positive 
    for (v_, _) in V
      v = matrix(ZZ, 1, nrows(bm), v_) * bm
      v = _rescale_primitive(v)
      if _crystallographic_condition(Q, v)
        push!(possible_vec, v)
      end
    end
  end
  is_empty(possible_vec) && return roots
  # we check if we can take the given direction vector or generate one if there is no one given
  if iszero(direction_vector)
    v1 = _generate_direction_vector(Q, possible_vec)
  else
    v1 = _check_direction_vector(Q, possible_vec, direction_vector)
  end
  # Now we want to sort the vectors by increase of the value $\frac{(v.v_1)^2}{v^2}$
  possible_vec_sort = Tuple{ZZMatrix,RationalUnion,RationalUnion}[]
  for v in possible_vec
    vQ = v * Q
    n = (vQ*v1)[1, 1]
    k = (vQ*transpose(v))[1, 1]
    push!(possible_vec_sort, (v, n, k))
  end
  sort!(possible_vec_sort; by=x -> (x[2]^2) // abs(x[3]))

  # We execute the algorithm with replacing `v0` by `v1`
  for (v, n) in possible_vec_sort
    if n < 0
      v = -v
    end
    vQ = v * Q
    if all((vQ*transpose(w))[1, 1] >= 0 for w in roots)
      push!(roots, v)
    end
  end
  return roots
end

@doc raw"""
    _generate_direction_vector(Q::ZZMatrix, possible_vec::Vector{ZZMatrix}) -> ZZMatrix

Return a `direction_vector` $v_1$ satisfying $(\~{v}, v_1) \neq 0$ 
for all `possible_vec` $\~{v}$ with $(v_0, \~{v}) = 0$
"""
function _generate_direction_vector(Q::ZZMatrix, possible_vec::Vector{ZZMatrix})
  l = ncols(Q)
  signal = 1 # initializing a stopping condition
  while signal < 100000000
    if l < 4
      v1_ = rand(-30:30, l)
    elseif l < 10
      v1_ = rand(-25:25, l)
    else
      v1_ = rand(-20:20, l)  # for higher `l` it is not necessary/efficient to choose a big random range
    end
    v1 = matrix(QQ, l, 1, v1_)
    @hassert :Vinberg 1 (transpose(v1)*Q*v1)[1, 1] < 0
    if all((v*Q*v1)[1, 1] != 0 for v in possible_vec)
      # To run the algorithm it suffices to require that v*Q*v1 != 0, but we get a random (right) solution. 
      # This is the case because it is a random choice in which direction the fundamental cone 
      # is built every time we use the algorithm.
      # To get always the same solution one should supply a `direction_vector` `v1`
      @vprintln :Vinberg 2 "direction vector v1 = $v1"
      return v1
    end
    signal += 1
  end
  @req signal < 10000000 "Choose another v0" # break the algorithm if no `v1` was found
end

@doc raw"""
    _check_direction_vector(Q::ZZMatrix, possible_vec::Vector{ZZMatrix}, direction_vector::ZZMatrix) -> ZZMatrix

We check if the given `direction_vector` $v_1$ satisfies $(\~{v}, v_1) \neq 0$ 
for all `possible_vec` $\~{v}$ with $(v_0, \~{v}) = 0$
"""
function _check_direction_vector(Q::ZZMatrix, possible_vec::Vector{ZZMatrix}, direction_vector::ZZMatrix)
  v1 = transpose(direction_vector)
  @req all((v*Q*v1)[1, 1] != 0 for v in possible_vec) "The direction vector cannot be orthogonal to one of the possible roots"
  return v1
end

@doc raw"""
    vinberg_algorithm(Q::ZZMatrix, upper_bound; v0::ZZMatrix, root_lengths::Vector{ZZRingElem}, direction_vector::ZZMatrix) -> Vector{ZZMatrix}

Return the fundamental roots `r` of a given hyperbolic reflection lattice with standard basis represented by its corresponding Gram matrix `Q` with 
squared length contained in `root_lengths` and by increasing order of the value $\frac{r.v_0)^2}{r^2}$, stopping at `upper_bound`.
If `root_lengths` is not defined it takes all possible values of $r^2$.
If `v0` lies on a root hyperplane and if there is no given `direction_vector` it is a random choice which reflection chamber
next to `v0` will be computed.

# Arguments
- `Q`: symmetric $\Z$ matrix of signature $(1, n)$ -- the corresponding Gram matrix
- `upper_bound`: the upper bound of the value $\frac{(r.v_0^2}{r^2}$
- `v0`: primitive row vector with $v_0^2 > 0$
- `root_lengths`: the possible integer values of $r^2$
- `direction_vector`: row vector `v1` with $v_0.v_1 = 0$ and $v.v_1 \neq 0$ for all possible roots `v` with $v.v_0 = 0$
"""
function vinberg_algorithm(Q::ZZMatrix, upper_bound; v0=ZZ[0;]::ZZMatrix, root_lengths=ZZRingElem[]::Vector{ZZRingElem}, direction_vector=ZZ[0;]::ZZMatrix, divisibilities::Union{Nothing,Dict{ZZRingElem,Vector{ZZRingElem}}}=nothing)
  @req is_symmetric(Q) "Matrix is not symmetric"
  
  v0 = _check_v0(Q, v0)
  if isempty(root_lengths)
    real_root_lengths = _all_root_lengths(Q)
  else
    real_root_lengths = _check_root_lengths(Q, root_lengths) # sort out impossible lengths
  end
  iteration = _distance_indices(upper_bound, real_root_lengths) # find the right order to iterate through $v.v_0$ and $v^2$
  roots = _distance_0(Q, v0, real_root_lengths, direction_vector) # special case $v.v_0 = 0$
  Qv = zero_matrix(ZZ, ncols(Q), 1)
  tmp2 = zero_matrix(ZZ, 1, 1)
  for (n, k) in iteration # search for vectors which solve $n = v.v_0$ and $k = v^2$
    @vprintln :Vinberg 1 "computing roots of squared length v^2=$(k) and v.v0 = $(n)"
    possible_Vec = short_vectors_affine(Q, v0, QQ(n), k)
    for v in possible_Vec
      if !isone(reduce(gcd, v))
        # v must be primitive.
        continue 
      end
      mul!(Qv, Q, transpose(v))
      #Qv = Q*transpose(v)
      if !(divisibilities isa Nothing)
        # filter for divisibilities
        reduce(gcd, Qv) in divisibilities[k] || continue
      end
      if _crystallographic_condition(Qv, k) && _has_non_obtuse_angles!(tmp2, Qv, roots)
        push!(roots, v)
      end
    end
  end
  return roots
end

@doc raw"""
    vinberg_algorithm(S::ZZLat, upper_bound; v0::ZZMatrix, root_lengths::Vector{ZZRingElem}, direction_vector::ZZMatrix) -> Vector{ZZMatrix}

Return the fundamental roots `r` of a given hyperbolic reflection lattice `S` with standard basis with squared length contained 
in `root_lengths` and by increasing order of the value $\frac{r.v_0)^2}{r^2}$, stopping at `upper_bound`.
If `root_lengths` is not defined it takes all possible values of $r^2$.
If `v0` lies on a root hyperplane and if there is no given `direction_vector`, 
then it is a random choice which reflection chamber next to `v0` will be computed.

# Arguments
- `S`: a hyperbolic $\Z$-lattice of signature $(1,0,n)$.
- `upper_bound`: the upper bound of the value $\frac{(r.v_0^2}{r^2}$
- `v0`: primitive row vector with $v_0^2 > 0$ given w.r.t. the ambient space
- `root_lengths`: the possible integer values of $r^2$
- `direction_vector`: row vector `v1` with $v_0.v_1 = 0$ and $v.v_1 \neq 0$ for all possible roots `v` with $v.v_0 = 0$, given w.r.t. the ambient space
- `divisibilities`: a dictionary; The keys are the root lengths and the values are the divisibilities for the given root length. If given requires that a fundamental root $r$ has one of the specified divisibilities.
"""
function vinberg_algorithm(S::ZZLat, upper_bound; v0=QQ[0;]::QQMatrix, root_lengths=ZZRingElem[]::Vector{ZZRingElem}, direction_vector=QQ[0;]::QQMatrix, divisibilities::Union{Nothing,Dict{ZZRingElem,Vector{ZZRingElem}}}=nothing)
  Q = gram_matrix(S)
  if !iszero(v0)
    v0 = solve(basis_matrix(S),v0; side=:left)
  end
  if !iszero(direction_vector)
    direction_vector = solve(basis_matrix(S),direction_vector; side=:left)
  end
  _v0 = ZZ.(v0)
  _check_direction_vector = ZZ.(direction_vector)
  roots = vinberg_algorithm(Q, upper_bound; v0=_v0, root_lengths=root_lengths, direction_vector=_direction_vector, divisibilities)
  B = basis_matrix(S)
  return [r * B for r in roots]
end

