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

Return all possible tuples $(n, l)$ with `n` a natural number and `l` contained in `root_lengths`,
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

  # The code for `short_vectors_affine` will, by default, turn `Q` and `v0`
  # into `QQMatrix`. This is due to the use of the function `_convert_type`
  # in Hecke, which requires this (in this context, since `v0` might not be a
  # vector in the lattice with Gram matrix `Q`, but only a vector of its
  # rational span). Since we call several times the function
  # `short_vectors_affine`, once per iteration, we do this conversion into
  # `QQMatrix` prior iterating on the list `iteration`.

  # The output of `short_vectors_affine` are `QQMatrix` with integer entries:
  # we thus recover honest `ZZMatrix` by taking `numerator(v)`.
  _Q = map_entries(QQ, Q)
  _v0 = map_entries(QQ, v0)
  for (n, k) in iteration # search for vectors which solve $n = v.v_0$ and $k = v^2$
    @vprintln :Vinberg 1 "computing roots of squared length v^2=$(k) and v.v0 = $(n)"
    possible_Vec = short_vectors_affine(_Q, _v0, QQ(n), QQ(k))
    for v in possible_Vec
      vZZ = numerator(v)
      if !isone(reduce(gcd, vZZ))
        # v must be primitive.
        continue 
      end
      mul!(Qv, Q, transpose(vZZ))
      #Qv = Q*transpose(v)
      if !(divisibilities isa Nothing)
        # filter for divisibilities
        reduce(gcd, Qv) in divisibilities[k] || continue
      end
      if _crystallographic_condition(Qv, k) && _has_non_obtuse_angles!(tmp2, Qv, roots)
        push!(roots, vZZ)
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


function _get_tau(f::ZZMatrix, Qb)
  tau = zero(Qb)
  for lambda in eigenvalues(Qb, f)
    if abs(lambda)>tau
      tau = abs(lambda)
    end
  end
  return tau
end

function _get_bilinearform(Lf::ZZLatWithIsom, Qb)
  gram_in_Qb = change_base_ring(Qb, gram_matrix(Lf))
  #return (a, b)->inner_product(ambient_space(L), a,b)[1,1] -> it doesn't go well as there is no convienient way to convert ZZLat in QQ field to the Qb field
  return function (a,b)
    tmp_a = change_base_ring(Qb, a)
    tmp_b = transpose(change_base_ring(Qb, b))
    mul!(tmp_b, mul!(tmp_a, tmp_a, gram_in_Qb), tmp_b)
    return tmp_b[1]
  end
end

function _get_C0(Lf::ZZLatWithIsom, tau::QQBarFieldElem)
  charPolyF = characteristic_polynomial(Lf)
  x = gen(parent(charPolyF))
  (n, remainder) = remove(charPolyF, x-1) # remove all x-1 factors
  return div(remainder, minpoly(parent(charPolyF), tau)) #remove Salem polynomial of salem number tau
end

function _get_Cfancy(Lf::ZZLatWithIsom, C0)
  return short_vectors_iterator(lattice(kernel_lattice(Lf, C0)), 2 , 2)
end


@doc raw"""
    _get_h(L::ZZLat, v, w, Qb, bi_form) -> ZZMatrix

Return a vector $h$ in the lattice `L`, such that $bi_form(h,h)>0$. 

Vectors 'v','w' are eigenvectors of isometry associated with L based on salem number $tau$ and $tau^(-1)$

'Qb' is a current field, where all calculations happen
"""
function _get_h(L::ZZLat, v, w, bi_form)
  l = number_of_rows(basis_matrix(L))
  r = rank(L)
  if r<4
    k = 10
  elseif  r<6
    k = 6
  elseif r<8
    k = 4
  elseif r<10
    k = 2
  else
    k = 1
  end
  h0 = v+w
  h0 = h0/bi_form(h0,h0)
  n = 1
  i = 1
  z = matrix(ZZ,1,l,rand(-k:k, l)) # random small vector in lattice basis
  n_h0 = map(x->floor(ZZRingElem, x) , h0*n)
  h = z+n_h0# in lattice basis
  h = map(x->floor(ZZRingElem, x) , h) # in lattice basis and rounded
  while bi_form(h,h) <= 0 || n == 10000  # block on 10000, because for a bigger n time consumption is way too big
    i+=1
    if i == 100
      n+=1
      i = 1
      n_h0 = map(x->floor(ZZRingElem, x) , h0*n) 
    end
    for x in 1:l
      z[1, x] = rand(-k:k)
    end
    add!(h, z, n_h0)
    #h = (z+n_h0) #in lattice basis
  end
  if n == 10000 throw(OverflowError("The upper limitation on h calculation is reached. Bigger h will not be able to process in reasonable time")) end
  return h::ZZMatrix
end

function _get_R(L, h::ZZMatrix)
  return short_vectors_affine(change_base_ring(ZZ, gram_matrix(L)),h,0,-2) #better to use iterator, but there is no equivalent function that returns iterator
end

@doc raw"""
    _process_finite_sets_of_h(h, f::QQMatrix, v, w, bi_form, L::ZZLat) -> Tuple{Bool, ZZMatrix}

The function calculates all $r$, that can be obstructing roots of `L`.
Calculation is made one by one and then the $r$ is checked
$r$ is based on pairs of integer $(a,b)$ of $-2x^2+2y^2+2aby>=x(a^2+b^2)$ with $a>0$, $b<0$, $x>0$,$y>0$
See Steps 7,8,9 of Algorithm 5.8 in [OY20](@cite)

Return a tuple of a boolean that represents if isometry f of the lattice is positive and ZZMatrix,
that represents an obstructing root
"""
function _process_finite_sets_of_h(h::ZZMatrix, f::ZZMatrix, v, w, bi_form, L::ZZLat)
  x = bi_form(h, h)
  y = bi_form(h, h*f)
  z = y^2-x^2
  if (z<0) return (true, zero(h)) end #then discriminant for a will be <0 and there is no root (a,b)=> no obstructing roots => positive
  b_min = -isqrt(round(ZZRingElem, 2*z/x, RoundDown))
  h_new = zero(h)
  tmp = zero(h)
  for b = b_min:-1
    a_max = round(ZZRingElem, (b*y+isqrt(round(ZZRingElem, (z*(b^2+2*x)), RoundDown)))/x, RoundDown)
    for a = 1:a_max
      add!(h_new,neg!(mul!(h_new,b,h),mul!(tmp,a,mul!(tmp,h,f))))
      #h_new = -b*h +a*h*f
      Rh = _get_R(L, h_new)
      for r in Rh
        if _check_R(r, v, w, bi_form) return false, r end
      end
    end
  end
  return true, zero(h)
end

function _check_R(r, v, w, bi_form)
  return bi_form(r, v)*bi_form(r, w) < 0
end

@doc raw"""
    isometry_is_positive(Lf::ZZLatWithIsom, h::Union{QQMatrix, Nothing} = nothing) -> Tuple{Bool, ZZMatrix}

Return whether the isometry of `Lf` is positive and an obstructing root if it exists. 

The isometry is called positive if it preserves a connected component of the set $\{x \in L \mid x^2>0, \wedge \forall r \in \{r \in L | r^2=-2\}: x.r\neq 0\}$. 

This implements see Algorithm 5.8 of [OY20](@cite). 

# Arguments
- `Lf` -- a hyperbolic lattice with an isometry of positive spectral radius. 
- `h` -- a vector of positive square given with respect to the ambient space of `Lf`. If no such vector is given, `h` will be calculated by function inside based on a random vector in `Lf`.

# Output
For positive isometries the second return value is a vector of zeros of the same degree as `h`.
"""

function isometry_is_positive(Lf::ZZLatWithIsom, h::Union{ZZMatrix, Nothing} = nothing)
  Qb = algebraic_closure(QQ);
  f = change_base_ring(ZZ, isometry(Lf))
  L = lattice(Lf)
  tau = _get_tau(f, Qb)
  bi_form = _get_bilinearform(Lf, Qb)

  # step 1
  C0 = _get_C0(Lf, tau)
  # step 2 - Check if C0 has obstructing roots => not positive
  if !isone(C0)
    Cfancy = _get_Cfancy(Lf, C0)
    first_element = iterate(Cfancy)
    if first_element !== nothing
      return (false, first_element[1][1])
    end
  end
  # step 3 - Prepare eigenvectors from tau and tau inverse 
  v = eigenspace(f, tau)
  w = eigenspace(f, tau^(-1))

  if bi_form(v,w)<0
    v = -v
  end
  # step 4 - Get the first h value if there is no in arguments
  if h === nothing
    h = _get_h(L,v,w, bi_form)
  end
  # step 5 - Get the R set based on current h value
  Rh = _get_R(L, h)
  # step 6 - Check all of the entries of R if there exists an obstructing root => positive
  for r in Rh
    if _check_R(r, v, w, bi_form) return false, r end
  end
  # step 6,7,8 are combined to process them iteratively
  return _process_finite_sets_of_h(h, f, v, w, bi_form, L)
end