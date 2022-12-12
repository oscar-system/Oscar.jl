# Generators of the orthogonal group of a torsion quadratic form.


#     r"""
#     Return the kernel of the orthogonal group under mod `p` reduction.
#
#     Let `L` be a `p`-adic lattice. For a positive integer `a` define
#     `O(L/p^aL)` as the image of `O(L)` in `GL(L/p^aL)`.
#     This function returns generators of the kernel of `O(L/p^aL)`--> `O(L/pL)`.
#
#     INPUT:
#
#     - ``G`` -- a non-singular, symmetric, integral, block diagonal
#       `p`-adic matrix with modular blocks of descending valuation
#     - ``a``-- an integer
#
#     OUTPUT:
#
#     - a list of `p`-adic matrices
#
#     EXAMPLES::
#
#         sage:
#     """
function _mod_p_to_a_kernel(G::Union{nmod_mat, fmpz_mod_mat}, a, p)
  n = ncols(G)
  R = base_ring(G)
  E = identity_matrix(R, n)
  ind, val = _block_indices_vals(G, p)
  # add a virtual even block at the end
  push!(ind, n+1)
  push!(val, val[end] - 1)

  # block diagonal contribution
  gens = typeof(G)[]
  for k in 1:length(ind)-1
    i1 = ind[k]
    i2 = ind[k+1]
    Gk = divexact(G[i1:i2-1,i1:i2-1], p^val[k])
    Gk_inv = inv(Gk)
    ni = i2 - i1
    Ek = identity_matrix(R, ni)
    Zk = zero_matrix(R, ni, ni)
    if p == 2 && a == 1
      genskmodp = _solve_X_ker(Zk, zeros(R, ni), diagonal(Gk_inv))
      gensk = [change_base_ring(R,lift(g)) for g in genskmodp]
    else
      # basis for the space of anti-symmetric matrices
      gensk = typeof(G)[]
      for i in 1:ni
        for j in 1:i-1
          gk = deepcopy(Zk)
          gk[i,j] = 1
          gk[j,i] = -1
          push!(gensk, gk)
        end
      end
    end
    for h in gensk
      g = deepcopy(E)
      g[i1:i2-1, i1:i2-1] = Ek + p^a*h*Gk_inv
      push!(gens, g)
    end
  end
  # generators below the block diagonal.
  for i in 1:n
    for j in 1:i
      g = deepcopy(E)
      g[i,j] = p^a
      flag = true
      for k in 1:(length(ind)-1)
        # exclude the block diagonal entries
        if ind[k] <= i < ind[k+1] && ind[k] <= j < ind[k+1]
          flag = false
          break
        end
      end
      if flag
        push!(gens, g)
      end
    end
  end
  return gens
end

#     r"""
#     Return the transformation to normal form.
#
#     INPUT:
#
#     - ``G`` -- `p`-adic symmetric matrix.
#
#     OUTPUT:
#
#     - ``D`` -- the normal form
#     - ``B`` -- a transformation matrix
#
function _normal(G::Union{fmpz_mod_mat, nmod_mat}, p)
  if p == 2
    D1, B1 = _jordan_2_adic(G)
    D2, B2 = _normalize(D1, p)
    D3, B3 = _two_adic_normal_forms(D2, p)
    B = B3 * B2 * B1
    return D3, B
  else
    D1, B1 = _jordan_odd_adic(G, p)
    D2, B2 = _normalize(D1, p)
    B = B2 * B1
    return D2, B
  end
end

#     Return generators of the orthogonal group modulo `p` for odd `p`.
#
#     `G` is the gram matrix of the bilinear form.
function _gens_form(G::Union{fmpz_mod_mat, nmod_mat}, form_constructor, p)
    R = base_ring(G)
    if ncols(G) == 0
        return typeof(G)[]
    end
    F = GF(p, 1)
    Gf = change_base_ring(F, lift(G))
    FF = GF(p)
    q = form_constructor(Gf)
    Oq = isometry_group(q)
    gensOq = gens(Oq)
    gensOq = [change_base_ring(FF,matrix(g)) for g in gensOq]
    gensOq = [change_base_ring(R,lift(g)) for g in gensOq]

    return gensOq
end

function _gens_qf(G::Union{fmpz_mod_mat, nmod_mat}, p)
  Gq = deepcopy(G)
  if p == 2
    r = ncols(G)
    for i in 1:r
      for j in 1:i-1
        Gq[i,j] = 0
      end
      Gq[i,i] = divexact(Gq[i,i],2)
    end
  end
  gensOq =  _gens_form(Gq, quadratic_form, p)
  return gensOq
end

_gens_af(G::Union{fmpz_mod_mat, nmod_mat}, p) = _gens_form(G, alternating_form, p)

_orthogonal_grp_gens_odd(G::Union{fmpz_mod_mat, nmod_mat},p) = _gens_form(G::Union{fmpz_mod_mat, nmod_mat}, quadratic_form, p)

#     Return generators of the orthogonal group of the bilinear form modulo `2`.
#
#     Here `G` is the gram matrix of the bilinear form.
#
#     INPUT:
#
#     - ``G`` -- a modular, symmetric `2`-adic matrix in homogeneous normal form
#
#     OUTPUT:
#
#     - a list of `2`-adic matrices well defined modulo `2`
#
function _orthogonal_gens_bilinear(G::Union{fmpz_mod_mat, nmod_mat})
  r = ncols(G)
  R = base_ring(G)
  # corner cases
  if r <= 1
    gens_1 = typeof(G)[]
  elseif r == 2 && mod(lift(G[r, r]), 2) == 1
    return [matrix(R, 2, 2, [0, 1, 1, 0])]
  # odd cases
  elseif r % 2 == 1
    # the space of points of 0 square is non degenerate && preserved
    # there the bilinear form is alternating && the group is a
    # sympletic group
    # the orthogonal complement has rank one && is invariant as well
    gens_1 = _gens_af(G[1:r-1, 1:r-1], 2)
    E1 = identity_matrix(R, 1)
    gens_1 = [diagonal_matrix([g, E1]) for g in gens_1]
  elseif _val(G[r, r], 2) == 0
    # the space of points of 0 square is degenerate with kernel of dim 1
    # the quotient by the kernel is again sympletic
    # but we obtain several possible lifts
    gens_1 = _gens_af(G[1:end-2,1:end-2], 2)
    gens_1 = [diagonal_matrix(g, identity_matrix(R, 2)) for g in gens_1]
    E = identity_matrix(R, r)
    for i in 1:r-2
      g = deepcopy(E)
      g[i, r-1] = 1
      g[i, r] = 1
      g[r-1:end, 1:r-1] = transpose(g[1:r-1, r-1:end]) * G[1:r-1, 1:r-1]
      push!(gens_1, g)
    end
    g = deepcopy(E)
    g[end-1:end, end-1:end] = matrix(R, 2, 2, [0, 1, 1, 0])
    push!(gens_1, g)
  else
    # even case
    # just a symplectic group nothing special
    gens_1 = _gens_af(G, 2)
  end
  # check that generators are isometries
  for g in gens_1
    err = g*G*transpose(g) - G
    @assert _min_val(err, 2) >= 1
  end
  return gens_1
end

#     Return generators of the orthogonal group of the quadratic form modulo `4`.
#
#     Here `G`represents a quadratic form on an $\FF_2$ vector space with
#     values in `\Zmod{4}` given by `q(x) = xGx^T \mod 4`.
#
#     INPUT:
#
#     - ``G`` -- a homogeneous, symmetric `2`-adic matrix in normal form
#
#     OUTPUT:
#
#     - a list of matrices. Generators of the orthogonal group modulo `2`.
function _orthogonal_grp_quadratic(G::Union{fmpz_mod_mat, nmod_mat})
  r = ncols(G)
  R = base_ring(G)
  # corner cases
  if r == 0
    return typeof(G)[]
  end
  v = _min_val(G, 2)
  G = divexact(G, 2^v)
  if r <= 1
    gens1 = typeof(G)[]
  elseif r == 2
    if _val(G[end,end], 2) == 0
      if mod(lift(G[end,end] + G[end-1,end-1]), 4) == 0
        gens1 = typeof(G)[]
      else
        gens1 = [matrix(R, 2, 2, [0, 1, 1, 0])]
      end
    elseif _val(G[end,end], 2) == 1
      gens1 = [matrix(R, 2, 2, [0, 1, 1, 0]),
              matrix(R,2,2, [0, 1, 1, 1])]
    else
      gens1 = [matrix(R, 2, 2, [0, 1, 1, 0])]
    end
  elseif r % 2 == 1   # usual cases
    # an odd case
    # the space of points of even square is preserved
    # there the quadratic form is classical
    # so is its orthogonal complement -> we get an invariant vector
    gens1 = _gens_qf(G[1:end-1,1:end-1], 2)
    # the invariant vector is the last one
    # continue the isometry as the identity there
    E1 = identity_matrix(R, 1)
    gens1 = [diagonal_matrix([g, E1]) for g in gens1]
  elseif _val(G[end,end],2) == 0  # now r % 2 == 0
    # odd case
    # the space of points of even square is preserved
    # but it is degenerate
    # modulo the degenerate part the form is
    if mod(lift(G[end,end] + G[end-1,end-1]), 4) == 0
      # the degenerate part v has q(v) = 0
      # thus q is preserved and we get
      # a classical form
      gens1 = _orthogonal_grp_quadratic(G[1:end-2,1:end-2])
      gens1 = [_lift(G, g, 0) for g in gens1]
      Id = identity_matrix(R, r - 2)
      z = zeros(R, r-2)
      A = typeof(z)[]
      for i in 1:r-2
        z = zeros(R, r-2)
        z[i] = R(1)
        push!(gens1,_lift(G, Id, z))
      end
    else
      # here the degenerate part v has q(v) = 2
      # thus only the bilinear form descends to the quotient
      # && we get a sympletic group
      @assert mod(lift(G[end,end] + G[end-1,end-1]), 4) == 2
      gens1 = _gens_af(G[1:end-2,1:end-2], 2)
      gens1 = [_lift(G, g, 0) for g in gens1]
      push!(gens1,_lift(G, identity_matrix(R, r-2), 1))
    end
  else
    # even case
    # the groups are classical we just need to bring them to the right basis
    gens1 = _gens_qf(G, 2)
  end
  # check that generators are isometries
  for g in gens1
    err = g*G*transpose(g)-G
    @assert _min_val(err, 2) >= 1
    @assert all(_val(d,2)>=2 for d in diagonal(err))
  end
  return gens1
end

#     Return a list of the lifts of f.
#
#     INPUT:
#
#     - ``q`` -- odd of scale `1` in homogeneous normal form
#     - ``f`` -- the `n-2 \times n-2` matrix to be lifted
#     - ``a`` -- ``0``, ``1`` or a row matrix depending on the case
#
#     OUTPUT:
#
#     - ``g`` -- the lift of ``f`` as determined by ``a``
#
#     EXAMPLES::
#
#         sage: from sage.groups.fqf_orthogonal_gens import _lift
#         sage: R = Zp(2,type='fixed-mod',prec=2,print_mode='terse', show_prec=false, print_pos=false)
#         sage: U = matrix(R,2,2,[0,1,1,0])
#         sage: W0 = matrix(R,2,2,[1,0,0,3])
#         sage: W1 = matrix(R,2,2,[1,0,0,1])
#         sage: q0 = matrix.block_diagonal([U,W0])
#         sage: g0 = matrix(R,2,2,[0,1,1,0])
#         sage: g0l = _lift(q0,g0,vector([1,1]))
#         sage: g0l
#         [0 1 1 1]
#         [1 0 1 1]
#         [3 3 0 3]
#         [1 1 1 2]
#
#     The essential property of the lifts is that is preserves the bilinear form
#     `\mod 2` && the quadratic `\mod 4`::
#
#         sage: (g0l'*q0*g0l - q0)
#         [0 0 0 0]
#         [0 0 0 0]
#         [0 0 0 0]
#         [0 0 0 0]
#
#     The parameter ``a`` is ``g0l1[-2,:-2]``::
#
#         sage: _lift(q0,g0,vector([0,1]))
#         [0 1 0 0]
#         [1 0 1 1]
#         [0 3 1 0]
#         [0 1 0 1]
#
#     In the second case one can lift any form preserving the bilinear form on the
#     small part. This is the whole symplectic group::
#
#         sage: q1 = matrix.block_diagonal([U,W1])
#         sage: g1 = matrix(R,2,2,[1,1,1,0])
#         sage: g1l = _lift(q1,g1,1)
#         sage: (g1l'*q1*g1l - q1)
#         [0 0 2 0]
#         [0 0 0 0]
#         [2 0 0 2]
#         [0 0 2 0]
#     """
function _lift(q::T, f::T, a) where T <: Union{fmpz_mod_mat, nmod_mat}
  R = base_ring(q)
  if mod(lift(q[end,end-1]), 2) != 0
    error("The form must be odd.")
  end
  # notation
  g = diagonal_matrix([f, identity_matrix(R,2)])
  b = identity_matrix(R,ncols(f)+2)
  b[end-1,end] = 1
  qb = b * q * transpose(b)
  G = qb[1:end-2,1:end-2]
  fG = f * G
  fGinv = inv(fG)

  if mod(lift(q[end-1,end-1] + q[end,end]), 4) == 2
    g[end, end-1] = a
    g[1:end-2,end-1] = diagonal(divexact(G - f*G*transpose(f),2))
    g[end,1:end-2] = transpose(fGinv * g[1:end-2,end-1])
  else
    if a == 0
      a = zeros(R, ncols(q)-2)
    end
    g[1:end-2,end-1] = a
    g[end,1:end-2] = transpose(fGinv * g[1:end-2,end-1])
    t = (g[end,1:end-2]*G*transpose(g[end,1:end-2]))[1,1]
    g[end,end-1] = divexact(t-mod(lift(t), 2), 2)
  end
  err = g*qb*transpose(g)-qb
  # check that lifting succeeded
  @assert _min_val(err, 2) >= 1
  @assert all(_val(d,2)>=2 for d in diagonal(err))
  return inv(b) * g * b
end



#     r"""
#     Return generators of the orthogonal groups of ``G`` modulo `p`.
#
#     Let `V = \Zp^n` && `b: V \times V \rightarrow \Zp` be the bilinear form
#     `b(x,y)= x^T G y`. This method computes generators of the image of
#     the orthogonal group `O(V,b)` under
#
#     ..MATH:
#
#         O(V,b) \rightarrow GL(V/pV)
#
#     INPUT::
#
#         -``G`` -- gram matrix of a non-degenerate, symmetric, bilinear
#           `p`-adic form.
#
#     OUTPUT::
#
#         - generators modulo `p`
#
#     EXAMPLES::
#
#         sage: from sage.groups.fqf_orthogonal_gens import _gens_mod_p
#         sage: R = Zp(3, type='fixed-mod', prec=10, print_mode='terse', show_prec=false, print_pos=false)
#         sage: G = matrix.diagonal(R, [3*1, 3*1])
#         sage: _gens_mod_p(G)
#         [
#         [0 1]  [0 1]
#         [2 0], [1 0]
#         ]
#         sage: G = matrix.diagonal(R, [1, 3, 3, 9, 2*27])
#         sage: _gens_mod_p(G)
#         [
#         [-1  0  0  0  0]  [1 0 0 0 0]  [1 0 0 0 0]  [ 1  0  0  0  0]
#         [ 0  1  0  0  0]  [0 0 1 0 0]  [0 0 1 0 0]  [ 0  1  0  0  0]
#         [ 0  0  1  0  0]  [0 2 0 0 0]  [0 1 0 0 0]  [ 0  0  1  0  0]
#         [ 0  0  0  1  0]  [0 0 0 1 0]  [0 0 0 1 0]  [ 0  0  0 -1  0]
#         [ 0  0  0  0  1], [0 0 0 0 1], [0 0 0 0 1], [ 0  0  0  0  1],
#         <BLANKLINE>
#         [ 1  0  0  0  0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]
#         [ 0  1  0  0  0]  [1 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]
#         [ 0  0  1  0  0]  [0 0 1 0 0]  [1 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]
#         [ 0  0  0  1  0]  [0 0 0 1 0]  [0 0 0 1 0]  [1 0 0 1 0]  [0 1 0 1 0]
#         [ 0  0  0  0 -1], [0 0 0 0 1], [0 0 0 0 1], [0 0 0 0 1], [0 0 0 0 1],
#         <BLANKLINE>
#         [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]
#         [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]
#         [0 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]
#         [0 0 1 1 0]  [0 0 0 1 0]  [0 0 0 1 0]  [0 0 0 1 0]  [0 0 0 1 0]
#         [0 0 0 0 1], [1 0 0 0 1], [0 1 0 0 1], [0 0 1 0 1], [0 0 0 1 1]
#         ]
function _gens_mod_p(G::Union{fmpz_mod_mat, nmod_mat}, p)
  n = ncols(G)
  R = base_ring(G)
  E = identity_matrix(R, n)
  indices, valuations = _block_indices_vals(G, p)
  push!(indices, n+1)
  gens1 = typeof(G)[]
  for k in 1:length(indices)-1
    i1 = indices[k]
    i2 = indices[k + 1]-1
    Gi = divexact(G[i1:i2, i1:i2], R(p)^valuations[k])
    gens_homog = _orthogonal_grp_gens_odd(Gi, p)
    for f in gens_homog
      g = deepcopy(E)
      g[i1:i2, i1:i2] = f
      push!(gens1, g)
    end
  end
  # generators below the block diagonal.
  for i in 1:n
    for j in 1:i-1
      g = deepcopy(E)
      g[i,j] = 1
      flag = true
      for k in 1:(length(indices)-1)
        if indices[k] <= i < indices[k+1] && indices[k] <= j < indices[k+1]
          g[i,j] = 0
          flag = false
          break
        end
      end
      if flag
        push!(gens1, g)
      end
    end
  end
  return gens1
end

# r"""
# Return the generators of the orthogonal groups of ``G`` modulo `2`.
#
# Let `V = \FF_2^n` && `b: V \times V \rightarrow \FF_2` be the bilinear form
# `b(x,y)= x^T G y`. Compute generators of `O(V,b)`.
#
# INPUT::
#
# -``G`` -- gram matrix of a non-degenerate, symmetric, bilinear
#   `2` form over `\FF_2` in normal form.
#
# """
function _gens_mod_2(G::Union{fmpz_mod_mat, nmod_mat})
    n = ncols(G)
    R = base_ring(G)
    E = identity_matrix(R, n)
    p = ZZ(2)
    ind0, val0 = _block_indices_vals(G, 2)
    par0 = Int[]
    push!(ind0, n+1)
    for k in 1:length(ind0)-1
        i = ind0[k+1] - 1    # last index of block i
        if _val(G[i,i], 2)>  val0[k]
          pa = 0
        else
          pa = 1
        end
        push!(par0, pa)
    end
    ind = Tuple{Int,Int}[]  # indices of the jordan blocks [(i1,i2),...]
    val = Int[]  # valuations of the jordan blocks
    par = Int[]  # parities of the jordan blocks
    k = 0
    for v in val0[1]+2:-1:val0[end]-1
        i = findfirst(x-> x==v, val0)
        if i isa Nothing
            push!(ind,(ind0[k+1],ind0[k+1]-1))
            push!(val,v)
            push!(par, 0)
        else
            k = i
            push!(ind,(ind0[k],ind0[k+1]-1))
            push!(val,v)
            push!(par,par0[k])
        end
    end
    val[end] = 0
    gens = typeof(G)[]
    for k in 3:length(val)-1
        if par[k + 1] == 1
            i1 = ind[k][1]
            i2 = ind[k][2]
            i3 = ind[k + 1][2]
            Gk = divexact(G[i1:i3,i1:i3],R(2)^val[k+1])
            gens_k = _gens_pair(Gk, 1 + i2-i1, false)
        elseif par[k-1] == 1
            i1 = ind[k-1][1]
            i2 = ind[k][1]
            i3 = ind[k][2]
            Gk = divexact(G[i1:i3,i1:i3],R(2)^val[k])
            gens_k = _gens_pair(Gk, i2-i1, true)
        else
            i1 = ind[k][1]
            i3 = ind[k][2]
            Gk = divexact(G[i1:i3,i1:i3],R(2)^val[k])
            gens_k = _orthogonal_grp_quadratic(Gk)
        end

        for h in gens_k
            g = deepcopy(E)
            g[i1:i3,i1:i3] = h
            push!(gens, g)
        end
    end
    # a change in convention
    trafo = deepcopy(E)
    for k in 2:length(ind)-1
        if par[k] == 1 && mod(ind[k][2]-ind[k][1],2) == 1
           i = ind[k][2]
           trafo[i-1:i,i-1:i] = matrix(R,2,2,[1,1,0,1])
        end
    end
    trafoinv = inv(trafo)
    Gt = trafo * G * transpose(trafo)
    # ker
    # row wise starting with the last row
    for k in length(ind)-1:-1:3
        pa = par[k - 2:k]
        i = ind[k][2]
        Gi = divexact(Gt[1:i,1:i], R(2)^val[k])
        gensK = _ker_gens(Gi, ind[k-1][1]-1, ind[k-1][2], pa)
        E = identity_matrix(R, n-i)
        gensK = [diagonal_matrix([g, E]) for g in gensK]
        gensK = [trafoinv * g * trafo for g in gensK]
        append!(gens, gensK)
    end
    return gens
end

function _gens_pair(G::Union{fmpz_mod_mat, nmod_mat}, k, on_second)
  gen = typeof(G)[]
  R = base_ring(G)
  n = ncols(G)
  G1 = G[1:k,1:k]  # 2^1 - modular
  G1inv = inv(divexact(G1, 2))
  G2 = G[k+1:end,k+1:end]  # 2^0 - modular
  E = identity_matrix(R, n)
  if on_second
    for f in _orthogonal_gens_bilinear(G2)
      a = diagonal(divexact(f*G2*transpose(f)-G2, 2))
      g = deepcopy(E)
      g[k+1:end, k+1:end] = f
      g[k+1:end, k] = a
      push!(gen, g)
    end
  else
    for f in _orthogonal_gens_bilinear(divexact(G1, 2))
      a = [divexact(x,2) for x in diagonal(f*G1*transpose(f) - G1)]
      g = deepcopy(E)
      g[1:k,1:k] = f
      g[1:k,end] = a
      g[k+1:end,1:k] = - G2 * transpose(divexact(g[1:k,k+1:end],2)) * transpose(inv(f)) * G1inv
      push!(gen,g)
    end
  end
  return gen
end

# Generators for the kernel of
# O(L/2^2L) ---> O(L/2 L)
# in the case of 3 consecutive jordan blocks
# with weird python indexing
function _ker_gens(G::Union{fmpz_mod_mat, nmod_mat}, i1, i2, parity)
  n = nrows(G)
  R = base_ring(G)
  E = identity_matrix(R, n)
  gens = typeof(G)[]
  e = n - 1
  if parity[3]==1 && mod(n - i2, 2)==0
    e = n - 2
  end
  for i in i2+1:n
    for j in 1:i2
      g = deepcopy(E)
      if parity == [0,0,0] || parity == [1,0,0]
        g[i,j] = 1
      elseif parity == [0,0,1]
        if (j-1 < i1) || (i-1 != e)
          g[i,j] = 1
        end
      elseif parity == [0,1,0] || parity == [1,1,0]
        if !(j-1 == i2 - 1)
          g[i,j] = 1
        end
      elseif parity == [0,1,1]
        if !((j-1 == i2 - 1) || (i-1 == e && j-1 >= i1))
          g[i,j] = 1
        end
      elseif parity == [1,0,1]
        if !(i-1 == e && j-1 == i1 - 1)
          g[i,j] = 1
        end
      elseif parity == [1,1,1]
        if !((i-1 == e && j-1 == i1 - 1) || (j-1 == i2 - 1))
          g[i,j] = 1
        end
      end
      if parity[1]==1 && parity[3]==1
        # compensate
        g[e+1, i1] = divexact((g[e+1,i1+1:i2]*G[i1+1:i2,i1+1:i2]*transpose(g[e+1,i1+1:i2]))[1,1],4)
        # the second row depends on the third
        g[i1+1:i2,i1] = - divexact((G[i1+1:i2,i1+1:i2]* transpose(g[i2+1:end,i1+1:i2]) * inv(G[i2+1:end,i2+1:end]))[:,end], 2)
      end
      if g[i,j] == 1   # no need to append the identity
        push!(gens, g)
      end
    end
  end
  return gens
end


#     r"""
#     Return generators.
#
#     INPUT:
#
#     - ``G`` -- a non-singular, integral, symmetric `p`-adic matrix in
#       descending normal form
#     - ``b`` -- a positive integer
#
#     EXAMPLES::
#
#         sage: from sage.groups.fqf_orthogonal_gens import _gens
#         sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=false, print_pos=false)
#         sage: U = matrix(R, 2, 2, [0, 1, 1, 0])
#         sage: V = matrix(R, 2, 2, [2, 1, 1, 2])
#         sage: W0 = matrix(R, 2, 2, [1, 0, 0, 3])
#         sage: W1 = matrix(R, 2, 2, [1, 0, 0, 1])
#         sage: G = matrix.block_diagonal([2*U, V])
#         sage: gens = _gens(G, 2)
#         sage: G = matrix.block_diagonal([2*U, W1])
#         sage: gens = _gens(G, 2)
#         sage: G = matrix.diagonal(R, [2, 1])
#         sage: gens = _gens(G, 2)
#         sage: G = matrix.block_diagonal([2*V, V])
#         sage: gens = _gens(G, 2)
#         sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
#         sage: A = AbelianGroupGap([2, 2, 4, 4])
#         sage: aut = A.aut()
#         sage: gens = [aut(g) for g in gens]
#         sage: Oq = aut.subgroup(gens)
#         sage: Oq.order()
#         1152
#
#     TESTS::
#
#         sage: R = Zp(3, type='fixed-mod', prec=10, print_mode='terse', show_prec=false, print_pos=false)
#         sage: G = matrix.diagonal(R,[3*1])
#         sage: gens = _gens(G,1)
#         sage: G = matrix.diagonal(R,[3*1,3*1])
#         sage: gens = _gens(G,2)
#         sage: G = matrix.diagonal(R,[2*27,9,3*1,3*1,1])
#         sage: gens = _gens(G,4)
#     """
function _gens(G::Union{fmpz_mod_mat, nmod_mat}, b, p)
    k = 1
    gens = typeof(G)[]
    while k <= b
        gen = _mod_p_to_a_kernel(G, k, p)
        gen = [hensel_qf(G, f, k+1, b+1, p) for f in gen]
        k *= 2
        append!(gens,gen)
    end
    if p == 2
        gen = _gens_mod_2(G)
    else
        gen = _gens_mod_p(G, p)
    end
    for g in gen
      push!(gens, hensel_qf(G,g,1,b,p))
    end
    return gens
end

# r"""
# Return generators of the orthogonal group of ``T``.
#
# INPUT:
#
# - ``T`` -- torsion orthogonal module in normal form.
# - ``deg`` -- avoids an infinite recursion
#
# OUTPUT:
#
# - a list of matrices -- the generators
#
# EXAMPLES::
#
#     sage: from sage.groups.fqf_orthogonal_gens import _compute_gens
#     sage: T = TorsionQuadraticForm(matrix.diagonal([2/3,2/3]))
#     sage: _compute_gens(T)
#     [
#     [  1 726]  [0 1]  [0 1]
#     [  3   1], [2 0], [1 0]
#     ]
# """
function _compute_gens(T::TorQuadMod)
  T.is_normal || error("T must be normal")

  # corner case
  invs = elementary_divisors(T)
  if length(invs) == 0
    return fmpz_mat[]
  end

  # normal form gens for the different primes
  blocks = []
  gensT_orders = [order(t) for t in gens(T)]
  n = length(gens(T))
  P = prime_divisors(exponent(T))
  for p in P
    indices = Int[]
    for k in 1:length(gensT_orders)
        if mod(gensT_orders[k], p) == 0
            push!(indices, k)
        end
    end
    push!(blocks, [p, indices])
  end

  # compute generators of the orthogonal groups
  gensG = fmpz_mat[]
  for (p, indices) in blocks
    # compute the generators of the p-primary part
    # the whole group is the direct product of the p-primary parts
    q_p = gram_matrix_quadratic(T)[indices, indices]
    b = valuation(invs[end], p)
    R = ResidueRing(ZZ,fmpz(p)^(b+5))
    G_p = change_base_ring(ZZ, q_p*p^b)
    G_p = change_base_ring(R, G_p)
    if p != 2
      # make sure each homogeneous block of G_p stays in normal form
      r = divexact(G_p, R(2))
    end
    # the generators in matrix form
    gens_mat = _gens(G_p, b, p)
    # extend as identity on the orthogonal complement
    E1 = identity_matrix(ZZ, indices[1]-1)
    E2 = identity_matrix(ZZ, n - indices[end])
    for g in gens_mat
      g = change_base_ring(ZZ, g)
      g = block_diagonal_matrix([E1, g, E2])
      push!(gensG, g)
    end
  end
  return gensG
end

# We carry the case when T is not semi-regular, but does split its
# radical quadratic.
#
# TODO: We have way too many generators, in both functions. There might be a way
# to reduce the number of generators along the way!
function _compute_gens_split_degenerate(T::TorQuadMod)
  # we need a splitting T = rd + N, where rd is the radical
  # quadratic and N is a normal form.
  #
  # If it is the case, since T splits as the sum of its primary parts,
  # then each primary part splits its radical quadratic too.
  # So, we split T into primary parts, since they don't "talk"
  # to each others via isometries, we compute the orthogonal group
  # of each primary part and we then "glue" the orthogonal groups
  # using diagonal matrices on an isometric module to T which is nice enough.
  @assert has_complement(radical_quadratic(T)[2])[1]

  # if T is "primary degenerate" then we can send it to the other function since
  # T has only one primary part.
  if is_prime_power_with_data(elementary_divisors(T)[end])[1]
    return _compute_gens_split_degenerate_primary(T)
  end

  # we first create some blocks corresponding to the primary parts of T.
  # We then compute the orthogonal sum of those blocks, which is a nice
  # TorQuadMod isometric to T: we get then an isometry between those two modules.
  gensOT = fmpz_mat[]
  pd = sort(prime_divisors(order(T)))
  blocks = TorQuadModMor[primary_part(T, pd[1])[2]]
  popfirst!(pd)
  while !isempty(pd)
    f = blocks[end]
    ok, j = has_complement(f)
    @assert ok
    _T = domain(j)
    _f = primary_part(_T, pd[1])[2]
    push!(blocks, compose(_f, j))
    popfirst!(pd)
  end
  Torth, _, _ = Hecke._orthogonal_sum_with_injections_and_projections(domain.(blocks))
  ok, phi = is_isometric_with_isometry(Torth, T)
  @assert ok # Same module with different basis, somehow

  # We compute the orthoognal group of each (primary) block. Since Torth is made
  # diagonally, we can generate its orthogonal group by taking some diagonal matrices
  # from taking the cartesian product on the generators of the blocks. We have nothing
  # to add since we can't map different primary parts to each other.
  orth_blocks = _compute_gens_split_degenerate_primary.(domain.(blocks))
  gensOTorth = fmpz_mat[]
  for x in Hecke.cartesian_product_iterator(orth_blocks, inplace=false)
    m = block_diagonal_matrix([f for f in x])
    push!(gensOTorth, m)
  end

  gensOTorth = TorQuadModMor[hom(Torth, Torth, g) for g in gensOTorth]
  gensOT = fmpz_mat[compose(compose(inv(phi), g), phi).map_ab.map for g in gensOTorth]
  unique!(gensOT)
  return gensOT
end

function _compute_gens_split_degenerate_primary(T::TorQuadMod)
  # we want only "primary" non semi regular modules that splits their
  # radical quadratic
  if is_semi_regular(T)
    N, i = normal_form(T)
    j = inv(i)
    gensOT = _compute_gens(N)
    gensOT = TorQuadModMor[hom(N, N, g) for g in gensOT]
    gensOT = fmpz_mat[compose(compose(i, g), j).map_ab.map for g in gensOT]
    return gensOT
  end

  @assert is_prime_power_with_data(elementary_divisors(T)[end])[1]
  rd, i = radical_quadratic(T)
  ok, j = has_complement(i)
  @assert ok

  # N now is isometric to the normal form of T, so it is in particular
  # semi-regular. We can already compute generators of the orth group of
  # its normal form, we then bring them back to N
  N = domain(j)
  @assert is_semi_regular(N)
  NN, NtoNN = normal_form(N)
  @assert is_bijective(NtoNN)
  NNtoN = inv(NtoNN)
  gensONN = TorQuadModMor[hom(NN, NN, m) for m in _compute_gens(NN)]
  gensON = fmpz_mat[compose(compose(NtoNN, g), NNtoN).map_ab.map for g in gensONN]
  n1 = nrows(gensON[1])

  # for the rd, since its quadratic form is trivial, automorphisms are just
  # automorphisms of the underlying abelian group.
  gensOrd = fmpz_mat[]
  OArd = automorphism_group(abelian_group(rd))
  for f in gens(OArd)
    push!(gensOrd, matrix(f))
  end
  n2 = nrows(gensOrd[1])

  # finally, we have to consider automorphism which maps N into rd: these are well
  # defined because N and rd are orthogonal and the quadratic form on rd is trivial.
  R, psi = hom(abelian_group(N), abelian_group(rd), task = :map)
  Ntord = fmpz_mat[matrix(psi(f)) for f in gens(R)]
  Torth, _, _ = Hecke._orthogonal_sum_with_injections_and_projections([rd, N])

  # Same module with different basis
  ok, phi = is_isometric_with_isometry(Torth, T)
  @assert ok

  # We have two kind of generators: block diagonal matrices coming from isometries
  # of N and rd, and those acting identically on N and rd but which send N to rd.
  # Combining all of them together, we have generators (maybe the set is too big
  # compared to what is needed) for the orthogonal group of Torth, and so of T.
  gensOTorth = fmpz_mat[]
  for x in gensOrd
    m = block_diagonal_matrix([x, identity_matrix(ZZ, n1)])
    push!(gensOTorth, m)
  end
  for x in gensON
    m = block_diagonal_matrix([identity_matrix(ZZ, n2), x])
    push!(gensOTorth, m)
  end
  for m in Ntord
    M = identity_matrix(ZZ, n1+n2)
    M[(n2+1):end, 1:n2] = m
    push!(gensOTorth, M)
  end
  gensOTorth = TorQuadModMor[hom(Torth, Torth, m) for m in gensOTorth]
  gensOT = fmpz_mat[compose(compose(inv(phi), g), phi).map_ab.map for g in gensOTorth]
  unique!(gensOT)
  return gensOT
end

