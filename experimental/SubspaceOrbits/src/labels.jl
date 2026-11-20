################################################################################
#
#  Labels for subspaces
#
#  A k-dimensional subspace of F_p^d is identified with its basis matrix in
#  reduced row echelon form.  The pivot columns K = {K[1] < ... < K[k]} form
#  its definition set; the free entries are those at (i, j) with j > K[i]
#  and j not in K.  Reading the free entries as base-p digits (row by row,
#  left to right, least significant first) numbers the matrices with
#  definition set K, and adding the number of matrices whose definition set
#  is lexicographically smaller gives a bijection onto 1:D, the label.
#
#  This is the labelling of [OBr90](@cite), Section 3.3.
#
################################################################################

# One visited bit per subspace: 2^36 subspaces are 8 GiB.
const MAX_NUMBER_OF_SUBSPACES = 2^36

struct LabelCodec
  p::Int
  d::Int
  k::Int
  D::Int                   # number of k-subspaces
  pivots::Matrix{Int}      # pivots[:, s] is the s-th definition set
  offsets::Vector{Int}     # labels of set s are offsets[s]+1 : offsets[s+1]
  # rankcum[i, j+1] = sum(binomial(d-t, k-i) for t in 1:j), used to rank
  # a definition set among all of them
  rankcum::Matrix{Int}
  # coarse[((lab-1) >> coarse_shift) + 1] is the definition set of the
  # first label in that bucket; replaces most of a binary search
  coarse::Vector{Int32}
  coarse_shift::Int
end

function number_of_subspaces(p::Int, d::Int, k::Int)
  D = ZZ(1)
  for i in 0:(k - 1)
    D = divexact(D * (ZZ(p)^(d - i) - 1), ZZ(p)^(i + 1) - 1)
  end
  return D
end

# Number of free entries in row i of a matrix whose i-th pivot is ki.
number_of_free_entries(C::LabelCodec, i::Int, ki::Int) = C.d - ki - (C.k - i)

function LabelCodec(p::Int, d::Int, k::Int)
  @req 0 < k < d "need 0 < k < d"
  nsubspaces = number_of_subspaces(p, d, k)
  @req nsubspaces <= MAX_NUMBER_OF_SUBSPACES "too many subspaces: $nsubspaces"
  D = Int(nsubspaces)

  nsets = binomial(d, k)
  pivots = Matrix{Int}(undef, k, nsets)
  offsets = Vector{Int}(undef, nsets)
  K = collect(1:k)
  offset = 0
  for s in 1:nsets
    pivots[:, s] = K
    offsets[s] = offset
    offset += p^sum(d - K[i] - (k - i) for i in 1:k)

    # lexicographically next k-subset of 1:d
    i = k
    while i >= 1 && K[i] == d - k + i
      i -= 1
    end
    i == 0 && break
    K[i] += 1
    for j in (i + 1):k
      K[j] = K[j - 1] + 1
    end
  end
  @assert offset == D

  rankcum = zeros(Int, k, d + 1)
  for i in 1:k, j in 1:d
    rankcum[i, j + 1] = rankcum[i, j] + binomial(d - j, k - i)
  end

  coarse_shift = max(0, ndigits(D; base=2) - 12)
  coarse = Vector{Int32}(undef, ((D - 1) >> coarse_shift) + 1)
  s = 1
  for b in eachindex(coarse)
    first_label = ((b - 1) << coarse_shift) + 1
    while s < nsets && offsets[s + 1] < first_label
      s += 1
    end
    coarse[b] = s
  end

  return LabelCodec(p, d, k, D, pivots, offsets, rankcum, coarse, coarse_shift)
end

# Index of the definition set containing the label `lab`.
@inline function definition_set_of_label(C::LabelCodec, lab::Int)
  offsets = C.offsets
  nsets = length(offsets)
  @inbounds begin
    s = Int(C.coarse[(((lab - 1) % UInt) >> (C.coarse_shift % UInt & 63)) % Int + 1])
    while s < nsets && offsets[s + 1] < lab
      s += 1
    end
  end
  return s
end

# Index of the definition set with the given increasing 1-based pivots.
@inline function definition_set_index(C::LabelCodec, pivs)
  rankcum = C.rankcum
  s = 1
  prev = 0
  @inbounds for i in 1:C.k
    ki = pivs[i]
    s += rankcum[i, ki] - rankcum[i, prev + 1]
    prev = ki
  end
  return s
end

################################################################################
#
#  Reference implementation on integer matrices
#
#  Used for input and output, and to test the label actions.  Matrices have
#  entries in 0:p-1.
#
################################################################################

function label_to_matrix(C::LabelCodec, lab::Int)
  @req 1 <= lab <= C.D "label out of range"
  s = definition_set_of_label(C, lab)
  K = view(C.pivots, :, s)
  rem = lab - 1 - C.offsets[s]
  M = zeros(Int, C.k, C.d)
  for i in 1:C.k
    M[i, K[i]] = 1
    for j in (K[i] + 1):C.d
      j in K && continue
      rem, M[i, j] = divrem(rem, C.p)
    end
  end
  return M
end

# `M` must be in reduced row echelon form, of rank k.
function matrix_to_label(C::LabelCodec, M::Matrix{Int})
  K = [findfirst(!is_zero, view(M, i, :))::Int for i in 1:C.k]
  s = definition_set_index(C, K)
  lab = C.offsets[s] + 1
  weight = 1
  for i in 1:C.k
    for j in (K[i] + 1):C.d
      j in K && continue
      lab += M[i, j] * weight
      weight *= C.p
    end
  end
  return lab
end

# Reduced row echelon form of a matrix of full row rank, in place.
function rref_mod_p!(M::Matrix{Int}, p::Int)
  k, d = size(M)
  for i in 1:k
    # pivot: the leftmost nonzero entry in rows i:k
    piv = (0, d + 1)
    for r in i:k
      c = findfirst(!is_zero, view(M, r, :))
      c !== nothing && c < piv[2] && (piv = (r, c))
    end
    r, c = piv
    M[[i, r], :] = M[[r, i], :]
    M[i, :] = mod.(M[i, :] .* invmod(M[i, c], p), p)
    for t in 1:k
      t == i && continue
      M[t, :] = mod.(M[t, :] .- M[t, c] .* M[i, :], p)
    end
  end
  return M
end

# Label of the image of the subspace with label `lab` under the matrix `g`.
function reference_image(C::LabelCodec, g::Matrix{Int}, lab::Int)
  return matrix_to_label(C, rref_mod_p!(mod.(label_to_matrix(C, lab) * g, C.p), C.p))
end

# Basis in reduced row echelon form of the orthogonal complement, with
# respect to the standard bilinear form, of the row space of `M`; `M` must
# be in reduced row echelon form.
function orthogonal_complement(M::Matrix{Int}, p::Int)
  k, d = size(M)
  K = [findfirst(!is_zero, view(M, i, :))::Int for i in 1:k]
  N = zeros(Int, d - k, d)
  r = 0
  for j in 1:d
    j in K && continue
    r += 1
    N[r, j] = 1
    for i in 1:k
      N[r, K[i]] = mod(-M[i, j], p)
    end
  end
  return rref_mod_p!(N, p)
end
