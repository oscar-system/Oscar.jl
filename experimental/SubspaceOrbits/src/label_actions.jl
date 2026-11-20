################################################################################
#
#  Label actions: the action of a matrix on labels
#
#  A label action maps a label to the label of the image of that subspace
#  under one of a fixed list of matrices.  Each row of a basis matrix is one
#  machine word, split into d fields of b bits, field j-1 holding column j.
#
#  The non-pivot columns of a definition set fall into k+1 runs: run r
#  consists of the columns between the r-th and the (r+1)-th pivot.  Deleting
#  the pivot columns from a row moves run r down by r fields, so packing and
#  unpacking the free entries takes one shift and one mask per run.
#
#  The dimension k of the subspaces, and for odd p also d, are type parameters
#  of the label actions, and the code for the image is generated without loops
#  over the rows: left to itself the compiler vectorizes these loops of at
#  most k steps, which on x86 costs more than it gains.
#
################################################################################

abstract type LabelAction end

# Shifts by 0 <= n < (number of bits of x).  The generic `<<` and `>>` also
# handle negative and too large n, which costs more than the shift itself.
@inline shl(x::T, n::Int) where {T<:Unsigned} = x << ((n % UInt) & (8 * sizeof(T) - 1))
@inline shr(x::T, n::Int) where {T<:Unsigned} = x >> ((n % UInt) & (8 * sizeof(T) - 1))

# runmasks[r+1, s] selects the fields of run r of definition set s.
function run_masks(::Type{T}, C::LabelCodec, b::Int) where {T<:Unsigned}
  nsets = size(C.pivots, 2)
  field = (one(T) << b) - one(T)
  masks = zeros(T, C.k + 1, nsets)
  for s in 1:nsets
    r = 0
    for j in 1:C.d
      if r < C.k && C.pivots[r + 1, s] == j
        r += 1
      else
        masks[r + 1, s] |= field << ((j - 1) * b)
      end
    end
  end
  return masks
end

# res[i, s] is the number of free entries in rows 1:i-1 of definition set s.
function free_entries_before(C::LabelCodec)
  nsets = size(C.pivots, 2)
  res = zeros(Int, C.k, nsets)
  for s in 1:nsets, i in 2:C.k
    res[i, s] = res[i - 1, s] + C.d - C.pivots[i - 1, s] - (C.k - i + 1)
  end
  return res
end

# Index of the definition set whose pivots are the bits of `pivmask`.
@inline function definition_set_index(C::LabelCodec, pivmask::UInt64, ::Val{k}) where {k}
  rankcum = C.rankcum
  s = 1
  prev = 0
  @inbounds for i in 1:k
    ki = trailing_zeros(pivmask) + 1
    pivmask &= pivmask - UInt64(1)
    s += rankcum[i, ki] - rankcum[i, prev + 1]
    prev = ki
  end
  return s
end

################################################################################
#
#  p = 2: one bit per entry
#
################################################################################

struct F2LabelAction{k} <: LabelAction
  codec::LabelCodec
  runmasks::Matrix{UInt64}
  rowshift::Matrix{Int}   # rowshift[i, s]: free entries before row i
  # A row vector is cut into chunks of `chunkbits` bits.  The entry
  # tables[g][(c << chunkbits) + v + 1] is the row vector with value v in
  # chunk c and zero elsewhere, times matrix g; so a row times g is the XOR
  # of one table entry per chunk.
  chunkbits::Int
  nchunks::Int
  tables::Vector{Vector{UInt64}}
end

# Up to this dimension a single table of 2^d words per matrix is used.
const F2_SINGLE_TABLE_MAX_DIM = 12

function F2LabelAction(C::LabelCodec, gens::Vector{Matrix{Int}})
  @req C.p == 2 "need p = 2"
  d = C.d
  chunkbits = d <= F2_SINGLE_TABLE_MAX_DIM ? d : 8
  nchunks = cld(d, chunkbits)

  tables = Vector{UInt64}[]
  for g in gens
    grows = [sum(UInt64(g[i, j]) << (j - 1) for j in 1:d) for i in 1:d]
    tab = zeros(UInt64, nchunks << chunkbits)
    for c in 0:(nchunks - 1), v in 1:((1 << chunkbits) - 1)
      i = c * chunkbits + trailing_zeros(v) + 1
      i > d && continue
      base = c << chunkbits
      tab[base + v + 1] = xor(tab[base + (v & (v - 1)) + 1], grows[i])
    end
    push!(tables, tab)
  end

  return F2LabelAction{C.k}(C, run_masks(UInt64, C, 1), free_entries_before(C), chunkbits,
                       nchunks, tables)
end

# A row vector times the matrix with table `tab`.
@inline function multiply(K::F2LabelAction, tab::Vector{UInt64}, row::UInt64)
  chunkbits = K.chunkbits
  chunkmask = shl(UInt64(1), chunkbits) - UInt64(1)
  @inbounds img = tab[(row & chunkmask) + 1]
  @inbounds for ch in 1:(K.nchunks - 1)
    img = xor(img, tab[(ch * (chunkmask % Int + 1)) + (shr(row, ch * chunkbits) & chunkmask) + 1])
  end
  return img
end

# The image rows row_1, ..., row_k are echelonized in the order in which
# they arise, with pivot bits piv_1, ..., piv_k (0-based); no branch depends
# on the data.
@generated function image(K::F2LabelAction{k}, tab::Vector{UInt64}, lab::Int) where {k}
  rows = [Symbol(:row_, i) for i in 1:k]
  pivs = [Symbol(:piv_, i) for i in 1:k]
  body = Expr(:block)

  for i in 1:k
    push!(body.args, quote
      # decode row i
      ki = C.pivots[$i, s]
      x = d - ki - $(k - i)
      c = shl(rem & (shl(UInt64(1), x) - UInt64(1)), ki - $i)
      rem = shr(rem, x)
      row = shl(UInt64(1), ki - 1)
      $([:(row |= shl(c, $r) & runmasks[$(r + 1), s]) for r in i:k]...)

      img = multiply(K, tab, row)

      # clear the pivot columns of the finished rows in the new row; a
      # finished row is zero in the other pivot columns, so the multipliers
      # can all be read off the new row
      acc = UInt64(0)
      $([:(acc = xor(acc, $(rows[t]) & -(shr(img, $(pivs[t])) & UInt64(1)))) for t in 1:(i - 1)]...)
      img = xor(img, acc)

      # clear the new pivot column in the finished rows
      $(pivs[i]) = trailing_zeros(img)
      $([:($(rows[t]) = xor($(rows[t]), img & -(shr($(rows[t]), $(pivs[i])) & UInt64(1)))) for t in 1:(i - 1)]...)
      $(rows[i]) = img
      pivmask |= shl(UInt64(1), $(pivs[i]))
    end)
  end

  # encode; a row is zero in the runs before its pivot
  push!(body.args, quote
    s2 = definition_set_index(C, pivmask, Val($k))
    digits = UInt64(0)
  end)
  for i in 1:k
    push!(body.args, quote
      rank = count_ones(pivmask & (shl(UInt64(1), $(pivs[i])) - UInt64(1)))
      c = UInt64(0)
      $([:(c |= shr($(rows[i]) & runmasks[$(r + 1), s2], $r)) for r in 1:k]...)
      digits |= shl(shr(c, $(pivs[i]) - rank), K.rowshift[rank + 1, s2])
    end)
  end

  return quote
    $(Expr(:meta, :inline))
    C = K.codec
    d = C.d
    runmasks = K.runmasks
    s = definition_set_of_label(C, lab)
    @inbounds begin
      rem = (lab - 1 - C.offsets[s]) % UInt64
      pivmask = UInt64(0)
      $body
      return C.offsets[s2] + 1 + (digits % Int)
    end
  end
end

################################################################################
#
#  Odd p: several bits per entry
#
#  With b bits per field, sums of several reduced rows fit into a word
#  without carries between fields, so a linear combination of rows is a few
#  word operations followed by one reduction of all fields modulo p.
#
#  There are two ways to reduce all fields at once.
#
#  - Binary long division.  One level subtracts p*2^(s-1) from every field
#    that is at least that large: adding 2^(b-1) - p*2^(s-1) to a field sets
#    its top bit exactly in that case.  Levels L, L-1, ..., 1 reduce fields
#    below p*2^L to fields below p, which needs p*2^(L-1) <= 2^(b-1).
#
#  - Division by multiplication.  If (v * m) >> s == div(v, p) for all field
#    values v that occur, and v * m fits into a field, then the quotients of
#    all fields are obtained with one multiplication of the word by m.  This
#    is faster but needs wider fields.
#
################################################################################

# How to reduce fields with values in 0:bound.
struct FieldReduction{T<:Unsigned}
  bound::Int
  by_multiplication::Bool
  multiplier::T
  shift::Int
  quotient_mask::T
  levels::Int
end

struct PackedLabelAction{T<:Unsigned,k,d} <: LabelAction
  codec::LabelCodec
  b::Int                    # bits per field
  binv::Int                 # (x * binv) >> 16 == div(x, b) for bit positions x
  fieldmask::T
  ones::T                   # 1 in every field
  red_add::Vector{T}        # red_add[s]: 2^(b-1) - p*2^(s-1) in every field
  red_sub::Vector{T}        # red_sub[s] = p*2^(s-1)
  nsum::Int                 # this many reduced rows can be summed
  nelim::Int                # this many row operations can be combined
  red_sum::FieldReduction{T}
  red_elim::FieldReduction{T}
  red_row::FieldReduction{T}   # for fields below p^2
  inverses::Vector{T}       # inverses[a] = a^-1 mod p
  powers::Vector{Int}       # powers[x+1] = p^x
  # expand[v+1] holds the base-p digits of v in consecutive fields, least
  # significant first; empty if that table would be too large
  expand::Vector{T}
  runmasks::Matrix{T}
  rowweight::Matrix{Int}    # rowweight[i, s] = p^(free entries before row i)
  # tables[g][(j-1)*p + a + 1] is a times row j of matrix g
  tables::Vector{Vector{T}}
end

const PACKED_EXPAND_TABLE_MAX_SIZE = 2^16
const PACKED_MAX_PRIME = 2^16
const PACKED_MAX_MULTIPLICATION_BOUND = 2^16

# Smallest L such that fields up to `maxval` are below p*2^L.
function reduction_levels(p::Int, maxval::Int)
  L = 0
  while p << L <= maxval
    L += 1
  end
  return L
end

function FieldReduction{T}(p::Int, b::Int, bound::Int) where {T<:Unsigned}
  levels = reduction_levels(p, bound)
  if bound <= PACKED_MAX_MULTIPLICATION_BOUND
    for shift in 1:min(b - 1, 30)
      m = cld(1 << shift, p)
      ZZ(bound) * m < ZZ(1) << b || continue
      all(v -> (v * m) >> shift == div(v, p), 0:bound) || continue
      ones = sum(one(T) << (j * b) for j in 0:(div(8 * sizeof(T), b) - 1))
      mask = ones * ((one(T) << (b - shift)) - one(T))
      return FieldReduction{T}(bound, true, T(m), shift, mask, levels)
    end
  end
  return FieldReduction{T}(bound, false, zero(T), 0, zero(T), levels)
end

function PackedLabelAction(C::LabelCodec, gens::Vector{Matrix{Int}})
  p = C.p
  d = C.d
  @req is_odd(p) "need p odd"
  @req p < PACKED_MAX_PRIME "prime too large: $p"
  levels_row = reduction_levels(p, p^2 - 1)
  for T in (UInt64, UInt128)
    b = fld(8 * sizeof(T), d)
    maxlevels = b - ndigits(p; base=2)
    maxlevels >= levels_row && return PackedLabelAction{T,C.k,d}(C, gens, b, maxlevels)
  end
  throw(ArgumentError("dimension $d and prime $p are too large"))
end

function PackedLabelAction{T,k,d}(C::LabelCodec, gens::Vector{Matrix{Int}}, b::Int,
                             maxlevels::Int) where {T<:Unsigned,k,d}
  p = C.p
  pack(v) = sum(T(v[j]) << ((j - 1) * b) for j in eachindex(v))

  fieldmask = (one(T) << b) - one(T)
  ones = pack(fill(1, d))
  red_add = [ones * ((one(T) << (b - 1)) - (T(p) << (s - 1))) for s in 1:maxlevels]
  red_sub = [T(p) << (s - 1) for s in 1:maxlevels]

  # fields up to `capacity` can be reduced by long division
  capacity = (ZZ(p) << maxlevels) - 1
  nsum = Int(min(fld(capacity, p - 1), d))
  nelim = Int(min(fld(capacity - (p - 1), p * (p - 1)), max(k - 1, 1)))
  red_sum = FieldReduction{T}(p, b, nsum * (p - 1))
  red_elim = FieldReduction{T}(p, b, (p - 1) + nelim * p * (p - 1))
  red_row = FieldReduction{T}(p, b, p^2 - 1)

  inverses = [T(invmod(a, p)) for a in 1:(p - 1)]
  powers = [p^x for x in 0:(d - k)]
  expand = T[]
  if powers[end] <= PACKED_EXPAND_TABLE_MAX_SIZE
    expand = [pack(digits(v; base=p, pad=d - k)) for v in 0:(powers[end] - 1)]
  end

  tables = Vector{T}[]
  for g in gens
    tab = zeros(T, d * p)
    for j in 1:d, a in 1:(p - 1)
      tab[(j - 1) * p + a + 1] = pack(mod.(a .* g[j, :], p))
    end
    push!(tables, tab)
  end

  return PackedLabelAction{T,k,d}(C, b, cld(1 << 16, b), fieldmask, ones, red_add, red_sub,
                             nsum, nelim, red_sum, red_elim, red_row, inverses,
                             powers, expand, run_masks(T, C, b),
                             p .^ free_entries_before(C), tables)
end

@inline function reduce_fields(K::PackedLabelAction{T}, x::T,
                               R::FieldReduction{T}) where {T}
  if R.by_multiplication
    q = shr(x * R.multiplier, R.shift) & R.quotient_mask
    return x - q * T(K.codec.p)
  end
  top = K.b - 1
  @inbounds for s in R.levels:-1:1
    m = shr(x + K.red_add[s], top) & K.ones
    x -= m * K.red_sub[s]
  end
  return x
end

# The x lowest base-p digits of v, in consecutive fields.
@inline function expand_digits(K::PackedLabelAction{T}, v::Int, x::Int) where {T}
  isempty(K.expand) || return @inbounds K.expand[v + 1]
  c = zero(T)
  for j in 0:(x - 1)
    v, a = divrem(v, K.codec.p)
    c |= shl(T(a), j * K.b)
  end
  return c
end

# A row vector times the matrix with table `tab`: the sum of one table entry
# per column.
@generated function multiply(K::PackedLabelAction{T,k,d}, tab::Vector{T},
                             row::T) where {T,k,d}
  entry(j) = :(tab[$j * p + ((shr(row, $j * b) & K.fieldmask) % Int) + 1])
  return quote
    $(Expr(:meta, :inline))
    p = K.codec.p
    b = K.b
    @inbounds if K.nsum >= $d
      img = +($([entry(j) for j in 0:(d - 1)]...))
    else
      # the sum of all entries does not fit into the fields
      img = zero(T)
      count = 0
      for j in 0:($d - 1)
        if count == K.nsum
          img = reduce_fields(K, img, K.red_sum)
          count = 1
        end
        img += $(entry(:j))
        count += 1
      end
    end
    return reduce_fields(K, img, K.red_sum)
  end
end

# The value of the fields of c, read as base-p digits; only the d-k lowest
# fields may be nonzero.  All of them are used, also where fewer would do,
# so that no branch depends on the data.
@generated function collect_digits(K::PackedLabelAction{T,k,d}, c::T) where {T,k,d}
  digit(j) = :((shr(c, $j * K.b) & K.fieldmask) % Int)
  return quote
    $(Expr(:meta, :inline))
    @inbounds return +($(digit(0)), $([:($(digit(j)) * K.powers[$(j + 1)]) for j in 1:(d - k - 1)]...))
  end
end

# The image rows row_1, ..., row_k are echelonized in the order in which
# they arise, with pivot fields piv_1, ..., piv_k (0-based).
@generated function image(K::PackedLabelAction{T,k,d}, tab::Vector{T},
                          lab::Int) where {T,k,d}
  rows = [Symbol(:row_, i) for i in 1:k]
  pivs = [Symbol(:piv_, i) for i in 1:k]
  body = Expr(:block)

  for i in 1:k
    # clear the pivot columns of the finished rows in the new row; a
    # finished row is zero in the other pivot columns, so the multipliers
    # can all be read off the new row
    eliminate = quote end
    if i > 1
      eliminate = quote
        acc = img
        count = 0
        $([quote
          if count == K.nelim
            acc = reduce_fields(K, acc, K.red_elim)
            count = 0
          end
          acc += (pT - (shr(img, $(pivs[t]) * b) & fieldmask)) * $(rows[t])
          count += 1
        end for t in 1:(i - 1)]...)
        img = reduce_fields(K, acc, K.red_elim)
      end
    end

    push!(body.args, quote
      # decode row i
      ki = C.pivots[$i, s]
      x = $d - ki - $(k - i)
      rem, v = divrem(rem, K.powers[x + 1])
      c = shl(expand_digits(K, v, x), (ki - $i) * b)
      row = shl(one(T), (ki - 1) * b)
      $([:(row |= shl(c, $r * b) & runmasks[$(r + 1), s]) for r in i:k]...)

      img = multiply(K, tab, row)
      $eliminate

      # normalize the new row
      $(pivs[i]) = (trailing_zeros(img) * K.binv) >> 16
      a = shr(img, $(pivs[i]) * b) & fieldmask
      img = reduce_fields(K, img * K.inverses[a % Int], K.red_row)

      # clear the new pivot column in the finished rows
      $([:($(rows[t]) = reduce_fields(K, $(rows[t]) + (pT - (shr($(rows[t]), $(pivs[i]) * b) & fieldmask)) * img, K.red_row)) for t in 1:(i - 1)]...)
      $(rows[i]) = img
      pivmask |= shl(UInt64(1), $(pivs[i]))
    end)
  end

  # encode; a row is zero in the runs before its pivot
  push!(body.args, quote
    s2 = definition_set_index(C, pivmask, Val($k))
    res = C.offsets[s2] + 1
  end)
  for i in 1:k
    push!(body.args, quote
      rank = count_ones(pivmask & (shl(UInt64(1), $(pivs[i])) - UInt64(1)))
      c = zero(T)
      $([:(c |= shr($(rows[i]) & runmasks[$(r + 1), s2], $r * b)) for r in 1:k]...)
      c = shr(c, ($(pivs[i]) - rank) * b)
      res += collect_digits(K, c) * K.rowweight[rank + 1, s2]
    end)
  end

  return quote
    $(Expr(:meta, :inline))
    C = K.codec
    b = K.b
    fieldmask = K.fieldmask
    runmasks = K.runmasks
    pT = T(C.p)
    s = definition_set_of_label(C, lab)
    @inbounds begin
      rem = lab - 1 - C.offsets[s]
      pivmask = UInt64(0)
      $body
      return res
    end
  end
end

################################################################################

function LabelAction(C::LabelCodec, gens::Vector{Matrix{Int}})
  return C.p == 2 ? F2LabelAction(C, gens) : PackedLabelAction(C, gens)
end

image(K::LabelAction, g::Int, lab::Int) = image(K, K.tables[g], lab)
