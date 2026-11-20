################################################################################
#
#  Orbits on labels
#
#  The acting group G is given by the matrices of a label action: first a
#  polycyclic generating sequence g_1, ..., g_n of a solvable normal
#  subgroup N, then further matrices h_1, ..., h_m which together with N
#  generate G.
#
#  1. Along the chain U_i = <g_i, ..., g_n>, the U_{i+1}-orbit B of a point
#     is a block, and the U_i-orbit is the union of B, B g_i, B g_i^2, ...
#     Each new point costs one image, and one more image per chain element
#     detects that the blocks are closed under g_i.  This is the orbit
#     algorithm of Leedham-Green, see [LNS84](@cite) and [OBr90](@cite),
#     Section 3.5.
#
#  2. The N-orbits are blocks for G.  A breadth-first search over blocks
#     computes the image of a whole block under h_j only if the image of its
#     first point is new.
#
#  With D points, B blocks and O orbits this takes D + B * (m - 1) + O * n
#  images.  For N = G step 2 is void; for N = 1 the blocks are points and
#  step 2 is the plain orbit algorithm with D * m images.
#
#  The points of the current orbit are kept in the order in which they are
#  found.  That order encodes a transversal: with f_i the number of blocks
#  joined by g_i, the point at position
#      q = sum(e_i * f_{i+1} * ... * f_n for i in 1:n)   (0 <= e_i < f_i)
#  of the N-orbit of the first point pt is  pt * g_n^e_n * ... * g_1^e_1,
#  and position q of block c is the image of position q of block 0 under
#  the product of the h_j on the path from block 0 to block c.
#
################################################################################

# The orbit that `label_orbits` is computing.  One object serves all orbits.
mutable struct LabelOrbit{L<:Unsigned}
  # Buffer for the points, of which orbit[1:npoints] is in use: it grows
  # geometrically and is not shrunk between orbits, so that orbits do not
  # allocate.  With `track = false` points whose images are all computed are
  # dropped from it during the search, see `TRIM_THRESHOLD`.
  orbit::Vector{L}
  npoints::Int
  rep::Int                 # first point, the smallest label of the orbit
  length::Int              # number of points of the orbit
  block_length::Int        # length of the N-orbits
  # factors[i]: number of blocks joined by g_i
  factors::Vector{Int}
  # closures[i]: 0-based position of pt * g_i^factors[i] in the U_{i+1}-orbit
  # of pt, or -1 if that was not determined
  closures::Vector{Int}
  # parents[c]: block c+1 is the image of block parents[c][1] (0-based)
  # under h_j, j = parents[c][2]
  parents::Vector{Tuple{Int,Int}}
end

function LabelOrbit{L}(npc::Int) where {L}
  return LabelOrbit{L}(Vector{L}(undef, 1024), 0, 0, 0, 0, zeros(Int, npc),
                       zeros(Int, npc), Tuple{Int,Int}[])
end

# Processed blocks are dropped once they take up this many entries, unless
# the caller asks for another threshold.
const TRIM_THRESHOLD = 2^22

# Make room for `len` more entries after the first `n`.
@inline function reserve!(orbit::Vector, n::Int, len::Int)
  n + len > length(orbit) && resize!(orbit, max(2 * length(orbit), n + len))
  return orbit
end

# Extend the orbit orbit[1:n], one U_{i+1}-orbit, to the U_i-orbit; `tab`
# belongs to g_i.  Return the new number of points, the number f of blocks,
# and the image of the first point under g_i^f.
function pc_step!(K::LabelAction, tab::Vector, orbit::Vector{L}, n::Int,
                  seen::BitVector) where {L}
  len0 = n
  start = 0
  f = 1
  @inbounds while true
    y = image(K, tab, Int(orbit[start + 1]))
    seen[y] && return n, f, y
    seen[y] = true
    reserve!(orbit, n, len0)
    orbit[n + 1] = y % L
    for idx in 2:len0
      z = image(K, tab, Int(orbit[start + idx]))
      seen[z] = true
      orbit[n + idx] = z % L
    end
    start = n
    n += len0
    f += 1
  end
end

# Extend the orbit orbit[1:len0], one N-orbit, to the G-orbit.  Return the
# number of points left in `orbit` and the number of blocks.  With
# `parents === nothing` processed blocks are dropped once they take up
# `trim_threshold` entries.
function block_search!(K::LabelAction, orbit::Vector{L}, len0::Int, seen::BitVector,
                       npc::Int, nextra::Int, parents, trim_threshold::Int) where {L}
  if len0 == 1
    return point_search!(K, orbit, seen, npc, nextra, parents, trim_threshold)
  end
  n = len0
  nblocks = 1
  dropped = 0
  b = 0
  @inbounds while b < nblocks
    base = b * len0
    for j in 1:nextra
      tab = K.tables[npc + j]
      y = image(K, tab, Int(orbit[base + 1]))
      seen[y] && continue
      seen[y] = true
      reserve!(orbit, n, len0)
      orbit[n + 1] = y % L
      for idx in 2:len0
        z = image(K, tab, Int(orbit[base + idx]))
        seen[z] = true
        orbit[n + idx] = z % L
      end
      n += len0
      nblocks += 1
      parents === nothing || push!(parents, (b, j))
    end
    b += 1
    if parents === nothing && b * len0 >= trim_threshold && 2 * b >= nblocks
      copyto!(orbit, 1, orbit, b * len0 + 1, n - b * len0)
      n -= b * len0
      dropped += b
      nblocks -= b
      b = 0
    end
  end
  return n, nblocks + dropped
end

# The case of `block_search!` where the blocks are single points.  Whether
# an image is new is hard to predict, and a branch on it would wait for the
# lookup in `seen`; so every image is appended to the orbit, and the number
# of points is advanced only if it is new.
function point_search!(K::LabelAction, orbit::Vector{L}, seen::BitVector, npc::Int,
                       nextra::Int, parents, trim_threshold::Int) where {L}
  tables = K.tables
  n = 1
  dropped = 0
  b = 0
  @inbounds while b < n
    u = Int(orbit[b + 1])
    reserve!(orbit, n, nextra)
    parents === nothing || reserve!(parents, n, nextra)
    for j in 1:nextra
      y = image(K, tables[npc + j], u)
      orbit[n + 1] = y % L
      parents === nothing || (parents[n] = (b, j))
      n += !seen[y]
      seen[y] = true
    end
    b += 1
    if parents === nothing && b >= trim_threshold && 2 * b >= n
      copyto!(orbit, 1, orbit, b + 1, n - b)
      n -= b
      dropped += b
      b = 0
    end
  end
  parents === nothing || resize!(parents, n - 1)
  return n, n + dropped
end

# Call `handle(state)` once for every orbit on the labels, in increasing
# order of the smallest label.  The first `npc` matrices of `K` form the
# polycyclic sequence, with relative orders `relorders`; the remaining ones
# are the further generators.
#
# With `track = true` the state holds the complete orbit and the data
# describing the transversal; otherwise `state.orbit` is not meaningful.
#
# The labels are stored in the smallest type that holds them.  That type and
# the threshold for dropping blocks can be set, so that small cases reach
# the code for large ones.
function label_orbits(handle, K::LabelAction, relorders::Vector{Int}; track::Bool,
                      trim_threshold::Int=TRIM_THRESHOLD,
                      label_type::Type{<:Unsigned}=default_label_type(K.codec))
  state = LabelOrbit{label_type}(length(relorders))
  return label_orbits(handle, K, relorders, track, trim_threshold, state)
end

default_label_type(C::LabelCodec) = C.D <= typemax(UInt32) ? UInt32 : UInt64

function label_orbits(handle::F, K::LabelAction, relorders::Vector{Int}, track::Bool,
                      trim_threshold::Int, state::LabelOrbit{L}) where {F,L}
  npc = length(relorders)
  nextra = length(K.tables) - npc
  seen = falses(K.codec.D)
  orbit = state.orbit

  l = findnext(!, seen, 1)
  while l !== nothing
    seen[l] = true
    orbit[1] = l % L
    n = 1
    state.rep = l

    for i in npc:-1:1
      len0 = n
      n, f, y = pc_step!(K, K.tables[i], orbit, n, seen)
      state.factors[i] = f
      # only then pt * g_i^f yields a new element of the stabilizer
      if track && f < relorders[i]
        state.closures[i] = findfirst(==(y % L), view(orbit, 1:len0))::Int - 1
      else
        state.closures[i] = -1
      end
    end

    state.block_length = n
    empty!(state.parents)
    n, nblocks = block_search!(K, orbit, n, seen, npc, nextra,
                               track ? state.parents : nothing, trim_threshold)
    state.npoints = n
    state.length = nblocks * state.block_length
    handle(state)

    l = findnext(!, seen, l + 1)
  end
  return nothing
end
