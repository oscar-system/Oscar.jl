################################################################################
# Tableaux
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
# by Tom Schmit and Ulrich Thiel; OSCAR-ified by Claudia He Yun and Matthias Zach.
################################################################################

################################################################################
#
#  Constructor and printing
#
################################################################################

@doc raw"""
    young_tableau([::Type{T}], v::Vector{Vector{<:IntegerUnion}}; check::Bool = true) where T <: IntegerUnion

Return the Young tableau given by `v` as an object of type `YoungTableau{T}`.

The element type `T` may be optionally specified, see also the examples below.

If `check` is `true` (default), it is checked whether `v` defines a tableau,
that is, whether the structure of `v` defines a partition.

# Examples
```jldoctest
julia> young_tableau([[1, 2, 3], [4, 5], [6]])
+---+---+---+
| 1 | 2 | 3 |
+---+---+---+
| 4 | 5 |
+---+---+
| 6 |
+---+

julia> young_tableau(Int8, [[1, 2, 3], [4, 5], [6]]) # save the elements in 8-bit integers
+---+---+---+
| 1 | 2 | 3 |
+---+---+---+
| 4 | 5 |
+---+---+
| 6 |
+---+
```
"""
young_tableau

function young_tableau(::Type{T}, v::Vector{Vector{TT}}; check::Bool = true) where {T <: IntegerUnion, TT <: IntegerUnion}
  if check
    @req _defines_partition(map(length, v)) "The input does not define a Young tableau: lengths of rows must be weakly decreasing"
  end
  return YoungTableau{T}(v)
end

young_tableau(v::Vector{Vector{T}}; check::Bool = true) where T <: IntegerUnion = young_tableau(T, v, check = check)

data(tab::YoungTableau) = tab.t

function Base.show(io::IO, tab::YoungTableau)
  print(io, data(tab))
end

function Base.show(io::IO, ::MIME"text/plain", tab::YoungTableau)
  if length(tab) == 0
    print(io, "Empty Young tableau")
    return nothing
  end
  # Translate the rows of `tab` into vectors of strings
  s = Vector{Vector{String}}(undef, length(tab))
  box_width = 0 # maximum length of an entry
  for i in 1:length(tab)
    r = tab[i]
    s[i] = Vector{String}(undef, length(r))
    for j in 1:length(r)
      x = string(r[j])
      box_width = max(length(x), box_width)
      s[i][j] = x
    end
  end
  # Pad the smaller boxes with whitespace (so that every box has the same width)
  for i in 1:length(s)
    for j in 1:length(s[i])
      x = s[i][j]
      while length(x) < box_width
        x = " "*x
      end
      s[i][j] = x
    end
  end
  # List of characters we use for the lines
  if is_unicode_allowed()
    pipe = "\u2502"
    horz_bar = "\u2500"
    cross = "\u253C"
    crossleft = "\u251C" # cross with left arm missing
    crossright = "\u2524" # cross with right arm missing
    crosstop = "\u252C" # cross with top arm missing
    crossbottom = "\u2534" # cross with bottom arm missing
    crosstopright = "\u2510" # cross with top and right arms missing
    crosstopleft = "\u250C" # cross with top and left arms missing
    crossbottomright = "\u2518" # cross with bottom and right arms missing
    crossbottomleft = "\u2514" # cross with bottom and left arms missing
  else
    pipe = "|"
    horz_bar = "-"
    # Without unicode, we just have one type of cross
    cross = "+"
    crossleft = "+"
    crossright = "+"
    crosstop = "+"
    crossbottom = "+"
    crosstopright = "+"
    crosstopleft = "+"
    crossbottomright = "+"
    crossbottomleft = "+"
  end
  hline_portion = ""
  while length(hline_portion) < box_width + 2
    hline_portion *= horz_bar
  end

  # Start actual printing

  # Print top line
  write(io, crosstopleft)
  for _ in 1:length(s[1]) - 1
    write(io, hline_portion, crosstop)
  end
  if length(s[1]) > 0
    write(io, hline_portion, crosstopright)
  end
  # Print rows and lines
  for i in 1:length(s)
    r = s[i]
    write(io, "\n")
    for b in r
      write(io, pipe, " ", b, " ")
    end
    write(io, pipe, "\n")
    if i == length(s)
      # Print bottom line
      write(io, crossbottomleft)
      for _ in 1:length(r) - 1
        write(io, hline_portion, crossbottom)
      end
      if length(r) > 0
        write(io, hline_portion, crossbottomright)
      end
    else
      # Print "normal" line
      # The next row should always be at most as long, but this is not checked
      # in the construction of a tableau, so let's make sure it doesn't mess up
      # printing.
      m = length(r)
      M = length(s[i + 1])
      if m == 0 && M == 0
        write(io, pipe)
        continue
      end
      write(io, crossleft)
      for _ in 1:min(m, M) - 1
        write(io, hline_portion, cross)
      end
      if m == M
        write(io, hline_portion, crossright)
      elseif m < M
        if m != 0
          write(io, hline_portion, cross)
        end
        for _ in m + 1:M - 1
          write(io, hline_portion, crosstop)
        end
        write(io, hline_portion, crosstopright)
      elseif m > M
        if M != 0
          write(io, hline_portion, cross)
        end
        for _ in M + 1:m - 1
          write(io, hline_portion, crossbottom)
        end
        write(io, hline_portion, crossbottomright)
      end
    end
  end
end

################################################################################
#
#  Array-like functionality
#
################################################################################

function Base.size(tab::YoungTableau)
  return size(data(tab))
end

function Base.length(tab::YoungTableau)
  return length(data(tab))
end

function Base.getindex(tab::YoungTableau, i::Int)
  return getindex(data(tab), i)
end

function Base.getindex(tab::YoungTableau, I::Vararg{Int, 2})
  return getindex(getindex(data(tab), I[1]), I[2])
end

push!(tab::YoungTableau{T}, r::Vector{T}) where T <: IntegerUnion = push!(tab.t, r)

################################################################################
#
#  Basic functionality
#
################################################################################

@doc raw"""
    shape(tab::YoungTableau)

Return the shape of the tableau `tab`, i.e. the partition given by the lengths
of the rows of the tableau.
"""
function shape(tab::YoungTableau{T}) where T
  # Line below DOES NOT CHECK that the lengths of the rows are weakly decreasing
  return partition(T[ length(tab[i]) for i = 1:length(tab) ]; check = false)
end

@doc raw"""
    weight(tab::YoungTableau)

Return the weight sequence of the tableau `tab` as an array whose `i`-th element
gives the number of times the integer `i` appears in the tableau.
"""
function weight(tab::YoungTableau)
  @req  is_semistandard(tab)  "Tableau must be (semi-)standard"
  
  if isempty(tab)
    return Int[]
  end

  # Computation of max must be changed if we want to permit non-semi-standard YT
  max = 0
  for i = 1:length(tab)
    if max < tab[i][end]
      max = tab[i][end]
    end
  end

  w = zeros(Int, max)
  for rows in tab
    for box in rows
      w[box] += 1
    end
  end
  return w
end

@doc raw"""
    reading_word(tab::YoungTableau)

Return the reading word of the tableau `tab` as an array, i.e. the word obtained
by concatenating the fillings of the rows, starting from the *bottom* row.

# Examples
```jldoctest
julia> reading_word(young_tableau([[1, 2, 3], [4, 5], [6]]))
6-element Vector{Int64}:
 6
 4
 5
 1
 2
 3
```
"""
function reading_word(tab::YoungTableau)
  w = zeros(Int, sum(shape(tab)))
  k = 0
  for i = length(tab):-1:1
    for j = 1:length(tab[i])
      k += 1
      w[k] = tab[i, j]
    end
  end
  return w
end

################################################################################
#
#  Semistandard tableaux
#
################################################################################

@doc raw"""
    is_semistandard(tab::YoungTableau)

Return `true` if the tableau `tab` is semistandard and `false` otherwise.

A tableau is called **semistandard** if the entries weakly increase along each
row and strictly increase down each column.

See also [`is_standard`](@ref).
"""
function is_semistandard(tab::YoungTableau)
  s = shape(tab)
  if isempty(s)
    return true
  end

  #correct shape
  for i = 1:length(s) - 1
    if s[i] < s[i + 1]
      return false
    end
  end

  #increasing first row
  for j = 2:s[1]
    if tab[1][j] < tab[1][j - 1]
      return false
    end
  end

  #increasing first column
  for i = 2:length(s)
    if tab[i][1] <= tab[i - 1][1]
      return false
    end
  end

  #increasing rows and columns
  for i = 2:length(tab)
    for j = 2:s[i]
      if tab[i][j] < tab[i][j - 1] || tab[i][j] <= tab[i - 1][j]
        return false
      end
    end
  end
  return true
end

@doc raw"""
    semistandard_tableaux(shape::Partition{T}, max_val::T = sum(shape)) where T <: IntegerUnion
    semistandard_tableaux(shape::Vector{T}, max_val::T = sum(shape)) where T <: IntegerUnion

Return an iterator over all semistandard Young tableaux of given shape `shape`
and filling elements bounded by `max_val`.

By default, `max_val` is equal to the sum of the shape partition (the number of
boxes in the Young diagram).

The list of tableaux is in lexicographic order from left to right and top
to bottom.
"""
function semistandard_tableaux(shape::Partition{T}, max_val::T = T(sum(shape))) where T <: IntegerUnion
  return SemiStandardTableaux(shape, max_val)
end

shape(S::SemiStandardTableaux) = S.shape
maximal_value(S::SemiStandardTableaux) = S.max_val

Base.eltype(::Type{SemiStandardTableaux{T}}) where T = YoungTableau{T}

function Base.show(io::IO, S::SemiStandardTableaux)
  print(io, "Iterator over semistandard Young tableaux of shape $(shape(S))")
end

Base.IteratorSize(::Type{<: SemiStandardTableaux}) = Base.SizeUnknown()

function iterate(S::SemiStandardTableaux{T}, state::Nothing = nothing) where T
  shape = Oscar.shape(S)
  max_val = maximal_value(S)

  if max_val < length(shape)
    return nothing
  elseif isempty(shape)
    return young_tableau(Vector{T}[], check = false), (Vector{T}[], 0, 0)
  end

  tab = [fill(T(i), Int(shape[i])) for i in 1:length(shape)]
  m = length(shape)
  n = Int(shape[m])

  return young_tableau([copy(row) for row in tab], check = false), (tab, m, n)
end

@inline function iterate(S::SemiStandardTableaux{T}, state::Tuple{Vector{Vector{T}}, Int, Int}) where T
  shape = Oscar.shape(S)
  max_val = maximal_value(S)
  tab, m, n = state

  if isempty(shape)
    return nothing
  end

  #raise one element by 1
  while !(tab[m][n] < max_val &&
         (n == shape[m] || tab[m][n] < tab[m][n + 1]) &&
         (m == length(shape) || shape[m + 1] < n || tab[m][n] + 1 < tab[m + 1][n]))
    if n > 1
      n -= 1
    elseif m > 1
      m -= 1
      n = Int(shape[m])
    else
      return nothing
    end
  end

  tab[m][n] += 1

  #minimize trailing elements
  if n < shape[m]
    i = m
    j = n + 1
  else
    i = m + 1
    j = 1
  end
  while (i <= length(shape) && j <= shape[i])
    if i == 1
      tab[1][j] = tab[1][j - 1]
    elseif j == 1
      tab[i][1] = tab[i - 1][1] + 1
    else
      tab[i][j] = max(tab[i][j - 1], tab[i - 1][j] + 1)
    end
    if j < shape[i]
      j += 1
    else
      j = 1
      i += 1
    end
  end
  m = length(shape)
  n = Int(shape[end])
  return young_tableau([copy(row) for row in tab], check = false), (tab, m, n)
end

function semistandard_tableaux(shape::Vector{T}, max_val::T = T(sum(shape))) where T <: IntegerUnion
  return semistandard_tableaux(partition(shape, check = false), max_val)
end

@doc raw"""
    semistandard_tableaux(box_num::T, max_val::T = box_num) where T <: Integer

Return an iterator over all semistandard Young tableaux consisting of `box_num`
boxes and filling elements bounded by `max_val`.
"""
function semistandard_tableaux(box_num::T, max_val::T = box_num) where {T <: IntegerUnion}
  # This basically does
  #   Iterators.flatten(semistandard_tableaux(p, max_val) for p in partitions(box_num))
  # but that would print horribly, so we do our own type.
  return SemiStandardTableauxFixedBoxNum(box_num, max_val)
end

box_num(S::SemiStandardTableauxFixedBoxNum) = S.box_num
maximal_value(S::SemiStandardTableauxFixedBoxNum) = S.max_val

Base.eltype(::Type{SemiStandardTableauxFixedBoxNum{T}}) where T = YoungTableau{T}

function Base.show(io::IO, S::SemiStandardTableauxFixedBoxNum)
  print(pretty(io), "Iterator over semistandard Young tableaux with ",
        ItemQuantity(box_num(S), "box", "boxes"))
end

Base.IteratorSize(::Type{<: SemiStandardTableauxFixedBoxNum}) = Base.SizeUnknown()

@inline function iterate(S::SemiStandardTableauxFixedBoxNum{T},
    state::Union{Nothing, Tuple{Partitions{T}, Tuple{Vector{T}, Int, Int}, SemiStandardTableaux{T}, Tuple{Vector{Vector{T}}, Int, Int}}} = nothing) where T
  if !isnothing(state)
    next_tab = iterate(state[3], state[4])
    next_tab !== nothing && return (next_tab[1], (state[1], state[2], state[3], next_tab[2]))
  end
  maximal_value(S) <= 0 && return nothing

  P = (state === nothing ? partitions(box_num(S)) : state[1])
  next_part = (state === nothing ? iterate(P) : iterate(P, state[2]))
  next_part === nothing && return nothing
  Sp = semistandard_tableaux(next_part[1], maximal_value(S))
  next_tab = iterate(Sp)
  while next_tab === nothing
    next_part = iterate(P, next_part[2])
    next_part === nothing && return nothing
    Sp = semistandard_tableaux(next_part[1], maximal_value(S))
    next_tab = iterate(Sp)
  end
  return next_tab[1], (P, next_part[2], Sp, next_tab[2])
end

@doc raw"""
    semistandard_tableaux(s::Partition{T}, weight::Vector{T}) where T <: Integer
    semistandard_tableaux(s::Vector{T}, weight::Vector{T}) where T <: Integer

Return an iterator over all semistandard Young tableaux with shape `s` and given
weight. This requires that `sum(s) = sum(weight)`.
"""
function semistandard_tableaux(s::Partition{T}, weight::Vector{T}) where {T <: IntegerUnion}
  return SemiStandardTableauxFixedShapeAndWeight(s, weight)
end

function semistandard_tableaux(s::Vector{T}, weight::Vector{T}) where {T <: IntegerUnion}
  return SemiStandardTableauxFixedShapeAndWeight(partition(s), weight)
end

shape(S::SemiStandardTableauxFixedShapeAndWeight) = S.shape
weight(S::SemiStandardTableauxFixedShapeAndWeight) = S.weight

Base.eltype(::Type{SemiStandardTableauxFixedShapeAndWeight{T}}) where T = YoungTableau{T}

function Base.show(io::IO, S::SemiStandardTableauxFixedShapeAndWeight)
  print(pretty(io), "Iterator over semistandard Young tableaux of shape ",
        shape(S), " and weight ", weight(S))
end

Base.IteratorSize(::Type{<: SemiStandardTableauxFixedShapeAndWeight}) = Base.SizeUnknown()

@inline function iterate(S::SemiStandardTableauxFixedShapeAndWeight{T}, state::Union{SemiStandardTableauxFixedShapeAndWeightState, Nothing} = nothing) where T
  s = shape(S)
  weight = Oscar.weight(S)

  if isnothing(state)
    state = SemiStandardTableauxFixedShapeAndWeightState{T}()
    if isempty(s)
      state.n = 0
      return young_tableau(Vector{T}[], check = false), state
    end

    state.n = 1
    state.increaseN = true
    state.tab = young_tableau([T[0 for j = 1:s[i]] for i = 1:length(s)], check = false)
    state.boxes_filled = zeros(Int, length(s))
    state.n_used_weight = zeros(Int, length(weight))
    state.row_pointer = zeros(Int, length(weight), maximum(weight))
  end

  # the integer that we currently try to write into the tableau
  n = state.n

  is_zero(n) && return nothing

  # whether we should increase or decrease n in the next iteration of the while loop
  increaseN = state.increaseN
  # the tableau filled in-place
  tab = state.tab
  # boxes_filled[i] is the number of boxes we successfully filled in row i
  boxes_filled = state.boxes_filled

  while n > 0
    if is_zero(weight[n])
      if increaseN
        n += 1
      else
        n -= 1
      end
      continue
    end

    # n is the last possible entry, so we try to fill the remaining boxes with n
    if n == length(weight)
      goodTableau = _is_good_tableau!(n, S, state)
      n -= 1
      increaseN = false
      if goodTableau
        # we discovered a tableau
        state.n = n
        state.increaseN = increaseN
        return young_tableau([copy(row) for row in tab], check = false), state
      end
      # we did not discover a tableau and have to go back to the step n - 1
      continue
    end

    # fill n into the boxes (modifies state in-place)
    state.n = n
    _fill_boxes!(n, S, state)
    n = state.n
    increaseN = state.increaseN
  end

  return nothing
end

function _fill_boxes!(n::Int, S::SemiStandardTableauxFixedShapeAndWeight{T}, state::SemiStandardTableauxFixedShapeAndWeightState{T}) where T
  s = shape(S)
  weight = Oscar.weight(S)

  # the tableau filled in-place
  tab = state.tab
  # boxes_filled[i] is the number of boxes we successfully filled in row i
  boxes_filled = state.boxes_filled
  # n_used_weight[i] is the number of entries equal to i
  n_used_weight = state.n_used_weight
  # row_pointer[k, l] is the row i where we put the lth entry k
  row_pointer = state.row_pointer

  m = n_used_weight[n]
  @assert m == 0 || m == weight[n] "Internal error; this should not happen"
  if m == 0
    # we did not fill a box with the entry n yet, so we start looking for
    # possible boxes in the "top left" corner
    i = 1
    while boxes_filled[i] == s[i]
      i += 1
    end
    j = boxes_filled[i] + 1
  elseif m == weight[n]
    # We already used the entry n as many times as possible, so we came back
    # here from an earlier iteration of the loop
    # We "forget" the last box filled with n and start filling one row below it
    i = row_pointer[n, m] + 1
    if i <= length(s)
      j = boxes_filled[i] + 1
    end
    m -= 1
    boxes_filled[i - 1] -= 1
  end

  while true
    if m == weight[n]
      # we filled n in as many boxes as possible, so we go to the step n + 1
      n_used_weight[n] = m
      state.n = n + 1
      state.increaseN = true
      return nothing
    end

    if i > length(s)
      # We arrived at the bottom of the tableau
      if m == 0
        # We did not fill any box with n yet, so the filling so far is already
        # illegal
        # Go back to step n - 1
        n_used_weight[n] = m
        state.n = n - 1
        state.increaseN = false
        return nothing
      else
        # The last box we filled with n must have been wrong, so we put i in
        # the next row after that and "forget" the box
        i = row_pointer[n, m] + 1
        if i <= length(s)
          j = boxes_filled[i] + 1
        end
        m -= 1
        boxes_filled[i - 1] -= 1
        continue
      end
    end

    if j <= s[i] && (i == 1 || (j <= boxes_filled[i - 1] && n > tab[i - 1][j]))
      # We can fill a box with n
      m += 1
      tab[i][j] = T(n)
      boxes_filled[i] += 1
      row_pointer[n, m] = i
      j += 1
    else
      # We cannot fill any box in row i with n, so we increase i
      i += 1
      if i <= length(s)
        j = boxes_filled[i] + 1
      end
    end
  end
  return nothing
end

# Check whether we get a semistandard tableau by filling n in the remaining boxes
# of state.tab
function _is_good_tableau!(n::Int, S::SemiStandardTableauxFixedShapeAndWeight{T}, state::SemiStandardTableauxFixedShapeAndWeightState{T}) where T
  s = shape(S)
  weight = Oscar.weight(S)

  # the tableau filled in-place
  tab = state.tab
  # boxes_filled[i] is the number of boxes we successfully filled in row i
  boxes_filled = state.boxes_filled

  for i = 1:length(s)
    for j = boxes_filled[i] + 1:s[i]
      tab[i][j] = T(n)
      if i != 1 && tab[i - 1][j] >= n
        # this filling does not give a semistandard tableau
        return false
      end
    end
  end
  return true
end

################################################################################
#
#  Standard tableaux
#
################################################################################

@doc raw"""
    is_standard(tab::YoungTableau)

Return `true` if the tableau `tab` is standard and `false` otherwise.

A tableau is called **standard** if it is semistandard and the entries
are in bijection with `1, ..., n`, where `n` is the number of boxes.

See also [`is_semistandard`](@ref).
"""
function is_standard(tab::YoungTableau)
  s = shape(tab)
  if isempty(s)
    return true
  end

  #correct shape
  for i = 1:length(s) - 1
    if s[i] < s[i + 1]
      return false
    end
  end

  #contains all numbers from 1 to n
  n = sum(s)
  numbs = falses(n)
  for i = 1:length(s)
    for j = 1:s[i]
      if tab[i][j] < 1 || tab[i][j] > n
        return false
      end
      numbs[tab[i][j]] = true
    end
  end
  if false in numbs
    return false
  end

  #increasing first row
  for j = 2:s[1]
    if tab[1][j] <= tab[1][j - 1]
      return false
    end
  end

  #increasing first column
  for i = 2:length(s)
    if tab[i][1] <= tab[i - 1][1]
      return false
    end
  end

  #increasing rows and columns
  for i = 2:length(s)
    for j = 2:s[i]
      if tab[i][j] <= tab[i][j - 1] || tab[i][j] <= tab[i - 1][j]
        return false
      end
    end
  end
  return true
end

@doc raw"""
    standard_tableaux(s::Partition)
    standard_tableaux(s::Vector{Integer})

Return an iterator over all standard Young tableaux of a given shape `s`.
"""
function standard_tableaux(s::Partition)
  return StandardTableaux(s)
end

shape(S::StandardTableaux) = S.shape

Base.eltype(::Type{StandardTableaux{T}}) where T = YoungTableau{T}

function Base.show(io::IO, S::StandardTableaux)
  print(io, "Iterator over standard Young tableaux of shape $(shape(S))")
end

Base.IteratorSize(::Type{<: StandardTableaux}) = Base.SizeUnknown()

@inline function iterate(S::StandardTableaux{T}, state::Union{Nothing, StandardTableauxState{T}} = nothing) where T
  s = shape(S)
  n_max = sum(s)

  if isnothing(state)
    if isempty(s)
      tab = young_tableau(Vector{T}[], check = false)
      return tab, StandardTableauxState{T}(-1, 0, 0, tab, Int[], Int[])
    end

    tab = young_tableau([ [T(0) for j in 1:s[i]] for i in 1:length(s)], check = false)
    sub_s = zeros(Int, length(s))
    tab[1][1] = T(1)
    sub_s[1] = 1
    tracker_row = zeros(Int, n_max)
    tracker_row[1] = 1
    n = 1
    i = 1
    j = 2

    state = StandardTableauxState{T}(n, i, j, tab, sub_s, tracker_row)
  end

  # Having this in a separate functions seems to make the compiler happy and
  # saves allocations
  fl, t = _discover_tableau!(S, state)
  !fl && return nothing
  return t, state
end

function _discover_tableau!(S::StandardTableaux{T}, state::StandardTableauxState{T}) where T
  s = shape(S)
  n_max = sum(s)

  n = state.n
  i = state.i
  j = state.j
  tab = state.tab
  sub_s = state.sub_s
  tracker_row = state.tracker_row

  tab2 = tab
  stop = false
  while n > 0 && !stop
    if n == n_max || i > length(s)
      if n == n_max
        tab2 = young_tableau([copy(row) for row in tab], check = false)
        stop = true
      end
      tab[tracker_row[n]][sub_s[tracker_row[n]]] = T(0)
      i = tracker_row[n] + 1
      if i <= length(s)
        j = sub_s[i] + 1
      end
      n -= 1
      sub_s[i - 1] -= 1
    elseif j <= s[i] && (i == 1 || j <= sub_s[i - 1])
      n += 1
      tab[i][j] = T(n)
      sub_s[i] += 1
      tracker_row[n] = i
      i = 1
      j = sub_s[1] + 1
    else
      i += 1
      if i <= length(s)
        j = sub_s[i] + 1
      end
    end
  end

  state.n = n
  state.i = i
  state.j = j
  state.tab = state.tab
  state.sub_s = sub_s
  state.tracker_row = tracker_row

  return stop, tab2
end

function standard_tableaux(s::Vector{T}) where T <: Integer
  return standard_tableaux(partition(s, check = false))
end

@doc raw"""
    standard_tableaux(n::IntegerUnion)

Return an iterator over all standard Young tableaux with `n` boxes.
"""
function standard_tableaux(n::IntegerUnion)
  # This basically does
  #   Iterators.flatten(standard_tableaux(p) for p in partitions(n))
  # but that would print horribly, so we do our own type.
  return StandardTableauxFixedBoxNum(n)
end

box_num(S::StandardTableauxFixedBoxNum) = S.box_num

Base.eltype(::Type{StandardTableauxFixedBoxNum{T}}) where T = YoungTableau{T}

function Base.show(io::IO, S::StandardTableauxFixedBoxNum)
  print(pretty(io), "Iterator over standard Young tableaux with ",
        ItemQuantity(box_num(S), "box", "boxes"))
end

Base.IteratorSize(::Type{<: StandardTableauxFixedBoxNum}) = Base.SizeUnknown()

@inline function iterate(S::StandardTableauxFixedBoxNum{T},
    state::Union{Nothing, Tuple{Partitions{T}, Tuple{Vector{T}, Int, Int}, StandardTableaux{T}, StandardTableauxState{T}}} = nothing) where T
  if !isnothing(state)
    next_tab = iterate(state[3], state[4])
    next_tab !== nothing && return (next_tab[1], (state[1], state[2], state[3], next_tab[2]))
  end

  P = (state === nothing ? partitions(box_num(S)) : state[1])
  next_part = (state === nothing ? iterate(P) : iterate(P, state[2]))
  next_part === nothing && return nothing
  Sp = standard_tableaux(next_part[1])
  next_tab = iterate(Sp)
  return next_tab[1], (P, next_part[2], Sp, next_tab[2])
end

################################################################################
#
#  Hook length
#
################################################################################

@doc raw"""
    hook_length(tab::YoungTableau, i::Integer, j::Integer)
    hook_length(lambda::Partition, i::Integer, j::Integer)

Return the hook length of the box with coordinates `(i, j)` in the Young tableau
`tab` respectively the Young diagram of shape `lambda`.

The **hook length** of a box is the number of boxes to the right in the same
row + the number of boxes below in the same column + 1.

See also [`hook_lengths`](@ref).
"""
hook_length

function hook_length(lambda::Partition, i::Integer, j::Integer)
  h = lambda[i] - j + 1
  k = i + 1
  while k <= length(lambda) && lambda[k] >= j
    k += 1
    h += 1
  end
  return h
end

function hook_length(tab::YoungTableau, i::Integer, j::Integer)
  return hook_length(shape(tab), i, j)
end

@doc raw"""
    hook_lengths(lambda::Partition)

Return the Young tableau of shape `lambda` in which the entry at position
`(i, j)` is equal to the hook length of the corresponding box.

See also [`hook_length`](@ref).
"""
function hook_lengths(lambda::Partition)
  if isempty(lambda)
    return young_tableau(Vector{Int}[], check = false)
  end
  tab = [ [hook_length(lambda, i, j) for j in 1:lambda[i]] for i in 1:length(lambda) ]
  return young_tableau(tab, check = false)
end

@doc raw"""
    number_of_standard_tableaux(lambda::Partition)

Return the number of standard Young tableaux of shape `lambda`.
"""
function number_of_standard_tableaux(lambda::Partition)
  n = sum(lambda)
  h = factorial(ZZ(n))
  for i = 1:length(lambda)
    for j = 1:lambda[i]
      h = div(h, ZZ(hook_length(lambda, i, j)))
    end
  end
  return h
end

################################################################################
#
#  Schensted insertion
#
################################################################################

@doc raw"""
    schensted(sigma::Vector{<:IntegerUnion})
    schensted(sigma::PermGroupElem)

Return the pair of standard Young tableaux (the insertion and the recording
tableau) corresponding to the permutation `sigma` under the Robinson-Schensted
correspondence.

# Examples
```jldoctest
julia> P, Q = schensted([3, 1, 6, 2, 5, 4]);

julia> P
+---+---+---+
| 1 | 2 | 4 |
+---+---+---+
| 3 | 5 |
+---+---+
| 6 |
+---+

julia> Q
+---+---+---+
| 1 | 3 | 5 |
+---+---+---+
| 2 | 4 |
+---+---+
| 6 |
+---+

```
"""
schensted

function schensted(sigma::Vector{T}) where T <: IntegerUnion
  if isempty(sigma)
    return young_tableau(Vector{T}[], check = false), young_tableau(Vector{T}[], check = false)
  end
  P = young_tableau([[sigma[1]]], check = false)
  Q = young_tableau([[1]], check = false)
  for i = 2:length(sigma)
    bump!(P, sigma[i], Q, i)
  end
  return P, Q
end

function schensted(::Type{T}, sigma::PermGroupElem) where T <: IntegerUnion
  return schensted(Vector{T}(sigma))
end

schensted(sigma::PermGroupElem) = schensted(Int, sigma)

@doc raw"""
    bump!(tab::YoungTableau, x::Int)

Insert the integer `x` into the tableau `tab` according to the bumping
algorithm by applying the Schensted insertion.
"""
function bump!(tab::YoungTableau, x::Integer)
  if isempty(tab)
    push!(tab, [x])
    return tab
  end

  i = 1
  while i <= length(tab)
    if tab[i, length(tab[i])] <= x
      push!(tab[i], x)
      return tab
    end
    j = 1
    while j <= length(tab[i])
      if tab[i, j] > x
        temp = x
        x = tab[i, j]
        tab[i][j] = temp
        i += 1
        break
      end
      j += 1
    end
  end
  push!(tab, [x])
  return tab
end

@doc raw"""
    bump!(tab::YoungTableau, x::Integer, Q::YoungTableau, y::Integer)

Insert the integer `x` into `tab` according to the bumping algorithm by applying
the Schensted insertion and insert the integer `y` into `Q` at the same position
as `x` in `tab`.
"""
function bump!(tab::YoungTableau, x::Integer, Q::YoungTableau, y::Integer)
  if isempty(tab)
    push!(tab, [x])
    push!(Q, [x])
    return tab, Q
  end
  i = 1
  while i <= length(tab)
    if tab[i, length(tab[i])] <= x
      push!(tab[i], x)
      push!(Q[i], y)
      return tab
    end
    j = 1
    while j <= length(tab[i])
      if tab[i, j] > x
        temp = x
        x = tab[i, j]
        tab[i][j] = temp
        i += 1
        break
      end
      j += 1
    end
  end
  push!(tab, [x])
  push!(Q, [y])

  return tab, Q
end
