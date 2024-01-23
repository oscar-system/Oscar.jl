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
#  Constructor and basic functionality
#
################################################################################

@doc raw"""
    young_tableau(v::Vector{Vector{T}}) where T <: IntegerUnion

Return the Young tableau given by `v` as an object of type `YoungTableau{T}`.

See[`YoungTableau`](@ref) for more details and examples.
"""
function young_tableau(v::Vector{Vector{T}}) where T <: IntegerUnion
  return YoungTableau(v)
end

data(tab::YoungTableau) = tab.t

# supercompact and one-line printing
function Base.show(io::IO, tab::YoungTableau)
  print(io, "Young tableau")
  # TODO: is there meaningful information to add in one-line mode?
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
  write(io, hline_portion, crosstopright)
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
      write(io, hline_portion, crossbottomright)
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
    shape(tab::YoungTableau{T})

Return the shape of a tableau, i.e. the partition given by the lengths of the
rows of the tableau.
"""
function shape(tab::YoungTableau{T}) where T
  return partition(T[ length(tab[i]) for i = 1:length(tab) ])
end

@doc raw"""
    weight(tab::YoungTableau)

The **weight** of a tableau is the number of times each number appears in the
tableau. The return value is an array whose `i`-th element gives the number of
times the integer `i` appears in the tableau.
"""
function weight(tab::YoungTableau)
  if isempty(tab)
    return Int[]
  end

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

The **reading word** of a tableau is the word obtained by concatenating the
fillings of the rows, starting from the *bottom* row. The word is here
returned as an array.

# Examples
```jldoctest
julia> reading_word(young_tableau([ [1, 2, 3] , [4, 5] , [6] ]))
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

A tableau is called **semistandard** if the entries weakly increase along each
row and strictly increase down each column.
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
    semistandard_tableaux(shape::Partition{T}, max_val::T = sum(shape)) where T <: Integer

Return a list of all semistandard tableaux of given shape and filling
elements bounded by `max_val`. By default, `max_val` is equal to the sum
of the shape partition (the number of boxes in the Young diagram). The
list of tableaux is in lexicographic order from left to right and top
to bottom.
"""
function semistandard_tableaux(shape::Partition{T}, max_val::T = sum(shape)) where T <: Integer
  SST = Vector{YoungTableau{T}}()
  len = length(shape)
  if max_val < len
    return SST
  elseif len == 0
    push!(SST, young_tableau(Vector{T}[]))
    return SST
  end
  tab = [Array{T}(fill(i, shape[i])) for i = 1:len]
  m = len
  n = shape[m]

  while true
    push!(SST, young_tableau([copy(row) for row in tab]))

    #raise one element by 1
    while !(tab[m][n] < max_val &&
           (n == shape[m] || tab[m][n] < tab[m][n + 1]) &&
           (m == len || shape[m + 1] < n || tab[m][n] + 1 < tab[m + 1][n]))
      if n > 1
        n -= 1
      elseif m > 1
        m -= 1
        n = shape[m]
      else
        return SST
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
    while (i <= len && j <= shape[i])
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
    m = len
    n = shape[len]
  end
end

@doc raw"""
    semistandard_tableaux(shape::Partition{T}, max_val::T = sum(shape)) where T <: Integer

Shortcut for `semistandard_tableaux(Partition(shape), max_val)`.
"""
function semistandard_tableaux(shape::Vector{T}, max_val::T = sum(shape)) where T <: Integer
  return semistandard_tableaux(partition(shape), max_val)
end

@doc raw"""
    semistandard_tableaux(box_num::T, max_val::T = box_num) where T <: Integer

Return a list of all semistandard tableaux consisting of `box_num`
boxes and filling elements bounded by `max_val`.
"""
function semistandard_tableaux(box_num::T, max_val::T = box_num) where T <: Integer
  @req box_num >= 0 "box_num >= 0 required"
  SST = Vector{YoungTableau{T}}()
  if max_val <= 0
    return SST
  end
  shapes = partitions(box_num)

  for s in shapes
    if max_val >= length(s)
      append!(SST, semistandard_tableaux(data(s), max_val))
    end
  end

  return SST
end


@doc raw"""
    semistandard_tableaux(s::Vector{T}, weight::Vector{T}) where T <: Integer

Return a list of all semistandard tableaux with shape `s` and given weight. This
requires that `sum(s) = sum(weight)`.
"""
function semistandard_tableaux(s::Vector{T}, weight::Vector{T}) where T <: Integer
  n_max = sum(s)
  @req n_max == sum(weight) "sum(s) == sum(weight) required"

  tabs = Vector{YoungTableau}()
  if isempty(s)
    push!(tabs, young_tableau(Vector{Int}[]))
    return tabs
  end
  ls = length(s)

  tab = young_tableau([ [0 for j = 1:s[i]] for i = 1:length(s)])
  sub_s = zeros(Integer, length(s))

  #tracker_row = zeros(Integer, n_max)

  function rec_sst!(n::Integer)

    #fill the remaining boxes if possible, else set them to 0
    if n == length(weight)
      for i = 1:ls
        for j = sub_s[i] + 1:s[i]
          tab[i][j] = n
          if i != 1 && tab[i - 1][j] == n
            for k = 1:i
              for l = sub_s[k] + 1:s[k]
                tab[i][j] = 0
              end
            end
            return
          end
        end
      end
      push!(tabs, young_tableau([copy(row) for row in tab]))

      return

      #skip to next step if weight[n] == 0
    elseif weight[n] == 0
      rec_sst!(n + 1)
      return
    end

    #here starts the main part of the function
    tracker_row = zeros(Integer, weight[n])
    i = 1
    while sub_s[i] == s[i]
      i += 1
    end
    j = sub_s[i] + 1

    m = 0
    while m >= 0
      if m == weight[n]   # jump to next recursive step
        rec_sst!(n + 1)
        tab[tracker_row[m]][sub_s[tracker_row[m]]] = 0
        i = tracker_row[m] + 1
        if i <= ls
          j = sub_s[i] + 1
        end
        m -= 1
        sub_s[i - 1] -= 1

      elseif i > ls
        if m == 0
          return
        else
          tab[tracker_row[m]][sub_s[tracker_row[m]]] = 0
          i = tracker_row[m] + 1
          if i <= ls
            j = sub_s[i] + 1
          end
          m -= 1
          sub_s[i - 1] -= 1
        end

      elseif j <= s[i] && (i == 1 || (j <= sub_s[i - 1] && n > tab[i - 1][j]))  #add an entry
        m += 1
        tab[i][j] = n
        sub_s[i] += 1
        tracker_row[m] = i
        j += 1

      else #move pointerhead
        i += 1
        if i <= ls
          j = sub_s[i] + 1
        end
      end
    end #while

  end #rec_sst!()

  rec_sst!(1)
  return tabs
end

function semistandard_tableaux(s::Partition{T}, weight::Partition{T}) where T <: Integer
  return semistandard_tableaux(Vector{T}(s), Vector{T}(weight))
end

################################################################################
#
#  Standard tableaux
#
################################################################################

@doc raw"""
    is_standard(tab::YoungTableau)

A tableau is called **standard** if it is semistandard and the entries
are in bijection with ``1, …, n``, where ``n`` is the number of boxes.
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
      if tab[i][j] > n
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

Return a list of all standard tableaux of a given shape.
"""
function standard_tableaux(s::Partition)
  tabs = Vector{YoungTableau}()
  if isempty(s)
    push!(tabs, young_tableau(Vector{Int}[]))
    return tabs
  end
  n_max = sum(s)
  ls = length(s)

  tab = young_tableau([ [0 for j = 1:s[i]] for i = 1:length(s)])
  sub_s = [0 for i = 1:length(s)]
  tab[1][1] = 1
  sub_s[1] = 1
  tracker_row = [0 for i = 1:n_max]
  tracker_row[1] = 1

  n = 1
  i = 1
  j = 2

  while n > 0
    if n == n_max || i > ls
      if n == n_max
        push!(tabs, young_tableau([copy(row) for row in tab]))
      end
      tab[tracker_row[n]][sub_s[tracker_row[n]]] = 0
      i = tracker_row[n] + 1
      if i <= ls
        j = sub_s[i] + 1
      end
      n -= 1
      sub_s[i - 1] -= 1
    elseif j <= s[i] && (i == 1 || j <= sub_s[i - 1])
      n += 1
      tab[i][j] = n
      sub_s[i] += 1
      tracker_row[n] = i
      i = 1
      j = sub_s[1] + 1
    else
      i += 1
      if i <= ls
        j = sub_s[i] + 1
      end
    end
  end

  return tabs
end

function standard_tableaux(s::Vector{T}) where T <: Integer
  return standard_tableaux(partition(s))
end

@doc raw"""
    standard_tableaux(n::Integer)

Return a list of all standard tableaux with n boxes.
"""
function standard_tableaux(n::Integer)
  @req n >= 0 "n >= 0 required"
  ST = Vector{YoungTableau}()
  for s in partitions(n)
    append!(ST, standard_tableaux(s))
  end
  return ST
end

################################################################################
#
#  Hook length
#
################################################################################

@doc raw"""
    hook_length(lambda::Partition, i::Integer, j::Integer)

Consider the Young diagram of a partition ``λ``. The **hook length** of a
box, is the number of boxes to the right in the same row + the number
of boxes below in the same column + 1. The function returns the hook
length of the box with coordinates `(i, j)`. The functions assumes that
the box exists.
"""
function hook_length(lambda::Partition, i::Integer, j::Integer)
  h = lambda[i] - j + 1
  k = i + 1
  while k <= length(lambda) && lambda[k] >= j
    k += 1
    h += 1
  end
  return h
end

@doc raw"""
    hook_length(tab::YoungTableau, i::Integer, j::Integer)

Shortcut for `hook_length(shape(tab), i, j)`.
"""
function hook_length(tab::YoungTableau, i::Integer, j::Integer)
  return hook_length(shape(tab), i, j)
end

@doc raw"""
    hook_lengths(lambda::Partition)

Return the tableau of shape ``λ`` in which the entry at position `(i, j)`
is equal to the hook length of the corresponding box.
"""
function hook_lengths(lambda::Partition)
  if isempty(lambda)
    return young_tableau(Vector{Int}[])
  end
  tab = [ [hook_length(lambda, i, j) for j in 1:lambda[i]] for i in 1:length(lambda) ]
  return young_tableau(tab)
end

@doc raw"""
    number_of_standard_tableaux(lambda::Partition)

Return the number $f^λ$ of standard tableaux of shape ``λ`` using the hook length formula

$$f^λ = \frac{n!}{\prod_{i, j} h_λ(i, j)}, $$

where the product is taken over all boxes in the Young diagram of ``λ`` and
``h_λ`` denotes the hook length of the box (i, j).

# References
1. Wikipedia, [Hook length formula](https://en.wikipedia.org/wiki/Hook_length_formula).
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
    schensted(sigma::Vector{Integer})
    schensted(sigma::Perm{T})

The Robinson–Schensted correspondence is a bijection between
permutations and pairs of standard Young tableaux of the same shape. For
a permutation sigma (given as an array), this function performs the
Schensted algorithm and returns the corresponding pair of standard
tableaux (the insertion and recording tableaux).

# Examples
```jldoctest
julia> P, Q = schensted([3, 1, 6, 2, 5, 4]);

julia> P
[[1, 2, 4], [3, 5], [6]]

julia> Q
[[1, 3, 5], [2, 4], [6]]
```

# References
1. Wikipedia, [Robinson–Schensted correspondence](https://en.wikipedia.org/wiki/Robinson–Schensted_correspondence)
"""
function schensted(sigma::Vector{T}) where T <: Integer
  if isempty(sigma)
    return young_tableau(Vector{T}[]), young_tableau(Vector{T}[])
  end
  P = young_tableau([[sigma[1]]])
  Q = young_tableau([[1]])
  for i = 2:length(sigma)
    bump!(P, sigma[i], Q, i)
  end
  return P, Q
end

function schensted(sigma::Perm{T}) where T <: Integer
  return schensted(sigma.d)
end

@doc raw"""
    bump!(tab::YoungTableau, x::Int)

Inserts the integer `x` into the tableau `tab` according to the bumping
algorithm by applying the Schensted insertion.

# References
1. Wolfram MathWorld, [Bumping Algorithm](https://mathworld.wolfram.com/BumpingAlgorithm.html)
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

Inserts `x` into tab according to the bumping algorithm by applying the
Schensted insertion. Traces the change with `Q` by inserting `y` at the same
position in `Q` as `x` in tab.
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
