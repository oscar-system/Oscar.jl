# TODO: correct typing

function cylinders(cylinder_diagram::CylinderDiagram)
  result = Vector{Vector{Vector{Int}}}()
  for i in 1:(cylinder_diagram.cycles_count)
    push!(result, [cylinder_diagram.bot[i], cylinder_diagram.top[i]])
  end
  return result
end

# computes the system of equations whose solutions are the realizable cylinder
# coordinates
function system_of_equations(cylinder_diagram::CylinderDiagram)
  M = zero_matrix(ZZ, cylinder_diagram.cycles_count, cylinder_diagram.separatrix_count)
  for (i, (bot, top)) in enumerate(cylinders(cylinder_diagram))
    for t in top
      M[i, t + 1] -= 1
    end
    for b in bot
      M[i, b + 1] += 1
    end
  end
  return M
end


# based on Donald E. Knuth The Art of Computer Programming Algorithm H p.300
# implementation from Sage
function product_gray_code(m::Vector{Int})
  n = 0
  k = 0
  # assumes each radix is >= 2
  new_m = Vector{Int}()
  mm = Vector{Int}()
  for i in m
    if i <= 0
      throw(ArgumentError("accepts only positive radices"))
    elseif i >= 2
      push!(new_m, i - 1)
      push!(mm, k)
      n += 1
    end
    k += 1
  end

  m = new_m
  f = collect(1:(n + 1))
  o = ones(Int, n)
  a = zeros(Int, n)
  results = Vector{Tuple{Int,Int}}()

  j = f[1]
  while j != n + 1
    f[1] = 1
    oo = o[j]
    a[j] += oo
    if a[j] == 0 || a[j] == m[j]
      f[j] = f[j + 1]
      f[j + 1] = j + 1
      o[j] = -oo
    end

    push!(results, (mm[j] + 1, oo))
    j = f[1]
  end

  return results
end

function horizontal_symmetry(cyl_diagram::CylinderDiagram)
  bot = [reverse(b) for b in cyl_diagram.bot]
  top = [reverse(t) for t in cyl_diagram.top]
  return CylinderDiagram(top, bot)
end

function vertical_symmetry(cyl_diagram::CylinderDiagram)
  bot = [reverse(b) for b in cyl_diagram.bot]
  top = [reverse(t) for t in cyl_diagram.top]
  return CylinderDiagram(bot, top)
end

function inverse(cyl_diagram::CylinderDiagram)
  return CylinderDiagram(cyl_diagram.top, cyl_diagram.bot)
end

# get origamis from cylinder diagram and realizable lengths
# from surface_dynamics by Vincent, Delecroix et al.
function origami_from_cylinder_coordinates(cyl_diagram::CylinderDiagram,
  lengths::Vector{Int}, heights::Vector{Int}, deg::Int)
  # the total width of each cylinder is the sum of the lengths of the separatrices
  widths = [sum(lengths[i + 1] for i in bot) for bot in cyl_diagram.bot]
  areas = [heights[i] * widths[i] for i in 1:(cyl_diagram.cycles_count)]
  S = symmetric_group(deg)

  v = [0]
  for a in areas
    push!(v, v[end] + a)
  end

  # prepare for twist
  sep_bottom_pos = zeros(Int, cyl_diagram.separatrix_count)
  for (i, bot) in enumerate(cyl_diagram.bot)
    w = 0
    for j in bot
      sep_bottom_pos[j + 1] = v[i] + w
      w += lengths[j + 1]
    end
  end

  # horizontal permutation
  lx = collect(1:v[end])
  for i in 1:(cyl_diagram.cycles_count)
    for j in v[i]:widths[i]:(v[i + 1] - 1)
      lx[j + widths[i]] = j
    end
  end
  lx = perm(S, [x + 1 for x in lx])

  ly = Int[]
  for i in 1:(cyl_diagram.cycles_count)
    append!(ly, (v[i] + widths[i]):(v[i + 1] - 1))
    append!(ly, zeros(Int64, widths[i]))
  end

  for (i, top_seps) in enumerate(cyl_diagram.top)
    top = Int[]
    rev_top_seps = reverse(top_seps)
    for k in rev_top_seps
      append!(top, sep_bottom_pos[k + 1]:(sep_bottom_pos[k + 1] + lengths[k + 1] - 1))
    end
    ly[(v[i + 1] - widths[i] + 1):(v[i + 1])] = top
  end

  no_twist = origami_disconnected(lx, perm(S, [x + 1 for x in ly]), deg)
  results = [no_twist]

  ly = [x + 1 for x in ly]

  for (i, o) in product_gray_code(widths)
    if o == 1
      insert!(ly, v[i + 1] - widths[i] + 1, popat!(ly, v[i + 1]))
    else
      insert!(ly, v[i + 1], popat!(ly, v[i + 1] - widths[i] + 1))
    end
    new_entry = origami_disconnected(lx, perm(S, ly), deg)
    push!(results, new_entry)
  end

  return results
end


"""
Enumerate all realizable separatrix lengths and cylinder heights
for square-tiled surfaces of total area `degree`.

This is a translation of surface_dynamics' `CylinderDiagram.widths_and_heights_iterator`.
https://flatsurf.github.io/surface-dynamics/, last changed on 07/08/2026.
"""
function possible_lengths_and_heights(cyl_diagram::CylinderDiagram, degree::Int)

  m = system_of_equations(cyl_diagram)

  nseps = cyl_diagram.separatrix_count
  ncyls = cyl_diagram.cycles_count

  # Compute lower bounds for individual separatrix lengths.
  min_lengths = ones(Int, nseps)

  for i in 1:ncyls
    pos_m = Int[]
    pos_p = Int[]

    for j in 1:nseps
      if m[i, j] == -1
        push!(pos_m, j)
      elseif m[i, j] == 1
        push!(pos_p, j)
      end
    end

    if length(pos_m) == 1
      j = pos_m[1]
      min_lengths[j] = max(min_lengths[j], length(pos_p))
    end

    if length(pos_p) == 1
      j = pos_p[1]
      min_lengths[j] = max(min_lengths[j], length(pos_m))
    end
  end

  # Minimum possible width of each cylinder.
  min_widths = Int[]

  for (bot, top) in cylinders(cyl_diagram)
    bot_min = sum(min_lengths[j + 1] for j in bot)
    top_min = sum(min_lengths[j + 1] for j in top)
    push!(min_widths, max(bot_min, top_min))
  end

  result = Vector{Tuple{Vector{Int},Vector{Int}}}()

  # Distribute total area among cylinders:
  #     a[1] + ... + a[ncyls] = degree.
  for a in compositions(degree, ncyls)

    if !all(a[i] >= min_widths[i] for i in 1:ncyls)
      continue
    end

    # For each cylinder: area = width * height.
    width_choices = Vector{Vector{Int}}()

    for i in 1:ncyls
      choices = [d for d in divisors(a[i]) if d >= min_widths[i]]
      push!(width_choices, choices)
    end

    # iterate over Cartesian product
    for widths in Iterators.product(width_choices...)

      # widths[i] runs over the divisors of a[i], so the division is exact
      heights = [div(a[i], widths[i]) for i in 1:ncyls]

      # For each cylinder, distribute its width among the
      # separatrices on its bottom boundary.
      bottom_seps = [bot for (bot, _) in cylinders(cyl_diagram)]

      length_choices = Vector{Vector{Vector{Int}}}()

      for i in 1:ncyls
        # all compositions of widths[i] into a sum of
        # length(bottom_seps[i]) many positive integers
        push!(length_choices, collect(compositions(widths[i], length(bottom_seps[i]))))
      end

      # Cartesian product over the composition choices.
      function recurse_lengths!(cyl_index::Int, lengths::Vector{Int})
        if cyl_index > ncyls
          # Check the cylinder-width relations.
          l_matrix = matrix(ZZ, nseps, 1, lengths)
          if is_zero(m * l_matrix)
            push!(result, (copy(lengths), copy(heights)))
          end
          return
        end

        for local_lengths in length_choices[cyl_index]
          new_lengths = copy(lengths)

          for j in eachindex(bottom_seps[cyl_index])
            sep = bottom_seps[cyl_index][j]
            # Separatrix labels are 0-based in CylinderDiagram.
            new_lengths[sep + 1] = local_lengths[j]
          end

          recurse_lengths!(cyl_index + 1, new_lengths)
        end
      end

      recurse_lengths!(1, zeros(Int, nseps))
    end
  end

  return result
end

# so far only genus 2 & 3 supported
# other strata may still be explored manually by using the read_cylinder_diagrams function
# genus 3 still contains errors and is thus not available right now, but will hopefully be fixed soon.
function origamis(stratum::Vector{Int}, degree::Int)::Set{Origami}
  file_name = ""
  if stratum == [1, 1]
    file_name = "h11.dat"
  elseif stratum == [2]
    file_name = "h2.dat"
  elseif stratum == [1, 1, 1, 1]
    file_name = "h1111.dat"
  elseif stratum == [4]
    file_name = "h4.dat"
  elseif stratum == [1, 1, 2] || stratum == [1, 2, 1] || stratum == [2, 1, 1]
    file_name = "h112.dat"
  elseif stratum == [2, 2]
    file_name = "h22.dat"
  else
    error("Stratum not supported!")
  end

  file_path = joinpath(@__DIR__, "..", "cylinder_diagrams", file_name)
  diagrams::Vector{CylinderDiagram} = read_cylinder_diagrams(file_path)


  representatives = Dict{CanonicalOrigamiKey, Origami}()
  for cc in diagrams
    for c in Set([cc, horizontal_symmetry(cc), vertical_symmetry(cc), inverse(cc)])
      lengths_heights = possible_lengths_and_heights(c, degree)
      for entry in lengths_heights
        oris = origami_from_cylinder_coordinates(c, entry[1], entry[2], degree)
        for o in oris
          key = canonical_origami_key(o)
          if !haskey(representatives, key)
            representatives[key] = o
          end
        end
      end
    end
  end
  return Set(values(representatives))
end

function read_cylinder_diagrams(filename::String)
  # Read the entire file content as a string
  content = open(filename, "r") do io
    read(io, String)
  end
  # Remove leading/trailing whitespace and outer brackets
  content = strip(content)
  if startswith(content, "[") && endswith(content, "]")
    content = content[2:(end - 1)]
  else
    error("File content must start with '[' and end with ']'")
  end
  # Split content into diagrams using commas not inside parentheses
  diagrams = split_diagrams(content)
  cylinders = Vector{CylinderDiagram}()
  for diagram in diagrams
    # Split diagram into permutations separated by spaces
    permutations = split(strip(diagram))
    bot_cycles = Vector{Vector{Int}}()
    top_cycles = Vector{Vector{Int}}()
    for permutation in permutations
      # Split each permutation into bot and top parts
      bot_str, top_str = split(permutation, "-"; limit=2)
      bot_cycle = parse_cycle(bot_str)
      top_cycle = parse_cycle(top_str)
      push!(bot_cycles, bot_cycle)
      push!(top_cycles, top_cycle)
    end
    cd = CylinderDiagram(bot_cycles, top_cycles)
    push!(cylinders, cd)
  end
  return cylinders
end

function split_diagrams(content::AbstractString)
  diagrams = String[]
  idx = 1
  len_content = length(content)
  while idx <= len_content
    diagram = ""
    depth = 0
    while idx <= len_content
      c = content[idx]
      if c == '('
        depth += 1
      elseif c == ')'
        depth -= 1
      elseif c == ',' && depth == 0
        idx += 1  # Skip the comma
        break
      end
      diagram *= c
      idx += 1
    end
    push!(diagrams, strip(diagram))
  end
  return diagrams
end

function parse_cycle(cycle_str::AbstractString)
  # Parse a cycle string of the form "(numbers)"
  cycle_str = strip(cycle_str)
  if startswith(cycle_str, "(") && endswith(cycle_str, ")")
    numbers_str = cycle_str[2:(end - 1)]
    numbers = split(numbers_str, ",")
    if isempty(numbers_str)
      return Int[]
    else
      return [parse(Int, strip(n)) for n in numbers]
    end
  else
    error("Invalid cycle format: $cycle_str")
  end
end
