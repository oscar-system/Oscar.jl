# TODO: correct typing

# the bijection between the cycles of the bot/top permutation
# is given implicitly by using lists and considering the order
# of the cycles
struct CylinderDiagram
  bot::Vector{Vector{Int}}
  top::Vector{Vector{Int}}
  cycles_count::Int
  separatrix_count::Int
  function CylinderDiagram(bot::Vector{Vector{Int}}, top::Vector{Vector{Int}})
    max = 1
    for cycle in bot
      for i in cycle
        if i > max
          max = i
        end
      end
    end
    new(bot, top, length(bot), max + 1)
  end
end

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
      M[i, t + 1] = 1
    end
    for b in bot
      M[i, b + 1] -= 1
    end
  end
  return M
end

# positive integer linear combinations of the computed rays
# give possible lengths for the separatrix of the cylinder diagram
# i.e. the resulting surface exists
function compute_rays(equations, separatrix_count::Int)::Vector{Vector{Int}}
  if is_zero(equations)
    n = separatrix_count
    return [Int[i == j for j in 1:n] for i in 1:n]
  end
  neg_identity = -identity_matrix(ZZ, separatrix_count)
  cone = cone_from_inequalities(neg_identity, equations)
  ray_list = rays(cone)
  return map(inner -> map(Int, inner), ray_list)
end

# generate all possible linear combinations
# m: number of parameters in linear combination, parameters have integer values from 1 to n
function all_combinations(n::Int, m::Int)::Vector{Vector{Int}}
  combinations = Vector{Vector{Int}}([[]])

  for _ in 1:m
    new_combinations = Vector{Vector{Int}}()
    for current in combinations
      for i in 1:n
        push!(new_combinations, [current..., i])
      end
    end
    combinations = new_combinations
  end
  return combinations
end

# lengths are realizable if they are solutions to the equations for the lengths
# arising from the cylinder diagram
function realizable_lengths(rays, n)::Vector{Vector{Int}}
  l = length(rays)
  combinations = partition_degree(l, n)

  result = Vector{Vector{Int}}()
  exceeded = false
  for combination in combinations
    linear_combination = combination[1] * rays[1]
    for i in 2:l
      linear_combination = linear_combination + combination[i] * rays[i]
      if sum(linear_combination) > n
        exceeded = true
        break
      end
    end
    if exceeded == false
      push!(result, linear_combination)
    end
    exceeded = false
  end
  return result
end

function realizable_lengths_of_cylinder_diagram(cyl::CylinderDiagram, n::Int)
  equations = system_of_equations(cyl)
  r = compute_rays(equations, cyl.separatrix_count)
  return realizable_lengths(r, n)
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

# get origamis from cylinder diagram and realizable lengths
# from surface_dynamics by Vincent, Delecroix et al.
function origami_from_cylinder_coordinates(
  cyl_diagram::CylinderDiagram, lengths::Vector{Int}, heights::Vector{Int}
)
  # the total width of each cylinder is the sum of the lengths of the separatrices
  widths = [sum(lengths[i + 1] for i in bot) for bot in cyl_diagram.bot]
  areas = [heights[i] * widths[i] for i in 1:(cyl_diagram.cycles_count)]

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
  lx = perm([x + 1 for x in lx])

  ly = Int[]
  for i in 1:(cyl_diagram.cycles_count)
    append!(ly, (v[i] + widths[i]):(v[i + 1] - 1))
    append!(ly, zeros(Int64, widths[i]))
  end

  for (i, top_seps) in enumerate(cyl_diagram.top)
    top = []
    rev_top_seps = reverse(top_seps)
    for k in rev_top_seps
      append!(top, sep_bottom_pos[k + 1]:(sep_bottom_pos[k + 1] + lengths[k + 1] - 1))
    end
    ly[(v[i + 1] - widths[i] + 1):(v[i + 1])] = top
  end

  no_twist = normal_form(origami(lx, perm([x + 1 for x in ly])))
  results = Set{Origami}([no_twist])

  ly = [x + 1 for x in ly]

  for (i, o) in product_gray_code(widths)
    if o == 1
      insert!(ly, v[i + 1] - widths[i] + 1, popat!(ly, v[i + 1]))
    else
      insert!(ly, v[i + 1], popat!(ly, v[i + 1] - widths[i] + 1))
    end
    new_entry = normal_form(origami(lx, perm(ly)))
    push!(results, new_entry)
  end

  return results
end

function cylinder_diagrams_h11()::Vector{CylinderDiagram}
  c1 = CylinderDiagram([[0, 3, 1, 2]], [[0, 3, 1, 2]])
  c2 = CylinderDiagram([[0], [1, 2, 3]], [[1], [0, 2, 3]])
  c3 = CylinderDiagram([[0, 3], [1, 2]], [[1, 3], [0, 2]])
  c4 = CylinderDiagram([[0, 1], [2], [3]], [[2, 3], [1], [0]])
  return [c1, c2, c3, c4]
end

# m: max length, n: max height
function origamis_in_h11(m::Int, n::Int)
  cylinder_diagrams = cylinder_diagrams_h11()
  origamis = Set{Origami}()
  for cyl in cylinder_diagrams
    equations = system_of_equations(cyl)
    rays = compute_rays(equations, cyl.separatrix_count)
    lengths = realizable_lengths(rays, m)
    heights = all_combinations(n, cyl.cycles_count)
    for l in lengths
      for h in heights
        o = origami_from_cylinder_coordinates(cyl, l, h)
        origamis = union!(origamis, o)
      end
    end
  end
  return origamis
end

# max: the number that will be partitioned
# n: the number of parts max will be divided into
function partition_degree(n::Int, max::Int, current::Vector{Int}=Vector{Int}(),
  combinations::Vector{Vector{Int}}=Vector{Vector{Int}}(), sumSoFar::Int=0)
  if sumSoFar > max
    return combinations  # Return early if the current sum exceeds the max
  end

  if length(current) == n
    if sumSoFar <= max
      push!(combinations, copy(current))  # Add the valid combination to the list
    end
    return combinations
  end

  for i in 1:(max - sumSoFar)
    push!(current, i)  # Append the current number to the combination
    partition_degree(n, max, current, combinations, sumSoFar + i)  # Recur with the updated combination and sum
    pop!(current)  # Remove the last element to backtrack
  end

  return combinations  # Return the list of all valid combinations
end

# kinda brute force
function possible_lengths_and_heights(cyl_diagram::CylinderDiagram, degree::Int)
  potential_heights = partition_degree(cyl_diagram.cycles_count, degree)
  potential_lengths = realizable_lengths_of_cylinder_diagram(cyl_diagram, degree)
  result = Vector{Vector{Vector{Int}}}()
  cyls = cylinders(cyl_diagram)
  for h in potential_heights
    for l in potential_lengths
      square_count = 0
      for i in 1:(cyl_diagram.cycles_count)
        length_sum = sum(l[j + 1] for j in cyls[i][1])
        square_count += length_sum * h[i]
      end

      if square_count == degree
        push!(result, [l, h])
      end
    end
  end

  return result
end

function origamis_in_h11(degree::Int)::Set{Origami}
  diagrams::Vector{CylinderDiagram} = cylinder_diagrams_h11()
  origamis = Set{Origami}()
  for c in diagrams
    lengths_heights = possible_lengths_and_heights(c, degree)
    for entry in lengths_heights
      o = origami_from_cylinder_coordinates(c, entry[1], entry[2])
      origamis = union!(origamis, o)
    end
  end
  return origamis
end

# so far only genus 2 & 3 supported
# other strata may still be explored manually by using the read_cylinder_diagrams function
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

  file_path = joinpath(dirname(@__FILE__), "..", "cylinder_diagrams", file_name)
  diagrams::Vector{CylinderDiagram} = read_cylinder_diagrams(file_path)

  origamis = Set{Origami}()
  for c in diagrams
    lengths_heights = possible_lengths_and_heights(c, degree)
    for entry in lengths_heights
      o = origami_from_cylinder_coordinates(c, entry[1], entry[2])
      origamis = union!(origamis, o)
    end
  end
  return origamis
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
