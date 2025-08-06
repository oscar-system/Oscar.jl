@doc raw"""
    show_dynkin_diagram(cartan_matrix::ZZMatrix) -> Nothing

Prints a string representation of the Dynkin diagram of the root system
with the given cartan matrix.
The labels of the nodes are the indices of the simple roots.

Currently, only cartan matrices of finite type are supported.
"""
function show_dynkin_diagram(cartan_matrix::ZZMatrix)
  type, label = cartan_type_with_ordering(cartan_matrix)
  return show_dynkin_diagram(type, label)
end

@doc raw"""
    show_dynkin_diagram(rs::RootSystem) -> Nothing

Prints a string representation of the Dynkin diagram of the given root system.
The labels of the nodes are the indices of the simple roots.

Currently, only root systems of finite type are supported.
"""
function show_dynkin_diagram(rs::RootSystem)
  type, label = root_system_type_with_ordering(rs)
  return show_dynkin_diagram(type, label)
end

@doc raw"""
    show_dynkin_diagram(type::Vector{Tuple{Symbol,Int}}) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type.
"""
function show_dynkin_diagram(type::Vector{Tuple{Symbol,Int}})
  return show_dynkin_diagram(type, 1:sum(t[2] for t in type; init=0))
end

@doc raw"""
    show_dynkin_diagram(type::Vector{Tuple{Symbol,Int}}, labels::AbstractVector{Int}) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type.
"""
function show_dynkin_diagram(type::Vector{Tuple{Symbol,Int}}, labels::AbstractVector{Int})
  @req length(labels) == sum(t[2] for t in type; init=0) "Invalid number of labels"
  offset = 0
  for (fam, rk) in type
    show_dynkin_diagram(fam, rk, labels[(offset + 1):(offset + rk)])
    offset += rk
    println()
    println()
  end
  return nothing
end

@doc raw"""
    show_dynkin_diagram(fam::Symbol, rk::Int) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type.
"""
function show_dynkin_diagram(fam::Symbol, rk::Int)
  return show_dynkin_diagram(fam, rk, 1:rk)
end

@doc raw"""
    show_dynkin_diagram(fam::Symbol, rk::Int, labels::Vector) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type, where the $i$-th node is labeled by `labels[i]`.
"""
function show_dynkin_diagram(fam::Symbol, rk::Int, labels::AbstractVector)
  @req is_cartan_type(fam, rk) "Invalid cartan type"
  @req length(labels) == rk "Invalid number of labels"
  D = ""
  if fam == :A
    D = join((labels[i] for i in 1:rk), " - ")
  elseif fam == :B
    D = join((labels[i] for i in 1:rk), " - ", " >=> ")
  elseif fam == :C
    D = join((labels[i] for i in 1:rk), " - ", " <=< ")
  elseif fam == :D
    main_path = join((labels[i] for i in 1:(rk - 2)), " - ")
    offset = length(main_path)
    # Documenter needs some character at the beginning of the first line (JuliaDocs/Documenter.jl#1159)
    D = '.' * repeat(' ', offset) * string(labels[rk - 1]) * '\n'
    D *= repeat(' ', offset) * "/\n"
    D *= main_path * '\n'
    D *= repeat(' ', offset) * "\\\n"
    D *= repeat(' ', offset + 1) * string(labels[rk])
  elseif fam == :E
    main_path = "$(labels[1]) - $(labels[3]) - "
    offset = length(main_path)
    main_path *= join((labels[i] for i in 4:rk), " - ")
    D = main_path * '\n'
    D *= repeat(' ', offset) * "|\n"
    D *= repeat(' ', offset) * string(labels[2])
  elseif fam == :F
    @assert rk == 4
    D = "$(labels[1]) - $(labels[2]) >=> $(labels[3]) - $(labels[4])"
  elseif fam == :G
    @assert rk == 2
    D = "$(labels[1]) <<< $(labels[2])"
  end
  isempty(D) && error("Unreachable")
  print(D)
  return nothing
end
