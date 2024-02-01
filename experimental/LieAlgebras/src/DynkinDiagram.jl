@doc raw"""
    show_dynkin_diagram(fam::Symbol, rk::Int) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type.
"""
function show_dynkin_diagram(fam::Symbol, rk::Int)
  show_dynkin_diagram(fam, rk, 1:rk)
end

@doc raw"""
    show_dynkin_diagram(fam::Symbol, rk::Int, labels::Vector{Int}) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type, where the $i$-th node is labeled by `labels[i]`.
"""
function show_dynkin_diagram(fam::Symbol, rk::Int, labels::AbstractVector{Int})
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
end
