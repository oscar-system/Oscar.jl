@doc raw"""
    show_dynkin_diagram(fam::Symbol, rk::Int) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type.
"""
function show_dynkin_diagram(fam::Symbol, rk::Int)
  @req is_cartan_type(fam, rk) "Invalid cartan type"
  D = ""
  if fam == :A
    D = join((string(i) for i in 1:rk), " - ")
  elseif fam == :B
    D = join((string(i) for i in 1:rk), " - ", " >=> ")
  elseif fam == :C
    D = join((string(i) for i in 1:rk), " - ", " <=< ")
  elseif fam == :D
    main_path = join((string(i) for i in 1:(rk - 2)), " - ")
    offset = length(main_path)
    # Documenter needs some character at the beginning of the first line (JuliaDocs/Documenter.jl#1159)
    D = '.' * repeat(' ', offset) * string(rk - 1) * '\n'
    D *= repeat(' ', offset) * "/\n"
    D *= main_path * '\n'
    D *= repeat(' ', offset) * "\\\n"
    D *= repeat(' ', offset + 1) * string(rk) * '\n'
  elseif fam == :E
    if rk == 6
      D = """
          1 - 3 - 4 - 5 - 6
                  |
                  2
          """
    elseif rk == 7
      D = """
          1 - 3 - 4 - 5 - 6 - 7
                  |
                  2
          """
    elseif rk == 8
      D = """
          1 - 3 - 4 - 5 - 6 - 7 - 8
                  |
                  2
          """
    end
  elseif fam == :F
    @assert rk == 4
    D = "1 - 2 >=> 3 - 4"
  elseif fam == :G
    @assert rk == 2
    D = "1 <<< 2"
  end
  isempty(D) && error("Unreachable")
  print(D)
end
