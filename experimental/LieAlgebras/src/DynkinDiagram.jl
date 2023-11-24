@doc raw"""
    show_dynkin_diagram(fam::Symbol, rk::Int) -> Nothing

Prints a string representation of the Dynkin diagram of the root system of
the given cartan type.
"""
function show_dynkin_diagram(fam::Symbol, rk::Int)
  @req is_cartan_type(fam, rk) "Invalid cartan type"
  D = ""

  if fam == :A
    for i in 1:(rk - 1)
      D = D * string(i) * " - "
    end
    D = D * string(rk)

  elseif fam == :B
    if rk == 1
      D = string(rk)
    else
      for i in 1:(rk - 2)
        D = D * string(i) * " - "
      end
      D = D * string(rk - 1) * " >=> " * string(rk)
    end

  elseif fam == :C
    if rk == 1
      D = string(rk)
    else
      for i in 1:(rk - 2)
        D = D * string(i) * " - "
      end
      D = D * string(rk - 1) * " <=< " * string(rk)
    end

  elseif fam == :D
    if rk >= 4
      for i in 1:(4 * rk - 10)
        D = D * " "
      end
      D = D * string(rk - 1) * "\n"
      for i in 1:(4 * rk - 11)
        D = D * " "
      end
      D = D * "/\n"
      for i in 1:(rk - 3)
        D = D * string(i) * " - "
      end
      D = D * string(rk - 2) * "\n"
      for i in 1:(4 * rk - 12)
        D = D * " "
      end
      D = D * " \\ \n"
      for i in 1:(4 * rk - 10)
        D = D * " "
      end
      D = D * string(rk)
    else
      error("This root system doesn't exist.")
    end

  elseif fam == :E
    if rk == 6
      D = "1 - 3 - 4 - 5 - 6\n        |\n        2"
    elseif rk == 7
      D = "1 - 3 - 4 - 5 - 6 - 7\n        |\n        2"
    elseif rk == 8
      D = "1 - 3 - 4 - 5 - 6 - 7 - 8\n        |\n        2"
    else
      error("This root system doesn't exist.")
    end

  elseif fam == :F
    if rk == 4
      D = "1 - 2 >=> 3 - 4"
    else
      error("This root system doesn't exist.")
    end
  elseif fam == :G
    if rk == 2
      D = "1 <<< 2"
    else
      error("This root system doesn't exist.")
    end
  else
    error("This root system doesn't exist.")
  end
  print(D)
end
