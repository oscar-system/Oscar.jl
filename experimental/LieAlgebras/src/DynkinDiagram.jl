@doc raw"""
    show_dynkin_diagram(fam::Symbol, rk::Int) -> String

Return a string representation of the Dynkin diagram of the root system of
the given type.
"""
function show_dynkin_diagram(fam::Symbol, rk::Int)
  @req _root_system_type_supported_by_GAP(fam, rk) "Unknown Dynkin type"
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
      D = "1 >>> 2"
    else
      error("This root system doesn't exist.")
    end
  else
    error("This root system doesn't exist.")
  end
  print(D)
end

function _root_system_type_supported_by_GAP(S, n)
  S in [:A, :B, :C, :D, :E, :F, :G] || return false
  n >= 1 || return false
  S == :D && n < 4 && return false
  S == :E && !(n in [6, 7, 8]) && return false
  S == :F && n != 4 && return false
  S == :G && n != 2 && return false
  return true
end
