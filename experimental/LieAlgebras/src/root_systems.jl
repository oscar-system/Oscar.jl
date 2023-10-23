@doc raw"""
    dynkin_diagram(S::Symbol, n::Int)

Return the Dynkin diagram of the root system of type `Sn`.
For the semantics of the arguments, refer to [root_system(S::Symbol, n::Int)` ref.
"""
function dynkin_diagram(S::Symbol, n::Int)
  @req _root_system_type_supported_by_GAP(S, n) "Unknown Dynkin type or not supported by GAP"
  D = ""

  if S == :A
    for i in 1:(n - 1)
      D = D * string(i) * " - "
    end
    D = D * string(n)

  elseif S == :B
    if n == 1
      D = string(n)
    else
      for i in 1:(n - 2)
        D = D * string(i) * " - "
      end
      D = D * string(n - 1) * " >=> " * string(n)
    end

  elseif S == :C
    if n == 1
      D = string(n)
    else
      for i in 1:(n - 2)
        D = D * string(i) * " - "
      end
      D = D * string(n - 1) * " <=< " * string(n)
    end

  elseif S == :D
    if n >= 4
      for i in 1:(4 * n - 10)
        D = D * " "
      end
      D = D * string(n - 1) * "\n"
      for i in 1:(4 * n - 11)
        D = D * " "
      end
      D = D * "/\n"
      for i in 1:(n - 3)
        D = D * string(i) * " - "
      end
      D = D * string(n - 2) * "\n"
      for i in 1:(4 * n - 12)
        D = D * " "
      end
      D = D * " \\ \n"
      for i in 1:(4 * n - 10)
        D = D * " "
      end
      D = D * string(n)
    else
      error("This root system doesn't exist.")
    end

  elseif S == :E
    if n == 6
      D = "1 - 3 - 4 - 5 - 6\n        |\n        2"
    elseif n == 7
      D = "1 - 3 - 4 - 5 - 6 - 7\n        |\n        2"
    elseif n == 8
      D = "1 - 3 - 4 - 5 - 6 - 7 - 8\n        |\n        2"
    else
      error("This root system doesn't exist.")
    end

  elseif S == :F
    if n == 4
      D = "1 - 2 >=> 3 - 4"
    else
      error("This root system doesn't exist.")
    end
  elseif S == :G
    if n == 2
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
