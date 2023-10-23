mutable struct RootSystem
  roots::Vector{Vector{Int}}
  simple_roots::Vector{Vector{Int}}
  positive_roots::Vector{Vector{Int}}
  root_system_type::Tuple{Symbol,Int}
  GAP_root_system::GAP.GapObj
  function RootSystem(S::Symbol, n::Int)
    # S is a symbol detailing the type of the indecomposable root system 
    # e.g. "A", "B", "C",... and n is an integer for the number of simple roots
    S1 = GAP.Obj(S)
    RS = GAP.Globals.RootSystem(S1, n)
    sR = Vector{Vector{Int}}(GAP.Globals.SimpleSystem(RS))
    Ro1 = Vector{Vector{Int}}(GAP.Globals.PositiveRoots(RS))
    Ro2 = Vector{Vector{Int}}(GAP.Globals.NegativeRoots(RS))
    Ro = vcat(Ro1, Ro2)
    t = (S, n)
    return new(Ro, sR, Ro1, t, RS)
  end
end

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
