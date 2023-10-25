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

###############################################################################
#
#   Basic manipulation
#
###############################################################################
@doc raw"""
    number_of_roots(R::RootSystem)

Return the numbers of roots in the root system `R`.
"""
number_of_roots(R::RootSystem) = size(R.roots)[1]

@doc raw"""
    number_of_roots(S::Symbol, n::Int)

Return the numbers of roots in the root system of type `S`
"""
number_of_roots(S::Symbol, n::Int) = number_of_roots(root_system(S, n))

@doc raw"""
    getindex(R::RootSystem, r::Int)

Return the `r`-th root of the root system `R`.
"""
getindex(R::RootSystem, r::Int) = getindex(R.roots, r)

@doc raw"""
    root_system_type(R::RootSystem)

Return the Dynkin type of the root system `R`.
"""
root_system_type(R::RootSystem) = R.root_system_type

root_system_type_string(R::RootSystem) =
  string(R.root_system_type[1]) * string(R.root_system_type[2])

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, R::RootSystem)
  print(io, "Root system of type $(root_system_type_string(R))")
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(R1::RootSystem, R2::RootSystem)
  return R1.root_system_type == R2.root_system_type && R1.roots == R2.roots
end

function Base.hash(R::RootSystem, h::UInt)
  b = 0x9d96557cb5f07773 % UInt
  h = hash(R.root_system_type, h)
  h = hash(R.roots, h)
  return xor(h, b)
end
###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    root_system(S::Symbol, n::Int) -> RootSystem
    
Return the root system of type `Sn` where `S` is a symbol consisting out of
a single letter `A`, `B`, `C`, `D`, `E`, `F`, `G`.
The allowed values for `n` depend on the choice of `S`.
"""
function root_system(S::Symbol, n::Int)
  @req _root_system_type_supported_by_GAP(S, n) "Unknown Dynkin type or not supported by GAP"
  return RootSystem(S, n)
end

###############################################################################
#
#   further functions
#
###############################################################################

function _root_system_type_supported_by_GAP(S, n)
  S in [:A, :B, :C, :D, :E, :F, :G] || return false
  n >= 1 || return false
  S == :D && n < 4 && return false
  S == :E && !(n in [6, 7, 8]) && return false
  S == :F && n != 4 && return false
  S == :G && n != 2 && return false
  return true
end

@doc raw"""
    cartan_matrix(S::Symbol, n::Int) -> Matrix{QQFieldElem}

Return the Cartan matrix of the root system of type `Sn`.
For the semantics of the arguments, refer to `root_system(S::Symbol, n::Int)` ref.
"""
function cartan_matrix(S::Symbol, n::Int)
  return cartan_matrix(root_system(S, n))
end

@doc raw"""
    cartan_matrix(R::RootSystem) -> Matrix{QQFieldElem}

Return the Cartan matrix of the root system `R`.
"""
function cartan_matrix(R::RootSystem)
  RS = R.GAP_root_system
  CG = GAP.Globals.CartanMatrix(RS)
  C = matrix(QQ, CG)
  return C
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

@doc raw"""
    dynkin_diagram(R::RootSystem)

Return the Dynkin diagram of the root system `R`
"""
function dynkin_diagram(R::RootSystem)
  return dynkin_diagram(R.root_system_type...)
end
