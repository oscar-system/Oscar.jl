<<<<<<< HEAD
mutable struct RootSystem
  roots::Vector{Vector{Int}}
  simple_roots::Vector{Vector{Int}}
  positive_roots::Vector{Vector{Int}}
<<<<<<< HEAD
<<<<<<< HEAD
  root_system_type::Tuple{Symbol, Int64}
  GAP_root_system::GAP.GapObj
=======
  root_system_type::Tuple{Symbol, Int64}

>>>>>>> changed to rebase to master
  function RootSystem(S::Symbol, n::Int64)
    # S is a symbol detailing the type of the indecomposable root system 
    # e.g. "A", "B", "C",... and n is an integer for the number of simple roots
    S1 = GAP.Obj(S)
    RS = GAP.Globals.RootSystem(S1, n)
    sR = Vector{Vector{Int}}(GAP.Globals.SimpleSystem(RS))
    Ro1 = Vector{Vector{Int}}(GAP.Globals.PositiveRoots(RS))
    Ro2 = Vector{Vector{Int}}(GAP.Globals.NegativeRoots(RS))
    Ro = vcat(Ro1, Ro2)
    t = (S, n)
<<<<<<< HEAD
    return new(Ro, sR, Ro1, t, RS)
=======
  root_system_type::Tuple{String, Int64}
  @doc raw"""
    RootSystem(S::String) -> RootSystem
  Return the root system of type `S` where `S` is a string consisting out of
  a letter `A`, `B`, `C`, `D`, `E`, `F`, `G` followed by an integer.
  For the exceptional root system, the integers are fixed, e.g. `G2`, `F4`, `E6`, ...
  """
  function RootSystem(S::String)
    # S is a string detailing the type of the indecomposable root system made out of a letter, e.g. "A", "B", "C",... and an integer for the number of simple roots
    S = parse_root_string(S)
    S1 = GAP.Obj(S[1])
    RS = GAP.Globals.RootSystem(S1, S[2])
    Ro_1 = GAP.Globals.PositiveRoots(RS)
    Ro_2 = GAP.Globals.NegativeRoots(RS)
    sR = GAP.Globals.SimpleSystem(RS)
    sR = [[sR[i][j] for j=1:length(sR[i])] for i=1:length(sR)]
    Ro1 = [[Ro_1[i][j] for j=1:length(Ro_1[i])] for i=1:length(Ro_1)]
    Ro2 = [[Ro_2[i][j] for j=1:length(Ro_2[i])] for i=1:length(Ro_2)]
    Ro = reduce(vcat, (Ro1, Ro2))
    return new(Ro,sR,Ro1,S)
>>>>>>> updated LieAlgebras.jl, root_systems.jl, simple_lie_algebra.jlm root_systems-test.jl, and runtests.jl
=======
    return new(Ro, sR, Ro1, t)
>>>>>>> changed to rebase to master
  end
=======
# LieAlg.jl : Implements root systems (by copying most structures from GAP)
#TODO: implement Weyl groups. Is it possible to give groups by generatrs and relations in Oscar?
using Oscar

GAP.Packages.load("sla")
 
import Base: show, size, getindex

mutable struct RootSystem <: AbstractVector{fmpq}
   roots::Vector
   simple_roots::Vector
   positive_roots::Vector
   root_system_type::String

   function RootSystem(S::String)
        l=length(S)
	n=parse(Int64,S[2:l])
	S1=GAP.Globals.String(S[1:1])[2:2]

	RS=GAP.Globals.RootSystem(S1,n)
	Ro_1=GAP.Globals.PositiveRoots(RS)
	Ro_2=GAP.Globals.NegativeRoots(RS)
	sR=GAP.Globals.SimpleSystem(RS)

        sR=[[sR[i][j] for j=1:length(sR[i])] for i=1:length(sR)]
	Ro1=[[Ro_1[i][j] for j=1:length(Ro_1[i])] for i=1:length(Ro_1)]
	Ro2=[[Ro_2[i][j] for j=1:length(Ro_2[i])] for i=1:length(Ro_2)]
	Ro=reduce(vcat, (Ro1, Ro2))
	new(Ro,sR,Ro1,S)
   end
>>>>>>> Add files via upload
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
@doc raw"""
    number_of_roots(R::RootSystem)

Return the numbers of roots in the root system `R`.
"""
number_of_roots(R::RootSystem) = size(R.roots)[1]

@doc raw"""
    number_of_roots(S::Symbol, n::Int64)
=======
 @doc raw"""
   number_of_roots(R::RootSystem)
=======
@doc raw"""
  number_of_roots(R::RootSystem)
>>>>>>> changed to rebase to master

Return the numbers of roots in the root system `R`
"""
number_of_roots(R::RootSystem) = size(R.roots)[1]

@doc raw"""
<<<<<<< HEAD
<<<<<<< HEAD
number_of_roots(S::String)
>>>>>>> updated root_systems.jl
=======
   number_of_roots(S::String)
>>>>>>> updated LieAlgebras.jl, root_systems.jl, simple_lie_algebra.jlm root_systems-test.jl, and runtests.jl
=======
  number_of_roots(S::String)
>>>>>>> changed to rebase to master

Return the numbers of roots in the root system of type `S`
"""
number_of_roots(S::Symbol, n::Int64) = number_of_roots(RootSystem(S, n))

@doc raw"""
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    getindex(R::RootSystem, r::Int)

Return the `r`-th root of the root system `R`.
=======
getindex(R::RootSystem, r::Int)
=======
   getindex(R::RootSystem, r::Int)
>>>>>>> updated LieAlgebras.jl, root_systems.jl, simple_lie_algebra.jlm root_systems-test.jl, and runtests.jl

Return the `r`=th root in the vector of roots in the root system `R`.
>>>>>>> updated root_systems.jl
=======
  getindex(R::RootSystem, r::Int)

Return the `r`-th root of the root system `R`.
>>>>>>> changed to rebase to master
"""
getindex(R::RootSystem, r::Int) = getindex(R.roots, r)

@doc raw"""
<<<<<<< HEAD
    root_system_type(R::RootSystem)
=======
  root_system_type(R::RootSystem)
>>>>>>> changed to rebase to master

Return the Dynkin type of the root system `R`.
"""
root_system_type(R::RootSystem) = R.root_system_type
=======

size(R::RootSystem,dim::Any)=size(R.roots,dim)

size(R::RootSystem)=size(R.roots)

getindex(R::RootSystem, r::Any) = getindex(R.roots, r)
>>>>>>> Add files via upload

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::RootSystem)
<<<<<<< HEAD
  print(io, "Root system of type ")
  show(io, string(R.root_system_type[1]) * string(R.root_system_type[2]))
<<<<<<< HEAD
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(R1::RootSystem, R2::RootSystem)
  return R1.roots == R2.roots &&
         R1.root_system_type != R2.root_system_type
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
    root_system(S::Symbol, n::Int64) -> root_system
    
Return the root system of type `S` where `S` is a string consisting out of
a letter `A`, `B`, `C`, `D`, `E`, `F`, `G` followed by an integer.
For the exceptional root system, the integers are fixed, e.g. `G2`, `F4`, `E6`, ...
"""
function root_system(S::Symbol, n::Int64)
  @req S in [:A, :B, :C, :D, :E, :F, :G] "Unknown Dynkin type, S needs to be a symbol in [:A, :B, :C, :D, :E, :F, :G]"
  @req n > 0 "n needs to be a positive integer"
  # S is a symbol detailing the type of the indecomposable root system 
  # e.g. "A", "B", "C",... and n is an integer for the number of simple roots
  return RootSystem(S, n)
=======
>>>>>>> changed to rebase to master
=======
   print(io, "Root system of type ")
   show(io, R.root_system_type)
>>>>>>> Add files via upload
end

###############################################################################
#
<<<<<<< HEAD
#   Comparison
#
###############################################################################

function Base.:(==)(R1::RootSystem, R2::RootSystem)
  return R1.roots == R2.roots &&
         R1.root_system_type != R2.root_system_type
end

function Base.hash(R::RootSystem, h::UInt)
  b = 0x9d96557cb5f07773 % UInt
  h = hash(R.positive_roots, h)
  h = hash(R.root_system_type, h)
  h = hash(R.roots, h)
  h = hash(R.simple_roots, h)
  return xor(h, b)
end
###############################################################################
#
#   Constructor
#
###############################################################################
@doc raw"""
  root_system(S::String) -> root_system
Return the root system of type `S` where `S` is a string consisting out of
a letter `A`, `B`, `C`, `D`, `E`, `F`, `G` followed by an integer.
For the exceptional root system, the integers are fixed, e.g. `G2`, `F4`, `E6`, ...
"""
function root_system(S::Symbol, n::Int64)
  @req S in [:A, :B, :C, :D, :E, :F, :G] "Unknown Dynkin type"
  # S is a symbol detailing the type of the indecomposable root system 
  # e.g. "A", "B", "C",... and n is an integer for the number of simple roots
  return RootSystem(S,n)
end
###############################################################################
#
#   further functions
#
###############################################################################
<<<<<<< HEAD
@doc raw"""
    cartan_matrix(S::Symbol, n::Int64) -> Matrix{QQFieldElem}

Return the Cartan matrix of the type of root system specified by `S`
  """
function cartan_matrix(S::Symbol, n::Int64)
  return cartan_matrix(root_system(S,n))
end

@doc raw"""
    cartan_matrix(R::RootSystem) -> Matrix{QQFieldElem}

Return the Cartan matrix of the type root system `R`
"""
function cartan_matrix(R::RootSystem)
  RS = R.GAP_root_system
  CG = GAP.Globals.CartanMatrix(RS)
  C = matrix(QQ, CG)
=======
 @doc raw"""
   cartan_matrix(S::String) -> Matrix{QQFieldElem}

 Return the Cartan matrix of the type of root system specified by `S`
  """
function cartan_matrix(S::Symbol, n::Int64)
  @req S in [:A, :B, :C, :D, :E, :F, :G] "Unknown Dynkin type"
  S1 = GAP.Obj(S)
  RS = GAP.Globals.RootSystem(S1, n)
  CG = GAP.Globals.CartanMatrix(RS)
  C = matrix(QQ, CG)
  return C
end

 @doc raw"""
   cartan_matrix(R::RootSystem) -> Matrix{QQFieldElem}

 Return the Cartan matrix of the type root system `R`
  """
function cartan_matrix(R::RootSystem)
  S = R.root_system_type
<<<<<<< HEAD
  S2 = S[1] * string(S[2])
  C = cartan_matrix(S2)
>>>>>>> updated LieAlgebras.jl, root_systems.jl, simple_lie_algebra.jlm root_systems-test.jl, and runtests.jl
=======
  S1 = GAP.Obj(S[1])
  RS = GAP.Globals.RootSystem(S1, S[2])
  CG = GAP.Globals.CartanMatrix(RS)
  C = matrix(QQ, CG)
>>>>>>> changed to rebase to master
  return C
end
 
@doc raw"""
<<<<<<< HEAD
    dynkin_diagram(S::String)

Return the Dynkin diagram of the type of root system specified by `S`
"""
function dynkin_diagram(S::Symbol, n::Int64)
  @req S in [:A, :B, :C, :D, :E, :F, :G] "Unknown Dynkin type"
  @req n >= 1 "We need a positive number of roots"
  D = ""
   
  if S == :A
=======
   dynkin_diagram(S::String)

 Return the Dynkin diagram of the type of root system specified by `S`
  """
<<<<<<< HEAD
function dynkin_diagram(S::String)
  S = parse_root_string(S)
  S1 = S[1]
	n = S[2]
	D = ""
	if S1 == "A"
>>>>>>> updated LieAlgebras.jl, root_systems.jl, simple_lie_algebra.jlm root_systems-test.jl, and runtests.jl
    for i = 1:(n-1)
      D = D * string(i) * " - "
    end
    D = D * string(n)
    
  elseif S == :B
    if n == 1
      D = string(n)
    else
      for i = 1:(n-2)
        D = D * string(i) * " - "
      end
      D = D * string(n-1) * " >=> " * string(n)
    end

  elseif S == :C
    if n == 1
      D = string(n)
    else
      for i = 1:(n-2)
        D = D * string(i) * " - "
      end
      D = D * string(n-1) * " <=< " * string(n)
    end
=======
function dynkin_diagram(S::Symbol, n::Int64)
  @req S in [:A, :B, :C, :D, :E, :F, :G] "Unknown Dynkin type"
   D = ""
  
  if S == :A
  
    for i = 1:(n-1)
      D = D * string(i) * " - "
    end
    D = D * string(n)
    
  elseif S == :B
  
    for i = 1:(n-2)
      D = D * string(i) * " - "
    end
    D = D * string(n-1) * " >=> " * string(n)
    
  elseif S == :C
  
    for i = 1:(n-2)
      D = D * string(i) * " - "
    end
    D = D * string(n-1) * " <=< " * string(n)
>>>>>>> changed to rebase to master
    
  elseif S == :D
    if n >= 4
      for i = 1:4*n-10
        D = D * " "
      end
      D = D * string(n-1) * "\n"
      for i = 1:4*n-11
        D = D * " "
      end
      D = D * "/\n"
      for i = 1:n-3
        D = D * string(i) * " - "
      end
      D = D * string(n-2) * "\n"
      for i = 1:4*n-12
        D = D * " "
      end
      D = D * " \\ \n"
      for i = 1:4*n-10
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
<<<<<<< HEAD
    dynkin_diagram(R::RootSystem)

Return the Dynkin diagram of the root system `R`
"""
=======
   dynkin_diagram(R::RootSystem)

 Return the Dynkin diagram of the root system `R`
  """
>>>>>>> updated LieAlgebras.jl, root_systems.jl, simple_lie_algebra.jlm root_systems-test.jl, and runtests.jl
function dynkin_diagram(R::RootSystem)
<<<<<<< HEAD
  S = R.root_system_type
  return dynkin_diagram(S[1],S[2])
=======
	S = R.root_system_type
	return dynkin_diagram(S[1],S[2])
>>>>>>> changed to rebase to master
=======
#   further functions
#
###############################################################################

function CartanMatrix(S::String)
		S1=S[1:1]
		l=length(S)
		n=parse(Int64,S[2:l])
		S1=GAP.Globals.String(S[1:1])[2:2]
		m=GAP.Globals.Int(n)

		RS=GAP.Globals.RootSystem(S1,n)
		CG=GAP.Globals.CartanMatrix(RS)
		C=Matrix{fmpq}(CG)
	return C
end

function CartanMatrix(R::RootSystem)
	S=R.root_system_type
	C=CartanMatrix(S)
	return C
end

function DynkinDiagram(S::String)
	S1=S[1:1]
	l=length(S)
	n=parse(Int64,S[2:l])
	D=""
	if S1 == "A"
		for i=1:(n-1)
			D=D* string(i) * " - "
		end
		D=D* string(n)
	elseif S1=="B"
		for i=1:(n-2)
			D=D* string(i) * " - "
		end
		D=D* string(n-1) * " >=> "*string(n)
	elseif S1=="C"
		for i=1:(n-2)
			D=D* string(i) * " - "
		end
		D=D* string(n-1) * " <=< "*string(n)
	elseif S1=="D"
		for i=1:4*n-10
			D=D*" "
		end
		D=D* string(n-1)*"\n"
		for i=1:4*n-11
			D=D*" "
		end
		D=D*"/\n"
		for i=1:n-3
			D=D* string(i) * " - "
		end
		D=D*string(n-2) *"\n"
		for i=1:4*n-12
			D=D*" "
		end
		D=D*" \\ \n"
		for i=1:4*n-10
			D=D*" "
		end
		D=D* string(n)
	elseif S1=="E"
		if n==6
			D="1 - 3 - 4 - 5 - 6\n        |\n        2"
		elseif n==7
			D="1 - 3 - 4 - 5 - 6 - 7\n        |\n        2"
		elseif n==8
			D="1 - 3 - 4 - 5 - 6 - 7 - 8\n        |\n        2"
		end
	elseif S1=="F"
		D="1 - 2 >=> 3 - 4"#this is changed from GAP
	elseif S1=="G"
		D="1 >>> 2"
	end
	print(D)
end

function DynkinDiagram(R::RootSystem)
	S=R.root_system_type
	DynkinDiagram(S)
>>>>>>> Add files via upload
end

