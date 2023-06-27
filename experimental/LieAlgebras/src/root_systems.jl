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
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

size(R::RootSystem,dim::Any)=size(R.roots,dim)

size(R::RootSystem)=size(R.roots)

getindex(R::RootSystem, r::Any) = getindex(R.roots, r)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::RootSystem)
   print(io, "Root system of type ")
   show(io, R.root_system_type)
end

###############################################################################
#
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
end

