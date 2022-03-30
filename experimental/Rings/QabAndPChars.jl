module QabModule

using ..Oscar

import Hecke: math_html
import Oscar: IJuliaMime

###############################################################################
#
#   Partial character functions
#
###############################################################################

mutable struct PartialCharacter{T}
  #A has generators of the lattice in rows
  A::fmpz_mat
	#images of the generators are saved in b
  b::Vector{T}
	#Delta are the indices of the cellular variables of the associated ideal
	#(the partial character is a partial character on Z^Delta)
  D::Set{Int64}
  function PartialCharacter{T}() where T 
    return new{T}()
  end

  function PartialCharacter{T}(mat::fmpz_mat, vals::Vector{T}) where T
    z = new{T}()
    z.A = mat
    z.b = vals
    return z
  end
end

function partial_character(A::fmpz_mat, vals::Vector{T}, variables::Set{Int} = Set{Int}()) where T <: FieldElem
  @assert nrows(A) == length(vals)
  z = PartialCharacter{T}(A, vals)
  if !isempty(variables)
    z.D = variables
  end
  return z
end

function (Chi::PartialCharacter)(b::fmpz_mat)
  @assert nrows(b) == 1
  @assert Nemo.ncols(b) == Nemo.ncols(Chi.A)
  s = can_solve_with_solution(Chi.A, b, side = :left)
  @assert s[1]
  return evaluate(FacElem(Dict([(Chi.b[i], s[2][1, i]) for i = 1:length(Chi.b)])))
end

function (Chi::PartialCharacter)(b::Vector{fmpz})
  return Chi(matrix(FlintZZ, 1, length(b), b))
end

function have_same_domain(P::PartialCharacter, Q::PartialCharacter)
  return have_same_span(P.A, Q.A)
end

function have_same_span(A::fmpz_mat, B::fmpz_mat)
  @assert ncols(A) == ncols(B)
  return hnf(A) == hnf(B)
end



function Base.:(==)(P::PartialCharacter{T}, Q::PartialCharacter{T}) where T <: FieldElem
  if P === Q
    return true
  end
  if !have_same_domain(P, Q)
  	return false
  end
  #now test if the values taken on the generators of the lattices are equal
	for i = 1:nrows(P.A)
		TestVec = view(P.A, i:i, 1:Nemo.ncols(P.A))
		if P(TestVec) != Q(TestVec)
			return false
		end
	end
  return true
end

function saturations(L::PartialCharacter{QabElem{T}}) where T
	#computes all saturations of the partial character L
  res = PartialCharacter{QabElem{T}}[]

  #first handle case wher the domain of the partial character is the zero lattice
  #in this case return L
  if iszero(L.A)
		push!(res, L)
		return res
  end

  #now not trivial case
  H = hnf(transpose(L.A))
  H = view(H, 1:ncols(H), 1:ncols(H))
  i, d = pseudo_inv(H)  #iH = d I_n
  #so, saturation is i' * H // d
  S = divexact(transpose(i)*L.A, d)

	B = Vector{Vector{QabElem{T}}}()
  for k = 1:nrows(H)
    c = i[1, k]
    for j = 2:ncols(H)
      c = gcd(c, i[j, k])
      if isone(c)
        break
      end
    end
    mu = evaluate(FacElem(Dict(Tuple{QabElem{T}, fmpz}[(L.b[j], div(i[j, k], c)) for j = 1:ncols(H)])))
		mu1 = roots(mu, Int(div(d, c)))
    push!(B,  mu1)
  end
  it = Hecke.cartesian_product_iterator(UnitRange{Int}[1:length(x) for x in B])
  vT = Vector{Vector{QabElem{T}}}()
  for I in it
    push!(vT, [B[i][I[i]] for i = 1:length(B)])
  end
  
  for k = 1:length(vT)
		#check if PChar(S,vT[k],L.D) puts on the right value on the lattice generators of L
		Pnew = partial_character(S, vT[k], L.D)
		flag = true	#flag if value on lattice generators is right
		for i = 1:Nemo.nrows(L.A)
			if Pnew(sub(L.A, i:i ,1:Nemo.ncols(L.A))) != L.b[i]
				flag = false
				println("found wrong saturation (for information), we delete it")
			end
		end
		if flag
			push!(res, Pnew)
		end
  end
  return res
end
end
