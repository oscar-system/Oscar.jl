using Oscar 

#constructs a "markedGB" consisting of two lists 
#gens, a vector of generators of a Gröbner basis
#markings, a vector "Markings" of G, corresponding to leading terms w.r.t some monomial ordering
#TODO: maybe initialize s.t. G is reduced/monic? 

struct markedGB
    gens::Vector{<:MPolyRingElem}
    markings::Vector{<:MPolyRingElem}
    function markedGB(gens::Vector{<:MPolyRingElem},markings::Vector{<:MPolyRingElem})
        if length(gens)!= length(markings)
            throw(ArgumentError("Inputs are of different length"))
        else
            if 0 in gens
                throw(ArgumentError("Gröbner basis contains the zero polynomial"))
            end
            new(gens, markings) 
        end
    end
      
end

#Given a polynomial p and a monomial lm
#if lm divides a term of p, return (true, q). where p = q*lm + "terms not divisible by lm"  
function new_divides_walk(p::MPolyRingElem, lm::MPolyRingElem)
  div = false
  newpoly = 0*p 
  for term in terms(p)
    (b, c) = divides(term, lm)
    if b
      newpoly += c
      div = true
    end
  end
  return div, newpoly
end

@doc raw"""
    normal_form(p::MPolyRingElem, MG::markedGB)

Computes the normal form of `p` with respect to the marked reduced Gröbner basis `MG`.
"""
function normal_form(p::MPolyRingElem, MG::markedGB)
  #queue = Set(terms(p))
  nf = zero(parent(p))

  markedGB = zip(MG.gens, MG.markings)

  while !iszero(p)
    m = first(terms(p))

    div = false
    for (g, lm) in markedGB
      (div, q) = divides(m, lm)
      if div
        p = p - q*g
        break
      end
    end

    if !div
      nf += m
      p -= m
    end
  end

  return nf
end

# Calculates whether the monomial x^g divides x^f
# divides(f::Vector{Int}, g::Vector{Int}) = all(g .<= f)

#Given a markedGB MG, reduce it by replacing each g with its normal form w.r.t G\{g} 
#NB: only works if MG is inclusion minimal
#In this case upon reduction, the markings are preserved
#Question: Can I eliminate the 'if' condition? probably
function reductionalg(MG::markedGB)
    newgens = Vector{typeof(first(MG.gens))}()
    newlm = Vector{typeof(first(MG.gens))}()
    for i in 1:length(MG.gens)
        if normal_form(MG.gens[i], markedGB(MG.gens[1:end .!= i], MG.markings[1:end .!= i])) != 0
        push!(newgens, normal_form(MG.gens[i], markedGB(MG.gens[1:end .!= i], MG.markings[1:end .!= i])))
        push!(newlm, MG.markings[i])
        end 
    end
    return markedGB(newgens, newlm) 
end
    
#=

tests

R, (x,y) = polynomial_ring(QQ, ["x","y"])

G = [x^3 + y^4 + y^6 + x*y^7, y^2 + 4*x^5 + x^2*y^7]
LM  = [x*y^7, 4*x^5]
MG = markedGB(G, LM)
reductionalg(MG)
normal_form(x^5, MG)

R, (x,y,z,w) = polynomial_ring(QQ, ["x","y","z","w"])
G = [x-2*y-z-w, z+3*w]
lm = [x, z]
MG = markedGB(G, lm)
p = w
normal_form(p, MG)
normal_form(x, MG)


KK = GF(19)
R, (x,y) = polynomial_ring(KK, ["x", "y"])
G = [x^2 + y^2, y]
Lm = [x^2, y]
MG = markedGB(G, Lm)
MG = reductionalg(MG)
normal_form(x^3 - 1, MG)
=#

