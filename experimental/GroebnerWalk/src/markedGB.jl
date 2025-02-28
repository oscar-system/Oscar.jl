using Oscar 
import Oscar: gens, base_ring
import Base: length, getindex, setindex!

#constructs a "markedGB" consisting of two lists 
#gens, a vector of generators of a Gröbner basis
#markings, a vector "Markings" of G, corresponding to leading terms w.r.t some monomial ordering
#TODO: maybe initialize s.t. G is reduced/monic? 

struct MarkedGroebnerBasis{S}
    gens::Vector{S}
    markings::Vector{S}

    function MarkedGroebnerBasis(gens::Vector{T}, markings::Vector{T}) where {T<:MPolyRingElem}
        @req (length(gens) == length(markings)) "Inputs are of different length"
        @req !(0 in gens) "Gröbner basis contains the zero polynomial"

        new{T}(gens, markings) 
    end
end

base_ring(G::MarkedGroebnerBasis) = parent(first(G.gens))
length(G::MarkedGroebnerBasis) = ngens(G)
ngens(G::MarkedGroebnerBasis) = length(gens(G))

gens(G::MarkedGroebnerBasis) = G.gens
markings(G::MarkedGroebnerBasis) = G.markings
gens_and_markings(G::MarkedGroebnerBasis) = zip(gens(G), markings(G))

#MarkedPolynomial = Tuple{<:MPolyRingElem, <:MPolyRingElem}
# getindex(G::MarkedGroebnerBasis, i::Int) = (G.gens[i], G.markings[i])
# function setindex!(G::MarkedGroebnerBasis, (f,m)::MarkedPolynomial, i::Int) 
#   G.gens[i] = f
#   G.markings[i] = m
# end

@doc raw"""
    some_gens_with_markings(G::MarkedGroebnerBasis, I::AbstractVector)

Return a subset of generators of G indexed by I as pairs of polynomials and leading terms.
"""
function some_gens_with_markings(G::MarkedGroebnerBasis, I::AbstractVector)
  gens_view = @view gens(G)[I]  
  markings_view = @view markings(G)[I]

  return gens_view, markings_view
end

function _normal_form(p::MPolyRingElem, G::AbstractVector{<:MPolyRingElem}, mark::AbstractVector{<:MPolyRingElem})
  nf = MPolyBuildCtx(parent(p))
  MG = zip(G, mark)

  while !iszero(p)
    m = first(AbstractAlgebra.terms(p))

    div = false
    for (g, lm) in MG
      (div, q) = divides(m, lm)
      if div
        p = p - q*g
        break
      end
    end

    if !div
      push_term!(nf, leading_coefficient_and_exponent(m)...)
      p -= m
    end
  end

  return finish(nf)
end

@doc raw"""
    normal_form(p::MPolyRingElem, MG::MarkedGroebnerBasis)

Compute the normal form of `p` with respect to the marked reduced Gröbner basis `MG`.
"""
normal_form(p::MPolyRingElem, MG::MarkedGroebnerBasis) = _normal_form(p, gens(MG), markings(MG))

@doc raw"""
    autoreduce(MG::MarkedGroebnerBasis)

Given a marked Gröbner basis $MG$, reduce it by replacing each $g\in MG$ with its normal form $g^G$. 
This method requires `MG` to be a inclusion-minimal Gröbner basis. 
"""
function autoreduce(MG::MarkedGroebnerBasis)
    newgens = Vector{MPolyRingElem}()
    newlm = Vector{MPolyRingElem}()
    n = length(MG)

    for (i, (g, mark)) in enumerate(gens_and_markings(MG))
        rest = some_gens_with_markings(MG, 1:n .!= i)
        nf = _normal_form(g, rest...)

        #Question: Can we eliminate the 'if' condition?
        # if !iszero(nf)
          push!(newgens, nf)
          push!(newlm, mark)
        # end 
    end
    return MarkedGroebnerBasis(newgens, newlm) 
end

@doc raw"""
    autoreduce!(MG::MarkedGroebnerBasis)

Given a marked Gröbner basis $MG$, reduce it inplace by replacing each $g\in MG$ with its normal form $g^G$. 
"""
function autoreduce!(MG::MarkedGroebnerBasis)
  n = length(MG)
  for (i, g) in enumerate(gens(MG))
    rest = some_gens_with_markings(MG, 1:n .!= i)
    nf = _normal_form(g, rest...)

    MG.gens[i] = nf
  end
end

