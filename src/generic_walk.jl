function generic_walk(G::Oscar.IdealGens, start::MonomialOrdering, target::MonomialOrdering)
  @vprintln :groebner_walk "Results for generic_walk"
  @vprintln :groebner_walk "Facets crossed for: "

  Lm = leading_term.(G, ordering = start)
  MG = MarkedGroebnerBasis(gens(G), Lm)
  v = next_gamma(MG, ZZ.([0]), start, target)

  @v_do :groebner_walk steps = 0
  while !isempty(v)
    MG = generic_step(MG, v, target)
    v = next_gamma(MG, v, start, target)

    @v_do :groebner_walk steps += 1
    @vprintln :groebner_walk v
    @vprintln :groebner_walk 2 G
  end

  @vprint :groebner_walk "Cones crossed: "
  @vprintln :groebner_walk steps
  return gens(MG)
end

@doc raw"""
    generic_step(MG::MarkedGroebnerBasis, v::Vector{ZZRingElem}, ord::MonomialOrdering)

Given the marked Gröbner basis `MG` and a facet normal vector `v`, compute the next marked Gröbner basis.
"""
function generic_step(MG::MarkedGroebnerBasis, v::Vector{ZZRingElem}, ord::MonomialOrdering)
  facet_Generators = facet_initials(MG, v)
  H = groebner_basis(
    ideal(facet_Generators); ordering=ord, complete_reduction=true, algorithm=:buchberger
  )

  H = MarkedGroebnerBasis(gens(H), leading_term.(H, ordering = ord))
  H = lift_generic(MG, H)

  autoreduce!(H)

  return H
end


#= 

----- thesis example (chap3) 
R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

o_s = lex(R)

o_t= weight_ordering([1,3,0], deglex(R))
S = canonical_matrix(o_s)
T = canonical_matrix(o_t)

I = ideal([x^2 + y*z, x*y + z^2])
G = groebner_basis(I, ordering = o_s, complete_reduction = true)

markedGB_generic_walk(G, o_s, o_t)


=#


#Auxiliary functions for the generic Gröbner walk, ordered by subroutine
#----------------------------------


#------next_gamma (Goal: Get next facet normal along generic path)

dropfirst(V::AbstractVector) = Iterators.drop(V, 1)
@doc raw"""
    difference_lead_tail(MG::MarkedGroebnerBasis)

Computes $a - b$ for $a$ a leading exponent and $b$ in the tail of some $g\in MG$. 
"""
function difference_lead_tail(MG::MarkedGroebnerBasis)
  (G,Lm) = gens(MG), markings(MG) 
  lead_exp = leading_exponent_vector.(Lm)
  
  v = zip(lead_exp, exponent_vectors.(G)) .|> splat((l, t) -> Ref(l).-t)
  
  return [ZZ.(v)./ZZ(gcd(v)) for v in unique!(reduce(vcat, v)) if !iszero(v)]
end

@doc raw"""
    next_gamma(
      MG::MarkedGroebnerBasis, w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering
    )

Given a marked Gröbner basis `MG`, a weight vector `w` and monomial orderings `start` and `target`,
returns the "next" facet normal, i.e. the bounding vector $v$ fulfilling $w<v$ and being minimal with respect to the facet preorder $<$.
"""
function next_gamma(
    MG::MarkedGroebnerBasis, w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering
  )
    V = filter_by_ordering(start, target, difference_lead_tail(MG))
    if w != ZZ.([0]) 
        V = filter_lf(w, start, target, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V 
        if facet_less_than(canonical_matrix(start), canonical_matrix(target),v, minV)
            minV = v
        end
    end
    return minV 
end

@doc raw"""
    facet_initials(MG::MarkedGroebnerBasis, v::Vector{ZZRingElem})

Given a marked Gröbner basis `MG` and a facet normal `v`, computes the Gröbner basis of initial forms 
of `MG` by truncating to all bounding vectors parallel to `v`.
"""
function facet_initials(MG::MarkedGroebnerBasis, v::Vector{ZZRingElem})
  R = base_ring(MG)
  inwG = Vector{MPolyRingElem}()

  ctx = MPolyBuildCtx(R)
  for (i, (g, m)) in gens_and_markings(MG) |> enumerate
    c, a = leading_coefficient_and_exponent(m)
    push_term!(ctx, c, a)

    for (d, b) in coefficients_and_exponents(g)
      if is_parallel(ZZ.(a - b), v)
        push_term!(ctx, d, b)
      end
    end

    push!(inwG, finish(ctx))
  end

  return inwG
end

@doc raw"""
    lift_generic(MG::MarkedGroebnerBasis, H::MarkedGroebnerBasis)

Given a marked Gröbner basis `MG` generating an ideal $I$ and a reduced marked Gröbner basis `H` of initial forms,
lift H to a marked Gröbner basis of I (with unknown ordering) by subtracting initial forms according to Fukuda, 2007.
"""
function lift_generic(MG::MarkedGroebnerBasis, H::MarkedGroebnerBasis)
  for i in 1:length(H.gens)
    H.gens[i] = H.gens[i] - normal_form(H.gens[i], MG)
    end
  return H
  end
  


#--------------------------




#helper functions 


#returns 'true' if Mv <_{lex} 0 , 'false' otherwise 
# <_{lex} is the lexicographic ordering on Q^n
#with the notation from my thesis, this is equivalent to v <_M 0
function new_less_than_zero(M::ZZMatrix, v::Vector{ZZRingElem})
  if is_zero(v)
    return false
  else
    i = 1
  end
  while dot(M[i, :], v) == 0
    i += 1
  end
  return dot(M[i, :], v) < 0
end
  
  
@doc raw"""
    filter_by_ordering(start::MonomialOrdering, target::MonomialOrdering, V::Vector{Vector{ZZRingElem}})

Computes all elements $v\in V$ with $0 <_{\texttt{target}} v$ and $v <_{\texttt{start}} 0$
"""
function filter_by_ordering(start::MonomialOrdering, target::MonomialOrdering, V::Vector{Vector{ZZRingElem}})
  pred = v -> (
    new_less_than_zero(canonical_matrix(target), ZZ.(v)) && !new_less_than_zero(canonical_matrix(start), ZZ.(v))
  )
  return unique!(filter(pred, V))
end

@doc raw"""
    matrix_less_than(M::ZZMatrix, v::Vector{ZZRingElem}, w::Vector{ZZRingElem})

Returns true if $Mv < Mw$ lexicographically, false otherwise.
"""
# matrix_lexicographic_less_than(M::ZZMatrix, v::Vector{ZZRingElem}, w::Vector{ZZRingElem}) = all(i -> dot(M[i, :], v) < dot(M[i, :], w), 1:size(M,1))
function matrix_lexicographic_less_than(M::ZZMatrix, v::Vector{ZZRingElem}, w::Vector{ZZRingElem})
    i = 1
    while dot(M[i, :], v) == dot(M[i, :], w) && i != number_of_rows(M) 
        i += 1 
    end
    return dot(M[i, :], v) < dot(M[i, :], w)
end

#Comment: It may make more sense to have the monomial orderings as inputs.
@doc raw"""
    facet_less_than(S::ZZMatrix, T::ZZMatrix, u::Vector{ZZRingElem}, v::Vector{ZZRingElem})

Returns true if $u < v$with respect to the facet preorder $<$. (cf. "The generic Gröbner walk" (Fukuda et a;. 2007), pg. 10)
"""
function facet_less_than(S::ZZMatrix, T::ZZMatrix, u::Vector{ZZRingElem}, v::Vector{ZZRingElem})
    i = 1
    while dot(T[i,:], u)v == dot(T[i,:], v)u && i != number_of_rows(T)
        i += 1
    end
    return matrix_lexicographic_less_than(S, dot(T[i,:], u)v, dot(T[i,:], v)u) 
end

#returns all elements of V smaller than w w.r.t the facet preorder

@doc raw"""
    filter_lf(w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering, V::Vector{Vector{ZZRingElem}})

Returns all elements of `V` smaller than `w` with respect to the facet preorder.
"""
function filter_lf(w::Vector{ZZRingElem}, start::MonomialOrdering, target::MonomialOrdering, V::Vector{Vector{ZZRingElem}})
    skip_indices = facet_less_than.(
      Ref(canonical_matrix(start)),
      Ref(canonical_matrix(target)),
      Ref(w),
      V
    )
    
    return unique!(V[skip_indices])
end


#returns true if u is a non-zero integer multiple of v

@doc raw"""
    is_parallel(u::Vector{ZZRingElem}, v::Vector{ZZRingElem})

Determines whether $u$ and $v$ are non-zero integer multiples of each other.
"""
function is_parallel(u::Vector{ZZRingElem}, v::Vector{ZZRingElem})
  return !iszero(v) && !iszero(u) && u./gcd(u) == v./gcd(v)
end

#TODO: add tests 