# Induced map of (homogeneous) Koszul complexes
#
# Let `S` be a graded ring and `a = (a₁,…,aₘ)` and `b = (b₁,…,bₙ)` two sequences 
# in `S` with `bⱼ = ∑ᵢcⱼᵢ ⋅ aᵢ` for some matrix `C = (cⱼᵢ) ∈ Sⁿˣᵐ`. 
# Then there is a natural map of Koszul complexes
#
#         b      b       b
#   ⋀⁰ Sⁿ ← ⋀¹Sⁿ ← ⋀² Sⁿ ← …
#     ↓id    ↓c     ↓c
#   ⋀⁰ Sᵐ ← ⋀¹Sᵐ ← ⋀² Sᵐ ← …
#         a      a       a
#
# This is to construct that morphism.
function (fac::InducedKoszulMorFactory)(self::AbsHyperComplexMorphism, I::Tuple)
  p = first(I)
  dom = fac.dom[p]
  cod = fac.cod[p]
  A = matrix(fac.A)
  n = ngens(fac.dom)
  m = ngens(fac.cod)
  img_gens = elem_type(cod)[]
  for (i, ind_dom) in enumerate(OrderedMultiIndexSet(p, n))
    ii = indices(ind_dom)::Vector{Int}
    img = zero(cod)
    for (j, ind_cod) in enumerate(OrderedMultiIndexSet(p, m))
      jj = indices(ind_cod)::Vector{Int}
      A_sub = A[ii, jj]
      c = det(A_sub)
      is_zero(c) && continue
      img += c*cod[j]
    end
    push!(img_gens, img)
  end
  return hom(dom, cod, img_gens)
end

function can_compute(fac::InducedKoszulMorFactory, self::AbsHyperComplexMorphism, i::Tuple)
  return can_compute_index(fac.dom, i) && can_compute_index(fac.cod, i)
end


# The return type for `coordinates` for ideals is not yet unified.
# This implementation is used to provide output as `SRow`
coordinates(::Type{OutputType}, f::T, I::Ideal{T}) where {OutputType<:SRow, T<:RingElem} = coordinates(f, I)

function coordinates(::Type{OutputType}, f::T, I::Ideal{T}) where {OutputType<:SRow, T<:MPolyRingElem}
  R = parent(f)
  c = coordinates(f, I)::MatrixElem
  sparse_row(R, [(i, c[1, i]) for i in 1:ngens(I) if !is_zero(c[1, i])])
end



underlying_morphism(phi::InducedKoszulMorphism) = phi.internal_morphism

