function (fac::StrandChainFactory)(c::AbsHyperComplex, i::Tuple)
  M = fac.orig[i]
  @assert is_graded(M) "module must be graded"
  R = base_ring(M)
  kk = coefficient_ring(R)
  return FreeMod(kk, length(all_exponents(fac.orig[i], fac.d)))
end

function can_compute(fac::StrandChainFactory, c::AbsHyperComplex, i::Tuple)
  return can_compute_index(fac.orig, i)
end


function (fac::StrandMorphismFactory)(c::AbsHyperComplex, p::Int, i::Tuple)
  I = collect(i)
  Next = I + (direction(c, p) == :chain ? -1 : 1)*[(k==p ? 1 : 0) for k in 1:dim(c)]
  next = Tuple(Next)
  orig_dom = fac.orig[i]
  orig_cod = fac.orig[next]
  dom = c[i]
  cod = c[next]

  orig_map = map(fac.orig, p, i)

  # Use a dictionary for fast mapping of the monomials to the 
  # generators of `cod`.
  cod_dict = Dict{Tuple{Vector{Int}, Int}, elem_type(cod)}(m=>cod[k] for (k, m) in enumerate(all_exponents(orig_cod, fac.d)))
  # Hashing of FreeModElem's can not be assumed to be non-trivial. Hence we use the exponents directly.
  img_gens_res = elem_type(cod)[]
  R = base_ring(orig_dom)
  vv = gens(R)
  for (e, i) in all_exponents(orig_dom, fac.d) # iterate through the generators of `dom`
    m = prod(x^k for (x, k) in zip(vv, e); init=one(R))*orig_dom[i]
    v = orig_map(m) # map the monomial
    # take preimage of the result using the previously built dictionary.
    # TODO: Iteration over the terms of v is VERY slow due to its suboptimal implementation.
    # We have to iterate manually. This saves us roughly 2/3 of the memory consumption and 
    # it also runs three times as fast. 
    w = zero(cod)
    for (i, b) in coordinates(v)
      #g = orig_cod[i]
      w += sum(c*cod_dict[(n, i)] for (c, n) in zip(AbstractAlgebra.coefficients(b), AbstractAlgebra.exponent_vectors(b)); init=zero(cod))
    end
    push!(img_gens_res, w)
  end
  return hom(dom, cod, img_gens_res)
end

function can_compute(fac::StrandMorphismFactory, c::AbsHyperComplex, p::Int, i::Tuple)
  return can_compute_map(fac.orig, p, i)
end

### User facing constructor
function strand(c::AbsHyperComplex{T}, d::Union{Int, FinGenAbGroupElem}) where {T<:ModuleFP}
  result = StrandComplex(c, d)
  inc = StrandInclusionMorphism(result)
  result.inclusion_map = inc
  return result, inc
end


# TODO: Code duplicated from `monomial_basis`. Clean this up!
function all_exponents(W::MPolyDecRing, d::FinGenAbGroupElem)
  D = W.D
  is_free(D) || error("Grading group must be free")
  h = hom(free_abelian_group(ngens(W)), W.d)
  fl, p = has_preimage_with_preimage(h, d)
  R = base_ring(W)
  B = Vector{Int}[]
  if fl
     k, im = kernel(h)
     #need the positive elements in there...
     #Ax = b, Cx >= 0
     C = identity_matrix(ZZ, ngens(W))
     A = reduce(vcat, [x.coeff for x = W.d])
     k = solve_mixed(transpose(A), transpose(d.coeff), C)
     B = Vector{Int}[k[ee, :] for ee in 1:nrows(k)]
  end
  return B
end

