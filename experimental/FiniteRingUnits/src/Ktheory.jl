function direct_product(R, projs::Vector, maps::Vector{<:RingMultMap}; simplify = true)
  _D, _fromD, _intoD = Oscar.biproduct(domain.(maps)...)
  if simplify
    D, StoD = snf(_D)
    Oscar.set_attribute!(D, :direct_product => nothing, :show => nothing)
    DtoS = inv(StoD)
    fromD = [StoD * d for d in _fromD]
    intoD = [d * DtoS for d in _intoD]
  else
    D = _D
    fromD = _fromD
    intoD = _intoD
  end
  # assemble the map :(
  return D, RingMultMap(R, D,
                    x -> sum(intoD[i](preimage(maps[i], projs[i](x))) for i in 1:length(maps)),
                    y -> sum(preimage(projs[i], image(maps[i], fromD[i](y))) for i in 1:length(maps))
                   )
end

function abelianization_of_unit_group(R::Union{FiniteRing, Oscar.Hecke.AbstractAssociativeAlgebra})
  rings, projs = decompose_into_indecomposable_rings(R)
  k1s = [_abelianization_of_unit_group(RR) for RR in rings]
  return direct_product(R, projs, k1s)
end

# non-decomposing version (only the map)
function _abelianization_of_unit_group(R::Union{FiniteRing, Oscar.Hecke.AbstractAssociativeAlgebra})
  u = unit_group_nonrecursive(R)
  A, GtoA = Oscar.maximal_abelian_quotient(Oscar.FinGenAbGroup, domain(u))
  return RingMultMap(R, A, x -> GtoA(u\x), y -> image(u, (preimage(GtoA, y))))
end

function k1_simple_ring(R::FiniteRing)
  @assert is_prime(characteristic(R))
  RtoA = Oscar.isomorphism(Oscar.MatAlgebra, R)
  A = codomain(RtoA)
  @assert Oscar.is_simple(A)
  # just need to understand if it is M_2(F_2) or not
  q = order(base_ring(A))
  if q != 2 || dim(A) != 4
    return _abelianization_of_unit_group(R)
  end
  Rab, RabtoR = abelianization_of_unit_group(R)
  Q, RabtoQ = quo(Rab, gens(Rab), false)
  S, StoQ = snf(Q)
  return RingMultMap(R, S, x -> StoQ\(RabtoQ(preimage(RabtoR, x))), y -> RabtoR(preimage(RabtoQ, StoQ(y))))
end

function k1_semisimple_pring(R::FiniteRing)
  rings, projs = decompose_semisimple_p_ring(R)
  k1s = [k1_simple_ring(R) for R in rings]
  D, f = direct_product(R, projs, k1s)
  return f
end

function k1_semisimple_ring(R::FiniteRing)
  rings, projs = decompose_into_p_rings(R)
  k1s = [k1_semisimple_pring(R) for R in rings]
  return direct_product(R, projs, k1s)
end

function k1(R::Union{FiniteRing, Oscar.Hecke.AbstractAssociativeAlgebra})
  rngs, projs = decompose_into_indecomposable_rings(R)
  k1s = [_k1_naive_nonrec(S) for S in rngs]
  return direct_product(R, projs, k1s; simplify = true)
end

function _k1_naive_nonrec(R::Union{FiniteRing, Oscar.Hecke.AbstractAssociativeAlgebra})
  M = Oscar.matrix_ring(R, 3)
  @vprintln :FiniteRings "Constructing matrix ring for $R"
  MR, MRtoM = finite_ring(M)
  @vprintln :FiniteRings "Done"
  @vprintln :FiniteRings "Compute abelianization of unit groups of matrix ring"
  MRtoMA = _abelianization_of_unit_group(MR)
  MA = domain(MRtoMA)
  @vprintln :FiniteRings "Compute abelianization of unit groups of ring"
  RAtoR = _abelianization_of_unit_group(R)
  RA = domain(RAtoR)
  imgs = elem_type(MA)[]
  for a in gens(RA)
    m = M(Oscar.diagonal_matrix(R, [RAtoR(a), one(R), one(R)]))
    b = MRtoMA\(MRtoM\(m))
    push!(imgs, b)
  end
  h = hom(RA, MA, imgs)
  K, KtoRA = kernel(h)
  #return quo(RA, KtoRA.(gens(K)))
  Q, RAtoQ = quo(RA, KtoRA.(gens(K)), false)
  return RingMultMap(R, Q, x -> RAtoQ(RAtoR\x), y -> RAtoR(preimage(RAtoQ, y)))
end

