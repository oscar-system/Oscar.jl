clear_denominators(A::MatrixType) where {T<:MPolyQuoLocRingElem, MatrixType<:MatrixElem{T}} = clear_denominators(map_entries(lift, A))

# use the lifts to do computations effectively over the base ring
# rather than the quotient ring.
function clear_denominators(v::FreeModElem{<:MPolyQuoLocRingElem})
  F = parent(v)
  L = base_ring(F)
  R = base_ring(L)
  d = lcm(lifted_denominator.(Vector(v)))
  u = elem_type(R)[]
  for a in Vector(v)
    push!(u, lifted_numerator(a)*div(d, lifted_denominator(a)))
  end
  Fb = base_ring_module(F)
  return sum([a*e for (a, e) in zip(u, gens(Fb))]), d
end

# constructs the matrix for the submodule I*F where 
# F is the free module of rank n over the base ring R 
# and I is the modulus.
function modulus_matrix(L::MPolyQuoLocRing, n::Int)
  I = modulus(underlying_quotient(L))
  m = ngens(I)
  R = base_ring(L)
  B = zero_matrix(R, 0, n)
  for j in 1:n
    A = zero_matrix(R, m, n)
    for i in 1:m
      A[i, j] = gen(I, i)
    end
    B = vcat(B, A)
  end
  return B
end

function syz(A::MatrixElem{<:MPolyQuoLocRingElem})
  B, D = clear_denominators(A)
  L = syz(vcat(B, modulus_matrix(base_ring(A), ncols(B))))
  return map_entries(base_ring(A), transpose(transpose(D) * transpose(L[:,1:nrows(D)])))
end

function ann(b::MatrixType, A::MatrixType) where {T<:MPolyQuoLocRingElem, MatrixType<:MatrixElem{T}}
  R = base_ring(A)
  R === base_ring(b) || error("matrices must be defined over the same ring")
  nrows(b) == 1 || error("only matrices with one row are allowed!")
  B = vcat(A, modulus_matrix(R, ncols(A)))
  m = nrows(B)
  n = ncols(B)
  Aext = vcat(b, B)
  L = syz(Aext)
  return ideal(R, vec(L[:, 1]))
end

function has_solution(A::MatrixType, b::MatrixType) where {T<:MPolyQuoLocRingElem, MatrixType<:MatrixElem{T}}
  S = base_ring(A)
  R = base_ring(S)
  S === base_ring(b) || error("matrices must be defined over the same ring")
  nrows(b) == 1 || error("only matrices with one row are allowed!")
  Aext = vcat(A, change_base_ring(S, modulus_matrix(S, ncols(A))))
  B, D = clear_denominators(Aext)
  c, u = clear_denominators(b)
  (success, y, v) = has_solution(B, c, inverted_set(S))
  success || return (false, zero_matrix(S, 1, ncols(b)))
  # We have B = D⋅Aext and c = u ⋅ b as matrices. 
  # Now [y z]⋅B = v⋅c ⇔ [y z]⋅D ⋅Aext = v ⋅ u ⋅ b ⇔ v⁻¹ ⋅ u⁻¹ ⋅ [y z] ⋅ D ⋅ Aext = b.
  # Take v⁻¹ ⋅ u⁻¹ ⋅ [y z] ⋅ D to be the solution x of x ⋅ Aext = b. 
  # Then for the first m components x' of x we have x' ⋅ A ≡ b mod I
  x = S(one(R), v*u[1,1])*map_entries(S, transpose(transpose(D) * transpose(y)))
  #x = S(one(R), v*u[1,1])*map_entries(S, y*D)
  xpart = zero_matrix(S, 1, nrows(A))
  for i in 1:nrows(A)
    xpart[1, i] = x[1, i]
  end
  return (success, xpart)
end

function pre_saturated_module(M::SubquoModule{T}) where {T<:MPolyQuoLocRingElem}
  has_attribute(M, :saturated_module) && return get_attribute(M, :saturated_module)::SubquoModule{elem_type(base_ring_type(T))}
  return get_attribute!(M, :pre_saturated_module) do
    S = base_ring(M)
    R = base_ring(S)
    (A, D) = clear_denominators(generator_matrix(M))
    relM = relations_matrix(M)
    (B, E) = clear_denominators(relM)
    mod_mat = modulus_matrix(S, ncols(B))
    B = vcat(B, mod_mat)
    #E = vcat(E, zero_matrix(R, nrows(mod_mat), ncols(E)))
    for i in 1:nrows(mod_mat)
      push!(E, sparse_row(base_ring(E)))
    end
    F = ambient_free_module(M)
    Fb = base_ring_module(F)
    Mb = SubquoModule(Fb, A, B)
    set_attribute!(M, :pre_saturation_data_gens, change_base_ring(S, D))
    set_attribute!(M, :pre_saturation_data_rels, change_base_ring(S, E))
    return Mb
  end::SubquoModule{elem_type(base_ring_type(T))}
end

# The kernel routine has to be overwritten since the base_ring_module of a
# free module does not have the modulus of the affine algebra under the 
# localization
function kernel(
    f::FreeModuleHom{DomType, CodType, Nothing}
  ) where {
    T<:MPolyQuoLocRingElem, 
    DomType<:FreeMod{T},
    CodType<:FreeMod{T}
  }
  S = base_ring(domain(f))
  A = representing_matrix(f)
  B, D = clear_denominators(A)
  Fb = base_ring_module(domain(f))
  Gb = base_ring_module(codomain(f))
  GGb, p = quo(Gb, elem_type(Gb)[g*e for g in gens(modulus(underlying_quotient(S))) for e in gens(Gb)])
  fb = hom(Fb, Gb, B)
  ffb = compose(fb, p)
  Kb, incb = kernel(ffb)
  Cb = representing_matrix(incb)
  C = change_base_ring(S, transpose(transpose(D) * transpose(Cb)))
  C = change_base_ring(S, transpose(transpose(D) * transpose(Cb)))
  #C = change_base_ring(S, Cb*D)
  K, inc = sub(domain(f), C)
  return K, inc
end


########################################################################
# lifting base_ring from MPolyQuoLocRing to MPolyLocRing               #
########################################################################

# MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}) where {BRT, BRET, RT, RET, MST}

#
### For a free module F = R^r over R = P/I, this returns a lifting map 
# to the module P^r/I*P^r. Note that this is an unnatural map since 
# the latter is an R-module only by accident.
@attr Any function _lifting_iso(F::FreeMod{T}) where {T<:MPolyQuoLocRingElem}
  M = _as_polyloc_module(F)
  function my_lift(v::FreeModElem{T}) where {T<:MPolyQuoLocRingElem}
    parent(v) === F || error("element does not have the right parent")
    w = elem_type(M)[lift(a)*M[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(M)
    return sum(w)
  end
  return my_lift
end

### For a free module F = R^r over R = P/I, this returns a lifting map 
# to the module P^r. Note that this is not a homomorphism of modules. 
@attr Any function _lifting_map(F::FreeMod{T}) where {T<:MPolyQuoLocRingElem}
  FP = _polyloc_module(F)
  function my_lift(v::FreeModElem{T}) where {T<:MPolyQuoLocRingElem}
    parent(v) === F || error("element does not have the right parent")
    w = [lift(a)*FP[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(FP)
    return sum(w)
  end
  return my_lift
end

### To a free module over R = P/I, return the free module over R 
# in the same number of generators
@attr FreeMod{<:MPolyLocRingElem{BRT, BRET, RT, RET, MST}} function _polyloc_module(F::FreeMod{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST}
  R = base_ring(F)
  P = localized_ring(R) # the polyloc ring
  r = rank(F)
  FP = FreeMod(P, r) 
  return FP
end

### Return the canonical projection FP -> F from the P-module FP to the 
# R-module F.
@attr FreeModuleHom{FreeMod{T}, FreeMod{MPolyQuoLocRingElem{T}}, MPolyQuoLocRing{T}} function _polyloc_module_restriction(F::FreeMod{MPolyQuoLocRingElem{T}}) where {T}
  R = base_ring(F)
  P = localized_ring(R)
  FP = _polyloc_module(F)
  return hom(FP, F, gens(F), R)
end

### Return the same module, but as a SubquoModule over the localized polynomial ring
@attr SubquoModule{<:MPolyLocRingElem{BRT, BRET, RT, RET, MST}} function _as_polyloc_module(F::FreeMod{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST}
  R = base_ring(F)
  P = localized_ring(R)
  I = modulus(R)
  FP = FreeMod(P, rank(F))
  IFP, inc = I*FP
  M, p = quo(FP, IFP)
  # # Manually set a groebner basis
  # # This is brittle! Please improve!
  # J = ideal(P, gens(groebner_basis(I))) # a groebner basis has already been computed.
  # JFP, _ = J*FP
  # G, _ = quo(FP, JFP)
  # gb = G.quo.gens
  # ord = default_ordering(M.quo)
  # gb.ordering = ord
  # gb.isGB = true
  # singular_generators(gb).isGB = true
  # M.quo.groebner_basis[ord] = gb
  return M
end

### Return an isomorphism with _as_polyloc_module(F)
@attr Any function _iso_with_polyloc_module(F::FreeMod{T}) where {T<:MPolyQuoLocRingElem}
  M = _as_polyloc_module(F)
  return hom(M, F, gens(F), base_ring(F))
end

### Return the preimage of M under the canonical projection P^r -> R^r 
# for R^r the ambient_free_module of M.
@attr Any function _polyloc_module(M::SubModuleOfFreeModule{T}) where {T<:MPolyQuoLocRingElem}
  F = ambient_free_module(M) 
  FP = _polyloc_module(F)
  v = elem_type(FP)[_lifting_map(F)(g) for g in gens(M)] 
  w = elem_type(FP)[f*e for e in gens(FP) for f in gens(modulus(base_ring(M)))]
  MP = SubModuleOfFreeModule(FP, vcat(v, w))
  return MP
end

@attr SubquoModule{<:MPolyLocRingElem{BRT, BRET, RT, RET, MST}} function _as_polyloc_module(M::SubquoModule{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST}
  F = ambient_free_module(M) 
  FP = _polyloc_module(F)
  v = [_lifting_map(F)(g) for g in ambient_representatives_generators(M)] 
  w = [f*e for e in gens(FP) for f in gens(modulus(base_ring(M)))]
  w_ext = vcat(w, elem_type(FP)[_lifting_map(F)(g) for g in relations(M)])
  MP = SubquoModule(FP, v, w_ext)
  return MP
end

@attr SubQuoHom{SubquoModule{<:MPolyLocRingElem{BRT, BRET, RT, RET, MST}}, SubquoModule{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}, MPolyQuoLocRing{BRT, BRET, RT, RET, MST}} function _iso_with_polyloc_module(F::SubquoModule{MPolyQuoLocRingElem{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST}
  M = _as_polyloc_module(F)
  return hom(M, F, gens(F), base_ring(F))
end

@attr Any function _lifting_iso(F::SubquoModule{T}) where {T<:MPolyQuoLocRingElem}
  M = _as_polyloc_module(F)
  function my_lift(v::SubquoModuleElem{T}) where {T<:MPolyQuoLocRingElem}
    parent(v) === F || error("element does not have the right parent")
    w = elem_type(M)[lift(a)*M[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(M)
    return sum(w)
  end
  return my_lift
end
