default_ordering(A::MonoidAlgebra) = MonomialOrdering(A, Orderings.SymbOrdering(:degrevlex, collect(1:ngens(A))))

singular_poly_ring(A::MonoidAlgebra, ord::Singular.sordering; keep_ordering=nothing) = singular_poly_ring(A.algebra) # ignores the ordering for now
singular_poly_ring(A::MonoidAlgebra, ord::MonomialOrdering; keep_ordering=nothing) = singular_poly_ring(A.algebra) 
singular_poly_ring(A::MonoidAlgebra; keep_ordering=nothing) = singular_poly_ring(A.algebra) 

(S::Singular.PolyRing)(a::MonoidAlgebraElem) = S(underlying_element(a))

function singular_module(F::FreeMod{<:MonoidAlgebraElem}, ordering::ModuleOrdering)
  A = base_ring(F)
  Sx = singular_poly_ring(A, induced_ring_ordering(ordering))
  return Singular.FreeModule(Sx, dim(F))
end

function standard_basis(F::ModuleGens{T}, reduced::Bool=false) where {T <: MonoidAlgebraElem}
  singular_assure(F)
  if reduced
    @assert Singular.has_global_ordering(base_ring(F.SF))
  end
  if singular_generators(F).isGB && !reduced
    return F
  end
  return ModuleGens(F.F, Singular.std(singular_generators(F), complete_reduction=reduced))
end

function (F::FreeMod{T})(svec::Singular.svector) where {T<:MonoidAlgebraElem}
  Rx = base_ring(F)
  row = _build_sparse_row(Rx, svec)
  return FreeModElem(row, F)
end

function _build_sparse_row(
    A::MonoidAlgebra{<:FieldElem, <:MPolyRing}, s::Singular.svector; 
    cast::Ring=A
  )
  # TODO: Avoid double allocation by overwriting this properly!
  return map_entries(A, _build_sparse_row(A.algebra, s; cast))
end

function _build_sparse_row(
    A::MonoidAlgebra{<:FieldElem, <:MPolyQuoRing}, s::Singular.svector; 
    cast::Ring=A
  )
  # TODO: Avoid double allocation by overwriting this properly!
  return map_entries(A, _build_sparse_row(base_ring(A.algebra), s; cast))
end

# copied and modified from ModuleGens.jl
function normal_form(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MonoidAlgebraElem}
  @assert M.F === GB.F
  @assert GB.isGB # TODO When Singular.jl can handle reduce with non-GB remove this

  P = isdefined(GB, :quo_GB) ? union(GB, GB.quo_GB) : GB

  singular_assure(P)
  singular_assure(M)

  red = _reduce(M.S, P.S)
  res = ModuleGens(M.F, red)
  oscar_assure(res)
  return res
end

function lift_std(M::ModuleGens{T}) where {T <: MonoidAlgebraElem}
  singular_assure(M)
  R = base_ring(M)
  G,Trans_mat = Singular.lift_std(singular_generators(M)) # When Singular supports reduction add it also here
  mg = ModuleGens(M.F, G)
  mg.isGB = true
  mg.S.isGB = true
  mg.ordering = default_ordering(M.F)
  mat = map_entries(R, transpose(Trans_mat))
  set_attribute!(mg, :transformation_matrix => mat)
  return mg, mat
end

function lift_std(M::ModuleGens{T}, ordering::ModuleOrdering) where {T <: MonoidAlgebraElem}
  M = ModuleGens(M.O, M.F, ordering)
  mg, mat = lift_std(M)
  mg.ordering = ordering
  return mg, mat
end

function sparse_row(
    A::MonoidAlgebra{<:FieldElem, <:MPolyRing}, 
    svec::Singular.svector, rng::AbstractUnitRange
  )
  pre_res = sparse_row(A.algebra, svec, rng)
  return map_entries(A, pre_res)
end

function sparse_row(
    A::MonoidAlgebra{<:FieldElem, <:MPolyQuoRing}, 
    svec::Singular.svector, rng::AbstractUnitRange
  )
  pre_res = sparse_row(base_ring(A.algebra), svec, rng)
  return map_entries(A, pre_res)
end

function syzygy_module(F::ModuleGens{T}; sub = FreeMod(base_ring(F.F), length(oscar_generators(F)))) where {T <: MonoidAlgebraElem}
  singular_assure(F)
  # TODO Obtain the Gröbner basis and cache it
  s = Singular.syz(singular_generators(F))
  return SubquoModule(sub, s)
end

function kernel(
    h::FreeModuleHom{<:FreeMod{T}, <:FreeMod{T}, Nothing}
  ) where {S<:FieldElem, T <: MonoidAlgebraElem{S}}
  is_zero(h) && return sub(domain(h), gens(domain(h)))
  is_graded(h) && return _graded_kernel(h)
  return _simple_kernel(h)
end

### Additional adjustments to get the graded aspects to run
function annihilator(N::SubquoModule{T}) where T <: MonoidAlgebraElem
  R = base_ring(N)
  N_quo = isdefined(N, :quo) ? N.quo : SubModuleOfFreeModule(ambient_free_module(N), Vector{elem_type(ambient_free_module(N))}())
  A = N.sub
  SA = singular_generators(A.gens) 
  B = N_quo
  SB = singular_generators(B.gens)
  res = Singular.quotient(SB, SA)
  return ideal(R, [R(f) for f in res])
end

function singular_generators(J::MonoidAlgebraIdeal)
  return singular_generators(underlying_ideal(J).gens)
end

function degree(v::FreeModElem{T}; check::Bool=true) where {T <: MonoidAlgebraElem}
  is_zero(coordinates(v)) && return zero(grading_group(parent(v)))
  @check is_homogeneous(v) "element is not homogeneous"
  i, c = first(coordinates(v))
  return parent(v).d[i] + degree(c)
end

### Stuff to get `MonoidAlgebraIdeal`s to work with the modules.
function _saturation(U::SubModuleOfFreeModule, J::MonoidAlgebraIdeal; iteration::Bool=false)
  F = ambient_free_module(U)
  SgU = singular_generators(U.gens)
  SgJ = singular_generators(J)
  SQ, _ = Singular.saturation(SgU, SgJ)
  MG = ModuleGens(F, SQ)
  return SubModuleOfFreeModule(F, MG)
end

