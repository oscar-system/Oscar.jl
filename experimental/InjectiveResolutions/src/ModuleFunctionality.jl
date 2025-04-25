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

function _degree_fast(a::MonoidAlgebraElem)
  return _degree_fast(underlying_element(a))
end

# `twist` stuff is only copy-pasted because the original signature was too narrow.
# This should be streamlined eventually.
function twist(M::ModuleFP{T}, g::FinGenAbGroupElem) where {T <: MonoidAlgebraElem}
 error("Not implemented for the given type")
end

function twist(M::SubquoModule{T}, g::FinGenAbGroupElem) where {T<:MonoidAlgebraElem}
 R = base_ring(M)
 @req parent(g) == grading_group(R) "Group element not contained in grading group of base ring"
 F = ambient_free_module(M)
 FN = twist(F, g)
 shift = hom(F, FN, gens(FN); check=false)
 return SubquoModule(FN, shift.(ambient_representatives_generators(M)), shift.(relations(M)))
 # The original code below was giving strange errors from somewhere in AA.
 GN = free_module(R, ngens(M))
 HN = free_module(R, length(relations(M)))
 a = hom(GN, F, ambient_representatives_generators(M); check=false)
 b = hom(HN, F, relations(M); check=false)
 A = matrix(a)
 B = matrix(b)
 N = subquotient(FN, A, B)
 return N
end

function twist(F::FreeMod{T}, g::FinGenAbGroupElem) where {T<:MonoidAlgebraElem}
 R = base_ring(F)
 @req parent(g) == grading_group(R) "Group element not contained in grading group of base ring"
 W = [x-g for x in F.d]
 G = graded_free_module(R, rank(F))
 G.d = W
 return G
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

### Speedup for free_resolution
# Copied from `src/Modules/UngradedModules/FreeResolutions.jl` and modified. 
# This can eventually be streamlined, but we did not want to touch the original 
# module files at the moment.

underlying_ring_type(A::MonoidAlgebra) = underlying_ring_type(typeof(A))
underlying_ring_type(::Type{T}) where {CT, RT, T<:MonoidAlgebra{CT, RT}} = RT

function free_resolution(M::SubquoModule{T}; 
    length::Int=0,
    algorithm::Union{Symbol, Nothing} = nothing
  ) where {T <: MonoidAlgebraElem}
  if isnothing(algorithm)
    algorithm = underlying_ring_type(base_ring(M)) <: MPolyQuoRing ? :sres : :fres
  end

  coefficient_ring(base_ring(M)) isa AbstractAlgebra.Field ||
      error("Must be defined over a field.")
  
  if underlying_ring_type(base_ring(M)) <: MPolyQuoRing
    !iszero(length) || error("Specify a length up to which a free resolution should be computed")
  end


  cc_complete = false

  #= Start with presentation =#
  pm = algorithm == :mres ? _presentation_minimal(M, minimal_kernel=false) : presentation(M)
  maps = [pm.maps[j] for j in 2:3]

  br = base_ring(M)
  kernel_entry          = image(pm.maps[1])[1]

  if ngens(kernel_entry) == 0
    cc = Hecke.ComplexOfMorphisms(Oscar.ModuleFP, pushfirst!(maps, pm.maps[1]), check = false, seed = -2)
    cc.fill     = _extend_free_resolution
    cc.complete = true
    return FreeResolution(cc)
  end

  singular_free_module  = singular_module(ambient_free_module(kernel_entry))
  singular_kernel_entry = Singular.Module(base_ring(singular_free_module),
                              [singular_free_module(repres(g)) for g in gens(kernel_entry)]...)

  #= This is the single computational hard part of this function =#
  if algorithm == :fres
    gbpres = Singular.std(singular_kernel_entry)
    res = Singular.fres(gbpres, length, "complete")
  elseif algorithm == :lres
    error("LaScala's method is not yet available in Oscar.")
    gbpres = singular_kernel_entry # or as appropriate, taking into account base changes
  elseif algorithm == :mres
    gbpres = singular_kernel_entry
    res = Singular.mres(gbpres, length)
  elseif algorithm == :nres
    gbpres = singular_kernel_entry
    res = Singular.nres(gbpres, length)
  elseif algorithm == :sres
    gbpres = Singular.std(singular_kernel_entry)
    res = Singular.sres(gbpres, length)
  else
    error("Unsupported algorithm $algorithm")
  end

  slen = iszero(res[Singular.length(res)+1]) ? Singular.length(res) : Singular.length(res)+1
  if length == 0 || slen < length
    cc_complete = true
  end

  if length != 0
    slen =  slen > length ? length : slen
  end

  #= Add maps from free resolution computation, start with second entry
   = due to inclusion of presentation(M) at the beginning. =#
  j   = 1
  while j <= slen
    if is_graded(M)
      codom = domain(maps[1])
      rk    = Singular.ngens(res[j])
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
      #ff = graded_map(codom, SM.matrix)
      ff = graded_map(codom, gens(SM); check=false)
      dom = domain(ff)
      insert!(maps, 1, ff)
      j += 1
    else
      codom = domain(maps[1])
      rk    = Singular.ngens(res[j])
      dom   = free_module(br, rk)
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
      insert!(maps, 1, hom(dom, codom, gens(SM); check=false))
      j += 1
    end
  end
  if cc_complete == true
    # Finalize maps.
    if is_graded(domain(maps[1]))
      Z = graded_free_module(br, 0)
    else
      Z = FreeMod(br, 0)
    end
    insert!(maps, 1, hom(Z, domain(maps[1]), Vector{elem_type(domain(maps[1]))}(); check=false))
  end

  cc = Hecke.ComplexOfMorphisms(Oscar.ModuleFP, maps, check = false, seed = -2)
  cc.fill     = _extend_free_resolution
  cc.complete = cc_complete
  set_attribute!(cc, :show => free_show, :free_res => M)
  set_attribute!(cc, :algorithm, algorithm)

  return FreeResolution(cc)
end

