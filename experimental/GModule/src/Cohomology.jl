module GrpCoh

using Oscar
import Oscar: action
import Oscar: induce
import Oscar: word
import Oscar: GAPWrap, pc_group, fp_group, direct_product, direct_sum, GAPGroup
import AbstractAlgebra: Group, Module, FPModule
import Base: parent

import Oscar: pretty, Lowercase, @show_name, @show_special

function __init__()
  Hecke.add_verbosity_scope(:GroupCohomology)
  Hecke.add_assertion_scope(:GroupCohomology)

  Hecke.add_verbosity_scope(:GaloisCohomology)
  Hecke.add_assertion_scope(:GaloisCohomology)
end

######################################################################
#
# to allow additive notation for multiplicative objects,
# eg
# M = MultGrp(K::AbsSimpleNumField) will result in
# M(a) + M(b) = M(ab)
#
# since the Gmodule stuff is strictly additive.
#
struct MultGrp{T} <: Oscar.Hecke.GrpAb
  data::Any #should be like parent_type{T}
  elem_rep::T

  function MultGrp(L::S) where {S}
    return MultGrp(L, one(L))
  end
  function MultGrp(L, a::typ) where typ
    return new{typ}(L, a)
  end
end

function Base.show(io::IO, M::MultGrp)
  print(io, "multiplicative group of $(M.data)")
end

struct MultGrpElem{T} <: Oscar.Hecke.GrpAbElem
  data::T
  parent::MultGrp{T}
end

function Base.show(io::IO, m::MultGrpElem)
  print(io, "$(m.data)")
end

(M::MultGrp{T})(a::T) where {T}  = MultGrpElem{T}(a, M)

Oscar.parent(a::MultGrpElem) = a.parent
Oscar.elem_type(::Type{MultGrp{T}}) where T = MultGrpElem{T}
Oscar.parent_type(::Type{MultGrpElem{T}}) where T = MultGrp{T}
Oscar.zero(a::MultGrpElem) = parent(a)(one(a.data))
Oscar.zero(a::MultGrp) = a(one(a.data))

import Base: ==, +, -, *

*(a::Integer, b::MultGrpElem{T}) where T = MultGrpElem{T}(b.data^a, parent(b))
*(a::ZZRingElem, b::MultGrpElem{T}) where T = MultGrpElem{T}(b.data^a, parent(b))
+(a::MultGrpElem{T}, b::MultGrpElem{T}) where T = MultGrpElem{T}(a.data*b.data, parent(a))
-(a::MultGrpElem{T}, b::MultGrpElem{T}) where T = MultGrpElem{T}(a.data//b.data, parent(a))
-(a::MultGrpElem{T}) where T = MultGrpElem{T}(inv(a.data), parent(a))
==(a::MultGrpElem{T}, b::MultGrpElem{T}) where T = a.data == b.data
Base.hash(a::MultGrpElem, u::UInt = UInt(1235)) = hash(a.data. u)


##############################################################
#
# basic G-Modules:
#
# a fin. gen group G acting on some module M
# the action is given via maps (should be automorphisms of M)
#
# very little assumptions in general.
#
@attributes mutable struct GModule{gT,mT}
  G::gT
  M::mT
  ac::Vector{Map} # automorphisms of M, one for each generator of G

  function GModule(M, G::T, ac::Vector{<:Map}) where {T}
    r = new{T,typeof(M)}()
    r.G = G
    r.ac = ac
    r.M = M
    @assert all(x -> domain(x) === codomain(x) === r.M, ac)
    if isa(G, Group)
      @assert length(ac) == ngens(G)
      @hassert :GroupCohomology 1 is_consistent(r)
    end
    return r
  end


  function GModule(G::T, ac::Vector{<:Map}) where {T}
    return GModule(domain(ac[1]), G, ac)
  end

  function GModule(ac::Vector{<:Map})
    T = typeof(domain(ac[1]))
    M = new{Nothing, T}()
    M.M = domain(ac[1])
    M.ac = ac
    return M
  end

  F::Group # G as an Fp-group (if set)
  mF::GAPGroupHomomorphism  # G -> F, maps G[i] to F[i]

  iac::Vector{Map} # the inverses of ac
end

function Base.:(==)(F::GModule{gT,mT}, E::GModule{gT,mT}) where {gT, mT}
  F.G == E.G || return false
  F.M == E.M || return false
  return F.ac == E.ac
end

function Base.hash(F::GModule{gT,mT}, h::UInt) where {gT, mT}
  b = 0x535bbdbb2bc54b46%UInt
  h = hash(F.G, h)
  h = hash(F.M, h)
  h = hash(F.ac, h)
  return xor(h, b)
end

function Base.show(io::IO, C::GModule)
  @show_name(io, C)
  @show_special(io, C)

  io = pretty(io)
  io = IOContext(io, :compact => true)
  print(io, "G-module for ", Lowercase(), C.G, " acting on ", Lowercase(), C.M)# , "\nvia: ", C.ac)
end

"""
Given an automorphism of some module for each generator of the
group `H`, return the `ZZ[H]` module.

Note: we do not check that this defined indeed a `ZZ[H]` module.
"""
function gmodule(H::Union{Nothing, Oscar.GAPGroup}, ac::Vector{<:Map})
  return GModule(H, ac)
end

#in case the group is trivial, (ngens == 0), then length(ac)=0
#and the modules cannot be inferred. Thus a version with the
#module...
function gmodule(M, H::Union{Nothing, Oscar.GAPGroup}, ac::Vector{<:Map})
  return GModule(M, H, ac)
end

"""
Check if the action maps satisfy the same relations
as the generators of `G`.
"""
function is_consistent(M::GModule)
  G, mG = fp_group_with_isomorphism(M) # mG: group(M) -> G
  V = Module(M)
  R = relators(G)
  for r = R
    w = word(r)
    a = action(M, preimage(mG, w[1]< 0 ? inv(gen(G, -w[1])) : gen(G, w[1])))
    for i=2:length(w)
      a = a* action(M, preimage(mG, w[i]< 0 ? inv(gen(G, -w[i])) : gen(G, w[i])))
    end
    all(x->a(x) == x, gens(V)) || (@show r; return false)
  end

  return true
end
##########################################################
#
# Basics for gmodules
#
# access and action
#
AbstractAlgebra.Group(C::GModule) = C.G
AbstractAlgebra.Module(C::GModule) = C.M
action(C::GModule) = C.ac

function inv_action(C::GModule)
  if !isdefined(C, :iac)
    C.iac = map(inv, C.ac)
  end
  return C.iac
end

function fp_group_with_isomorphism(C::GModule)
  if !isdefined(C, :F)
    iso = isomorphism(FPGroup, group(C), on_gens=true)
    C.F, C.mF = codomain(iso), iso
  end
  return C.F, C.mF
end

#TODO? have a GModuleElem and action via ^?
"""
For an array of objects in the module, compute the image under the
action of `g`, ie. an array where each entry is mapped.
"""
function action(C::GModule, g, v::Array)
  @assert parent(g) == Group(C)
  if has_attribute(C, :is_trivial)
    return v
  end

  if has_attribute(C, :action)
    return get_attribute(C, :action)(C, g, v)
  end

  ac = action(C)
  f = findfirst(isequal(g), gens(Group(C)))
  if f !== nothing
    return map(ac[f], v)
  end

  iac = inv_action(C)
  f = findfirst(isequal(inv(g)), gens(Group(C)))
  if f !== nothing
    return map(iac[f], v)
  end

  F, mF = fp_group_with_isomorphism(C)
  for i = word(mF(g))
    if i > 0
      v = map(ac[i], v)
    else
      v = map(iac[-i], v)
    end
  end
  return v
end

function Base.deepcopy_internal(C::GModule, ::IdDict)
  return C
end

function Base.deepcopy_internal(G::FinGenAbGroup, ::IdDict)
  return G
end


"""
The image of `v` under `g`
"""
function action(C::GModule, g, v)
  return action(C, g, [v])[1]
end

import Base.:^
function ^(f::AbstractAlgebra.Generic.ModuleHomomorphism, n::Int64)
  if n < 0
    n = -n
    f = pseudo_inv(f)
  end
  a = id_hom(domain(f))
  b = f
  i = 1
  while i < n
    if n & i != 0
      a *= b
    end
    b = b*b
    i *= 2
  end
  return a
end

"""
The operation of `g` on the module as an automorphism.
"""
function action(C::GModule, g)
  @assert parent(g) == Group(C)
  if has_attribute(C, :is_trivial)
    return C.ac[1]
  end

  ac = action(C)
  G = Group(C)
  f = findfirst(isequal(g), gens(G))
  if f !== nothing
    return ac[f]
  end
  iac = inv_action(C)
  f = findfirst(isequal(inv(g)), gens(G))
  if f !== nothing
    return iac[f]
  end

  h = id_hom(C.M)
#  if isa(C.G, PcGroup)
#    return map_word(g, ac, genimgs_inv = iac, init = h)
#  end
  F, mF = fp_group_with_isomorphism(C)
  return map_word(mF(g), ac, genimgs_inv = iac, init = h)

  for i = word(mF(g))
    if i > 0
      h = h*ac[i]
#      v = map(ac[i], v)
    else
      h = h*iac[-i]
#      v = map(iac[-i], v)
    end
  end
  return h
end

function AbstractAlgebra.Generic.add_direct_sum_injection!(a::FinGenAbGroupElem, i::Int, b::FinGenAbGroupElem)
  C = get_attribute(parent(a), :direct_product)
  n = 0
  for j=1:i-1
    n += ngens(C[j])
  end
  _a = Hecke.FinGenAbGroupElem(C[i], _view_window(a.coeff, 1, n+1, 1, n+ngens(C[i])))
  add!(_a, _a, b)
  return a
end

function _coeffs(a::FinGenAbGroupElem)
  return a.coeff
end
function _coeffs(a::AbstractAlgebra.FPModuleElem)
  return a.v
end
function Oscar.zero!(a::AbstractAlgebra.FPModuleElem)
  zero!(a.v)
end

#TODO: write an action function that does not use create the matrix...
"""
For a Z[U]-Module C and a map U->G, compute the induced module:
    ind_U^G(C) = C otimes Z[G]
where the tensor product is over Z[U].
The induced module is returned as a product of copies of C. it also returns
  - the transversal used
  - the projections
  - the injections

  If D and mDC are given then mDC: D -> C.M has to be a Z[U] linear
homomorphism. I this case a Z[G] linear map to the induced module
is returned.
"""
function induce(C::GModule{<:Oscar.GAPGroup}, h::Map, D = nothing, mDC = nothing)
  U = domain(h)
  G = codomain(h)
  @assert U == C.G
  @assert D === nothing || mDC !== nothing
  @assert D === nothing || (domain(mDC) == D.M && codomain(mDC) == C.M &&
                            D.G == codomain(h))
  iU = image(h)[1]

# ra = right_coset_action(G, image(h)[1]) # will not always match
# the transversal, so cannot use.
# See https://github.com/gap-system/gap/issues/5337
# for a discussion whether to return both transversal and action on it.
  g = right_transversal(G, iU)
  S = symmetric_group(length(g))
  ra = hom(G, S, [S([findfirst(x->x*inv(z*y) in iU, g) for z = g]) for y in gens(G)])

  #= C is Z[U] module, we need
    C otimes Z[G]

    any pure tensor c otimes g can be "normalized" g = u*g_i for one of the
    reps fixed above, so c otimes g = c otimes u g_i == cu otimes g_i

    For the G-action we thus get
    (c otimes g_i)g = c otimes g_i g = c otimes u_i g_j (where the j comes
                                                         from the coset action)
                    = cu_i otimes g_j
  =#

#  @assert isdefined(C.M, :hnf)
  indC, pro, inj = direct_product([C.M for i=1:length(g)]..., task = :both)
#  @assert isdefined(indC, :hnf)
  if hasfield(typeof(indC), :__attrs)
    set_attribute!(indC, :induce => (h, g))
  end
  ac = []
  iac = []
  o = one(base_ring(C))
  for s = gens(G)
    sigma = ra(s)
    u = [ preimage(h, g[i]*s*g[i^sigma]^-1) for i=1:length(g)]
    au = [action(C, x) for x = u]
    @assert all(in(iU), u)
    im_q = []
    q = zero(indC)
    for _q = 1:ngens(indC)
      zero!(q)
      _coeffs(q)[_q] = o

      X = zero(indC)
      for i=1:length(g)
#        add!(X, inj[i^sigma](au[i](pro[i](q))))
        p = pro[i](q)
        if !iszero(p)
          AbstractAlgebra.Generic.add_direct_sum_injection!(X, i^sigma, au[i](p))
        end
      end
      @assert !iszero(X)
      push!(im_q, X)
#      push!(im_q, sum(inj[i^sigma](action(C, preimage(h, u[i]), pro[i](q))) for i=1:length(g)))
    end
    push!(ac, hom(indC, indC, [x for x = im_q]))

    s = inv(s)
    sigma = ra(s)
    u = [ preimage(h, g[i]*s*g[i^sigma]^-1) for i=1:length(g)]
    au = [action(C, x) for x = u]
    @assert all(in(iU), u)
    im_q = []
    q = zero(indC)
    for _q = 1:ngens(indC)
      zero!(q)
      _coeffs(q)[_q] = o
      X = zero(indC)
      for i=1:length(g)
        p = pro[i](q)
        if !iszero(p)
          AbstractAlgebra.Generic.add_direct_sum_injection!(X, i^sigma, au[i](p))
        end
      end
      @assert !iszero(X)
      push!(im_q, X)
#      push!(im_q, sum(inj[i^sigma](action(C, preimage(h, u[i]), pro[i](q))) for i=1:length(g)))
    end
    push!(iac, hom(indC, indC, [x for x = im_q]))

  end
  iC = GModule(G, [x for x = ac])
  iC.iac = [x for x = iac]
  if D === nothing
    return iC, g, pro, inj
  end
  #= for a Z[G]-modul D s.th. D has a Z[U]-lin embedding into C,
    compute the Z[G]-lin embedding into the induced module.
    a -> sum a g_i^-1 otimes g_i
    works (direct computation with reps and cosets)
  =#
  h = hom(D.M, iC.M, [sum(inj[i](mDC(action(D, inv(g[i]), h))) for i=1:length(g)) for h = gens(D.M)])
  return iC, h
end

function Oscar.quo(C::GModule{<:Any, <:AbstractAlgebra.FPModule}, mDC::Generic.ModuleHomomorphism)
  q, mq = Oscar.quo(C.M, image(mDC)[1])
  S = GModule(C.G, [hom(q, q, [mq(x(preimage(mq, t))) for t = gens(q)]) for x = C.ac])
  return S, mq
end

function Oscar.quo(C::GModule, mDC::Map{FinGenAbGroup, FinGenAbGroup}, add_to_lattice::Bool = true)
  q, mq = Oscar.quo(C.M, image(mDC)[1], add_to_lattice)
  S = GModule(C.G, [FinGenAbGroupHom(pseudo_inv(mq)*x*mq) for x = C.ac])
  if isdefined(C, :iac)
    S.iac = [FinGenAbGroupHom(pseudo_inv(mq)*x*mq) for x = C.iac]
  end
  return S, mq
end

function Oscar.direct_product(C::GModule...; task::Symbol = :none)
  @assert task in [:sum, :prod, :both, :none]
  G = C[1].G
  @assert all(x->x.G == G, C)
  mM, pro, inj = direct_product([x.M for x = C]..., task = :both)

  mC = gmodule(G, [hom_direct_sum(mM, mM, [action(C[i], g) for i=1:length(C)]) for g = gens(G)])
  mC.iac = [hom_direct_sum(mM, mM, [action(C[i], inv(g)) for i=1:length(C)]) for g = gens(G)]

  if task == :none
    return mC
  elseif task == :sum
    return mC, inj
  elseif task == :prod
    return mC, pro
  else
    return mC, pro, inj
  end
end

function Oscar.one(f::FinGenAbGroupHom)
  A = domain(f)
  @assert A === codomain(f)
  return hom(A, A, gens(A))
end

function Oscar.one(f::Generic.ModuleHomomorphism)
  A = domain(f)
  @assert A === codomain(f)
  return hom(A, A, gens(A))
end

Oscar.direct_sum(C::GModule...; task::Symbol = :none) = Oscar.direct_product(C...; task = task)

import Hecke.⊕
⊕(C::GModule...) = Oscar.direct_sum(C...; task = :none)

function Oscar.tensor_product(C::GModule{<:Any, FinGenAbGroup}...; task::Symbol = :map)
  @assert all(x->x.G == C[1].G, C)

  T, mT = Oscar.tensor_product([x.M for x = C]...; task = :map)
  TT = gmodule(T, C[1].G, [hom_tensor(T, T, [action(C[i], g) for i=1:length(C)]) for g = gens(C[1].G)])
  if task == :map
    return TT, mT
  else
    return TT
  end
end

function Oscar.tensor_product(C::GModule{T, <:FPModule}, Cs::GModule{T, <:FPModule}...; task::Symbol = :map) where {T <: GAPGroup}
  return Oscar.tensor_product(GModule{T, <:FPModule}[C, Cs...]; task)
end

function Oscar.tensor_product(C::Vector{<:GModule{<:GAPGroup, <:FPModule}}; task::Symbol = :map)
  @assert all(x->x.G == C[1].G, C)
  @assert all(x->base_ring(x) == base_ring(C[1]), C)

  T, mT = Oscar.tensor_product([x.M for x = C]...; task = :map)
  #this is slow. It should be done like in abelian groups using
  #the tensor product of the maps
  TT = gmodule(T, C[1].G, [hom(T, T, [mT(Tuple(action(C[i], g, preimage(mT, t)[i]) for i=1:length(C))) for t = gens(T)] ) for g = gens(C[1].G)])
  if task == :map
    return TT, mT
  else
    return TT
  end
end

#TODO: cannot extend as it is not yet defined...
#function Oscar.tensor_power(C::GModule, n::Int; task = :none)
#  return tensor_product([C for i=1:n]; task)
#end

#TODO: do not use the tensor product, directly write down
#      the module and the action (both symmetric and alternating)
#TODO: higher powers?
function symmetric_square(C::GModule)
  t, mt = tensor_product(C, C; task = :map)
  g = gens(C.M)
  return sub(t, sub(t.M, [mt((g[i], g[j])) + mt((g[j], g[i])) for i=1:length(g) for j=1:i])[2])[1]
end

export symmetric_square

function alternating_square(C::GModule)
  t, mt = tensor_product(C, C; task = :map)
  g = gens(C.M)
  return sub(t, sub(t.M, [mt((g[i], g[j])) - mt((g[j], g[i])) for i=2:length(g) for j=1:i-1])[2])[1]
end

export alternating_square


import Hecke.⊗
⊗(C::GModule...) = Oscar.tensor_product(C...; task = :none)


function Oscar.tensor_product(F::FPModule{T}, Fs::FPModule{T}...; task = :none) where {T}
  return Oscar.tensor_product([F, Fs...]; task)
end
function Oscar.tensor_product(F::Vector{<:FPModule{T}}; task = :none) where {T}
  @assert all(x->base_ring(x) == base_ring(F[1]), F)
  d = prod(dim(x) for x = F)
  G = free_module(base_ring(F[1]), d)
  if task == :none
    return G
  end

  g = vec(collect(Base.Iterators.ProductIterator(Tuple(gens(g) for g = reverse(F)))))

  function pure(g...)
    @assert length(g) == length(F)
    @assert all(i-> parent(g[i]) == F[i], 1:length(F))

    return G(vec(collect(prod(x) for x = Base.Iterators.product([h.v for h = reverse(g)]...))))
  end
  function pure(T::Tuple)
    return pure(T...)
  end

  function inv_pure(t::Generic.FreeModuleElem)
    p = Base.findall(i -> !iszero(t[i]), 1:ngens(G))
    if length(p) == 0
      return Tuple(collect(zero(g) for g = F))
    end
    @assert length(p) == 1
    @assert t[p[1]] == 1
    return reverse(g[p[1]])
  end

  return G, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = F])), G, pure, inv_pure)
end

⊗(F::Generic.FreeModule...) = Oscar.tensor_product(F...; task = :none)

function Oscar.restrict(C::GModule, U::Oscar.GAPGroup)
  fl, m = is_subgroup(U, C.G)
  @assert fl
  return gmodule(U, [action(C, m(g)) for g = gens(U)])
end

function Oscar.restrict(C::GModule, m::Map)
  U = domain(m)
  return gmodule(U, [action(C, m(g)) for g = gens(U)])
end

function Oscar.inflate(C::GModule, h)
  G = domain(h)
  U = codomain(h)
  @assert U == group(C)
  return gmodule(G, [action(C, h(g)) for g = gens(G)])
end

export GModule, gmodule, word, confluent_fp_group
export action, cohomology_group, extension
export induce, is_consistent, istwo_cocycle, all_extensions
export split_extension, extension_with_abelian_kernel
export CoChain, MultGrp, MultGrpElem
export is_stem_extension, is_central

_rank(M::FinGenAbGroup) = torsion_free_rank(M)
_rank(M) = rank(M)
_rank(M::AbstractAlgebra.FPModule{<:FieldElem}) = dim(M)

Oscar.dim(C::GModule) = _rank(C.M)
Oscar.base_ring(C::GModule) = base_ring(C.M)
Oscar.base_ring_type(::Type{GModule{gT, mT}}) where {gT, mT} = base_ring_type(mT)
Oscar.group(C::GModule) = C.G

###########################################################
#
# Supporting group theory
# To be moved and revised eventually
###########################################################

"""
For an element of an fp-group, return a corresponding word as a sequence
of integers. A positive integers indicates the corresponding generator,
a negative one the inverse.
"""
function word(y::FPGroupElem)
  # TODO: get rid of this
  return letters(y)
end

"""
The relations defining 'G' as an array of pairs.
"""
function Oscar.relations(G::Oscar.GAPGroup)
  rels = relators(G)
  if length(rels) == 0
    T = eltype(rels)
    return Tuple{T,T}[]
  end
  z = one(parent(rels[1]))
  return [(x, z) for x in rels]
end


######################################################
#
#
# Main goal: cohomology computations.
#
# So "empty" structure for parent of co-chains and a co-chain.
# currently co-chains are dumb: the values need to all be known
# on creation. They should support lazy filling.
# Also possibly they should be "modules", ie. inherit addition
# and scalar multiplication (and possibly? the G-operation)
#

struct AllCoChains{N, G, M} #Int (dim), Group(elem), Module(elem)
end

struct CoChain{N, G, M}
  C::GModule
  d::Dict{NTuple{N, G}, M}
  D::Function

  function CoChain{N, G, M}(C::GModule, d::Dict{NTuple{N, G}, M}) where {N, G, M}
    return new{N, G, M}(C, d)
  end

  function CoChain{N, G, M}(C::GModule, d::Dict{NTuple{N, G}, M}, D) where {N, G, M}
    return new{N, G, M}(C, d, D::Function)
  end
end

function Base.show(io::IO, C::CoChain{N}) where {N}
  print(io, "$N-cochain with values in ", C.C.M)
end

function Base.keys(C::CoChain{N}) where {N}
  return Iterators.ProductIterator(Tuple([C.C.G for i=1:N]))
end

Oscar.Nemo.elem_type(::Type{AllCoChains{N,G,M}}) where {N,G,M} = CoChain{N,G,M}
Oscar.Nemo.parent_type(::Type{CoChain{N,G,M}}) where {N,G,M} = AllCoChains{N,G,M}
Oscar.parent(::CoChain{N,G,M}) where {N, G, M} = AllCoChains{N, G, M}()

function differential(C::CoChain{N, G, M}) where {N, G, M}
  n = N+1
  d = Pair{NTuple{n, G}, M}[]
  for g = Iterators.ProductIterator(Tuple([C.C.G for i=1:n]))
    v = action(C.C, g[n], C(Tuple(g[i] for i=1:n-1))) +
        sum((-1)^i* C(Tuple(vcat(G[g[j] for j=1:i-1], G[g[i]*g[i+1]], G[g[j] for j=i+2:n]))) for i=1:n-1; init = zero(C.C.M)) +
        (-1)^n*C(Tuple(g[i] for i=2:n))
    push!(d, g=>v)
  end
  return CoChain{n, G, M}(C.C, Dict(d))
end

function action(C::CoChain, g::PermGroupElem)
  CC = deepcopy(C)
  if isdefined(C, :D)
    for x = keys(CC.d)
      CC.d[x] = action(C.C, g, C.d[x])
    end
    CC.D = x->action(C.C, g, C.D(x))
  else
    for x = keys(C)
      CC.d[x] = action(C.C, g, C[x])
    end
  end

  return CC
end

function Oscar.restrict(c::CoChain{N, G, M}, m::Map) where {N, G, M}
  if isdefined(c, :D)
    return CoChain{N, G, M}(restrict(c.C, m), c.d, c.D)
  else
    return CoChain{N, G, M}(restrict(c.C, m), c.d)
  end 
end

#TODO: careful with lazy chains!!!
function ==(c::CoChain{N, G, M}, d::CoChain{N, G, M}) where {N, G, M}
  @assert c.C === d.C
  @assert !isdefined(c, :D)
  return all(c(x) == d(x) for x = keys(c))
end

function Base.hash(c::CoChain{N, G, M}, h::UInt) where {N, G, M}
  # this is a very bad hash, but it is correct
  return hash(c.C, h)
end

function +(c::CoChain{N, G, M}, d::CoChain{N, G, M}) where {N, G, M}
  @assert c.C === d.C
  if isdefined(c, :D) && isdefined(d, :D)
    return CoChain{N, G, M}(c.C, Dict( x=> c(x) + d(x) for x = keys(c.d)), x->c.D(x) + d.D(x))
  else
    return CoChain{N, G, M}(c.C, Dict( x=> c(x) + d(x) for x = keys(c)))
  end
end


function -(c::CoChain{N, G, M}, d::CoChain{N, G, M}) where {N, G, M}
  @assert c.C === d.C
  if isdefined(c, :D) && isdefined(d, :D)
    return CoChain{N, G, M}(c.C, Dict( x=> c(x) - d(x) for x = keys(c.d)), x->c.D(x) + d.D(x))
  else
    return CoChain{N, G, M}(c.C, Dict( x=> c(x) - d(x) for x = keys(c)))
  end
end

"""
Evaluate a 0-cochain
"""
(C::CoChain{0})() = first(values(C.d))
(C::CoChain{0})(::Tuple{}) = first(values(C.d))

#TODO: should this rather be a map from a 1-tuple of group elements?
"""
Evaluate a 1-cochain, a 1-cochain is a map from the group into the
module
"""
function (C::CoChain{1})(g::Oscar.BasicGAPGroupElem)
  if haskey(C.d, (g,))
    return C.d[(g,)]
  end
  F, mF = fp_group_with_isomorphism(C.C)
  G = parent(g)
  @assert G == group(C.C)
  @assert ngens(F) == ngens(G)
  @assert all(i-> preimage(mF, gen(F, i)) == gen(G, i), 1:ngens(G))
  w = word(mF(g))
  t = zero(Module(C.C))
  ac = action(C.C)
  iac = inv_action(C.C)
  G = Group(C.C)
  #TODO: build up the group element step by step
  #      and store the values: (use Dimino code in Hecke)
  #XXX: this is wrong!, compare to the H_one code below
  # problem is that F and G might have different gens
  # needs different Gap code: write g as a word in the
  # generators of G and use this.
  # also inverses are more complicated.
  for i = w
    if i > 0
      t = ac[i](t)+C.d[(gen(G, i),)]
    else
      t = iac[-i](t-C.d[(gen(G, -i),)])
    end
  end
  C.d[(g,)] = t
  return t
end
(C::CoChain{1})(g::NTuple{1, <:Oscar.BasicGAPGroupElem}) = C(g[1])

#should support lazy via call-back.
#which should be used by default...
"""
Evaluate a 2-cochain, a 2-cochain is a map from pairs of group elements
into the module
"""
function (C::CoChain{2})(g::Oscar.GAPGroupElem, h::Oscar.GAPGroupElem)
  if haskey(C.d, (g,h))
    return C.d[(g,h)]
  end
  @assert isdefined(C, :D)
  return C.d[(g,h)] = C.D((g, h))
end
(C::CoChain{2})(g::NTuple{2, <:Oscar.GAPGroupElem}) = C(g[1], g[2])

#TODO: re-write to get the maps! To support Q/Z as well
"""
  H^0(G, M)

Return a module (same type as M) that abstractly represent
the 0-cohomology group as well as a map realizing this via explicit
co-chains
"""
function H_zero(C::GModule)
  z = get_attribute(C, :H_zero)
  if z !== nothing
    return domain(z), z
  end
  G = Group(C)
  M = Module(C)
  id = hom(M, M, gens(M))
  ac = action(C)
  k = kernel(id - ac[1])[1]
  for i=2:length(ac)
    k = intersect(k, kernel(id - ac[i])[1])
    if isa(k, Tuple) #GrpAb: intersect yield ONLY intersection,
      k = k[1]
    end
  end
  fl, inj = is_sub_with_data(k, M)
  @assert fl

  #this is fix, now it "should" be mod by norm? Only for Tate, see below
  z = MapFromFunc(k, AllCoChains{0,elem_type(G),elem_type(M)}(), x->CoChain{0,elem_type(G),elem_type(M)}(C, Dict(() => inj(x))), y->preimage(inj, y()))
  set_attribute!(C, :H_zero => z)
  return k, z
end

function H_zero_tate(C::GModule)
  z = get_attribute(C, :H_zero_tate)
  if z !== nothing
    return domain(z), z
  end
  G = Group(C)
  M = Module(C)
  #fix under action modulo norm (trace) = sum over all elem in group

  id = hom(M, M, gens(M))
  ac = action(C)
  k = kernel(id - ac[1])[1]
  for i=2:length(ac)
    k = intersect(k, kernel(id - ac[i])[1])
  end
  N = sum(action(C, g) for g = group(C))

  i = image(N)[1]
  fl, inj = is_sub_with_data(i, k)
  q, mq = quo(k, image(inj)[1])
  fl, Inj = is_sub_with_data(k, M)

  z = MapFromFunc(q, AllCoChains{0,elem_type(G),elem_type(M)}(), x->CoChain{0,elem_type(G),elem_type(M)}(C, Dict(() => Inj(preimage(mq, x)))), y->mq(preimage(Inj, y())))
  set_attribute!(C, :H_zero_tate => z)

  if isfinite(G) && isa(q, FinGenAbGroup)
    q.exponent = order(G)
  end

  return q, z
end

Oscar.is_sub_with_data(G::FinGenAbGroup, H::FinGenAbGroup) = is_subgroup(G, H)
#= TODO
 - break out coboundaries and cochains
 - depending on the module type:
   - intersect yields an embedding (Z-module) or not GrpAb
   - make sure that image/ kernel are consistent
   - preimage
   - issubset yields (for GrpAb) only true/ false, not the map
   - is_subgroup cannot apply to modules
   - quo does ONLY work if B is a direct submodule of A (Z-modules)
   - matrix is used to get "the matrix" from a hom
   - zero_hom/ zero_obj/ identity_hom is missing
   - Janko-Module-Homs have different types, they probably need to
     come under a common abstract type or be more selective
=#


"""
Code of the H^1(G, M) computation:
returns homomorphisms A and B s.th.

   M_1 -A-> M_2 -B-> M_3

satisfies

  H^1 = kern(B)/image(A)

Or, kern(B) are the 1-co-chains, image(A) the 1-co-boundaries.

If M is a free abelian group (Z^n), then this is used in the solvable
quotient to compute the H^1 of Q^n/Z^n via duality.
"""
function H_one_maps(C::GModule; task::Symbol = :maps)
  @assert task in [:maps, :all]
  #= idea, after Holt:
  H^1 = crossed homs. due to action on the right
  f(ab) = f(a)^b + f(b)
  if G=<g_1, ..., g_r | r_1, ..., r_l>
  then X in H^1 iff X(r_i) = 0 for all i
  X:G->M is given as X in M^r, where X(g_i) = X[i]
  X(r_i) corresponds to some map phi_i : M^r -> M
  phi_i = oplus h_j M for some homs h_j coming from the word in r
  so, a kernel computation again
  =#

  G = Group(C)
  n = ngens(G)
  M = Module(C)
  D, pro, inj = direct_product([M for i=1:n]..., task = :both)

  F, mF = fp_group_with_isomorphism(C)
  @assert ngens(F) == ngens(G)
  @hassert :GroupCohomology 1 all(i->mF(gen(F, i)) == gen(G, i), 1:ngens(G))

  R = relators(F)
#  @assert G == F

  K = D
  ac = action(C)
  iac = inv_action(C)
  idM = hom(M, M, gens(M)) #identity map to start with
                           #TODO: require an identity_hom constructor

  Kr, pKr, iKr = direct_product([M for i=R]..., task = :both)
  gg = nothing
  i = 1
  for r = R
    W = word(r)
    g = idM
    P = hom(D, M, [zero(M) for i=1:ngens(D)])
    for w in W
      if w < 0
        #by above: f(ab) = f(a)^b + f(b)
        #thus      0 = f(1) = f(a a^-1) = f(a)^(a^-1) + f(a^-1)
        P = P*iac[-w]-pro[-w]*iac[-w]
        g = g*iac[-w]
      else
        P = P*ac[w]+pro[w]
        g = g*ac[w]
      end
    end
    @hassert :GroupCohomology 1 all(x -> x == g(x), gens(M))
    if gg === nothing
      gg = P*iKr[i]
    else
      gg += P*iKr[i]
    end
    i += 1
  end
  #K is Z[1]  - the co-cycles
  #TODO: is kernel(g) directly faster than the method above (H_zero)
  #      where kernel(g) is computed slice by slice?
  #TODO: cache the expensive objects!!!

  @vtime :GroupCohomology 1 g = sum((ac[i] - idM)*inj[i] for i=1:n)
  if task == :all
    return g, gg, pro, inj, mF
  else
    return g, gg
  end
end

"""
  H^1(G, M)

Return an abstract module (of the same type as M) describing the
first co-homology group. Furthermore, the second return value
is a map realising elements of H^1 as explicit co-cycles.
"""
function H_one(C::GModule)
  z = get_attribute(C, :H_one)
  if z !== nothing
    return domain(z), z
  end

  g, gg, pro, inj, _ = H_one_maps(C, task = :all)

  K = kernel(gg)[1]
  D = domain(gg)
  lf, lft = is_sub_with_data(K, D)
  @assert lf

  Q, mQ = quo(K, image(g)[1])

  M = Module(C)
  G = group(C)
  n = ngens(G)

  z = MapFromFunc(
    Q, AllCoChains{1, elem_type(G), elem_type(M)}(),
    x->CoChain{1,elem_type(G),elem_type(M)}(C, Dict([(gen(G, i),) => pro[i](lft(preimage(mQ, x))) for i=1:ngens(G)])),
    y->mQ(preimage(lft, sum(inj[i](y(gen(G, i))) for i=1:n))))

  set_attribute!(C, :H_one => z)
  return Q, z
  #need to ALSO return the coboundary(s)
end


function confluent_fp_group_pc(G::Oscar.GAPGroup)
   @vtime :GroupCohomology 1 g = isomorphism(PcGroup, G)
   P = codomain(g)
   f = GAPWrap.IsomorphismFpGroupByPcgs(GAP.Globals.FamilyPcgs(GapObj(P)), GAP.Obj("g"))
   @req f != GAP.Globals.fail "Could not convert group into a group of type FPGroup"
   H = FPGroup(GAPWrap.Image(f))
   R = relations(H)
   ru = Vector{Tuple{Vector{Int}, Vector{Int}}}()
   for r = R
     push!(ru, (map(Int, GAP.Globals.LetterRepAssocWord(GapObj(r[1]))),
                map(Int, GAP.Globals.LetterRepAssocWord(GapObj(r[2])))))
  end
  i = 0
  ex = []
  for r = ru
    i += 1
    @assert length(r[2]) == 0
    if r[1][1] == r[1][2] #order relation!
      j = 2
      while j <= length(r[1]) && r[1][j] == r[1][1]
        j += 1
      end
      ru[i] = (r[1][1:j-1], -1 .* reverse(r[1][j:end]))
      r = ru[i]
      push!(ex, ([-r[1][1]], vcat([r[1][i] for i=2:length(r[1])], -r[2])))
    else #conjugator rel
      @assert r[1][1] < 0 && -r[1][1] == r[1][3]
      @assert r[1][2] < 0 && -r[1][2] == r[1][4]
      ru[i] = ([r[1][3], r[1][4]], vcat([r[1][4], r[1][3]], -1 .* reverse(r[1][5:end])))
    end
  end
  append!(ru, ex)

  @vtime :GroupCohomology 1 return H, GAPGroupHomomorphism(H, P, GAP.Globals.InverseGeneralMapping(f))*inv(g), ru
end


"""
Compute an isomorphic fp-group and a confluent system of
relations given as pairs of words.

Return the new group, the isomorphism and the confluent relations.
"""
function confluent_fp_group(G::Oscar.GAPGroup)
  @vtime :GroupCohomology 1 C = GAP.Globals.ConfluentMonoidPresentationForGroup(GapObj(G))
  #has different generators than G! So the action will have to
  #be adjusted to those words. I do not know if a RWS (Confluent) can
  #just be changed...
  k = C.monhom #[2] #hopefully the monhom entry in 4.12 it will be the name
  M = GAPWrap.Range(k)
  g = [GAP.Globals.PreImageElm(k, x) for x = GAP.Globals.GeneratorsOfMonoid(M)]
  g = map(GAPWrap.UnderlyingElement, g)
  g = map(GAP.Globals.LetterRepAssocWord, g)
  @assert all(x->length(x) == 1, g)
  g = map(x->Int(x[1]), g)
  R = GAP.Globals.RelationsOfFpMonoid(M)

  ru = Vector{Tuple{Vector{Int}, Vector{Int}}}()
  for r = R
    push!(ru, (map(x->g[Int(x)], GAP.Globals.LetterRepAssocWord(r[1])),
               map(x->g[Int(x)], GAP.Globals.LetterRepAssocWord(r[2]))))
  end

  #now to express the new gens as words in the old ones:

  Fp = FPGroup(GAPWrap.Range(C.fphom))
  return Fp, GAPGroupHomomorphism(Fp, G, GAP.Globals.InverseGeneralMapping(C.fphom)), ru
end


#############################
#
# H^2
#
# Thought of as group extensions given via a rewriting system.
# so we need collection...
#
mutable struct CollectCtx
  r::Vector{Tuple{Vector{Int}, Vector{Int}}} #the rules, RWS

  d1::Dict{Int, Int} #rules where lhs has length 1

  d2::Dict{Tuple{Int, Int}, Vector{Int}} # length 2 prefixes

  f::Function #(w::Vector{Int}, r::Int, p::Int)
              #to be called in addition (to play with the tail(s))
              #w the word, to be "reduced" using rule no r at pos p

  T::Any
  function CollectCtx(R::Vector{Tuple{Vector{Int}, Vector{Int}}})
    n = new()
    n.r = R
    n.d1 = Dict{Int, Int}()
    n.d2 = Dict{Tuple{Int, Int}, Vector{Int}}()
    for i = 1:length(R)
      r = R[i]
      if length(r[1]) == 1
        #still confused about that one...
        # but I have a rule [-1] -> [1, 2]
#        @assert length(r[2]) == 1
        n.d1[r[1][1]] = i
        continue
      end
      @assert length(r[1]) > 1
      p = (r[1][1], r[1][2])
      if Base.haskey(n.d2, p)
        push!(n.d2[p], i)
      else
        n.d2[p] = [i]
      end
    end
    for p = keys(n.d2)
      sort!(n.d2[p], lt = (a,b) -> isless(R[a], R[b]))
    end
    return n
  end
end

function Base.collect(w::Vector{Int}, C::CollectCtx)
  d1 = C.d1
  d2 = C.d2
  R = C.r
  do_f = isdefined(C, :f)

  nc = 0
  i = 1
  while true
    nc += 1
    if i>length(w)
      return w
    end
    if haskey(d1, w[i])
      if do_f
        C.f(C, w, d1[w[i]], i)
      end
      w = vcat(view(w, 1:i-1), R[d1[w[i]]][2], view(w, i+1:length(w)))
      i = 1
      continue
    end

    if i>=length(w)
      return w
    end


    if haskey(d2, (w[i], w[i+1]))
      for r = d2[(w[i], w[i+1])]
        if length(R[r][1]) + i-1 <= length(w) &&
           R[r][1] == view(w, i:i+length(R[r][1])-1)
          if do_f
            C.f(C, w, r, i)
          end
          w = vcat(view(w, 1:i-1), R[r][2], view(w, i+length(R[r][1]):length(w)))
          i = 0
          break
        end
      end
    end
    i += 1
  end
  return w
end

function H_two_maps(C::GModule; force_rws::Bool = false, redo::Bool = false)
  H_two(C; force_rws, redo, maps_only = true)
  return get_attribute(C, :H_two_maps)
end

#= Hulpke-Dietrich:
UNIVERSAL COVERS OF FINITE GROUPS
https://arxiv.org/pdf/1910.11453.pdf
almost the same as Holt
=#
#TODO: lazy = true, or even remove it
function H_two(C::GModule; force_rws::Bool = false, redo::Bool = false, lazy::Bool = !false, maps_only::Bool = false)
  z = get_attribute(C, :H_two)
  if !redo && z !== nothing
    return domain(z[1]), z[1], z[2]
  end

  G = Group(C)
  M = Module(C)

  @vprint :GroupCohomology 1 "starting H^2 for group of size $(order(G)) and module with $(ngens(M)) gens\n"

  if is_finite(M) && gcd(order(G), order(M)) == 1
    @show :all_trivial
    H2 = quo(M, M)[1]
    z = MapFromFunc(H2, AllCoChains{2,elem_type(G),elem_type(M)}(),
    x->CoChain{2, elem_type(G),elem_type(M)}(C, Dict{NTuple{2, elem_type(G)}, elem_type(M)}(), t->zero(M)), 
    y->zero(H2))
    #TODO: check if this is indeed correct (it is not!!!)
    #      re-organize to do minimal work per task
    #      add the other returned stuff in the trivial case as well
    zz = function(x)
      error("broken")
      return true, CoChain{1,elem_type(G),elem_type(M)}(C, Dict([(gen(G, i),) =>zero(M) for i=1:ngens(G)]))
    end
    set_attribute!(C, :H_two => (z, zz))

    return H2, z, zz
  end

  id = hom(M, M, gens(M), check = false)
  @vtime :GroupCohomology 1 F, mF = fp_group_with_isomorphism(C) #mF: F -> G

  if !force_rws && (isa(G, PcGroup) || is_solvable(G))
    @vprint :GroupCohomology 2 "using pc-presentation ...\n"
    FF, mFF, R = confluent_fp_group_pc(G) #mFF: FF -> G
    use_pc = true
  else
    @vprint :GroupCohomology 2 "using generic rws ...\n"
    FF, mFF, R = confluent_fp_group(G) #mFF: FF -> G
    use_pc = false
  end
  #now map the action generators (for gens(G)) to the gens for the RWS
  ac = []
  iac = []
  @vprint :GroupCohomology 2 "computing action for the gens in the rws..\n"
  @vtime :GroupCohomology 2 for g = gens(FF)
    f = action(C, mFF(g))
    push!(ac, f)
    #should we inv(f) or build inv f as a product as above???
    @vtime :GroupCohomology 3 push!(iac, inv(f))
  end

  c = CollectCtx(R)

  #rules with length(LHS) == 1 and rules of the form
  # [a a^-1] -> [], [a^-1 a] -> [] do not get tails
  pos = Vector{Int}()
  n = 0
  for i = 1:length(R)
    r = R[i]
    if length(r[1]) == 1
      push!(pos, 0)
      continue
    end
    if length(r[1]) == 2 && length(r[2]) == 0 && r[1][1] == -r[1][2]
      push!(pos, 0)
      continue
    end
    n += 1
    push!(pos, n)
  end

  @vprint :GroupCohomology 1 "need $n tails\n"

  if n == 0
    D = sub(M, elem_type(M)[])[1]
    pro = []
    inj = []
  else
    D, pro, inj = direct_product([M for i=1:n]..., task = :both)
  end


  #when collecting (i.e. applying the RWS we need to also
  #use the tails:  g v h -> gh h(v)
  #and if [gh] -> [x] with tail t, then
  #       gh -> x t, so
  #       g v h -> gh h(v) -> x t+h(v)
  # Hulpke calls this the normal version: reduced group word
  # at the beginning, module at the end, the tail.
  # collect will call the extra function c.f if set in the
  # CollectCtx
  # if use_rws we investigate all overlaps
  # otherwise, we know it a pc-presentation and thus fewer tests
  # are needed.
  is_triv = has_attribute(C, :is_trivial)

  function symbolic_collect(C::CollectCtx, w::Vector{Int}, r::Int, p::Int)
    #w = ABC and B == r[1], B -> r[2] * tail[r]
    # -> A r[2] C C(tail)
    # C = c1 c2 ... C(tail):
#    @assert w[p:p+length(R[r][1])-1] == R[r][1]

    if pos[r] == 0
      return
    end
#    T = pro[pos[r]]
    T = nothing
    if is_triv
      C.T += pro[pos[r]]
    else
      for i=view(w, p+length(R[r][1]):length(w))
        if i < 0
          T = (T === nothing) ? iac[-i] : T*iac[-i]
  #        T = T*iac[-i]
        else
          T = (T === nothing) ? ac[i] : T * ac[i]
  #        T = T*ac[i]
        end
      end
      C.T += T === nothing ? pro[pos[r]] : pro[pos[r]] * T
    end
  end
  c.f = symbolic_collect

  E = D
  #need zero hom this is too slow
  Z = hom(D, M, [zero(M) for i=1:ngens(D)], check = false)
  all_T = typeof(Z)[]

  @vprint :GroupCohomology 2 "building relations...\n"
  for i = 1:length(R)
    r = R[i]
    for j=1:length(R)
#      i == j && continue
      s = R[j]
      #we want overlaps, all of them:
      #r[1] = AB, s[1] = BC this is what we need to find...
      #(then we call collect on r[2]C and As[2] they should agree)
      #if use_pc, then l can only be 1:
      #From Holt: those are the r[1] we need
      # (i i .. i)
      #   (i .. i i)
      # (i .. i)
      #      (i j)
      # (i j)
      #  ( j ... j)
      # (i j)
      #   (j i)
      if use_pc
        l_max = 1
      else
        l_max = min(length(s[1]), length(r[1]))
      end

      for l=1:l_max
        if view(r[1], length(r[1])-l+1:length(r[1])) == view(s[1], 1:l)
          #TODO  AB    -> Ss  s,t are tails
          #       BC   -> Tt
          #      (AB)C -> SsC -> SC C(s)
          #      A(BC) -> ATt -> AT t
          if pos[i] > 0
            if is_triv
              c.T = pro[pos[i]]
            else
              first = true
              for h = view(s[1], l+1:length(s[1]))
                if h < 0
                  if first
                    c.T = iac[-h] 
                    first = false
                  else
                    c.T = c.T * iac[-h]
                  end
                else
                  if first
                    c.T = ac[h] 
                    first = false
                  else
                    c.T = c.T * ac[h]
                  end
                end
              end
              if first
                c.T = pro[pos[i]]
              else
                c.T = pro[pos[i]] * c.T
              end
            end
          else
            c.T = Z
          end
          z1 = collect(vcat(r[2], view(s[1], l+1:length(s[1]))), c)
          T = c.T
          c.T = Z
          z2 = collect(vcat(view(r[1], 1:length(r[1])-l), s[2]), c)
          if pos[j] > 0
            c.T += pro[pos[j]]
          end
          @hassert :GroupCohomology 1 z1 == z2
          push!(all_T, T-c.T)
        end
      end
    end
  end

  @vprint :GroupCohomology 1 "found $(length(all_T)) relations\n"

  if length(all_T) == 0
    Q = sub(M, elem_type(M)[])[1]
    mm = hom(D, Q, elem_type(Q)[zero(Q) for i=1:ngens(D)], check = false)
  else
    mm = Oscar.direct_sum(all_T)
    Q = codomain(mm)
  end
  if length(ac) == 0
    B = sub(M, elem_type(M)[])[1]
    B_pro = []
    B_inj = []
  else
    B, B_pro, B_inj = direct_product([M for i=1:length(ac)]..., task = :both)
  end
  CC = hom(B, D, elem_type(D)[zero(D) for i=1:ngens(B)], check = false)
  for i=1:length(R)
    if pos[i] == 0
      continue
    end
    r = R[i]
    if false && length(r[1]) == 1
      continue
    end
    #we have words r[1] and r[2] of shape g_1 g_2 ....
    #they need to be replaced by g_1 pro[1] g_2 pro[2]
    #and then sorted: g_1 pro[1] g_2 pro[2] ... ->
    #                 g_1 g_2 (pro[1] * g_2 + pro[2]) ...
    if r[1][1] < 0
      T = -B_pro[-r[1][1]]*iac[-r[1][1]]
    else
      T = B_pro[r[1][1]]
    end
    for j=2:length(r[1])
      if r[1][j] < 0
        T = (T-B_pro[-r[1][j]])*iac[-r[1][j]]
      else
        T = T*ac[r[1][j]] + B_pro[r[1][j]]
      end
    end

    if length(r[2]) == 0
      S = hom(B, M, [zero(M) for g = gens(B)], check = false)
    elseif r[2][1] < 0
      S = -B_pro[-r[2][1]]*iac[-r[2][1]]
    else
      S = B_pro[r[2][1]]
    end
    for j=2:length(r[2])
      if r[2][j] < 0
        S = (S-B_pro[-r[2][j]])*iac[-r[2][j]]
      else
        S = S*ac[r[2][j]] + B_pro[r[2][j]]
      end
    end

#    @assert issubset(image((T-S)*inj[pos[i]])[1], E)

    CC += (T-S)*inj[pos[i]]
  end

  set_attribute!(C, :H_two_maps => (CC, mm))
  if maps_only
    return
  end

  @vprint :GroupCohomology 2 "computing 2-cycles...\n"
  @vtime :GroupCohomology 2 E, mE = kernel(mm)
  @hassert :GroupCohomology 1 all(x->all(y->iszero(y(mE(x))), all_T), gens(E))
  @hassert :GroupCohomology 1 all(x->iszero(mm(mE(x))), gens(E))

  @vprint :GroupCohomology 2 "now the 2-boundaries...\n"
  @vtime :GroupCohomology 2 i, mi = image(CC)
  @vprint :GroupCohomology 2 "and the quotient...\n"
  @vtime :GroupCohomology 2 H2, mH2 = quo(E, i)
  if isfinite(G) && isa(H2, FinGenAbGroup)
    H2.exponent = order(G)
  end
  #we know |G| is an exponent - this might help

  function TailFromCoChain(cc::CoChain{2})
    #for all tails, ie. rules with pos[r]>0, we need to use
    #the 2-chain to compute inv(r[2])*r[1]
    #rule with tail: r[1] = r[2]*t, so t = r[2]^-1*r[1]
    T = zero(D)
    for r=1:length(pos)
      if pos[r] == 0
        continue
      end
      t1 = zero(M)
      g1 = one(G)
      w = R[r][1]
      for i=1:length(w)
        if w[i] > 0
          t1 = ac[w[i]](t1)+cc(g1, mFF(gen(FF, w[i])))
          g1 = g1*mFF(gen(FF, w[i]))
        else
          #need to mult by (w, 0)^-1 = (w^-1, -cc(w, w^-1))
          #so (g, t) (w, 0)^-1 = (g w^-1, t^w^-1 - cc(w, w^-1) + cc(g, w^-1)
          t1 = iac[-w[i]](t1)-cc(mFF(gen(FF, -w[i])), inv(mFF(gen(FF, -w[i]))))+cc(g1, inv(mFF(gen(FF, -w[i]))))
          g1 = g1*mFF(inv(gen(FF, -w[i])))
        end
      end

      t2 = zero(M)
      g2 = one(G)
      w = R[r][2]
      for i=1:length(w)
        if w[i] > 0
          t2 = ac[w[i]](t2)+cc(g2, mFF(gen(FF, w[i])))
          g2 = g2*mFF(gen(FF, w[i]))
        else
          #need to mult by (w, 0)^-1 = (w^-1, -cc(w, w^-1))
          #so (g, t) (w, 0)^-1 = (g w^-1, t^w^-1 - cc(w, w^-1) + cc(g, w^-1)
          t2 = iac[-w[i]](t2)-cc(mFF(gen(FF, -w[i])), inv(mFF(gen(FF, -w[i]))))+cc(g2, inv(mFF(gen(FF, -w[i]))))
          g2 = g2*mFF(inv(gen(FF, -w[i])))
        end
      end
      @assert g1 == g2

      #=
      w = [-x for x = reverse(R[r][2])]
      append!(w, R[r][1])
      #w is inv(r[2])*r[1]

      t = zero(M)
      g = one(G)
      for i=1:length(w)
        if w[i] > 0
          t = ac[w[i]](t)+cc(g, mFF(gen(FF, w[i])))
          g = g*mFF(gen(FF, w[i]))
        else
          #need to mult by (w, 0)^-1 = (w^-1, -cc(w, w^-1))
          #so (g, t) (w, 0)^-1 = (g w^-1, t^w^-1 - cc(w, w^-1) + cc(g, w^-1)
          t = iac[-w[i]](t)-cc(mFF(gen(FF, -w[i])), inv(mFF(gen(FF, -w[i]))))+cc(g, inv(mFF(gen(FF, -w[i]))))
          g = g*inv(gen(G, -w[i]))
        end
      end
      #= Maybe? Not clear what I actually want/ need here, darn
      wrong currently...
      if length(R[r][2]) > 0
        r2 = R[r][2][1] < 0 ? inv(gen(FF, -R[r][2][1])) : gen(FF, R[r][2][1])
        for i=2:length(R[r][2])
          r2 *= R[r][2][i] < 0 ? inv(gen(FF, -R[r][2][i])) : gen(FF, R[r][2][i])
        end
      else
        r2 = one(FF)
      end
      @show r2, mFF(r2)
      t = t + cc(mFF(inv(r2)), mFF(r2))
      =#
      =#
      T += inj[pos[r]](t1-t2)
    end
    return T
  end

  function TailToCoChain(t)
    c.f = function(C::CollectCtx, w::Vector{Int}, r::Int, p::Int)
      #w = ABC and B == r[1], B -> r[2] * tail[r]
      # -> A r[2] C C(tail)
      # C = c1 c2 ... C(tail):
      @assert w[p:p+length(R[r][1])-1] == R[r][1]

      if pos[r] == 0
        return
      end
      T = pro[pos[r]](t)
      for i=w[p+length(R[r][1]):end]
        if i < 0
          T = iac[-i](T)
        else
          T = ac[i](T)
        end
      end
      C.T += T
    end

    di = Dict{NTuple{2, elem_type(G)}, elem_type(M)}()
    #= if I figure out how to extend from generators
    w = [word(order(G) == 1 ? one(domain(mFF)) : preimage(mFF, g)) for g = gens(G)]
    for i=1:ngens(G)
      for j=1:ngens(G)
        c.T = zero(M)
        collect(vcat(w[i], w[j]), c)
        di[(gen(G, i), gen(G, j))] = c.T
      end
    end
    =#
    for g = G
      for h = G
        c.T = zero(M)
        if order(G) > 1
          gg = collect(word(preimage(mFF, g)), c)
          hh = collect(word(preimage(mFF, h)), c)
          c.T = zero(M)
          d = collect(vcat(gg, hh), c)
        end
        di[(g, h)] = c.T
      end
    end
    return CoChain{2,elem_type(G),elem_type(M)}(C, di)
  end

  function TailToCoChainLazy(t, g, h)
    c.f = function(C::CollectCtx, w::Vector{Int}, r::Int, p::Int)
      #w = ABC and B == r[1], B -> r[2] * tail[r]
      # -> A r[2] C C(tail)
      # C = c1 c2 ... C(tail):
      @assert w[p:p+length(R[r][1])-1] == R[r][1]

      if pos[r] == 0
        return
      end
      T = pro[pos[r]](t)
      for i=w[p+length(R[r][1]):end]
        if i < 0
          T = iac[-i](T)
        else
          T = ac[i](T)
        end
      end
      C.T += T
    end

    c.T = zero(M)
    if order(G) > 1
      gg = collect(word(preimage(mFF, g)), c)
      hh = collect(word(preimage(mFF, h)), c)
      c.T = zero(M)
      d = collect(vcat(gg, hh), c)
    end
    return c.T
  end

  symbolic_chain = function(g, h)
    c.f = symbolic_collect
    if order(G) == 1
      w = word(preimage(mFF, one(G)))
    else
      c.T = Z
      wg = collect(word(preimage(mFF, g)), c)
      wh = collect(word(preimage(mFF, h)), c)
      w = vcat(wg, wh)
    end
    c.T = Z
    @assert is_zero(Z)
    w = collect(w, c)
    return mE*c.T, w
  end

  set_attribute!(C, :H_two_symbolic_chain => (symbolic_chain, mH2))

  function is_coboundary(cc::CoChain{2}; reduce::Bool = false)
    t = TailFromCoChain(cc)
    fl, b = has_preimage_with_preimage(CC, t)
    if !fl
      return false, nothing
    end
    if reduce
      k, mk = kernel(CC)
      if !is_trivial(k)
        m = vcat(ZZMatrix[mk(x).coeff for x = gens(k) if !iszero(x)]...)
        m = lll(m)
        b = parent(b)(Hecke.MultDep.size_reduce(m, b.coeff))
      end
    end

    d = Dict{Tuple{elem_type(G), }, elem_type(M)}()
    # t gives, directly, the images of the generators (of FF)
    im_g = [B_pro[i](b) for i=1:ngens(FF)]
    # otherwise: sigma(g, h) + sigma(gh) = sigma(g)^h + sigma(h)
    # this gives the images for the inverses, and then for everything
    im_gi = [cc((mFF(gen(FF, i)), mFF(inv(gen(FF, i))))) - iac[i](im_g[i]) for i=1:ngens(FF)]
    @assert domain(mFF) == FF
    @assert codomain(mFF) == G == group(cc.C)
    for g = G
      m = zero(M)
      h = one(G)
      w = word(preimage(mFF, g))

      for i=1:length(w)
        if w[i] < 0
          m = iac[-w[i]](m)+im_gi[-w[i]]-cc((h, mFF(inv(gen(FF, -w[i])))))
          h = h*mFF(inv(gen(FF, -w[i])))
        else
          m = ac[w[i]](m)+im_g[w[i]]-cc((h, mFF(gen(FF, w[i]))))
          h = h*mFF(gen(FF, w[i]))
        end
      end
      d[(g,)] = m
      @assert g == h
    end
    return true, CoChain{1,elem_type(G),elem_type(M)}(C, d)
  end

  z2 = function(y)
    T = TailFromCoChain(y)
    return mH2(preimage(mE, T))
  end

  if !lazy
    z = (MapFromFunc(H2, AllCoChains{2,elem_type(G),elem_type(M)}(),
                   x->TailToCoChain(mE(preimage(mH2, x))), z2),
             is_coboundary)
           else
    z = (MapFromFunc(H2, AllCoChains{2,elem_type(G),elem_type(M)}(),
          x->CoChain{2, elem_type(G),elem_type(M)}(C, Dict{NTuple{2, elem_type(G)}, elem_type(M)}(), t->TailToCoChainLazy(mE(preimage(mH2, x)), t[1], t[2])),
                   z2),
#                   x->TailToCoChain(mE(preimage(mH2, x))), z2),
#                         y->TailFromCoChain(y), D, AllCoChains{2,elem_type(G),elem_type(M)}()),
             is_coboundary)
  end

  set_attribute!(C, :H_two => z)
  return H2, z[1], z[2]
  #now the rest...
  #(g, m)*(h, n) = (gh, m^h+n+gamma(g, h)) where gamma is "the" 2-cocycle
  #using tails:
  # gmhn -> gh h(m)+n -> x t+h(m) + n where x is the reduced
  #                                   word under collection and t is the
  #                                   "tail"
  # so gamma(g, h) = t
  # given gamma need the tails:
  # need to implement the group operation for the extension
  # (g, v)(h, u) -> (gh, v^h + u + gamma(g, h))
  # then the rules with tails need to be evaluated at
  # the group generators (g_i, 0)
  # r -> s gives a relation r s^-1 which should evaluate, using gamma
  # to (0, t) where t is the tail for this rule
end

"""
For a vector of maps
  f_i : M -> N_i
compute the map
  f : M -> prod N_i : m -> (f_i(m))_i
"""
function Oscar.direct_sum(a::Vector{<:Union{<:Generic.ModuleHomomorphism{<:RingElement}, FinGenAbGroupHom}})
  D = direct_product([codomain(x) for x = a]...; task = :none)
  return hom(domain(a[1]), D, hcat([matrix(x) for x = a]...))
end

function Base.sum(a::Vector{FqMatrix})
  c = deepcopy(a[1])
  for i=2:length(a)
    add!(c, c, a[i])
  end
  return c
end

function Base.iszero(x::AbstractAlgebra.Generic.ModuleHomomorphism)
  return iszero(x.matrix)
end

function istwo_cocycle(c::CoChain{2})
  C = c.C
  G = C.G
  for g = G
    for h = G
      for k = G
        #= if (g*h)(x) = h(g(x)), then the cocycle should be
             X[(g*h, k)] X[(g, h)] == mA(g)(X[(h, k)]) X[(g, hk)]
           if (g*h)(x) = h(g(x)) then we should get
             X[(g, hk)] X[(h, k)]  == mA(k)(X[(g, h)]) X[(gh, k)]

             (Debeerst, PhD, (1.1) & (1.2))

             However, if we mix the conventions, all bets are off...
        =#
        a = c(g, h*k) + c(h, k) - action(C, k, c(g, h))- c(g*h, k)
#        @show a, iszero(a) || valuation(a)
iszero(a) || (@show g, h, k, a ; return false)
        @assert iszero(a) # only for local stuff...|| valuation(a) > 20
      end
    end
  end
  return true
end

# create a free module with type "compatible" with that of `M`
_similar_free_module(M::FinGenAbGroup, n::Int) = free_abelian_group(n)
_similar_free_module(M::AbstractAlgebra.FPModule, n::Int) = free_module(base_ring(M), n)
_similar_free_module(M::Oscar.ModuleFP, n::Int) = FreeMod(base_ring(M), n)

"""
Compute
  0 -> C -I-> hom(Z[G], C) -q-> B -> 0
To allow "dimension shifting": H^(n+1)(G, C) - H^n(G, q)
returns (I, q), (hom(Z[G], C), B)
"""
function dimension_shift(C::GModule)
  zg, _, _ = regular_gmodule(C)
  Z = _similar_free_module(zg.M, 0)
  @assert is_consistent(zg)
  H, mH = Oscar.GModuleFromGap.ghom(zg, C)
  @assert is_consistent(H)
  #from https://math.mit.edu/classes/18.786/2018/LectureNotes29.pdf
  #around 29.8
  #Drew's notes.
  #the augmentation map on the (canonical) generators is 1
  inj = hom(C.M, H.M, [preimage(mH, hom(zg.M, C.M, [c for g in 1:ngens(zg.M)])) for c in gens(C.M)])
  @assert is_G_hom(C, H, inj)
  B, q = quo(H, image(inj)[2])

  return (inj, q), (H, B)

  #would be nice, but looses the G-operation.
  #XXX: we don't have homs for GModules
  #   : sice we also don't have elements
  #   : do we need elements for homs?
  Z1 = hom(Z, C.M, elem_type(C.M)[zero(C.M) for i in 1:ngens(Z)])
  Z2 = hom(q.M, Z, [zero(Z) for x = gens(q.M)])
  return cochain_complex([Z1, inj, mq, Z2])
end

function fixed_module(H::GModule)
  K = H.M
  id = hom(K, K, gens(K))
  mK = hom(K, K, gens(K))
  for x = H.ac
    k = intersect(K, kernel(x-id)[1])
    fl, mk = is_sub_with_data(k, K)
    K = k
    mK = mk*mK
  end
  return gmodule(H.G, [FinGenAbGroupHom(mK*g*pseudo_inv(mK)) for g = H.ac]), mK
end

function dimension_shift_left(C::GModule)
  G = C.G
  if isa(C.M, FinGenAbGroup)
    zg, ac, em = regular_gmodule(FinGenAbGroup, G, ZZ)
    Z = Hecke.zero_obj(zg.M)
  elseif isa(C.M, AbstractAlgebra.FPModule{<:FieldElem})
    zg, ac, em = regular_gmodule(G, base_ring(C))
    Z = free_module(base_ring(C), 0)
  else
    error("unsupported module")
  end
  @assert is_consistent(zg)
  H, mH = tensor_product(zg, C)
  @assert is_consistent(H)
  #XXX broken from here onwards
  pro = hom(H.M, C.M, [pseudo_inv(mH)(x)[2] for x = gens(H.M)])
  @assert is_G_hom(H, C, pro)
  k, mk = kernel(pro)
  if isa(C.M, FinGenAbGroup)
    K = gmodule(C.G, [FinGenAbGroupHom(mk*g*pseudo_inv(mk)) for g = H.ac])
  else
    K = gmodule(C.G, [mk*g*pseudo_inv(mk) for g = H.ac])
  end

  return (mk, pro), (K, H)

  #would be nice, but looses the G-operation.
  #XXX: we don't have homs for GModules
  #   : sice we also don't have elements
  #   : do we need elements for homs?
  Z1 = hom(Z, C.M, elem_type(C.M)[zero(C.M) for i = gens(Z)])
  Z2 = hom(q.M, Z, [zero(Z) for x = gens(q.M)])
  return cochain_complex([Z1, inj, mq, Z2])
end


"""
Compute H^3 via dimension-shifting:
There is a short exact sequence
  1 -> A -> Hom(Z[G], A) -> B -> 1
thus
  H^3(G, A) = H^2(G, B)
as Hom(Z[G], A) is induced hence has trivial cohomology.
Currently only the group is returned
"""
function H_three(C::GModule{<:Oscar.GAPGroup, <:Any})
  (inj, mq), (H, q) = dimension_shift(C)

#  return q, mq, inj, H
  #possibly, to get 3-chains:
  # 2 chain in q
  # preimage mq 2 chain in H
  # differential 3 chain in H
  # preimage inj 3 chain in C
  H, mH, _ = H_two(q)
  function chain(a)
    @assert parent(a) == H
    c = mH(a)
    d = map_entries(pseudo_inv(mq), c, parent = H)
    e = differential(d)
    f = map_entries(pseudo_inv(inj), e, parent = C)
  end
  #preimage under the differential is hard...
  return H, chain
end

function is_right_G_module(C::GModule)
  #tests if the action is right-linear
  G = C.G
  return all(action(C, g)*action(C, h) == action(C, g*h) for g in gens(G), h in gens(G))
end

function is_left_G_module(C::GModule)
  #tests if the action is left-linear
  G = C.G
  return all(action(C, h)*action(C, g) == action(C, g*h) for g = gens(G) for h = gens(G))
end

"""
For a gmodule `C` compute the `i`-th cohomology group
where `i` can be `0`, `1` or `2`. (or `3` ...)
Together with the abstract module, a map is provided that will
produce explicit cochains.
"""
function cohomology_group(C::GModule, i::Int; Tate::Bool = false)
  #should also allow modules...
  if Tate
    @assert is_finite(group(C))
  end
  if i==0
    if Tate
      return H_zero_tate(C)
    else
      return H_zero(C)
    end
  elseif i==1
    return H_one(C)
  elseif i==2
    return H_two(C)
  elseif i==3
    return H_three(C)
  end
  error("only H^0, H^1 and H^2 are supported")
end

# better use `isomorphism` directly
function fp_group_with_isomorphism(M::AbstractAlgebra.FPModule{<:FinFieldElem})
  iso = isomorphism(FPGroup, M, on_gens=true)
  return codomain(iso), iso
end


#########################################################
function Oscar.matrix(M::FreeModuleHom{FreeMod{QQAbFieldElem}, FreeMod{QQAbFieldElem}})
  return M.matrix
end

###########################################################

#=
function get_collector(G::GapObj)
  @show G
  return GAP.evalstr("x -> FamilyObj(x.1)!.rewritingSystem")(G)
end
=#

"""
Compute an isomorphic pc-group `G` (and the isomorphism from `M` to `G`).
If `refine` is true,
the pc-generators will all have prime relative order, thus the
group should be safe to use.
If `refine` is false, then the relative orders are just used from the hnf
of the relation matrix.
"""
function pc_group_with_isomorphism(M::FinGenAbGroup; refine::Bool = true)
  @assert is_finite(M)
  R = rels(M)
  h = hnf(R)
  if nrows(h) > ncols(h)
    h = view(h, 1:ncols(h), :)
  end
  @assert nrows(h) == ncols(h)
  if refine
    hm = elem_type(M)[]
    for i=1:nrows(h)
      lf = collect(factor(h[i,i]).fac)
      for (p,k) = lf
        v = divexact(h[i,i], p^k)*M[i]
        for j=1:k-1
          push!(hm, v)
          v *= p
        end
        push!(hm, v)
      end
    end
    M, mM = sub(M, hm) #without simplify is guaranteed to keep the gens!
  else
    mM = hom(M, M, gens(M))
  end

  G = free_group(ngens(M))
  h = rels(M)
  @assert !any(x->h[x,x] == 1, 1:ncols(h))

  C = collector(nrows(h))
  set_relative_orders!(C, [h[i,i] for i in 1:nrows(h)])
  for i in 1:ngens(M)
    r = Pair{Int, ZZRingElem}[]
    for j in i+1:ngens(M)
      push!(r, j => -h[i, j])
    end
    set_power!(C, i, r)
  end
  B = pc_group(C)
  FB = GAP.Globals.FamilyObj(GAP.Globals.Identity(GapObj(B)))

  function Julia_to_gap(a::FinGenAbGroupElem)
    r = ZZRingElem[]
    for i=1:ngens(M)
      if !iszero(a[i])
        push!(r, i)
        push!(r, a[i])
      end
    end
    return GAP.Globals.ObjByExtRep(FB, GAP.Obj(r; recursive = true))
  end

  function Gap_to_julia(a::GapObj)
    e = GAPWrap.ExtRepOfObj(a)
    z = zeros(ZZRingElem, ngens(M))
    for i=1:2:length(e)
      if !iszero(e[i+1])
        z[e[i]] = e[i+1]
      end
    end
    return M(z)
  end

#  @assert is_isomorphic(B, fp_group(M))

  return B, MapFromFunc(
    codomain(mM), B,
    y->PcGroupElem(B, Julia_to_gap(preimage(mM, y))),
    x->image(mM, Gap_to_julia(GapObj(x))))
end

# `refine` is irrelevant because `M` is elementary abelian.
function pc_group_with_isomorphism(M::AbstractAlgebra.FPModule{<:FinFieldElem}; refine::Bool = true)
  iso = isomorphism(PcGroup, M, on_gens=true)
  return codomain(iso), iso
end


"""
Given a 2-cocycle, return the corresponding group extension, ie. the large
group, the injection of the abelian group and the quotient as well as a map
that given a tuple of elements in the group and the abelian group returns
the corresponding elt in the extension.

If the gmodule is defined via a pc-group and the 1st argument is the
`Type{PcGroup}`, the resulting group is also pc.
"""
function extension(::Type{FPGroup}, c::CoChain{2,<:Oscar.GAPGroupElem})
  C = c.C
  G = Group(C)
  F = codomain(isomorphism(FPGroup, G, on_gens=true))
  M = Module(C)
  ac = action(C)
  iac = inv_action(C)
  mfM = inv(isomorphism(FPGroup, M))
  fM = domain(mfM)

  N = free_group(ngens(G) + ngens(fM))
  function fMtoN(g)
    return reduce(*, [gen(N, ngens(G)+abs(w))^sign(w) for w = word(g)], init = one(N))
  end
  #TODO: this "loop" has been implemented several times....
  s = map(fMtoN, relators(fM))
  for R = relators(F)
    t = zero(M)
    g = one(G)
    r = one(N)
    for w = word(R)
      if w > 0
        t = ac[w](t) + c(g, gen(G, w))
        g = g*gen(G, w)
        r = r*gen(N, w)
      else
        t = iac[-w](t) + c(g, inv(gen(G, -w))) - c(gen(G, -w), inv(gen(G, -w)))
        g = g*inv(gen(G, -w))
        r = r*inv(gen(N, -w))
      end
    end
    push!(s, r*inv(fMtoN(preimage(mfM, t))))
  end
  for i=1:ngens(G)
    for j=1:ngens(fM)
      #m[j]*g[i] = g[i] m[j]^g[i]
      t = preimage(mfM, ac[i](gen(M, j)))
      push!(s, gen(N, ngens(G)+j)*gen(N, i)*inv(fMtoN(t)) * inv(gen(N, i)))
    end
  end
  Q, mQ = quo(N, s)
  @assert ngens(Q) == ngens(N)
  MtoQ = hom(fM, Q, gens(fM), gens(Q)[ngens(G)+1:end]; check = false)
  QtoG = hom(Q, G, gens(Q), vcat(gens(G), [one(G) for i=1:ngens(fM)]); check = false)
  @assert domain(mfM) ==fM
  @assert codomain(mfM) == M

  function GMtoQ(g::GAPGroupElem, m)
    @assert parent(m) == M
    @assert parent(g) == G
    h1 = hom(free_group(G), N, gens(free_group(G)), [N[i] for i=1:ngens(G)]; check = false)
    h2 = hom(free_group(fM), N, gens(free_group(fM)), [N[i+ngens(G)] for i=1:ngens(fM)]; check = false)
    return mQ(h1(underlying_word(g))*h2(underlying_word(preimage(mfM, m))))
  end

  return Q, pseudo_inv(mfM)*MtoQ, QtoG, GMtoQ
end

function extension(::Type{PcGroup}, c::CoChain{2,<:Oscar.PcGroupElem})
  C = c.C
  G = Group(C)
  @assert isa(G, PcGroup)
  M = Module(C)
  ac = action(C)
  iac = inv_action(C)
  fM, mfM = pc_group_with_isomorphism(M)

  N = free_group(ngens(G) + ngens(fM))
  Gp = GAP.Globals.Pcgs(GapObj(G))
  @assert length(Gp) == ngens(G)
#  @assert all(x->Gp[x] == GapObj(gen(G, x)), 1:ngens(G))
  Go = GAP.Globals.RelativeOrders(Gp)

  Mp = GAP.Globals.Pcgs(GapObj(fM))
  @assert length(Mp) == ngens(fM) == ngens(M)
#  @assert all(x->Mp[x] == GapObj(gen(fM, x)), 1:ngens(M))
  #problem/ TODO: Z/100Z has a useful GAP-pc-group has 4 gens (of
  #order 2, 2, 5, 5
  #so need to switch GAP to the other Pc-Groups and/or drop this 
  #assert
  Mo = GAP.Globals.RelativeOrders(Mp)

  CN = GAP.Globals.SingleCollector(GapObj(N), GAP.Globals.Concatenation(Go, Mo))
  FN = GAP.Globals.FamilyObj(GapObj(N[1]))

  for i=1:ngens(fM)
    lp = deepcopy(GAPWrap.ExtRepOfObj(Mp[i]^Mo[i]))
    for k=1:2:length(lp)
      lp[k] += ngens(G)
    end
    m = GAP.Globals.ObjByExtRep(FN, lp)
    GAP.Globals.SetPower(CN, i+ngens(G), m)
    for j=i+1:ngens(fM)
      p = Mp[j]^Mp[i]
      @assert p == Mp[j]
      lp = deepcopy(GAPWrap.ExtRepOfObj(p))
      for k=1:2:length(lp)
        lp[k] += ngens(G)
      end
      GAP.Globals.SetConjugate(CN, j+ngens(G), i+ngens(G), GAP.Globals.ObjByExtRep(FN, lp))
    end
  end

  fMtoN = function(x)
    lp = deepcopy(GAPWrap.ExtRepOfObj(GapObj(x)))
    for k=1:2:length(lp)
      @assert lp[k] > 0
      lp[k] += ngens(G)
    end
    return GAP.Globals.ObjByExtRep(FN, lp)
  end

  word = function(y)
    z = GAPWrap.UnderlyingElement(y)
    return map(Int, GAP.Globals.LetterRepAssocWord(z))
  end

  #for W = (w1, ... w_n) compute ((w1, 0), ..., (wn, 0))
  #and return the tail only.
  word_to_elem = function(W)
    t = zero(M)
    g = one(G)
    r = one(N)
    for w = W
      if w > 0
        t = ac[w](t) + c(g, gen(G, w))
        g = g*gen(G, w)
        r = r*gen(N, w)
      else
        t = iac[-w](t) + c(g, inv(gen(G, -w))) - c(gen(G, -w), inv(gen(G, -w)))
        g = g*inv(gen(G, -w))
        r = r*inv(gen(N, -w))
      end
    end
    return t
    return fMtoN(mfM(t))
  end

  #to lift the pc-relations:
  # F^p = w (order relation)
  #  compute (F, 0)^p = (?, t) = (?, 0)(1, t)
  #  compute (w, 0)   = (?, s) = (?, 0)(1, s)
  #  so (?, 0) = (w, 0)(1,s)^-1= (w, 0)(1,-s) if chain is normalized
  #  thus (F, 0)^p = (?, 0)(1, t) = (w, 0)(1,-s)(1, t)
  #  the ? should be identical, namely the collected version of w
  #  then (F, 0)^p = (w, t-s) might be the answer
  # F^G = w (conjugate relation): same
  #  (F, 0)^(G, 0) = (?, t) = (?, 0)(1, t)
  #  (w, 0)        = (?, s) = (?, 0)(1, s)
  #  thus (F, 0)^(G, 0) = (w, t-s)
  for i=1:ngens(G)
    p = Gp[i]^Go[i]
    pp = GAP.Globals.ObjByExtRep(FN, GAPWrap.ExtRepOfObj(p))
    m = fMtoN(mfM(word_to_elem([i for k=1:Go[i]])-word_to_elem(word(p))))
    GAP.Globals.SetPower(CN, i, pp*m)
    for j=i+1:ngens(G)
      p = Gp[j]^Gp[i]
      m = fMtoN(mfM(word_to_elem([-i, j, i])-word_to_elem(word(p))))
      pp = GAP.Globals.ObjByExtRep(FN, GAPWrap.ExtRepOfObj(p))
      GAP.Globals.SetConjugate(CN, j, i, pp*m)
    end
    for j=1:ngens(fM)
      m = fMtoN(mfM(action(C, gen(G, i), preimage(mfM, gen(fM, j)))))
      GAP.Globals.SetConjugate(CN, j+ngens(G), i, m)
    end
  end

#  l = GAP.Obj([])
#  GAP.Globals.FinitePolycyclicCollector_IsConfluent(CN, l)
#  @show l

#  z = GAP.Globals.GroupByRwsNC(CN)
#  s = GAP.Globals.GapInputPcGroup(z, GAP.Obj("Z"))
  Q = PcGroup(GAP.Globals.GroupByRws(CN))
  fQ = GAP.Globals.FamilyObj(GapObj(one(Q)))
  mQ = hom(N, Q, gens(N), gens(Q); check = false)

  @assert ngens(Q) == ngens(N)
  MtoQ = hom(fM, Q, gens(fM), gens(Q)[ngens(G)+1:end]; check = false)
  QtoG = hom(Q, G, gens(Q), vcat(gens(G), [one(G) for i=1:ngens(fM)]); check = false)
  @assert domain(mfM) == M
  @assert codomain(mfM) == fM
#  @assert is_surjective(QtoG)
#  @assert is_injective(MtoQ)

  mfG = epimorphism_from_free_group(G)
  mffM = epimorphism_from_free_group(fM)

  function GMtoQ(wg, m)
    wm = syllables(preimage(mffM, mfM(m)))
    for (ind, exp) in wm
      push!(wg, ind+ngens(G))
      push!(wg, exp)
    end
    return mQ(FPGroupElem(N, GAP.Globals.ObjByExtRep(FN, GAP.Obj(wg))))
  end

  return Q, mfM*MtoQ, QtoG, GMtoQ
end

"""
    extension_with_abelian_kernel(X::Oscar.GAPGroup, M::Oscar.GAPGroup)

For a group `K` and an abelian normal subgroup `M`, construct `M` as a
`X/M`-modul amd find the 2-cochain representing `X`.
"""
function extension_with_abelian_kernel(X::Oscar.GAPGroup, M::Oscar.GAPGroup)
  g, pi = quo(X, M)
  fl, i = is_subgroup(M, X)
  phi = isomorphism(FinGenAbGroup, M)
  A = codomain(phi)
  gg = isomorphism(PermGroup, g)
  C = gmodule(codomain(gg), [hom(A, A, [phi(preimage(i, i(preimage(phi, a))^preimage(pi, preimage(gg, x)))) for a = gens(A)]) for x = gens(codomain(gg))])

  H2, mH2, _ = cohomology_group(C, 2)

  c = Dict( (gg(a), gg(b)) => phi(preimage(i, preimage(pi, a)*preimage(pi, b)*inv(preimage(pi, a*b)))) for a = g for b = g)
  return C, CoChain{2, PermGroupElem, FinGenAbGroupElem}(C, c)
end

function Oscar.automorphism_group(F::AbstractAlgebra.Generic.FreeModule{<:FinFieldElem})
  G = GL(dim(F), base_ring(F))
  set_attribute!(G, :aut_group=>F)
  return G, MapFromFunc(G, Hecke.MapParent(F, F, "homomorphisms"),
                         x->hom(F, F, matrix(x)),
                         y->G(matrix(y)))
end

function (G::MatrixGroup{T})(h::AbstractAlgebra.Generic.ModuleHomomorphism{T}) where T
  return G(matrix(h))
end

function (G::MatrixGroupElem{T})(h::AbstractAlgebra.FPModuleElem{T}) where T
  return h*G
end

function Oscar.hom(g::MatrixGroupElem)
  G = parent(g)
  p = get_attribute(G, :aut_group)
  p === nothing && error("Matrix group must be the automorphism group of some module")
  return hom(p, p, matrix(g))
end

"""
Let C be a G-module with G action on M. Then this function find the
subgroup of 'Aut(M) x Aut(G)' that is compatible with the G-module
structure:

Let 'h: G to Aut(M)'  from C, then '(alpha, gamma) in Aut(M) x Aut(G)'
are compatible iff
  for all 'a in M' and 'g in G' we have
  'h(gamma(g))(a) == alpha(inv(h(g))(alpha^-1(a)))'

The group and the action on 2-cochains is returned. Cochains in the
same orbit parametrize the same group.
"""
function compatible_pairs(C::GModule)
  #not exported
  G = C.G
  M = C.M

  autG = automorphism_group(G)
  autM = automorphism_group(M)
  if isa(autM, Tuple)
    autM = autM[1]
  end

  D, emb, pro = direct_product(autM, autG, morphisms = true)

  function action_on_chain(g::GAPGroupElem, c::CoChain{2, S, T}) where {S, T}
    al = pro[1](func(g)::elem_type(D))::elem_type(autM)
    ga = pro[2](func(g)::elem_type(D))::elem_type(autG)
    if isdefined(c, :D)
      d = Dict{Tuple{S, S}, T}(ab=> al(c((inv(ga)(ab[1]), inv(ga)(ab[2])))) for ab = keys(c.d))
      return CoChain{2, S, T}(C, d, ab->al(c((inv(ga)(ab[1]), inv(ga)(ab[2])))))
    else
      d = Dict{Tuple{S, S}, T}(ab=> al(c((inv(ga)(ab[1]), inv(ga)(ab[2])))) for ab = keys(c))
      return CoChain{2, S, T}(C, d)
    end
  end

  h = hom(G, autM, [autM(x) for x = C.ac])
  im = image(h)[1]
  ke = kernel(h)[1]

  if order(im) == 1
    func = Base.identity
    return D, action_on_chain
  end

  N, mN = normalizer(autM, im)
  S, mS = stabilizer(autG, ke, (x,y)->image(hom(y), x)[1])

  NS, em, pr = direct_product(N, S, morphisms = true)
  #from Holt/ Eick, ... Handbook of Computational Group Theory, P319
  S, func = stabilizer(NS, h, (x, y) -> hom(G, autM,
            [inv(pr[1](y))*x(inv(pr[2](y))(g))*pr[1](y) for g = gens(G)]))
  #the more direct naive (and slow) approach...
  #C = sub(D, [t for t in preimage(pro[1], N)[1] if all(ag -> action(C, pro[2](t)(ag[2]), ag[1]) == pro[1](t)(action(C, ag[2], inv(pro[1](t))(ag[1]))), Iterators.product(gens(M), gens(G)))])

  return S, action_on_chain
end

function split_extension(C::GModule)
  #bypasses the Cohomolgy computation, hopefully
  c = Dict((g, h) => zero(C.M) for g = C.G for h = C.G)
  S = elem_type(C.G)
  T = elem_type(C.M)
  return extension(FPGroup, CoChain{2, S, T}(C, c))
end

function split_extension(::Type{PcGroup}, C::GModule{<:PcGroupElem})
  #bypasses the Cohomolgy computation, hopefully
  c = Dict((g, h) => zero(C.M) for g = C.G for h = C.G)
  S = elem_type(C.G)
  T = elem_type(C.M)
  return extension(PcGroup, CoChain{2, S, T}(C, c))
end


function all_extensions(C::GModule)
  @assert isfinite(C.M)
  if gcd(order(C.M), order(C.G)) == 1
    return [split_extension(C)]
  end
  H2, mH2, _ = H_two(C, lazy = true)
  if order(H2) == 1
    return [extension(mH2(zero(H2)))]
  end
  T, mT = compatible_pairs(C)
  G = gset(T, (a, g) -> preimage(mH2, mT(g, mH2(a))), collect(H2), closed = true)
  O = orbits(G)
  all_G = []
  for o = O
    E, NtoE, EtoG, EGtoM = extension(mH2(representative(o)))
    push!(all_G, E)
  end
  return all_G
end

function Oscar.image(M::Map{FinGenAbGroup, <:Oscar.GAPGroup})
  A = domain(M)
  G = codomain(M)
  s, ms = sub(G, map(M, gens(A)))
  return s, ms
end

@doc """
    is_central(NtoE::Map) -> Bool

Tests if the domain is in the center of the codomain.    
"""
function Oscar.is_central(NtoE::Map{<:Union{Group, FinGenAbGroup}, <:Group})
  E = codomain(NtoE)
  n = map(NtoE, gens(domain(NtoE)))
  return all(x->all(y->y*x == x*y, n), gens(E))
  return is_subset(N, center(E)[1])
end

@doc """
    is_stem_extension(NtoE::Map) -> Bool

Tests if the domain is in the center and the derived subgroup of the codomain.
"""
function is_stem_extension(NtoE::Map{<:Union{Group, FinGenAbGroup}, <:Group}; is_central_known::Bool = false)
  E = codomain(NtoE)
  N = image(NtoE)[1]
  E = codomain(NtoE)
  if !is_central_known
    is_central(NtoE) || return false
  end
  return is_subset(N, derived_subgroup(E)[1])
end
#= from Max: thinking about 2-co-chains and group extensions
 N -> E -> G via sigma in H^2(G, N)
 central iff operation of G on N is trivial
 stem:
 E' is generated by commutators [a, b] for a, b in E
 [(g, a), (h, b)] = (g, a)^-1 (h, b)^-1 (g, a) (h, b)
    = (g^-1, - a^(g^-1) - sigma(g, g^-1)) *
      (h^-1, - b^(h^-1) - sigma(h, h^-1)) *
      (gh, a^h + b + sigma(g, h))
    = (g^-1, -a - sigma(g, g^-1)) * (h^-1, -b - sigma(h, h^-1)) *
      (gh, a+b+sigma(g, h)) #central: a^g = a!
    = ([g, h], -a -sigma(g, g^-1) -b - sigma(h, h^-1) + a+b+sigma(g, h))
    = ([g, h], sigma(g, h) - sigma(g, g^-1) - sigma(h, h^-1)) 

    This shows that the 2.nd component of a commutator in E does not
    depend on the values in N, only on the 1st component (and the 2-co-chain)

  (g, a)^(h, b) = (h, b)^-1 * (g, a) * (h, b) 
  = (h^-1, -b - sigma(h, h^-1)) (gh, a+b+sigma(g, h)) 
  = (g^h, a + sigma(g, h) - sigma(h, h^-1))

  We want to test if N in E`, N in E looks like (1, n), so we
  need words in G that are in G' and evaluate to 1, so
  relators for G`

  Let G = <gi|i>
  Let psi: F -> G xi -> gi the map from a free group
  Let phi: G' -> G the embedding of the derived subgroup

  Let r be a relator for G` as an element in G, then a preimage in F
  will write it as a word in the gi
  Instead of decomposing this as a word in [gi,gj] we just supplement the
  gi to (gi, a) for any a (e.g.). In any commutator relator evaluation
  the choice of a is irrelevant, it will cancel.

  So: N in E' iff <pi(rj((gi, 0)_i)) | j> = M where 
    pi:E -> N: (1, n) -> n is the preimage under the injection
=#    

    
(G::MatrixGroup{FqFieldElem, FqMatrix})(a::GAP.GapObj) = Oscar.group_element(G, a)

@doc raw"""
    gmodule_class_reps(M::Union{<:AbstractAlgebra.FPModule, FinGenAbGroup}, G::Oscar.GAPGroup) -> Vector{GModule}

Find representaives for all possible operations of `G` on `M`.
"""
function gmodule_class_reps(M::Union{<:AbstractAlgebra.FPModule, FinGenAbGroup}, G::Oscar.GAPGroup) #the cohomology wants it
  A = automorphism_group(M)
  if isa(A, Tuple)
    A = A[1]
  end
  l = GAP.Globals.AllHomomorphismClasses(GapObj(G), GapObj(A))
  all_G = []
  for i = l
    C = gmodule(G, [hom(A(i(GapObj(g)))) for g = gens(G)])
    push!(all_G, C)
  end
  return all_G
end

function all_extensions(M::Union{<:AbstractAlgebra.FPModule, FinGenAbGroup}, G::Oscar.GAPGroup) #the cohomology wants it
  all_G = []
  for M = gmodule_class_reps(M, G)
    append!(all_G, all_extensions(M))
  end
  return all_G
end

function fp_group(c::CoChain{2})
  return extension(c)[1]
end

function pc_group(c::CoChain{2, <:Oscar.PcGroupElem})
  return extension(PcGroup, c)[1]
end

function Oscar.permutation_group(c::CoChain{2})
  g = extension(FPGroup, c)[1]
  return permutation_group(g)
end

end #module

using .GrpCoh

export gmodule, fp_group, pc_group, induce, cohomology_group, extension
export permutation_group, is_consistent, istwo_cocycle, GModule
export split_extension, all_extensions, extension_with_abelian_kernel
export is_stem_extension, is_central
export alternating_square, symmetric_square
