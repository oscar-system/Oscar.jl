###############################################################################
#
#  Helpers
#
###############################################################################

### Gluing factory

# Return the two torsion quadratic modules from the initial problem
ambient_modules(Fac::ZZLatGluingFactory) = Fac.ambient_modules

# Return the groups used to do the computations of local gluings
local_classifying_groups(Fac::ZZLatGluingFactor) = Fac.local_classifying_groups

# Return the modules computed from the initial conditions, where we have to
# for anti-isometric glue groups
conditions_modules(Fac::ZZLatGluingFactory) = Fac.conditions_modules

# Return the signature of a primitive extension as setup in the inital problem
overlattice_signature(Fac::ZZLatGluingFactory) = Fac.sign

# Return the parity type of a primitive extension as setup in the initial
# problem. It can be :even for even overlattices, :odd for odd lattices and
# :both for any integral overlattices
overlattice_parity(Fac::ZZLatGluingFactory) = Fac.par

# Return the list of possible order for the glue groups as setup by the
# conditions of the initial problem
possible_glue_order(Fac::ZZLatGluingFactory) = Fac.glue_order

# Test whether a given lattice has wanted genus, as decided by the initial
# conditions of the problem
function test_overlattice(Fac::ZZLatGluingFactory, L::ZZLat)
  if isdefined(Fac, :genus_over)
    genus(L) in Fac.genus_over || return false
  end
  return true
end

### Gluing

function gluing_data(x::ZZLatGluing)
  return x.glue_map, x.inv_glue_map, x.glue_group_left, x.stabilizer_left, x.glue_group_right, x.stabilizer_right
end

# The type ZZLatGluing is symmetrical, we just reverse the order of the fields
# to see it from right to left
function Oscar.inv(x::ZZLatGluing)
  phi, iphi, i1, j1, i2, j2 = gluing_data(x)
  return ZZLatGluing(iphi, phi, i2, j2, i1, j1)
end

function glue_map(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.glue_map
  elseif i == 2
    return x.inv_glue_map
  end
end

function glue_group(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.glue_group_left
  elseif i == 2
    return x.glue_group_right
  end
end

function stabilizer_glue_group(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.stabilizer_left
  elseif i == 2
    return x.stabilizer_right
  end
end

### Gluing ambient

function _gluing_ambient(
    q1::TorQuadModule,
    q2::TorQuadModule,
    par::Symbol,
  )
  as_bilinear_module = par === :even ? false : true
  return ZZLatGluingAmbient(_direct_sum_with_embeddings_orthogonal_groups(q1, q2; as_bilinear_module)
end


### More on genera

# Return all genera with given form, signature and parity conditions
# There are at most 2 of them
# par is either :even, or :odd, or :both (actually anything :odd or :even
# would do here)
function _integer_genera(
  q::TorQuadModule,
  sign::Tuple{Int, Int},
  par::Symbol,
)
  p, n = sign
  GKs = Set{ZZGenus}()
  if par != :odd
    ok, _G = is_genus_with_genus(q, (p, n); parity=2)
    ok && push!(GKs, _G)
  end

  if par != :even
    ok, _G = is_genus_with_genus(q, (p, n); parity=1)
    ok && push!(GKs, _G)
  end
  return GKs
end

function __representatives(g::ZZGenus)
  return get_attribute(g, :representatives, representatives(g))
end

### Handle parity

@doc raw"""
    orthogonal_group_bilinear(T::TorQuadModule) -> AutomorphismGroup{TorQuadModule}

Return the orthogonal group for the bilinear form on ``T``, i.e. seeing ``T``
as a torsion bilinear form. If the modulus of the quadratic form on ``T`` is
equal to the modulus of its bilinear form, then this is the same as calling
`orthogonal_group(T)`.
"""
@attr AutomorphismGroup{TorQuadModule} function orthogonal_group_bilinear(T::TorQuadModule)
  return __orthogonal_group(q; as_bilinear_module=true)
end

function __orthogonal_group(
  q::TorQuadModule;
  as_bilinear_module::Bool=false,
)
  if !as_bilinear_module || modulus_bilinear_form(q) == modulus_quadratic_form(q)
    return orthogonal_group(q)
  end
  qb = _as_finite_bilinear_module(q)
  return _orthogonal_group(q, orthogonal_group_gens(qb); check=false)
end

### To be moved to Hecke eventually

# Return the elementary divisors of the largest finite abelian group
# which embed in both G1 and G2
function _maximal_common_subgroup_snf(
  G1::FinGenAbGroup,
  G2::FinGenAbGroup,
)
  @assert is_finite(G1) && is_finite(G2)
  s1 = elementary_divisors(G1)
  s2 = elementary_divisors(G2)
  s = Array{ZZRingElem}(undef, min(length(s1), length(s2)))
  for i in 0:length(s)-1
    s[end-i] = gcd(s1[end-i], s2[end-i])
  end
  return s
end

# Given the elementary divisors elG of a finite abelian group G, return the
# set of all possible lists of valuations for the elementary divisors of a
# subgroup of G whose p-Sylow subgroup has order p^v. The lists in output
# contain only positive valuations and they are sorted in decreasing order.
function _psubgroups_types(
  elG::Vector{ZZRingElem},
  p::ZZRingElem,
  v::Int,
)
  @assert v > 0
  testv = isequal(v)∘sum
  p_subtypes = Set{Vector{Int}}()
  if isempty(elG)
    return p_subtypes
  end
  Gtype = Int[valuation(a, p) for a in elG]
  reverse!(Gtype)
  filter!(!=(0), Gtype)
  types = _subpartitions(Gtype)
  for t in types
    testv(t) || continue
    filter!(!=(0), t)
    push!(p_subtypes, t)
  end
  return p_subtypes
end

@doc raw"""
    torsion_subgroup(
      G::FinGenAbGroup,
      n::IntegerUnion,
      add_to_lattice:Bool = true,
    )

Return the subgroup of ``G`` consisting of elements of order dividing ``n``.
"""
function torsion_subgroup(
  G::FinGenAbGroup,
  n::IntegerUnion,
  add_to_lattice::Bool = true,
)
  f = FinGenAbGroupHom(G, G, ZZ(n)*identity_matrix(ZZ, ngens(G)))
  return kernel(f, add_to_lattice)
end

@doc raw"""
    torsion_subgroup(
      T::TorQuadModule,
      n::IntegerUnion,
    )

Return the subgroup of ``T`` consisting of elements of order dividing ``n``.
"""
function torsion_subgroup(
  T::TorQuadModule,
  n::IntegerUnion,
)
  Sab, j = torsion_subgroup(abelian_group(T), n)
  return sub(T, TorQuadModuleElem[T(i(s)) for s in gens(Sab)])
end


#function elementary_divisors(s::ZZLocalGenus)
#  p = prime(s)
#  els = Array{ZZRingElem}(undef, rank(s))
#  i = Int(1)
#  for sym in symbol(s)
#    for _ in 1:sym[2]
#      els[i] = p^first(sym)
#      i += 1
#    end
#  end
#  return els
#end
#
#function elementary_divisors(G::ZZGenus)
#  lgs = local_symbols(G)
#  elG = elementary_divisors(first(lgs))
#  for i in 2:length(lgs)
#    mul!.(elG, elementary_divisors(lgs[i]))
#  end
#  filter!(!=(1), elG)
#  return elG
#end
#
#function _ptype_from_elementary_divisors(
#  v::Vector{ZZRingElem},
#  p::ZZRingElem,
#)
#  s = reverse!(Int[valuation(a, p) for a in v])
#  filter!(!=(0), s)
#  return s
#end

###############################################################################
#
#  Init functions
#
###############################################################################

# Initialize the gluing factory with the conditions from the problem
# By construction `Fac` should already know about the ambient modules, the
# signature and the wanted parity for the primitive extensions. If the local
# classifying groups have not been set, the function initialize them to be the
# orthogonal groups (maybe as finite bilinear modules, depending on the value
# of Fac.par) of the ambient modules.
function init_gluing_factory!(
  Fac::ZZLatGluingFactory;
  left_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  right_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_exponent::Union{nothing, IntegerUnion}=nothing,
  glue_order::AbstractVector{T}=Int[],
  glue_elementary_divisors::Vector{Vector{ZZRingElem}}=Vector{ZZRingElem}[],
  form_over::Vector{TorQuadModule}=TorQuadModule[],
  genus_over::Vector{ZZGenus}=ZZGenus[],
) where T <: IntegerUnion
  init_gluing_conditions!(
    Fac,
    left_glue_annihilator,
    right_glue_annihilator,
    glue_exponent,
    glue_order,
    glue_elementary_divisors,
    form_over,
  )
  assert_has_local_classifying_groups!(Fac)
end

# Using the glue annihilators and glue_exponent, the function sets up the
# so-called conditions modules where to look for glue groups. From this, the
# fonction filters through the remaining initial conditions to determine the
# actual conditions for the gluing computed in the factory. This list of
# actual conditions will consist of:
# - a set of genera for the primitive extensions (optional: only if genus_over
#   or form_over is non empty)
# - a set of elementary divisors for the glue groups (optional: only if
#   glue_elementary_divisors is not empty)
# - a set of orders for the glue groups.
#
# If part of the initial conditions, the function runs through the list
# form_over and genus_over to compute the possible genera for an overlattice
# within the context setup by the conditions modules, Fac.par and Fac.sign.
# If none works, Fac.genus_over is the empty set and the following functions
# will abort (because of impossible conditions). Similarly with the list
# glue_elementary_divisors.
#
# For Fac.glue_order, several things can happen:
# - if the initial conditions genus_over, form_over or glue_elementary_divisors
#   are not empty, then the function infers lists of possible orders from the
#   non-empty ones and intersect them. If moreover glue_order is not empty,
#   the previous list of infered orders is intersected with glue_order.
# - otherwise, if glue_order is non-empty, the functions filters which orders
#   could actually work in the context set up by the conditions modules.
# - otherwise, Fac.glue_order consists of all the divisors of the order of the
#   maximal abelian group embedding in both of the conditions modules.
# In the two first cases, if the final list of orders infered from the initial
# conditions is empty, Fac.glue_order is the empty and the rest of the
# algorithm will abort too.
function init_gluing_conditions!(
  Fac::ZZLatGluingFactory,
  left_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  right_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_exponent::Union{nothing, IntegerUnion}=nothing,
  glue_order::AbstractVector{T}=ZZRingElem[],
  glue_elementary_divisors::Vector{Vector{ZZRingElem}}=Vector{ZZRingElem}[],
  form_over::Vector{TorQuadModule}=TorQuadModule[],
  genus_over::Vector{ZZGenus}=ZZGenus[],
) where T <: IntegerUnion

  # Remember if some conditions were initially set
  ignore_genus_over = isempty(genus_over) && isempty(form_over)
  ignore_elem_divs = isempty(glue_elementary_divisors)
  ignore_glue_order = ignore_genus_over && ignore_elem_divs && isempty(glue_order)
  empty_order_cond = isempty(glue_order)

  qM, qN = ambient_modules(Fac)

  # Left glue group annihilated by the endomorphism left_glue_annihilator
  if !isnothing(left_glue_annihilator)
    @assert domain(left_glue_annihilator) === codomain(left_glue_annihilator) === qM
    VM, VMinqM = kernel(left_glue_annihilator)
  else
    VMinqM = id_hom(qM)
    VM = qM
  end

  # Right glue group annihilated by the endomorphism right_glue_annihilator
  if !isnothing(right_glue_annihilator)
    @assert domain(right_glue_annihilator) === codomain(right_glue_annihilator) === qN
    VN, VNinqN = kernel(right_glue_annihilator)
  else
    VNinqN = id_hom(qN)
    VN = qN
  end

  # If glue_exponent is a positive integer, the exponent of the glue groups
  # should divide it
  if !isnothing(glue_exponent) && glue_exponent > 0
    filter!(Base.Fix1(is_divisible_by, glue_exponent)∘last, elementary_divisors)
    VM, _VMinVM = torsion_subgroup(domain(VMinqM), e)
    VMinqM = compose(_VMinVM, VMinqM)
    VN, _VNinVN = torsion_subgroup(domain(VNinqN), e)
    VNinqN = compose(_VNinVN, VNinqN)
  end

  # This fixes the conditions modules where we do everything
  Fac.conditions_modules = (VMinqM, VNinqN)

  elG = _maximal_common_subgroup_snf(VM, VN)
  pds = isempty(elG) ? ZZRingElem[] : prime_divisors(last(elG))
  ged = Dict{ZZRingElem, Dict{Vector{Int}, Vector{ZZLatGluing}}}(p => Dict{Vector{Int}, Vector{ZZLatGluing}}())
  Fac.glue_group_parent_snf = elG
  Fac.primes_of_interest = Set(pds)
  Fac.local_gluings_primary = ged

  if !isempty(genus_over)
    filter!(isequal(Fac.sign)∘signature_pair, genus_over)
  end
  _genus_over = Set(genus_over)

  k1 = prod(elG; init=ZZ(1))
  _glue_order = Set(filter!(>(0), glue_order))
  for o in _glue_order
    !is_divisible_by(k1, o) && delete!(_glue_order, o)
  end

  _glue_elem_divs = Set{Vector{ZZRingElem}}()
  for v in glue_elementary_divisors
    if isempty(v)
      push!(_glue_elem_divs, v)
      continue
    end
    any(<=(0), v) && continue
    length(v) > length(elG) && continue
    flag = is_divisible_by(elG[end-length(v)+1], first(v))
    !flag && continue
    for i in 1:length(v)-1
      flag &= is_divisible_by(v[i+1], v[i])
      flag &= is_divisible_by(elG[end-i+1], v[end-i+1])
      !flag && break
    end
    if flag
      push!(_glue_elem_divs, v)
    end
  end

  for q in form_over
    !isone(modulus_bilinear_form(q)) && continue
    union!(_genus_over, _integer_genera(q, Fac.sign, Fac.par))
  end

  k2 = order(qM)*order(qN)
  if !ignore_genus_over
    _tmp = Set{ZZRingElem}()
    for G in _genus_over
      bool, _k = divides(k2, numerator(abs(det(G))))
      if !bool
        delete!(_genus_over, G)
        continue
      end
      bool, _k = is_square_with_sqrt(_k)
      if !bool
        delete!(_genus_over, G)
        continue
      end
      if empty_order_cond || (_k in _glue_order)
        push!(_tmp, _k)
      else
        delete!(_genus_over, G)
        continue
      end
    end
    if empty_order_cond
      union!(_glue_order, _tmp)
    else
      intersect!(_glue_order, _tmp)
    end
  end

  if !ignore_elem_divs
    _tmp = Set{ZZRingElem}()
    for v in _glue_elem_divs
      _k = prod(v; init=ZZ(1))
      if (ignore_genus_over && empty_order_cond) || (_k in _glue_order)
        push!(_tmp, _k)
      else
        delete!(_glue_elem_divs, v)
        continue
      end
    end
    if ignore_genus_over && empty_order_cond
      union!(_glue_order, _tmp)
    else
      intersect!(_glue_order, _tmp)
    end
  end

  if !ignore_glue_order
    Fac.glue_order = _glue_order
  else
    Fac.glue_order = Set(prime_divisors(k1))
  end

  if !ignore_genus_over
    Fac.genus_over = _genus_over
  end

  if !ignore_elem_divs
    Fac.glue_elementary_divisors = _glue_elem_divs
  end
  return nothing
end

function assert_has_local_classifying_groups!(Fac::ZZLatGluingFactory)
  isdefined(Fac, :local_classifying_groups) && return nothing
  q1, q2 = ambient_modules(Fac)
  as_bilinear_module = (Fac.par != :even)
  Oq1 = __orthogonal_group(q1; as_bilinear_module)
  Oq2 = __orthogonal_group(q2; as_bilinear_module)
  Fac.local_classifying_groups = (Oq1, Oq2)
  return nothing
end

###############################################################################
#
#  Local gluings
#
###############################################################################

function __subgroups_orbit_representatives_and_stabilizers_primary_subtype(
  Fac::ZZLatGluingFactory,
  Winq::TorQuadModuleMap,
  O::AutomorphismGroup{TorQuadModule},
  p::ZZRingElem,
  subtype::Vector{Int},
  i::Int = -1,
)
  W = domain(Winq)
  flag = isdefined(Fac, :Ctx) && (0 < i <= length(Fac.Ctx.modules) && (Fac.Ctx.modules[i] === codomain(Winq))
  if flag && haskey(Fac.Ctx.orb_and_stab, (i, p, subtype))
    subs = Fac.Ctx.orb_and_stab[(i, p, subtype)]
  elseif first(subtype) == 1
    W = domain(Winq)
    _, WWinW = elementary_part(W, p)
    WWinq = compose(WWinW, Winq)
    subs = Oscar._subgroups_orbit_representatives_and_stabilizers_elementary(WWinq, O, p^(length(subtype)), p)
    if flag
      Fac.Ctx.orb_and_stab[(i, p, subtype)] = subs
    end
  else
    subs = _subgroups_orbit_representatives_and_stabilizers_primary_subtype(Winq, O, p, subtype)
    if flag
      Fac.Ctx.orb_and_stab[(i, p, subtype)] = subs
    end
  end
  return subs
end

function _local_gluings_primary!(
  Fac::ZZLatGluingFactory,
  p::ZZRingElem,
  subtype::Vector{Int},
)
  as_bilinear_module = (Fac.par != :even)
  ged = Fac.local_gluings_primary
  @assert haskey(ged, p)
  if haskey(ged[p], subtype)
    return ged[p][subtype]
  end
 
  loc_glue_p = ZZLatGluing[]
  flag1, flag2 = false, false
  i1, i2 = -1, -1
  if isdefined(Fac, :Ctx)
    vi = Fac.vertex_identification
    if vi[1] > 0
      flag1 = true
      i1 = vi[1]
    end
    if vi[2] > 0
      flag2 = true
      i2 = vi[2]
    end
  end

  O1, O2 = local_classifying_groups(Fac)
  V1inq1, V2inq2 = conditions_modules(Fac)
  V1 = domain(V1inq1)
  V2 = domain(V2inq2)

  W1, W1inV1 = primary_part(V1, p)
  W1inq1 = compose(W1inV1, V1inq1)
  W2, W2inV2 = primary_part(V2, p)
  W2inq2 = compose(W2inV2, V2inq2)

  subs1 = __subgroups_orbit_representatives_and_stabilizers_primary_subtype(Fac, W1inq1, O1, p, subtype, i1)
  subs2 = __subgroups_orbit_representatives_and_stabilizers_primary_subtype(Fac, W2inq2, O2, p, subtype, i2)

  for (H1inq1, stab1) in subs1
    H1 = domain(H1inq1)
    for (H2inq2, stab2) in subs2
      H2 = domain(H2inq2)
      ok, phi = is_anti_isometric_with_anti_isometry(H1, H2; as_bilinear_module)
      !ok && continue
      x = ZZLatGluing(phi, inv(phi), H1inq1, stab1, H2inq2, stab2)
      push!(loc_glue_p, x)
    end
  end
  ged[p][subtype] = loc_glue_p
  return loc_glue_p
end

function _trivial_gluing(
  Fac::ZZLatGluingFactory,
)
  q1, q2 = ambient_modules(Fac)
  s1, s2 = id_hom.(local_classifying_groups(Fac))
  z1, z1inq1 = sub(q1, TorQuadModuleElem[])
  z2, z2inq2 = sub(q2, TorQuadModuleElem[])
  phi = hom(z1, z2, zero_matrix(ZZ, 0, 0))
  iphi = hom(z2, z1, zero_matrix(ZZ, 0, 0))
  return ZZLatGluing(phi, iphi, z1inq1, s1, z2inq2, s2)
end

function _local_glue_maps_ord(Fac::ZZLatGluingFactory, o::ZZRingElem)
  isone(o) && return ZZLatGluing[_trivial_gluing(Fac))]

  elG = Fac.glue_group_parent_snf
  _loc = Vector{ZZLatGluing}[]
  for (p, v) in factor(o)
    ptypes = _psubgroups_types(elG, p, v)
    loc_p = ZZLatGluing[]
    for subtype in ptypes
      append!(loc_p, _local_gluings_primary!(Fac, p, subtype))
    end
    isempty(loc_p) && return ZZLatGluing[]
    push!(_loc, loc_p)
  end

  return merge_glue_maps(Fac, _loc)
end

function _local_glue_maps_eldiv(Fac::ZZLatGluingFactory, v::Vector{ZZRingElem})
  isempty(v) && return ZZLatGluing[_trivial_gluing(Fac))]

  pds = Fac.primes_of_interest
  _loc = Vector{ZZLatGluing}[]
  for p in pds
    subtype = reverse!(Int[valuation(a, p) for a in v])
    filter!(!=(0), subtype)
    isempty(subtype) && continue
    loc_p = _local_gluings_primary!(Fac, p, subtype)
    isempty(loc_p) && return ZZLatGluing[]
    push!(_loc, loc_p)
  end

  return merge_glue_maps(Fac, _loc)
end

function merge_glue_maps(
  Fac::ZZLatGluingFactory,
  _loc::Vector{Vector{ZZLatGluing}},
)
  if isempty(_loc)
    return ZZLatGluing[]
  elseif isone(length(_loc))
    return only(_loc)
  end

  res = ZZLatGluing[]
  q1, q2 = ambient_modules(Fac)
  for x in Hecke.cartesian_product_iterator(_loc; inplace=true)
    gens1 = TorQuadModuleElem[]
    gens2 = TorQuadModuleElem[]
    for y in x
      H1 = domain(glue_group(y, 1))
      H2 = domain(glue_group(y, 2))
      for a in gens(H1)
        push!(gens1, q1(lift(a)))
      end
      for a in gens(H2)
        push!(gens2, q2(lift(a)))
      end
    end
    H1, H1inq1 = sub(q1, gens1)
    H2, H2inq2 = sub(q2, gens2)
    j1 = stabilizer_glue_group(first(x), 1)
    j2 = stabilizer_glue_group(first(x), 2)
    for i in 2:length(_loc)
      _j1 = stabilizer_glue_group(x[i], 1)
      _j2 = stabilizer_glue_group(x[i], 2)
      stab1, jj, _ = intersect(domain(j1), domain(_j1))
      j1 = compose(jj, j1)
      stab2, jj, _ = intersect(domain(j2), domain(_j2))
      j2 = compose(jj, j2)
    end
    @hassert :ZZLatWithIsom 1 domain(j1) == first(stabilizer(codomain(j1), H1inq1))
    @hassert :ZZLatWithIsom 1 domain(j2) == first(stabilizer(codomain(j2), H2inq2))
    phi = hom(H1, H2, gens(H1), gens(H2))
    @assert is_anti_isometry(phi; as_bilinear_module=(Fac.par != :even))
    push!(res, ZZLatGluing(phi, inv(phi), H1inq1, j1, H2inq2, j2))
  end
  return res
end

function local_glue_maps(
  Fac::ZZLatGluingFactory,
)
  res = ZZLatGluing[]
  if isdefined(Fac, :glue_elementary_divisors)
    eldiv = Fac.glue_elementary_divisors
    for v in eldiv
      append!(res, _local_glue_maps_eldiv(Fac, v))
    end
  else
    l = possible_glue_order(Fac)
    for o in l
      append!(res, _local_glue_maps_ord(Fac, o))
    end
  end
  return res
end

###############################################################################
#
#  Operations on gluings
#
###############################################################################

function _split_orbit_left_group(
  x::ZZLatGluing,
  O1::AutomorphismGroup{TorQuadModule},
)
  res = ZZLatGluing[]
  phi, _, H1inq1, s1, H2inq2, s2 = gluing_data(x)
  @assert domain(O1) === codomain(H1inq1)

  H1, H2 = domain(H1inq1), domain(H2inq2)
  q1 = codomain(H1inq1)
  stab1 = domain(s1)
  Oq1 = codomain(s1)

  elq1 = elementary_divisors(q1)
  if isone(order(q1)) || elq1[1] == elq1[end]
    iso1 = isomorphism(PermGroup, Oq1)
  else
    iso1 = id_hom(Oq1)
  end

  @vprintln :ZZLatWithIsom 1 "Split orbit of left glue group"
  splits = double_cosets(codomain(iso), first(iso(stab1)), first(iso(O1)))
  @vprintln :ZZLatWithIsom 1 "    done: $(length(splits))"
  for _h in splits
    h = hom(iso\(representative(_h)))
    ih = inv(h)
    _H1, _H1inq1 = sub(q1, elem_type(q1)[h(H1inq1(a)) for a in gens(H1)])
    _phi = hom(_H1, H2, elem_type(H2)[phi(H1(lift(ih(_H1inq1(a))))) for a in gens(_H1)])
    _iphi = inv(_phi)
    _stab1, _ = Oscar._as_subgroup(Oq1, Oscar.GAPWrap.ConjugateSubgroup(GapObj(stab1), GapObj(Oq1(h))))
    _stab1, _, j1 = intersect(_stab1, O1)
    push!(res, ZZLatGluing(_phi, _iphi, _H1inq1, j1, H2inq2, s2))
  end
  return res
end

function _split_orbit_right_group(
  x::ZZLatGluing,
  O2::AutomorphismGroup{TorQuadModule},
)
  return inv.(_split_orbit_left_group(inv(x), O2))
end

function _pullback_left(
  x::ZZLatGluing,
  f::TorQuadModuleMap,
)
  phi, iphi, i1, s1, i2, s2 = gluing_data(x)
  @assert codomain(gamma) === codomain(i1)
  invf = inv(f)
  q0 = domain(f)
  Oq0 = _orthogonal_group(q0, TorQuadModuleMap[f * hom(g) * invf for g in gens(codomain(s1))]; check=false)
  _, s0 = sub(Oq0, elem_type(Oq0)[Oq0(f * hom(s1(g)) * invf; check=false) for g in gens(domain(s1))])
  H0, i0 = sub(q0, TorQuadModuleElem[f\(i1(a)) for a in gens(domain(i1))])
  H0toH1 = hom(H0, domain(i1), TorQuadModuleElem[i1\(f(i0(a))) for a in gens(H0)])
  psi = H0toH1 * phi
  ipsi = inv(psi)
  return ZZLatGluing(psi, ipsi, i0, s0, i2, s2)
end

function _all_glue_maps(
  x::ZZLatGluing,
)
  res = ZZLatGluing[]
  phi, iphi, H1inq1, s1, H2inq2, s2 = gluing_data(x)
  q1 = codomain(H1inq1)
  q2 = codomain(H2inq2)
  H1 = domain(H1inq1)
  H2 = domain(H2inq2)
  stab1 = domain(s1)
  stab2 = domain(s2)

  OH1 = orthogonal_group(H1)
  imOH1 = elem_type(OH1)[OH1(restrict_automorphism(x, H1inq1); check=false) for x in gens(stab1)]
  act1 = hom(stab1, OH1, imOH1)
  im1, _ = image(act1)
  OH2 = orthogonal_group(H2)
  imOH2 = elem_type(OH2)[OH2(restrict_automorphism(x, H2inq2); check=false) for x in gens(stab2)]
  act2 = hom(stab2, OH2, imOH2)
  im2, _ = image(act2)

  elH2 = elementary_divisors(H2)
  if isone(order(H2)) || elH2[1] == elH2[end]
    iso2 = isomorphism(PermGroup, OH2)
  else
    iso2 = id_hom(OH2)
  end
  OH22 = first(iso2(OH2))
  im22 = first(iso2(im2))
  @vprintln :ZZLatWithIsom 1 "Compute right transversals"
  orb2 = collect(right_transversal(OH22, im22))
  @vprintln :ZZLatWithIsom 1 "    done: $(length(orb2))"
  
  stab1gamma = elem_type(OH2)[OH2(iphi * hom(g) * phi; check=false) for g in gens(im1)]
  S1, _ = sub(OH2, stab1phi)
  S2 = first(iso2(S1))
  omega = gset(S2, (x,g) -> x*OH22(g), orb2)
  for _g in orbits(omega)
    g = hom(iso1\representative(_g))
    phig = compose(phi, g)
    iphig = compose(inv(g), iphi)
    push!(res, ZZLatGluing(phig, iphig, H1inq1, s1, H2inq2, s2))
  end
  return res
end

function _overlattice(
  x::ZZLatGluing,
  D::ZZLatGluingAmbient,
)
  phi, iphi, i1, s1, i2, s2 = expand(x)
  H1, H2 = domain(i1), domain(i2)
  H1inD = i1 * D.j1
  H2inD = i2 * D.j2
  O1, O2 = orthogonal_group(H1), orthogonal_group(H2)
  im1 = elem_type(O1)[O1(restrict_automorphism(g, i1); check=false) for g in gens(domain(s1))]
  act1 = hom(domain(s1), O1, im1)
  im2 = elem_type(O2)[O2(restrict_automorphism(g, i2); check=false) for g in gens(domain(s2))
  act2 = hom(domain(s2), O2, im2)
  V, extinD = _overlattice_with_graph(phi, i1, i2)
  _disc, _stab = _glue_stabilizers(phi, act1, act2, D.k1, D.k2, extinD)
  qV = discriminant_group(V)
  OqV = orthogonal_group(V)
  phi2 = hom(qV, _disc, elem_type(_disc)[_disc(lift(x)) for x in gens(qV)])
  iphi2 = inv(phi2)
  GV, _ = sub(OqV, elem_type(OqV)[OqV(phi2 * g * iphi2; check=false) for g in _stab])
  M, T = _as_sublattice(V, relations(q1), relations(q2))
  return V, M, T, GV
end

###############################################################################
#
#  Primitive extensions
#
###############################################################################

### Local glue maps

function _local_glue_maps(
  q1::TorQuadModule,
  q2::TorQuadModule,
  sign_over::Tuple{Int, Int},
  parity::Symbol;
  Ctx=nothing,
  vi::Tuple{Int, Int}=(-1, -1),
  kwargs...,
)
  Fac = ZZLatGluingFactory(q1, q2, sign_over, parity)
  if !isnothing(Ctx)
    Fac.Ctx = Ctx
    Fac.vertex_identification = vi
  end
  init_gluing_factory!(Fac; kwargs...)
  res = local_glue_maps(Fac)
  return res
end

### Unimodular case

function unimodular_primitive_extensions(
  M::ZZLat,
  N::ZZLat;
  right_action::Union{MatGroup{QQFieldElem, QQMatrix}, Nothing}=nothing,
  right_discriminant_action::Union{AutomorphismGroup{TorQuadModule}, Nothing}=nothing,
  left_action::Union{MatGroup{QQFieldElem, QQMatrix}, Nothing}=nothing,
  left_discriminant_action::Union{AutomorphismGroup{TorQuadModule}, Nothing}=nothing,
  first::Bool=false,
  exist_only::Bool=false,
  even::Bool=(is_even(M) && is_even(N)), #low priority
  parity::Symbol=(even ? :even : :both),
) where T <: Hecke.IntegerUnion

  @req is_integral(M) && is_integral(N) "Only available for integral lattices"

  if !isnothing(right_discriminant_action)
    GMbar = right_discriminant_action
  elseif !isnothing(right_action)
    GMbar, _ = image(discriminant_representation(M, right_action; full=false))
  elseif first || exist_only
    qM = discriminant_group(M)
    GMbar = Oscar._orthogonal_group(qM, TorQuadModuleMap[id_hom(qM)]; check=false)
  else
    GMbar, _ = image_in_Oq(M)
  end

  if !isnothing(left_discriminant_action)
    GNbar = left_discriminant_action
  elseif !isnothing(left_action)
    GNbar, _ = image(discriminant_representation(N, left_action; full=false))
  elseif first || exist_only
    qN = discriminant_group(N)
    GNbar = Oscar._orthogonal_group(qN, TorQuadModuleMap[id_hom(qN)]; check=false)
  else
    GNbar, _ = image_in_Oq(N)
  end

  if !is_even(M) || !is_even(N)
    (parity == :even) && Tuple{ZZLat, ZZLat, ZZLat}[]
  end
  return _unimodular_primitive_extensions(M, N, parity, GMbar, GNbar; first, exist_only)
end

function _unimodular_primitive_extensions(
  M::ZZLat,
  N::ZZLat,
  parity::Symbol,
  GM::AutomorphismGroup{TorQuadModule},
  GN::AutomorphismGroup{TorQuadModule};
  first::Bool=false,
  exist_only::Bool=false,
)
  results = Tuple{ZZLat, ZZLat, ZZLat}[]
  qM = discriminant_group(M)
  @assert qM === domain(GM)
  qN = discriminant_group(N)
  @assert qN === domain(GN)

  parity == :even && (!is_even(M) || !is_even(N)) && return false, results

  as_bilinear_module = (parity != :even)
  ok, phi = is_anti_isometric_with_anti_isometry(qM, qN; as_bilinear_module)
  !ok && return false, results
  if exist_only
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      return true, results
    end
  end

  D, inj = direct_sum(qM, qN; cached=false, as_bilinear_module)
  qMinD, qNinD = inj
  if first
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      L, _ = _overlattice_with_graph(phi, qMinD, qNinD)
      M2, N2 = _as_sublattices(L, M, N)
      push!(results, (L, M2, N2))
      return true, results
    end
  end

  iphi = inv(phi)
  if parity == :even
    OqN = orthogonal_group(qN)
  else
    OqN = orthogonal_group_bilinear(qN)
  end

  genC = elem_type(OqN)[OqN(iphi * hom(g) * phi; check=false) for g in gens(GM)]
  GMphi, _ = sub(OqN, genC)
  elqN = elementary_divisors(qN)
  if isone(order(qN)) || elqN[1] == elqN[end]
    iso = isomorphism(PermGroup, OqN)
  else
    iso = id_hom(OqN)
  end
  reps = double_cosets(codomain(iso), first(iso(GN)), first(iso(GMphi)))
  for _g in reps
    g = iso\(representative(_g))
    phig = compose(phi, hom(g))
    L, _ = _overlattice_with_graph(phig, qMinD, qNinD)
    if parity == :even
      is_even(L) || continue
    elseif parity == :odd
      !is_even(L) || continue
    end
    exist_only && return true, results
    M2, N2 = _as_sublattices(L, M, N)
    push!(results, (L, M2, N2))
    first && return true, results
  end
  return length(res) > 0, results
end

function _primitive_extensions_coprime_left(
  M::ZZLat,
  N::ZZLat,
  parity::Symbol,
  GM::AutomorphismGroup{TorQuadModule},
  GN::AutomorphismGroup{TorQuadModule};
  first::Bool=false,
  exist_only::Bool=false,
)
  results = Tuple{ZZLat, ZZLat, ZZLat}[]
  qM = discriminant_group(M)
  @assert qM === domain(GM)
  qN = discriminant_group(N)
  @assert qN === domain(GN)

  pds = prime_divisors(order(qM))
  gensHN = TorQuadModuleElem[]
  for p in pds
    _H, j = primary_part(qN, p)
    append!(gensHN, j.(gens(_H)))
  end
  HN, HNinqN = sub(qN, gensHN)

  as_bilinear_module = (partiy != :even)
  ok, phi = is_anti_isometric_with_anti_isometry(qM, HN; as_bilinear_module)
  !ok && return false, results
  if exist_only
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      return true, results
    end
  end

  D, inj = direct_sum(qM, qN; cached=false, as_bilinear_module)
  qMinD, qNinD = inj
  HNinD = compose(HNinqN, qNinD)
  if first
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      L, _ = _overlattice_with_graph(phi, qMinD, HNinD)
      M2, N2 = _as_sublattices(L, M, N)
      push!(results, (L, M2, N2))
      return true, results
    end
  end

  iphi = inv(phi)
  if parity == :even
    is_even(M) && is_even(N) || return results
    OHN = orthogonal_group(HN)
  else
    OHN = orthogonal_group_bilinear(HN)
  end
  
  resHN = restrict_automorphism_group(GN, HNinqN)
  imHN, _ = image(resHN)
  genC = elem_type(OHN)[OHN(iphi * hom(g) * phi; check=false) for g in gens(GM)]
  GMphi, _ = sub(OHN, genC)
  elHN = elementary_divisors(HN)
  if isone(order(HN)) || elHN[1] == elHN[end]
    iso = isomorphism(PermGroup, OHN)
  else
    iso = id_hom(OHN)
  end
  reps = double_cosets(codomain(iso), first(iso(imHN)), first(iso(GMphi)))
  for _g in reps
    g = iso\(representative(_g))
    phig = compose(phi, hom(g))
    L, _ = _overlattice_with_graph(phig, qMinD, HNinD)
    if parity == :even
      is_even(L) || continue
    elseif parity == :odd
      !is_even(L) || continue
    end
    exist_only && return true, results
    M2, N2 = _as_sublattices(L, M, N)
    push!(results, (L, M2, N2))
    first && return true, results
  end
  return length(res) > 0, results
end

function _primitive_extensions_coprime_right(
  M::ZZLat,
  N::ZZLat,
  parity::Symbol,
  GM::AutomorphismGroup{TorQuadModule},
  GN::AutomorphismGroup{TorQuadModule};
  kwargs...
)
  r = _primitive_extensions_coprime_left(N, M, parity, GN, GM; kwargs...)
  t = [1,3,2]
  return eltype(r)[x[t] for x in r]
end

###############################################################################
#
#  Primitive embeddings
#
###############################################################################

### Unimodular top lattice

function _primitive_embeddings_in_unimodular(
  G1::ZZGenus,
  G2s::Vector{ZZGenus};
  check::Bool=true,
)
  @assert is_unimodular(G1)
  results = Vector{ZZLat, ZZLat, ZZLat}[]

  G2s = filter(<(rank(G1))∘rank, G2s)
  filter!(<(signature_pair(G1))∘signature_pair, G2s)
  if is_even(G1)
    filter!(is_even, G2s)
  end

  for G2 in G2s
    append!(results, _primitive_embeddings_in_unimodular_safe(G1, G2))
  end
  return results
end

function _primitive_embeddings_in_unimodular(
  G1::ZZGenus,
  G2::ZZGenus,
)
  @assert is_unimodular(G1)
  if is_even(G1) && !is_even(G2)
    return results
  elseif rank(G1) <= rank(G2)
    return results
  elseif !(signature(G2) < signature(G1))
    return results
  end
  return _primitive_embeddings_in_unimodular_safe(G1, G2)
end

function _primitive_embeddings_in_unimodular_safe(
  G1::ZZGenus,
  G2::ZZGenus,
)
  if is_even(G1)
    parity = :even
    par = :even
  else
    parity = :odd
    par = :both
  end
  sign = signature_pair(G1) .- signature_pair(G2)
  q = rescale(discriminant_group(G2), -1; cached=false)
  GKs = _integer_genera(q, sign, par)
  isempty(GKs) && continue
 
  Ns = ZZLat[]
  for GK in GKs
    append!(Ns, __representatives(GK))
  end
  Ms = __representatives(G2)
  for M in Ms, N in Ns
    append!(results, unimodular_primitive_extensions(M, N; parity))
  end
  return results
end

### Coprime det

function _primitive_embeddings_coprime_det_safe(
  G1::ZZGenus,
  G2::ZZGenus,
)
 if is_even(G1)
    parity = :even
    par = :even
  else
    parity = :odd
    par = :both
  end
  sign = signature_pair(G1) .- signature_pair(G2)
  q1 = discriminant_group(G1)
  q2 = discriminant_group(G2)
  q, _ = direct_sum(q1, rescale(q2, -1; cached=false); cached=false)
  GKs = _integer_genera(q, sign, par)
  isempty(GKs) && continue

  Ns = ZZLat[]
  for GK in GKs
    append!(Ns, __representatives(GK))
  end
  Ms = __representatives(G2)
  for M in Ms
    GM, _ = image_in_Oq(M)
    for N in Ns
      GN, _ = image_in_Oq(N)
      tmp = _primitive_extensions_coprime_left(M, N, parity, GM, GN)
      append!(results, tmp)
    end
  end
  return results
end

### General case

function _primitive_embeddings(
  G1s::Vector{ZZGenus},
  G2::ZZGenus,
)
  results = Tuple{ZZLat, ZZLat, ZZLat}[]
  G1s = filter(>(rank(G2))∘rank, G1s)
  filter!(>(signature_pair(G2))∘signature_pair, G1s)
  if !is_even(G2)
    filter!(!is_even, G1s)
  end
  q2 = discriminant_group(G2)
  Ctx = ZZLatGluingCtx()
  push!(Ctx.modules, q2)

  for G1 in G1s
    if is_unimodular(G1)
      append!(results, _primitive_embeddings_in_unimodular_safe(G1, G2))
    elseif isone(gcd(numerator(det(G1)), numerator(det(G2))))
      append!(results, _primitive_embeddings_coprime_det_safe(G1, G2))
    else
      append!(results, _primitive_embeddings_generic_safe(G1, G2; Ctx, vi=(-1, 1), q2)
    end
  end
  return results
end

function _primitive_embeddings_generic(
  G1::ZZGenus,
  G2::ZZGenus,
)
  if is_even(G1) && !is_even(G2)
    return results
  elseif rank(G1) <= rank(G2)
    return results
  elseif !(signature(G2) < signature(G1))
    return results
  end
  return _primitive_embeddings_generic_safe(G1, G2)
end

function _primitive_embeddings_generic_safe(
  G1::ZZGenus,
  G2::ZZGenus;
  Ctx=nothing,
  vi::Tuple{Int, Int}=(-1, -1),
  q2::TorQuadModule=discriminant_group(G2),
)
  R = rescale(representative(G1), -1; cached=false)
  U = hyperbolic_plane_lattice()
  if is_even(G1)
    parity = :even
    T, _ = direct_sum(R, U; cached=false)
  else
    parity = :both
    T, _ = direct_sum(R, U, U; cached=false)
  end
  q1n = discriminant_group(T)
  signT = signature_tuple(T)[[1,3]]
  signV = signT .+ signature_pair(G2)
  signK = signature_pair(G1) .- signature_pair(G2)

  lgm = _local_glue_maps(q2, q1n, signV, parity; Ctx, vi=reverse(vi))
  isempty(lgm) && continue

  Ms = __representatives(G2)
  results = Tuple{ZZLat, ZZLat, ZZLat}[]
  for M in Ms
    qM = discriminant_group(M)
    D = _gluing_ambient(qM, q1n, parity)
    GM, _ = image_in_Oq(M)
    _, phiM = is_isometric_with_isometry(qM, q2)
    for x in lgm
      xM = _pullback_left(x, phiM)
      xMs = _split_orbit_left_group(xM, GM)
      for _x in xMs, y in _all_glue_maps(_x)
        V, M2, T2, GV = _overlattice(y, D)
        resV = Tuple{ZZLat, ZZLat, ZZLat}[]
        qV = domain(GV)
        qK = rescale(qV, -1; cached=false)
        GKs = _integer_genera(qK, signK, parity)
        isempty(GKs) && continue
        Ks = reduce(vcat, Vector{ZZLat}[__representatives(GK) for GK in GKs])
        for K in Ks
          pe = unimodular_primitive_extensions(V, K; right_discriminant_action=GV, parity)
          append!(resV, pe)
        end
        for (S, V2, W2) in resV
          T3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(T2), zero_matrix(QQ, rank(T2), degree(W2)-degree(T2))))
          L = lll(orthogonal_submodule(S, T3))
          genus(L) == G1 || continue
          M3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(M2), zero_matrix(QQ, rank(M2), degree(W2)-degree(M2))))
          @assert is_sublattice(L, M3)
          @assert is_primitive(L, M3)
          N = orthogonal_submodule(L, M3)
          bM = coordinates(basis_matrix(M3), L)
          bN = coordinates(basis_matrix(N3), L)
          L = Hecke.integer_lattice(; gram=gram_matrix(L))
          M3 = lattice_in_same_ambient_space(L, bM)
          N = lattice_in_same_ambient_space(L, bN)
          push!(results, (L, M3, N))
        end
      end
    end
  end
  return results
end
