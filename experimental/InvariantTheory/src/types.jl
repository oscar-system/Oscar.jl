# These are objects that carry information about a linearly reductive group.
# As of now it is only implemented for SLm, its direct products and tensors.
mutable struct LinearlyReductiveGroup
  field::Field
  group::Tuple{Symbol,Int}
  defining_ideal::MPolyIdeal
  natural_representation::MatElem

  function LinearlyReductiveGroup(sym::Symbol, m::Int, fld::Field) #have not decided the representation yet
    R, _ = polynomial_ring(fld, :z => (1:m, 1:m); cached=false)
    return LinearlyReductiveGroup(sym, m, R)
  end

  function LinearlyReductiveGroup(sym::Symbol, m::Int, pring::MPolyRing) #the ring input is the group ring
    G = new()
    fld = base_ring(pring)
    @req sym === :SL "only :SL is supported"
    # By Nagata, a group is linearly reductive in characteristic p > 0 only if
    # its identity component is a torus.
    @req characteristic(fld) == 0 "the field must have characteristic 0, as SLm is not linearly reductive in positive characteristic for m >= 2"
    G.field = fld
    @req m^2 == ngens(pring) "ring not compatible"
    G.group = (sym, m)
    M = transpose(matrix(pring, m, m, gens(pring)))
    G.natural_representation = M
    G.defining_ideal = ideal([det(M) - 1])
    #base ring of M has to be the same as the representation matrix when that is created later.
    return G
  end

  # Whether I defines a linearly reductive group cannot be checked.
  function LinearlyReductiveGroup(I::MPolyIdeal)
    fld = coefficient_ring(base_ring(I))
    @req fld isa Field "the coefficients of the defining ideal must lie in a field"
    G = new()
    G.defining_ideal = I
    G.field = fld
    return G
  end
end

# Objects of type LinearlyReductiveGroup can be embedded in GLn (for some n) via
# a representation. This defines the action on a vector space.
# We set up an object RepresentationLinearlyReductiveGroup that carries information
# about this representation.
mutable struct RepresentationLinearlyReductiveGroup
  group::LinearlyReductiveGroup
  rep_mat::MatElem

  #stores if the representation is on symmetric forms, and of which degree.
  sym_deg::Tuple{Bool,Int}

  #representation of group G over symmetric degree d
  function RepresentationLinearlyReductiveGroup(G::LinearlyReductiveGroup, d::Int)
    R = new()
    R.group = G
    R.rep_mat = _rep_mat(G, d)
    R.sym_deg = (true, d)
    return R
  end

  #matrix M is the representation matrix. does not check M.
  function RepresentationLinearlyReductiveGroup(G::LinearlyReductiveGroup, M::MatElem)
    @req base_ring(M) == base_ring(defining_ideal(G)) "Defining ideal and representation matrix must have same parent ring"
    R = new()
    R.group = G
    R.rep_mat = M
    R.sym_deg = (false, 0)
    return R
  end
end

# Invariant rings of linearly reductive groups
@attributes mutable struct RedGroupInvarRing{FldT,PolyRingElemT,PolyRingT} <: InvarRing
  field::FldT
  poly_ring::PolyRingT # graded

  group::LinearlyReductiveGroup
  representation::RepresentationLinearlyReductiveGroup

  fundamental::Vector{PolyRingElemT}
  presentation::MPolyAnyMap{MPolyQuoRing{PolyRingElemT},PolyRingT,Nothing,PolyRingElemT}

  #Invariant ring of reductive group G (in representation R), no other input.
  function RedGroupInvarRing(R::RepresentationLinearlyReductiveGroup) #here G already contains information n and rep_mat
    G = group(R)
    K = base_ring(G)
    n = ncols(R.rep_mat)
    poly_ring, _ = graded_polynomial_ring(K, :X => 1:n; cached=false)
    z = new{typeof(K),elem_type(poly_ring),typeof(poly_ring)}()
    z.representation = R
    z.group = G
    z.field = K
    z.poly_ring = poly_ring
    return z
  end

  #to compute invariant ring ring^G where G is the reductive group of R.
  function RedGroupInvarRing(R::RepresentationLinearlyReductiveGroup, ring::MPolyDecRing)
    n = ncols(R.rep_mat)
    n == ngens(ring) || error("The given polynomial ring is not compatible.")
    G = group(R)
    K = base_ring(G)
    z = new{typeof(K),elem_type(ring),typeof(ring)}()
    z.representation = R
    z.group = G
    z.field = K
    z.poly_ring = ring
    return z
  end
end

struct TorusGroup
  field::Field
  rank::Int
end

struct RepresentationTorusGroup
  group::TorusGroup
  weights::Vector{Vector{ZZRingElem}}
end

@attributes mutable struct TorGroupInvarRing{FldT,PolyRingElemT,PolyRingT} <: InvarRing
  field::FldT
  poly_ring::PolyRingT #graded

  group::TorusGroup
  representation::RepresentationTorusGroup

  fundamental::Vector{PolyRingElemT}
  presentation::MPolyAnyMap{MPolyQuoRing{PolyRingElemT},PolyRingT,Nothing,PolyRingElemT}

  #Invariant ring of reductive group G (in representation R), no other input.
  function TorGroupInvarRing(R::RepresentationTorusGroup) #here G already contains information n and rep_mat
    n = length(weights(R))
    super_ring, _ = graded_polynomial_ring(base_ring(group(R)), :X => 1:n; cached=false)
    return TorGroupInvarRing(R, super_ring)
  end

  #to compute invariant ring ring^G where G is the reductive group of R.
  function TorGroupInvarRing(R::RepresentationTorusGroup, ring_::MPolyDecRing)
    K = base_ring(group(R))
    z = new{typeof(K),elem_type(ring_),typeof(ring_)}()
    n = length(weights(R))
    z.field = K
    z.poly_ring = ring_
    z.representation = R
    z.group = group(R)
    return z
  end
end
