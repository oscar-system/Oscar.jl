################################################################################
# Computation of Automorphism groups of Enriques surfaces following
# https://doi.org/10.1307/mmj/20195769
# and
# http://www.arxiv.org/abs/1909.10813
# Initial Version by Simon Brandhorst
################################################################################

@doc raw"""
    EnriquesBorcherdsCtx

Holds data for computing invariants of families of Enriques surfaces with fixed root 
invariants, in particular the automorphism group of a general element. 

The setting is the following:
Let ``Y`` be an Enriques surface over an algebraically closed field of characteristic not ``2``
and ``π : X → Y`` its universal covering K3 surface. We denote by ``S_Y`` and ``S_X`` the 
numerical lattices of ``Y`` and ``X``.

Most importantly the type stores a chain of lattices
```math
S_Y(2) \subseteq S_X \subseteq L_{1,25}
```
where ``L_{1,15}`` is an even unimodular lattice of signature ``(1,25)``. 

An easy way to construct examples is to call [`generic_enriques_surface(n::Int)`](@ref).
Here the generic 1-nodal Enriques surface:
```jldoctest
julia> Y = generic_enriques_surface(1)
Enriques Borcherds context
  with det(SX) = 1024
  with root invariant [(:A, 1)]

julia> numerical_lattice(Y)
Integer lattice of rank 10 and degree 10
with gram matrix
[-2    0    0    1    0    0    0    0    0    0]
[ 0   -2    1    0    0    0    0    0    0    0]
[ 0    1   -2    1    0    0    0    0    0    0]
[ 1    0    1   -2    1    0    0    0    0    0]
[ 0    0    0    1   -2    1    0    0    0    0]
[ 0    0    0    0    1   -2    1    0    0    0]
[ 0    0    0    0    0    1   -2    1    0    0]
[ 0    0    0    0    0    0    1   -2    1    0]
[ 0    0    0    0    0    0    0    1   -2    1]
[ 0    0    0    0    0    0    0    0    1   -2]

```
"""
mutable struct EnriquesBorcherdsCtx
  # SY < SX < L26
  L26::ZZLat
  SX::ZZLat # numerical lattice of X
  SY::ZZLat # numerical lattice of Y with form rescaled by 2
  numerical_lattice::ZZLat
  # the following are given with respect to the basis of SY
  initial_walls::Vector{ZZMatrix}
  initial_rays::Vector{ZZMatrix}
  initial_isotropic_rays::Vector{ZZMatrix}
  # automorphisms of the initial chamber
  # written w.r.t. the basis of SY
  initial_automorphisms::Vector{ZZMatrix}
  initial_automorphisms_mod2::Vector{FqMatrix}  # same order as above but mod 2
  initial_chamber::K3Chamber   # a non-degenerate L26/SX chamber
  # whether a -2 vector defines an outer wall is a mod 2 condition
  roots_mod2::Set{FqMatrix}  # mod 2 images of splitting roots
  gramSY::ZZMatrix  # to avoid coercions
  gramSX::ZZMatrix  # to avoid coercions
  membership_test
  # needed for the membership test
  # the following are in SY/2SY w.r.t the basis of SY mod 2
  Dplus_perp::FqMatrix
  Dplus::FqMatrix
  imgs_mod2::Set{FqMatrix}
  # w.r.t the basis given by Dplus
  Gplus
  Gplus_mat::MatrixGroup{FqFieldElem, FqMatrix}
  volume_index::ZZRingElem
  orderGbar::ZZRingElem

  @doc raw"""
      EnriquesBorcherdsCtx(SY::ZZLat, SX::ZZLat, L26::ZZLat, weyl::ZZMatrix)

  # Input:
  -`SY` `SX` and `L26` must be an ascending chain of lattices in the same quadratic space. 
  - `weyl` -- a Weyl vector of `L26` given with respect to the basis of the lattice `L26`. 
  """
  function EnriquesBorcherdsCtx(SY::ZZLat, SX::ZZLat, L26::ZZLat, weyl::ZZMatrix; check::Bool=true)
    # X K3 ---> Y Enriques
    ECtx = new(L26, SX, SY)
    ECtx.gramSY = change_base_ring(ZZ, gram_matrix(SY))
    ECtx.gramSX = change_base_ring(ZZ, gram_matrix(SX))

    @vprintln :EnriquesAuto 2 "computing Borcherds context"
    dataY,_ = BorcherdsCtx(L26, SY, weyl; compute_OR=false, check)
    dataY.membership_test = (x -> true)
    ECtx.initial_chamber = chamber(dataY, dataY.weyl_vector; check)
    # SY + Sm < SX is a primitive extension with glue map phi: D(Sm) -> D(SY)
    Sm = orthogonal_submodule(SX, SY)
    phi, inc_Dminus, inc_Dplus = glue_map(SX, Sm, SY;check=false)
    # H_Sm = pi_Sm(SX) note that H_Sm/Sm is
    H_Sm = cover(domain(phi))
    sv2 = [1//2*(i[1]*basis_matrix(Sm)) for i in short_vectors(Sm, 4) if i[2] == 4]
    sv2 = [domain(phi)(i) for i in sv2 if i in H_Sm]
    ECtx.roots_mod2 = Set(change_base_ring(GF(2), solve(basis_matrix(SY), matrix(QQ, 1, 26, 2*lift(phi(i))); side=:left)) for i in sv2)

    # Cook up the membership test
    @vprint :EnriquesAuto 2 "computing orthogonal group\n"
    OSm = orthogonal_group(Sm)
    GSm,_= image_in_Oq(Sm)
    phiSm = hom(OSm, GSm, GSm.(gens(OSm)); check=false)
    DSm = discriminant_group(Sm)
    Dminus = domain(phi)
    @vprint :EnriquesAuto 2 "computing stabilizer "
    stab_Dminus, inc_stab = stabilizer(GSm, inc_Dminus)
    @vprint :EnriquesAuto 2 "done\n"
    Dminus_perp, inc_Dminus_perp = orthogonal_submodule(DSm, Dminus)
    res, inc = restrict_automorphism_group(stab_Dminus, inc_Dminus_perp)
    @vprint :EnriquesAuto 2 "computing kernel "
    Gminus, inc_Gminus = kernel(inc)
    @vprintln :EnriquesAuto 2 "done"
    ECtx.orderGbar = order(Gminus)

    # compute Gplus
    Q = orthogonal_submodule(L26,SY)
    DQ = discriminant_group(Q)
    glue_SY_Q,inc1,inc2 = glue_map(L26,SY, Q)
    DO = codomain(glue_SY_Q)

    @vprint :EnriquesAuto 2 "preparing membership test\n"
    DSY = discriminant_group(SY)
    ODSY = orthogonal_group(DSY)
    tmp = [preimage(phiSm, i) for i in small_generating_set(Gminus)]
    tmp_hom = [hom(DO,DO, [DO(lift(j)*i) for j in gens(DO)]) for i in tmp]
    gens_Gplus = elem_type(ODSY)[ODSY(inv(inc1)*glue_SY_Q*i*inv(glue_SY_Q)*inc1) for i in tmp_hom]
    Gplus,_ = sub(ODSY, gens_Gplus)
    ECtx.Gplus = Gplus
    @assert order(Gplus) == order(Gminus)
    # iso_Gplus_minus = hom(Gplus, Gminus, gens(Gplus), small_generating_set(Gminus))

    Dplus = codomain(phi)
    ODplus = orthogonal_group(Dplus)
    # turn this into a mod 2 matrix group
    # with respect to the basis of SY mod 2
    B = 1//2*basis_matrix(SY)
    gens_Gplus_mat = FqMatrix[]
    for g in gens(Gplus)
      g2 = reduce(vcat, [change_base_ring(GF(2), solve(B, matrix(QQ, 1, 26, lift(g(DSY(vec(B[i, :]))))); side=:left))  for i in 1:10])
      push!(gens_Gplus_mat, g2)
    end
    Gplus_mat = matrix_group(GF(2), 10, gens_Gplus_mat)
    ECtx.Gplus_mat = Gplus_mat

    snf_Dplus, i_snf = snf(Dplus)
    gens_Dplus = [change_base_ring(GF(2), solve(basis_matrix(SY), matrix(QQ, 1, 26, 2*lift(inc_Dplus(i_snf(i)))); side=:left)) for i in gens(snf_Dplus)]

    gens_Dplus = reduce(vcat, gens_Dplus)
    ECtx.Dplus = gens_Dplus


    Dplus_perp, _ = orthogonal_submodule(DSY, Dplus)
    Dplus_perp,_ = snf(Dplus_perp)

    # needed for hash
    ECtx.Dplus_perp = reduce(vcat,[change_base_ring(GF(2), 2*solve(basis_matrix(SY), matrix(QQ, 1, 26, lift(i)); side=:left)) for i in gens(Dplus_perp)])


    @vprint :EnriquesAuto 2 "done\n"
    ECtx.membership_test = membership_test_as_group
    #upper_bound = mass(ECtx)*length(initial_automorphisms(ECtx))
    return ECtx
  end
end
    
@doc raw"""
    invariant_lattice_of_K3_cover(Y::EnriquesBorcherdsCtx) -> ZZLat
    
Return the invariant lattice of the universal covering K3 surface of `Y` under the Enriques involution. 
"""
invariant_lattice_of_K3_cover(Y::EnriquesBorcherdsCtx) = Y.SY

@doc raw"""
    numerical_lattice_of_K3_cover(Y::EnriquesBorcherdsCtx) -> ZZLat
    
Return the numerical lattice of the universal covering K3 surface of ``Y``.
"""
numerical_lattice_of_K3_cover(Y::EnriquesBorcherdsCtx) = Y.SX

@doc raw"""
    numerical_lattice(Y::EnriquesBorcherdsCtx) -> ZZLat

Return the numerical lattice of the Enriques surface ``Y``.
"""
function numerical_lattice(Y::EnriquesBorcherdsCtx) 
  if !isdefined(Y,:numerical_lattice) 
    Y.numerical_lattice = lattice(rescale(rational_span(Y.SY), 1//2))
  end 
  return Y.numerical_lattice 
end
  
function Base.show(io::IO, ::MIME"text/plain", Y::EnriquesBorcherdsCtx)
  io = pretty(io)
  println(io, "Enriques Borcherds context")
  println(io, Indent(), "with det(SX) = $(det(Y.SX))")
  print(io, "with root invariant ")
  print(io, root_invariant(Y), Dedent())
end

    
function Base.show(io::IO, dat::EnriquesBorcherdsCtx)
  io = pretty(io)
  print(io, "Enriques Borcherds context with det(SX) = $(det(dat.SX))")
end
  
@doc raw"""
    EnriquesBorcherdsCtx(SY2::ZZLat, SX::ZZLat) -> EnriquesBorcherdsCtx 
  
Return a context object for Borcherds' method for Enriques surfaces. 
  
Let ``\pi: X \to Y`` be the universal cover of an Enriques surface and ``\epsilon`` the Enriques involution. 
# Input: 
- `SY2` -- the invariant lattice of the Enriques involution in the numerical lattice `SX` of ``X``. 
"""
function EnriquesBorcherdsCtx(SY2::ZZLat, SX::ZZLat; ample=nothing)
  @req is_subset(SY2, SX) "SY2 must be contained in SX"
  tmp = rescale(SY2, 1//2)
  @req is_unimodular(tmp) && is_even(tmp) && signature_tuple(tmp) == (1, 0, 9) "SY2 must be isomorphic to E10(2)"
  B = solve(basis_matrix(SX), basis_matrix(SY2); side=:left)
  L26, SX1, iSX1 = embed_in_unimodular(lattice(rational_span(SX)), 1, 25, primitive=true,even=true)
  SY1 = lattice(ambient_space(SX1),solve(basis_matrix(SX),basis_matrix(SY2);side=:left)*basis_matrix(SX1))
  if ample === nothing 
    ample = ample_class(SX1)
    # orthogonal projection of ample to SY1
    B = basis_matrix(SY1)
    GB = gram_matrix(ambient_space(SY1))*transpose(B)
    ample = solve(B*GB, 2*ample*GB; side=:left)
  end
  @assert only(ample*gram_matrix(SY2)*transpose(ample))>0
  weyl, u =  borcherds_method_preprocessing(L26, SY1; ample)
  return EnriquesBorcherdsCtx(SY1, SX1, L26, weyl)
end

@doc raw"""
    enriques_surface_automorphism_group(SY2::ZZLat, SX::ZZLat; ample::Union{ZZMatrix,Nothing}=nothing)
    
Compute the automorphism group of an Enriques surface, its nef cone and its smooth rational curves.
  
Let ``\pi: X \to Y`` be the universal cover of an Enriques surface and ``\epsilon`` the Enriques involution. Let ``S_Y`` be the numerical lattice of ``Y``. 
This function computes the image of the natural map 
```math
\varphi_Y \colon \mathrm{Aut}_{s}(Y) \to O(S_Y \otimes \mathbb{F}_2)
```
where ``\mathrm{Aut}_{s}(Y) \leq \mathrm{Aut}(Y)`` denotes the subgroup of semi-symplectic automorphisms.
Note that the kernel of ``\varphi_Y`` is finite.

See [BS22](@cite), [BS22*1](@cite), and [BRS23](@cite) for background and algorithms. 

# Input: 
- `SY2` -- the invariant lattice of the Enriques involution in the numerical lattice `SX` of ``X``. 
- `ample` -- optionally an ample class as a ``1\times x 10`` matrix representing an ample class w.r.t. the basis of `SY2`; if not given, some arbitrary Weyl-chamber is picked.
  
# Output
See [`borcherds_method(::EnriquesBorcherdsCtx)`](@ref) for a description of the output. 
"""
function enriques_surface_automorphism_group(SY2::ZZLat, SX::ZZLat; ample::Union{ZZMatrix,Nothing}=nothing)
  return borcherds_method(EnriquesBorcherdsCtx(SY2::ZZLat, SX::ZZLat; ample))
end 
  
@doc raw"""
    splitting_roots_mod2(Y::EnriquesBorcherdsCtx)
    
Return the image of the splitting roots of ``Y`` in ``S_Y \otimes \mathbb{F}_2``.
"""
function splitting_roots_mod2(Y::EnriquesBorcherdsCtx)
  return Y.roots_mod2
end

@doc raw"""
    root_invariant(Y::EnriquesBorcherdsCtx)
    
Return the root invariant of ``Y``.
"""
function root_invariant(Y::EnriquesBorcherdsCtx)
  Sm = orthogonal_submodule(Y.SX, Y.SY)
  phi, inc_Dminus, inc_Dplus = glue_map(Y.SX, Sm, Y.SY)
  # H_Sm = pi_Sm(SX)
  H_Sm = cover(domain(phi))
  R = rescale(2*H_Sm, 1//2)
  return root_lattice_recognition(R)[1]
end

@doc raw"""
    mass(ECtx::EnriquesBorcherdsCtx)
    
Return the mass of ``Y`` 

Let ``V_0`` be a complete set of representatives of the ``\mathrm{Aut}(Y)``-orbits of ``L_{26}|S_Y``-chambers.
The mass satisfies
```math
\mathrm{mass}(Y)=\sum_{D \in V_0} \frac{1}{\sharp \mathrm{Aut}_G(D)} = \frac{\mathrm{ind}(D_0)}{\sharp \bar G_{X-}}
```
where ``D_0`` is the initial chamber and ``\mathrm{ind}(D_0)`` is its volume index. 
"""
function mass(ECtx::EnriquesBorcherdsCtx)
  return chamber_invariants(ECtx)[1]//ECtx.orderGbar
end

function membership_test_as_group(data::EnriquesBorcherdsCtx, f::FqMatrix)
  # test if f lies in Gplus
  # f must
  # act as identity on Dplus^perp
  return f in data.Gplus_mat
end

function initial_rays(ECtx::EnriquesBorcherdsCtx)
  if !isdefined(ECtx, :initial_rays)
    ECtx.initial_rays = rays(ECtx.initial_chamber)
  end
  return ECtx.initial_rays
end

function initial_isotropic_rays(ECtx::EnriquesBorcherdsCtx)
  if !isdefined(ECtx, :initial_isotropic_rays)
    iso_rays = filter(r->0==r*ECtx.gramSY*transpose(r), rays(ECtx.initial_chamber))
    ECtx.initial_isotropic_rays = iso_rays
  end
  return ECtx.initial_isotropic_rays
end


function _assure_membership_test_as_set(ECtx::EnriquesBorcherdsCtx)
  if isdefined(ECtx,:imgs_mod2)
    return
  end
  @vprint :EnriquesAuto 2 "computing $(order(ECtx.Gplus_mat)) images mod 2\n"
  imgs_mod2 = Set{FqMatrix}()
  for g in ECtx.Gplus_mat
    h = matrix(g)
    push!(imgs_mod2, ECtx.Dplus*h)
  end
  ECtx.imgs_mod2 = imgs_mod2
end

function membership_test_set(data::EnriquesBorcherdsCtx, f::FqMatrix)
  # test if f lies in Gplus
  # f must
  # act as identity on Dplus^perp
  data.Dplus_perp*f == data.Dplus_perp || return false
  # its restriction to Dplus
  # must lie in in Gplus_res
  f_of_Dplus = data.Dplus*f
  return f_of_Dplus in data.imgs_mod2
end

# in principle it would be enough to just store
# tau and parent wall
# then tau can be recomputed via the spanning tree
# that would reduce storage by a factor of about 10
# ... but increase computation time
@doc raw""" 
    EnriquesChamber 

An ``S_Y(2)|L_{26}`` chamber.

The chamber is represented as ``D_0^\tau`` for ``D_0`` 
the initial chamber and ``\tau \in O(S_Y)``. 
"""
struct EnriquesChamber  # no need to make this mutable.
  data::EnriquesBorcherdsCtx
  tau::ZZMatrix
  parent_wall::ZZMatrix # for the spanning tree, the wall via which this chamber was computed
end

function ==(x::EnriquesChamber, y::EnriquesChamber)
  @req x.data === y.data "Chambers must be associated to the same Enriques surface"
  return x.tau == y.tau
end

function Base.hash(D::EnriquesChamber, h::UInt)
  return hash((D.data, D.tau), h)
end

# x given in the basis of the numerical lattice SY
function Base.in(x::ZZMatrix, D::EnriquesChamber)
  G = D.data.gramSY
  Gx = G*transpose(x)
  return all((v*Gx)[1,1]>=0 for v in walls(D))
end

@doc raw"""
    fingerprint(D::EnriquesChamber)
   
Return an invariant of the chamber under ``G``-congruence.

The computation of this invariant is a bit expensive
"""
function fingerprint(D::EnriquesChamber)
  # Let ``D_1 = D_0^g`` and ``D_2 = D_0^h`` 
  # We want inv(g)*h in G => fingerprint(D_1) == fingerprint(D_2)
  # t.f.a.e.
  # D0*g G-cong D0*g
  # g^-1 f h in G for some f in aut(D0)
  # B g^-1 f h cong B mod 2 for some f in aut(D0)
  # B g^-1 f cong B h^-1 for some f in aut(D0)
  # B g^-1 aut(D0) = Bh^-1 aut(D0)
  # We can hash this set, or something canonically computed from it
  B = D.data.Dplus_perp
  T = B*inv(D.tau)
  function _lt(x::FqMatrix,y::FqMatrix)
    n = nrows(x)
    m = ncols(x)
    for i in 1:n
      for j in 1:m
        x0 = iszero(x[i,j])
        y0 = iszero(y[i,j])
        if x0 && !y0
          return true
        elseif !x0 && y0
          return false
        end
      end
    end
    # now x == y
    return false
  end
  S = sort!([T*f for f in initial_automorphisms(D.data)],lt=_lt)
  return S
  # this fingerprint will be cheaper if we precompute the sum
  S = sum(change_base_ring(GF(2),i) for i in initial_automorphisms(D.data))
  return B*inv(D.tau)*S
end

@doc raw"""
    initial_automorphisms(Y::EnriquesBorcherdsCtx) 
    
Compute the stabilizer of ``D_0`` in ``O(S_Y)``.
"""
function initial_automorphisms(Y::EnriquesBorcherdsCtx)
  if !isdefined(Y, :initial_automorphisms)
    @vprint :EnriquesAuto 2 "computing automorphisms"
    Y.initial_automorphisms = aut(Y.initial_chamber)
    @vprintln :EnriquesAuto 2 " done found $(length(Y.initial_automorphisms)) automorphisms"
  end
  return Y.initial_automorphisms
end

@doc raw"""
    initial_automorphisms(Y::EnriquesBorcherdsCtx)
    
Compute the image of ``\mathrm{Aut}(D_0)`` in ``O(S_Y \otimes \mathbb{F}_2)`` for the initial chamber ``D_0`` 
and store it in ``Y``.
"""
function initial_automorphisms_mod2(Y::EnriquesBorcherdsCtx)
  if !isdefined(Y, :initial_automorphisms_mod2)
    Y.initial_automorphisms_mod2 = [change_base_ring(GF(2), i) for i in initial_automorphisms(Y)]
  end
  return Y.initial_automorphisms_mod2
end

@doc raw"""
    initial_walls(Y::EnriquesBorcherdsCtx)

Compute and store the walls of the initial chamber ``D_0``.
"""
function initial_walls(Y::EnriquesBorcherdsCtx)
  if !isdefined(Y, :initial_walls)
    @vprint :EnriquesAuto 2 "computing walls\n"
    Y.initial_walls = walls(Y.initial_chamber)
  end
  return Y.initial_walls
end

function Base.show(io::IO, ::MIME"text/plain", D::EnriquesChamber)
  io = pretty(io)
  println(io, LowercaseOff(), "EnriquesChamber")
  print(io, Indent())
  show(io, MIME("text/plain"), change_base_ring(GF(2), D.tau))
  print(io, Dedent())
end

function Base.show(io::IO, D::EnriquesChamber)
  io = pretty(io)
  print(io, LowercaseOff(), "EnriquesChamber")
end

@doc raw"""
    adjacent_chamber(D::EnriquesChamber, v::ZZMatrix) -> EnriquesChamber

Return return the chamber adjacent to `D` via the wall defined by `v`.
"""
function adjacent_chamber(D::EnriquesChamber, v::ZZMatrix)
  gramS = D.data.gramSY
  ref = _reflection(gramS, v)
  tau2 = D.tau*ref
  Dnew = EnriquesChamber(D.data, tau2, v)
  return Dnew
end

@doc raw"""
    rays(D::EnriquesChamber)
    
Return the list of primitive ray generators of ``D``.
"""
function rays(D::EnriquesChamber)
  r = initial_rays(D.data)
  return [i*D.tau for i in r]
end

@doc raw"""
    isotropic_rays(D::EnriquesChamber)
    
Return the list of primitive isotropic ray generators of ``D``.
"""
function isotropic_rays(D::EnriquesChamber)
  r = initial_isotropic_rays(D.data)
  return [i*D.tau for i in r]
end

@doc raw"""
    walls(D::EnriquesChamber)
    
Return the list of walls of ``D``.
"""
function walls(D::EnriquesChamber)
  gramS = D.data.gramSY
  walls0 = initial_walls(D.data)
  return [r*D.tau for r in walls0]
end

@doc raw"""
    walls_defined_by_rational_curves(D::EnriquesChamber)
    
Return the subset of walls of ``D`` defined by classes of smooth rational curves on ``Y``.
"""
function walls_defined_by_rational_curves(D::EnriquesChamber)
  return [r for r in walls(D) if GF(2).(r) in D.data.roots_mod2] 
end

@doc raw"""
    hom(D1::EnriquesChamber, D2::EnriquesChamber)
    
Return the set of elements of ``G_Y^0`` mapping ``D_1`` to ``D_2`` where 

```math
G_Y^0 = \mathrm{Aut}^*(Y)\mathrm{W}(Y) = \{f \in O(S_Y) \mid f \mbox{ extends to an isometry of }S_X\mbox{ acting trivially on the discriminant group of }S_X\}
```
"""
function hom(D1::EnriquesChamber, D2::EnriquesChamber)
  result = ZZMatrix[]
  tau1inv = inv(D1.tau)
  t1 = change_base_ring(GF(2), tau1inv)
  t2 = change_base_ring(GF(2), D2.tau)
  for (i,g2) in enumerate(initial_automorphisms_mod2(D1.data))
    h2 = t1*g2*t2
    if D1.data.membership_test(D1.data, h2)
      g = initial_automorphisms(D1.data)[i]
      h = tau1inv*g*D2.tau
      push!(result, h)
    end
  end
  return result
end

@doc raw"""
    is_congruent_with_data(D1::EnriquesChamber, D2::EnriquesChamber) -> Bool, ZZMatrix
    
Return whether `D1` and `D2` are congruent under ``G_Y^0`` and if yes
an element ``g \in G_Y^0`` with ``D1^g=D2``.
"""
function is_congruent_with_data(D1::EnriquesChamber, D2::EnriquesChamber)
  tau1inv = inv(D1.tau)
  t1 = change_base_ring(GF(2),tau1inv)
  t2 = change_base_ring(GF(2), D2.tau)
  for (i,g2) in enumerate(initial_automorphisms_mod2(D1.data))
    h2 = t1*g2*t2
    if D1.data.membership_test(D1.data, h2)
      g = initial_automorphisms(D1.data)[i]
      h = tau1inv*g*D2.tau
      return true, h
    end
  end
  return false, tau1inv
end

aut(D::EnriquesChamber) = hom(D, D)

@doc raw"""
    generic_enriques_surface(n::Int)  -> EnriquesBorcherdsCtx
 
Return the ``(\tau,\overline{\tau})``-generic Enriques surface of number ``n`` as in Table 1.1 of [BS22](@cite).

#Input:
``n`` -- an integer between ``1`` and ``184`` different from ``88`` and ``146``.
"""
function generic_enriques_surface(n::Int)
  @req 1<=n<=184 "n must be a number between 1 and 184"
  @req n!=88 && n!=146  "Entries 88 and 146 cannot be constructed. See Remark 1.16 of [BS22](@cite)"
  SY, SX, L26, w, u = load(joinpath(oscardir, "data/TauTaubarGenericEnriquesSurfaces/TauTaubarGenericEnriquesSurfaceNo$(n).mrdi"))
  Y = EnriquesBorcherdsCtx(SY, SX, L26, w; check=false)
end

@doc raw"""
    borcherds_method(Y::EnriquesBorcherdsCtx; max_nchambers=-1)
    
Compute ``\mathrm{Aut}_{s}(Y)`` of the Enriques surface ``Y`` using Borcherds method. 

Here ``\mathrm{Aut}_{s}(Y)`` denotes the group consisting of semi-symplectic automorphisms of ``Y``. 
The quotient ``\mathrm{Aut}(Y)/\mathrm{Aut}_{s}(Y)`` is a finite group and known to be cyclic in most cases [BG24](@cite).

Let ``\pi \colon X \to Y`` be the K3 cover of ``Y``. 

# Input:
- `Y` -- represents an Enriques surface in terms of ``\pi^*\colon S_Y(2) \hookrightarrow S_X``.
- `max_nchambers` -- abort the computation after `max_nchambers` chambers have been computed; return the generators, chambers and curves computed so far. They may not generate the full automorhism group, and not cover all orbits of chambers or rational curves. 

# Output:
1. A matrix group, the image of ``\mathrm{Aut}_{s}(Y) \to O(S_Y)`` 
2. A complete list of representatives of the ``G:=G_Y^0``-congruence classes of ``S_Y|S_X``-chambers.
3. A list of ``(-2)``-vectors representing smooth rational curves on ``Y`` 
   such that any rational curve of ``Y`` is in the same ``\mathrm{Aut}(Y)``-orbit as at least one class in the list. 

All vectors and matrices are represented with respect to the basis of the [`numerical_lattice(::EnriquesBorcherdsCtx)`](@ref)

# Examples 
Let ``Y`` be an Enriques surface with finite automorphism group of type ``I``. Assume that it is very general so that the covering K3 surface ``X`` has picard rank ``19`` and no extra non-symplectic automorphisms.
In the terminology of [BS22](@cite) it is an ``(E_8+A_1,E_8+A_1)``-generic Enriques surface, namely the one of Number 172 in Table 1.1. Following Section 7.4 of loc. cit. we compute some invariants. 

The group ``\mathrm{Aut}(Y)`` is a dihedral group of order ``8`` and its image ``\mathrm{Aut}^*(Y)`` in ``O(S_Y)`` is a group of order ``4``. It is this group we compute.
```jldoctest EnriquesAut
julia> Y = generic_enriques_surface(172)
Enriques Borcherds context
  with det(SX) = 4
  with root invariant [(:A, 1), (:E, 8)]

julia> autY, chambersY, rational_curves_reprY =  borcherds_method(Y);

julia> order(autY)
4
```
The nef-cone of ``Y`` is is the union of the ``\mathrm{Aut}(Y)``-orbits of the chambers in `chambersY`. 
Since the automorphism group of this particular Enriques surface ``Y`` is finite, 
the nef-cone is a finite rational polyhedral cone. 
In our case there is a single chamber ``D`` and it fixed by ``\mathrm{Aut}(Y)``. 
Hence ``D`` is the nef cone of ``Y``. It has ``12`` walls 
and each of them is defined by a smooth rational curve.
```jldoctest EnriquesAut
julia> D = only(chambersY);

julia> length(aut(D))
4

julia> length(walls_defined_by_rational_curves(D))
12

julia> length(walls(D))
12
```

The surface has ``12`` rational curves and ``9`` elliptic fibrations.
```jldoctest EnriquesAut
julia> rational_curvesY = gset(autY, (x,g) -> x*matrix(g), rational_curves_reprY);

julia> length(rational_curvesY)
12

julia> length(isotropic_rays(D))
9
```

The ``\mathrm{Aut}(Y)``-orbits of the rational curves split as ``2+2+2+2+4``.
```jldoctest EnriquesAut; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> [length(i) for i in orbits(rational_curvesY)]
5-element Vector{Int64}:
 2
 2
 4
 2
 2
```

The elliptic fibrations split into ``4`` orbits, as we shall see by two different methods:
```jldoctest EnriquesAut
julia> elliptic_classes = gset(autY,(x,g)->x*matrix(g), isotropic_rays(D));

julia> [length(i) for i in orbits(elliptic_classes)]
4-element Vector{Int64}:
 2
 1
 2
 4

julia> isomorphism_classes_elliptic_fibrations(Y)
4-element Vector{Tuple{Int64, Tuple{Vector{Tuple{Symbol, Int64}}, Vector{Tuple{Symbol, Int64}}}, TorQuadModuleElem}}:
 (2, ([(:E, 8)], []), [0 1 0 0 0 0 0 0 0 0])
 (135, ([(:A, 1)], [(:A, 7)]), [1 0 1 0 0 0 0 0 0 0])
 (270, ([(:D, 8)], []), [0 1 1 0 0 0 0 0 0 0])
 (120, ([(:E, 7)], [(:A, 1)]), [1 1 0 1 0 0 0 0 0 0])
```

Finally, we check the mass formula in this rather trivial case.
```jldoctest EnriquesAut
julia> mass(Y) == sum(1//length(aut(D)) for D in chambersY)
true

```
"""
function borcherds_method(Y::EnriquesBorcherdsCtx; max_nchambers=-1)
  S = Y.SY
  # for G-sets
  F = FreeModule(ZZ,rank(S), cached=false)
  # initialization
  n = rank(S)
  D = EnriquesChamber(Y, identity_matrix(ZZ, n), zero_matrix(ZZ, 1, n))
  waiting_list = EnriquesChamber[D]

  chambers = Dict{UInt64,Vector{EnriquesChamber}}()
  chambers[hash(fingerprint(D))] = EnriquesChamber[D]
  waiting_list = [D]

  automorphisms = Set{ZZMatrix}()
  rational_curves = Set{ZZMatrix}()

  # apparently computing the small generating set is very slow
  autD = aut(D)
  autD_grp = matrix_group(autD)
  @vprint :EnriquesAuto 4 "computing small generating set "
  autD_mod2 = matrix_group([change_base_ring(GF(2),i) for i in autD])
  iso = hom(autD_grp, autD_mod2, gens(autD_mod2),check=false)
  autD = [matrix(preimage(iso,i)) for i in small_generating_set(autD_mod2)]
  if order(autD_mod2) != length(autD)
    K,i = kernel(iso)
    append!(autD, matrix.(gens(K)))
  end
  @vprintln :EnriquesAuto 4 "done"
  # the following was too slow
  #order(autD_grp)  # somehow computing the order makes it faster
  #autD = [matrix(g) for g in small_generating_set(autD_grp)]
  for f in autD
    isone(f) && continue
    push!(automorphisms, f)
  end
  massY = mass(Y)
  mass_explored = QQ(0)
  ntry = 0
  nchambers = 1
  while length(waiting_list) > 0
    ntry = ntry + 1
    if mod(ntry, 5)==0
      @vprint :EnriquesAuto 2 "largest bucket: $(maximum(length(i) for i in values(chambers))) "
      @vprint :EnriquesAuto 1 "buckets: $(length(chambers)) explored: $(nchambers) unexplored: $(length(waiting_list)) gens: $(length(automorphisms)); rat. curv. $(length(rational_curves))\n"
      @vprint :EnriquesAuto 1 "mass left $(massY - mass_explored)"
    end
    D = popfirst!(waiting_list)

    autD = aut(D)
    mass_explored = mass_explored + inv(QQ(length(autD)))
    # we need the orbits of the walls only
    if length(autD) > 1
      # compute the orbits
      @vprint :EnriquesAuto 3 "computing orbits"
      Omega = [F(v) for v in walls(D)]
      W = gset(matrix_group(autD),Omega)
      vv = F(D.parent_wall)
      wallsDmodAutD = [representative(w).v for w in orbits(W) if !(vv in w)]
      @vprintln :EnriquesAuto 3 " done"
    else
      # the minus shouldn't be necessary ... but who knows?
      wallsDmodAutD = (v for v in walls(D) if !(v==D.parent_wall || -v==D.parent_wall))
    end
    # compute the adjacent chambers to be explored
    for v in wallsDmodAutD
      if change_base_ring(GF(2), v) in D.data.roots_mod2 # is v an outer wall of D
        # v comes from a rational curve
        push!(rational_curves, v)
        continue
      end
      Dv = adjacent_chamber(D, v)
      # check G-congruence
      @vprint :EnriquesAuto 5 "checking G-congruence "
      fp = hash(fingerprint(Dv))
      if !haskey(chambers, fp)
        chambers[fp] = EnriquesChamber[]
      end
      is_explored = false
      for E in chambers[fp]
        b, g = is_congruent_with_data(Dv, E)
        if b
          # enough to add a single homomorphism
          push!(automorphisms, g)
          is_explored = true
          break
        end
      end
      @vprintln :EnriquesAuto 5 "done"
      if is_explored
        continue
      end
      push!(chambers[fp], Dv)
      push!(waiting_list, Dv)
      nchambers = nchambers+1
    end
    if max_nchambers != -1 && ntry > max_nchambers
      return collect(automorphisms), reduce(append!,values(chambers), init=EnriquesChamber[]), collect(rational_curves), false
    end
  end
  @vprint :EnriquesAuto 1 "$(length(automorphisms)) automorphism group generators\n"
  @vprint :EnriquesAuto 1 "$(nchambers) congruence classes of chambers \n"
  @vprint :EnriquesAuto 1 "$(length(rational_curves)) orbits of rational curves\n"
  @vprint :EnriquesAuto 1 "mass $(massY)\n"
  @vprint :EnriquesAuto 1 "mass explored $(mass_explored)\n"
  @assert massY == mass_explored
  return matrix_group(collect(automorphisms)), reduce(append!,values(chambers), init=EnriquesChamber[]), collect(rational_curves), true
end

@doc raw"""
  frame_lattice(L::ZZLat, v::QQMatrix)
  
Return a sublattice ``S`` of ``v^\perp`` surjecting onto ``v^\perp/\mathbb{Z}v``.
"""  
function frame_lattice(L::ZZLat, v::QQMatrix)
  @req v*gram_matrix(ambient_space(L))*transpose(v)==0 "vector must be isotropic"
  F = orthogonal_submodule(L, lattice(ambient_space(L), v))
  vF = change_base_ring(ZZ, solve(basis_matrix(F), v; side=:left))
  @assert gcd(vF[1,:])==1
  b = Hecke._complete_to_basis(vF)
  return lll(lattice(ambient_space(F), b[1:end-1,:]*basis_matrix(F)))
end

@doc raw"""
    chamber_invariants(Y::EnriquesBorcherdsCtx)

Return a triple `[volindex, name, root_type]` of invariants of the induced chambers.

The invariants are taken from Table 1.1 and Table 1.2 of [BS22*1](@cite).

# Output:
- The number of Vinberg chambers contained in ``D_0`` is ``O(S_Y \otimes \mathbb{F}_2)/\mathrm{volindex}``.
- The name of the embedding ``S_Y(2) \to L_{1,25}`` stored in ``Y``. 
- The type of the root sublattice of the orthognal complement of ``S_Y(2)`` in ``L_{1,25}``.
"""
function chamber_invariants(Y::EnriquesBorcherdsCtx)
  Sm = orthogonal_submodule(Y.L26, Y.SY)
  rt = root_lattice_recognition(Sm)[1]
  sort!(rt)
  d = _emb_types()
  j = findfirst(i->rt == i[3], d)
  return d[j]
end

@doc raw"""
    _emb_types()
    
Invariants of embeddings of ``E_{10}(2) \to L_{1,15}``.
"""
function _emb_types()
  d = [
  [2^12 * 3^5 * 5^2 * 7, "12A", [(:D, 8)]],
  [2^12 * 3^3 * 5 * 7,"12B", [(:A, 7)]],
  [2^8 * 3^4 * 5 * 7,"20A", [(:D, 4), (:D, 5)]],
  [2^10 * 3^2 * 5 * 7, "20B", [(:D, 4), (:D, 4)]],
  [2^6 * 3^3 * 5 * 7, "20C", vcat([(:A, 1) for i in 1:10], [(:D, 6)])],
  [2^6 * 3^3 * 5 * 7, "20D", [(:A, 3), (:A, 4)]],
  [2^7 * 3^4 * 5, "20E", vcat([(:A, 1) for i in 1:5], [(:A, 5)])],
  [2^9 * 3^2 * 5, "20F", [(:A, 3), (:A, 3)]],
  [2^7 * 3^2 * 5, "40A", vcat([(:A, 1) for i in 1:4], [(:A, 3),(:A, 3)])],
  [2^3 * 3^2 * 5 * 7, "40B", vcat([(:A, 1) for i in 1:8], [(:D, 4),(:D, 4)])],
  [2^3 * 3^2 * 5 * 7, "40C", vcat([(:A, 1) for i in 1:6], [(:A, 3)])],
  [2^5 * 3^2 * 5, "40D", vcat([(:A, 1) for i in 1:12], [(:D, 4)])],
  [2^5 * 3^2 * 5, "40E", vcat([(:A, 1) for i in 1:2], [(:A, 2),(:A, 2)])],
  [2^5 * 3^2,"96A", [(:A, 1) for i in 1:8]],
  [2^3 * 3^2,"96B", [(:A, 1) for i in 1:16]],
  [2^3 * 3^2,"96C", [(:A, 1) for i in 1:4]]]
  for i in d
    sort!(i[3])
  end
  return d
end


@doc raw"""
    isomorphism_classes_elliptic_fibrations(Y::EnriquesBorcherdsCtx)

Return the list of isomorphism classes of elliptic fibrations on the Enriques surface ``Y``.
Let 
```math 
\mathcal{P}^e \colon \mathcal{M}_{En,e} \to \mathcal{M}_{En}, \quad  (Y,f) \mapsto Y
```
be the forgetful map from the moduli space of elliptic Enriques surfaces to the moduli space of Enriques surfaces. Here ``f`` is the numerical class of a half-fiber.

# Output:
A list of triples with the following entries:
- the ramification index of ``\mathcal{P}^e`` in ``(V,f)``
- the simple fibers, the double fibers
- representative of the half fiber ``f`` modulo ``2``
"""
function isomorphism_classes_elliptic_fibrations(Y::EnriquesBorcherdsCtx)
  isotropic_Y = [i for i in discriminant_group(Y.SY) if quadratic_product(i)==0 && !is_zero(i)]
  X = gset(Y.Gplus, (x,g)->g(x), isotropic_Y)
  return  [(length(fbar), reducible_fibers(Y, representative(fbar)), representative(fbar)) for fbar in orbits(X)]
end 
  

@doc raw"""
    reducible_fibers(Y::EnriquesBorcherdsCtx, fbar::TorQuadModuleElem) -> simple fibers, multiple fibers

Return the ADE types and multiplicity of the reducible singular fibers of the genus-1-fibration induced by ``f``.

# Input:
- `fbar` -- the class of ``f/2`` in the discriminant group of ``S_Y(2)`` where ``f`` is the class of a half fiber 
  of an elliptic fibration on ``Y``

# Output:
- The first return value is a list of the ADE-types of the simple fibers.
- The second return value is a list of the ADE-types of the double fibers.
"""
function reducible_fibers(Y::EnriquesBorcherdsCtx, fbar::TorQuadModuleElem)
  SY2 = Y.SY 
  SX = Y.SX
  L26 = Y.L26
  DY = discriminant_group(SY2)
  Sm = orthogonal_submodule(SX, SY2)
  phi, inc_Dplus, _ = glue_map(SX, SY2, Sm)
  Dplus = domain(phi)
  sv2 = [1//2*(i[1]*basis_matrix(Sm)) for i in short_vectors(Sm, 4) if i[2]==4]
  sv2 = [codomain(phi)(i) for i in sv2 if i in cover(codomain(phi))]
  DeltabarY = Set([inc_Dplus(inv(phi)(i)) for i in sv2])
  ebar = 0*DY[1]
  for x in DY
    if inner_product(x,fbar)!=0
      ebar = x
      break
    end
  end
  @assert !iszero(ebar)

  Delta_fbar = [x for x in DeltabarY if inner_product(fbar,x) ==0]
  Delta1 = [ x for x in Delta_fbar if !(x+fbar in Delta_fbar)]
  Delta2 = [x for x in Delta_fbar if x+fbar in Delta_fbar && inner_product(x,ebar)==0]
  return _identify_fibers(Delta1,fbar), _identify_fibers(Delta2,fbar)
end

function _identify_fibers(D::Vector{TorQuadModuleElem},fbar)
  n = length(D)
  g = zero_matrix(GF(2),n,n)
  for i in 1:n
    for j in i+1:n
      if !iszero(inner_product(D[i],D[j]))
        g[i,j] = 1
        g[j,i] = 1
      end
    end
  end
  return _identify_ADE(g)
end

function _identify_ADE(g::MatElem)
  G = graph_from_adjacency_matrix(Undirected,g)
  # We identify an ADE lattice by
  # ((number of roots) // 2, rank of gram mod 2)
  An = [(1, 0), (3, 2), (6, 2), (10, 4), (15, 4), (21, 6), (28, 6), (36, 8), (45, 8), (55, 10)]
  Dn = [(12, 2), (20, 4), (30, 4), (42, 6), (56, 6), (72, 8), (90, 8)]
  En = [(36, 6), (63, 6), (120, 8)]
  fiber_types = Tuple{Symbol,Int}[]
  for c in connected_components(G)
    l = (length(c),rank(g[c,c]))
    if l in An
      i = findfirst(==(l),An)
      t = (:A,i)
    elseif l in Dn
      i = 3+findfirst(==(l),Dn)
      t = (:D,i)
    elseif l in En
      i = 5+findfirst(==(l),En)
      t = (:E,i)
    end
    push!(fiber_types,t)
  end
  sort!(fiber_types)
  return fiber_types
end

function contraction_type(enr::EnriquesBorcherdsCtx, hbar)
  phi = enr.roots_mod2
  k = base_ring(hbar)
  SY = change_base_ring(k, 1//2*gram_matrix(enr.SY))
  SYh = SY*transpose(hbar)
  phi_h = reduce(vcat,[r for r in phi if r*SYh ==0], init=zero_matrix(k,0,10))
  g = phi_h*SY*transpose(phi_h)
  return _identify_ADE(g)
end 

@doc raw"""
    isomorphism_classes_polarizations(Y::EnriquesBorcherdsCtx, h::ZZMatrix)

Return the isomorphism classes of numerical quasi-polarizations of ``Y`` which are of 
the same type as ``h`` along with some invariants. 

See [BG24](@cite).
# Output 
The first return value is the degree of the the forgetful map
```math
\mathcal{P}^h\colon \mathcal{M}_{En,h} \to \mathcal{M}_{En}, \quad (Y,h') \mapsto Y
```
of the moduli space of numerically ``h``-quasi-polarized Enriques complex surfaces.

The second return value is a list of triples representing the fiber ``(\mathcal{P}^h)^{-1}(Y)``
- the ramificaton index of ``\mathcal{P}^h`` at ``(Y,h')``
- the ADE types of the smooth rational curves orthogonal to ``h'``
- a matrix ``\bar g`` mod 2 such that ``h' = h^{gw}`` for ``w`` the unique element of the Weyl group of ``Y`` such that ``h^gw`` is nef. 

# Input:
- `Y` -- representing the Enriques surface in question 
- `h` -- row matrix representing an element in ``S_Y`` with respect to its basis matrix.
"""
function isomorphism_classes_polarizations(Y::EnriquesBorcherdsCtx, h::ZZMatrix)
  @req size(h)==(1,10) "not a row vector of length 10"
  k = GF(2)
  SY = lattice(rational_span(Y.SY))
  h = QQ.(h)
  GYbar = Y.Gplus_mat
  @req (h*gram_matrix(SY)*transpose(h))[1,1]>0 "h must have positive square"
  @vprintln :EnriquesAuto 1 "computing O(SY,h)"
  Oh = stabilizer_in_orthogonal_group(SY, h)
  gensOhbar = [change_base_ring(k, matrix(g)) for g in gens(Oh)]
  Ohbar = matrix_group(gensOhbar)
  
  S = matrix_group(vcat([matrix(g) for g in gens(GYbar)], gensOhbar))
  one_vinberg = 46998591897600
  degPh = one_vinberg//order(Ohbar)
  @vprintln :EnriquesAuto 1 "degree of forgetful map P^h: $(degPh)"
  qSY = discriminant_group(SY)
  SGS = small_generating_set(orthogonal_group(qSY))
  E = 1//2*identity_matrix(QQ, 10)
  OE10bar = matrix_group([change_base_ring(k,reduce(vcat,[matrix(ZZ,1,10,2*lift(qSY(E[i,:])^g)) for i in 1:10])) for g in SGS])
  phi = isomorphism(PermGroup, OE10bar)
  oE10bar, psi = smaller_degree_permutation_representation(codomain(phi))
  f = compose(phi, psi)
  ohbar,_ = f(Ohbar)
  ohbar,_ = sub(oE10bar, small_generating_set(ohbar))
  gYbar,_ = f(GYbar)
  gYbar,_ = sub(oE10bar, small_generating_set(gYbar))
  @vprintln :EnriquesAuto 1 "computing double cosets"
  DC = double_cosets(oE10bar, ohbar, gYbar)
  dc = []
  i = 0
  hbar = change_base_ring(k, h)
  for gG in DC
    i = i+1
    if mod(i,10000)==0
      @vprintln :EnriquesAuto 1 i
    end
    g = representative(gG)
    n = order(gG)
    ramification = n//order(ohbar)
    m = matrix(inv(f)(g))
    representative_mod2 = hbar*m 
    c = Oscar.contraction_type(Y, representative_mod2)
    push!(dc, [ramification, c, m])
  end
  @assert sum(i[1] for i in dc) == degPh
  return degPh, dc
end
