export action_on_ideal,
       action_on_polynomial,
       is_invariant_ideal,
       is_semi_invariant_polynomial,
       projective_group_action,
       pullback,
       push_forward,
       symmetric_intersections,
       underlying_moduli_space_of_modules

### Action on homogeneous ideals/polynomials

@doc Markdown.doc"""
    is_semi_invariant_polynomial(rep::LinRep, f::S) where S <: MPolyElem_dec -> Bool

Given a linear representation `rep` on a finite vector space `V` and a
homogeneous polynomial `f` of degree `d` in the ($\mathbb Z$-graded) polynomial
algebra `R` associated to `V`, return whether `f` is semi-invariant under the induced
action on `R_d`.
""" 
function is_semi_invariant_polynomial(rep::LinRep, f::S) where S <: MPolyElem_dec
  @req is_homogeneous(f) "f must be a homogeneous polynomial"
  R = parent(f)
  @req base_ring(R) === splitting_field(representation_ring(rep)) "The coefficient of f are in the wrong field"
  @req ngens(R) == dimension_representation(rep) "There is no induced action on the polynomial ring of f"
  R1, R1toR = homogeneous_component(R, 1)
  repd = dual_representation(rep)
  for m in matrix_representation(repd)
    poly_m = [R1toR(R1(reverse(vec(collect(m[i,:]))))) for i in 1:nrows(m)]
    fd = evaluate(f, poly_m)
    ok, k = divides(fd, f)
    if !ok || !is_constant(k)
      return false
    end
  end
  return true
end

@doc Markdown.doc"""
    action_on_polynomial(rep::LinRep, f::S) where S <: MPolyElem_dec -> LinRep

Given a linear representation `rep` on a finite vector space `V` and a semi-invariant
homogeneous polynomial `f` of degree `d` in the ($\mathbb Z$-graded) polynomial
algebra `R` associated to `V`, return the 1-dimensional linear representation
corresponding to the induced action on `f`.
""" 
function action_on_polynomial(rep::LinRep, f::S) where S <: MPolyElem_dec
  @req is_homogeneous(f) "f must be a homogeneous polynomial"
  R = parent(f)
  @req base_ring(R) === splitting_field(representation_ring(rep)) "The coefficient of f are in the wrong field"
  @req ngens(R) == dimension_representation(rep) "There is no induced action on the polynomial ring of f"
  R1, R1toR = homogeneous_component(R, 1)
  F = base_ring(R)
  repd = dual_representation(rep)
  coll = eltype(matrix_representation(rep))[]
  for m in matrix_representation(rep)
    poly_m = [R1toR(R1(reverse((vec(collect(m[i,:])))))) for i in 1:nrows(m)]
    fd = evaluate(f, poly_m)
    ok, k = divides(fd, f)
    @req ok && is_constant(k) "Polynomial is not semi-invariant"
    push!(coll, matrix(base_ring(R),1,1,[collect(coefficients(k))[1]])::eltype(matrix_representation(rep)))
  end
  return linear_representation(representation_ring(rep), coll)
end

@doc Markdown.doc"""
    is_invariant_ideal(rep::LinRep, I::S)
                                where S <: MPolyIdeal{<: MPolyElem_dec} -> Bool

Given a linear representation `rep` on a finite vector space `V` and a
homogeneous polynomial `I` in the ($\mathbb Z$-graded) polynomial algebra `R`
associated to `V`, return whether `I` is invariant under the induce action on
`R`.

A homogeneous ideal `I` is invariant under a linear action of a group on its
associated polynomial algebra `R` if and only if the `d`-homogeneous component
`I_d` of `I` is invariant under the corresponding action on `R_d`.
""" 
function is_invariant_ideal(rep::LinRep, I::S) where S <: MPolyIdeal{<: MPolyElem_dec}
  R = base_ring(I)
  @req base_ring(R) === splitting_field(representation_ring(rep)) "The coefficient of f are in the wrong field"
  @req ngens(R) == dimension_representation(rep) "There is no induced action on the polynomial ring of I"
  gene = minimal_generating_set(I)
  R1, R1toR = homogeneous_component(R, 1)
  repd = dual_representation(rep)
  for f in gene for m in matrix_representation(repd)
    poly_m = [R1toR(R1(reverse(vec(collect(m[i,:]))))) for i in 1:nrows(m)]
    fd = evaluate(f, poly_m)
    fd in I || return false
  end
  end
  return true
end

@doc Markdown.doc"""
    action_on_ideal(rep::LinRep, I::S)
                              where S <: MPolyIdeal{<: MPolyElem_dec} -> LinRep

Given a linear representation `rep` on a finite vector space `V` and an invariant
homogeneous polynomial `I` in the ($\mathbb Z$-graded) polynomial algebra `R`
associated to `V`, such that `I = I_d` for some posiitve integer `d, return the
linear representation corresponding to the induced action on `I`.
""" 
function action_on_ideal(rep::LinRep, I::S) where S <: MPolyIdeal{<: MPolyElem_dec}
  R = base_ring(I)
  @req base_ring(R) === splitting_field(representation_ring(rep)) "The coefficient of f are in the wrong field"
  @req ngens(R) == dimension_representation(rep) "There is no induced action on the polynomial ring of I"
  gene = minimal_generating_set(I)
  d = Int(degree(gene[1]).coeff[1])
  @req all(f -> degree(f) == degree(gene[1]), gene) "All generators of I should have same total degree"
  @req is_invariant_ideal(rep, I) "I is not invariant"
  rd = homogeneous_polynomial_representation(rep, d)
  Rd, RdtoR = homogeneous_component(R, d)
  M = [transpose(matrix(reverse(vec(collect((RdtoR\f).v))))) for f in gene]
  M = reduce(vcat, M)
  return action_on_submodule(rd, M)
end

###############################################################################

### Parametrising certain symmetric intersections

@doc Markdown.doc"""
    underlying_moduli_space_of_modules(symci::SymmetricIntersections{S, T, U, V})
                                        where {S, T, U, V} -> CharacterGrassmannian{S, T, U}

Given a parametrizing space `symci` for ideals defining symmetric intersections,
return the underlying symmetric Grassmannian parametrizing submodule of the corresponding
homogeneous part of a polynomial algebra generating the ideals in `symci`.
"""
underlying_moduli_space_of_modules(symci::SymmetricIntersections) = symci.para

@doc Markdown.doc"""
    projective_group_action(symci::SymmetricIntersections{S, T, U, V})
                                               where {S, T, U, V} -> ProjRep{S, T, U, V}

Given a parametrizing space `symci` for ideals defining symmetric intersections,
return the projective representation of the group on the ambient projective space which
encodes the symmetry of the complete intersections defined by the ideals in `symci`.
"""
projective_group_action(symci::SymmetricIntersections) = symci.prep

@doc Markdown.doc"""
    parametrization_data(symci::SymmetricIntersections)
                                        -> Vector{Tuple{Vector{MPolyElem_dec}, Int}}

Given a parametrizing space `symci` for ideals defining symmetric intersections,
return a list of tuples `(B, n)` where `B`'s' are bases for the factor spaces of `symci`
and `n` is the size of a sample in the factor space spanned by `B` to be taken.
"""
function parametrization_data(symci::SymmetricIntersections)
  pd = parametrization_data(symci.para)
  j = symci.j
  pd2 = Tuple{Vector{Vector{MPolyElem_dec}}, Int}[]
  for (B, n) in pd
    B2 = Vector{MPolyElem_dec}[]
    for b in B
      vv = reverse.(vec.(collect.([b[i,:] for i in 1:nrows(b)])))
      vv = domain(j).(vv)
      push!(B2, j.(vv))
    end
    push!(pd2, (B2, n))
  end
  return pd2
end

@doc Markdown.doc"""
    standard_element(symci::SymmetricIntersections)
                                             -> MPolyIdeal_dec

Given a parametrizing space `symci` for ideal defining symmetric intersections,
return an specific ideal in `symci` made up of all possible samples.
"""
function standard_element(symci::SymmetricIntersections)
  std_el = standard_element(symci.para)
  std_el2 = MPolyElem_dec[]
  j = symci.j
  S = codomain(j)
  for B in std_el
    for b in B
      vv = reverse.(vec.(collect.([b[i,:] for i in 1:nrows(b)])))
      vv = j.(domain(j).(vv))
      append!(std_el2, vv)
    end
  end
  return ideal(S, std_el2)
end

@doc Markdown.doc"""
    symmetric_complete_intersections(prep::ProjRep{S, T, U, V}, d::Int, t::Int;
                                     j::MapFromFunc = nothing,
                                     check::Bool = true) -> Vector{SymmetricIntersections{S, T, U, V}}

Given a projective representation `prep` of a finite group `G` on a `F`-vector space `V`, where `F`
is a field of characteristic zero, and two integers `d` and `t`, return the parametrizing space
of ideals defining intersections of `t` hypersurfaces of common degree `d` in $\mathbb P(V)$
which are invariant under the action of `G` given by `prep`.

The optional argument `j` should be an injection from the `d`-th homogeneous component of a
$\mathbb Z$-graded polynomial algebra `S` associated to `V`, seen as an abstract vector space, to `S`.

If `check === true` and `j` is provided, then the function checks `j` satisfies the previous
requirements.

The output is given under the form of a dictionnary whose keys are degree `t` `F`-characters `chi`
of a fixed Schur cover `E` of `G` and the corresponding values are parametrizing space for ideals
for which the action of `E` on a set of generators is given by `chi`.
"""
function symmetric_intersections(prep::ProjRep, d::Int, t::Int; j = nothing, check::Bool = true)
  RR = representation_ring_linear_lift(prep)
  F = splitting_field(RR)
  if j === nothing
    S, _ = grade(PolynomialRing(F, "x" => 0:dimension_linear_lift(prep)-1)[1])
    _, j = homogeneous_component(S, d)
  elseif check
    V = domain(j)
    @req base_ring(V) === F "Incompatible map j"
    @req dim(V) == binomial(dimension_linear_lift(prep)+d-1, d) "Incompatible map j"
    S = codomain(j)
    @req base_ring(S) === F "Incompatible map j"
    @req ngens(S) === dimension_linear_lift(prep) "Incompatible map j"
  end
  rd = homogeneous_polynomial_representation(prep, d)
  M = invariant_grassmannian(rd, t)
  chis = submodule_character.(irreducible_components(M))
  reps = SymmetricIntersections[]
  for chi in chis
    N = irreducible_component(M, chi)
    symci = SymmetricIntersections(prep, N, j)
    push!(reps, symci)
  end
  return reps
end

function Base.show(io::IO, ::MIME"text/plain", symci::SymmetricIntersections)
  M = symci.para
  j = symci.j
  RR = representation_ring_linear_lift(symci.prep)
  F = splitting_field(RR)
  G = underlying_group(symci.prep)
  n = ngens(codomain(j))
  t = submodule_dimension(M)
  d = total_degree(j(gens(domain(j))[1]))
  ty = "($d, "
  for i in 2:t-1
    ty *= "$d, "
  end
  ty *= "$d)"
  println(io, "Parameter space for intersections of type")
  println(io, ty)
  println(io, "in the $(n-1)-dimensional projective space over")
  println(io, F)
  println(io, "and invariant under the action of")
  print(io, G)
end

function Base.show(io::IO, symci::SymmetricIntersections)
  print(io, "Parameter space for symmetric intersections")
end

@doc Markdown.doc"""
    symmetric_intersections(id::Tuple{Int, Int}, n::Int, d::Int, t::Int)
    symmetric_intersections(G::Oscar.GAPGroup, n::Int, d::Int, t::Int)
                              -> Vector{Tuple{ProjRep, Vector{SymmetricIntersections}}}
                              
Given a small group `G` of ID `id` and three integers `n`, `d` and `t`, return a
parametrisation for the ideals defining intersections of `t` hypersurfaces
of the same degree `d` inside $\mathbb P^{n-1}$, which are invariant under the action of
`G` on the latter projective space.

The output is given by a list of dictionnaries, one for each class of faithful
linear action of `G` on $\mathbb P^{n-1}$
"""
function symmetric_intersections(G::Oscar.GAPGroup, n::Int, d::Int, t::Int)
  pfr = faithful_projective_representations(G, n)
  is_empty(pfr) && return InvariantGrassmannian, grade(PolynomialRing(QQ, "x" => 0:n-1)[1])
  F = splitting_field(representation_ring_linear_lift(pfr[1]))
  S, _ = grade(PolynomialRing(F, "x" => 0:n-1)[1])
  _, j = homogeneous_component(S, d)
  res = Tuple{ProjRep, Vector{SymmetricIntersections}}[]
  for prep in pfr
    D = symmetric_intersections(prep, d, t, j = j, check = false)
    push!(res, (prep, D))
  end
  return res
end

symmetric_intersections(id::Tuple{Int, Int}, n::Int, d::Int, t::Int) = symmetric_intersections(small_group(id[1], id[2]), n, d, t)

