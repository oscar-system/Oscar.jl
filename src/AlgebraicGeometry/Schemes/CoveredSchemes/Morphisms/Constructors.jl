########################################################################
# Constructors for CoveredSchemeMorphism                               #
########################################################################

### This type has no external constructors.

@doc raw"""
    identity_map(X::AbsCoveredScheme) -> AbsCoveredSchemeMorphism

Given a covered scheme `X`, return the identity map of `X` seen as a
covered scheme morphism.

# Examples
```jldoctest
julia> P, x = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> X = variety(ideal(x[1:2]))
Projective variety
  in projective 2-space over QQ with coordinates [x, y, z]
defined by ideal (y, x)

julia> Xcov = covered_scheme(X)
Scheme
  over rational field
with default covering
  described by patches
    1: scheme((y//z), (x//z))
  in the coordinate(s)
    1: [(x//z), (y//z)]

julia> identity_map(Xcov)
Covered scheme morphism
  from scheme over QQ covered with 1 patch
    1a: [(x//z), (y//z)]   scheme((y//z), (x//z))
  to scheme over QQ covered with 1 patch
    1b: [(x//z), (y//z)]   scheme((y//z), (x//z))
given by the pullback function
  1a -> 1b
    (x//z) -> 0
    (y//z) -> 0
```
"""
@attr AbsCoveredSchemeMorphism function identity_map(X::AbsCoveredScheme)
  C = default_covering(X)
  id_cov = identity_map(C)
  result = CoveredSchemeMorphism(X, X, id_cov)
  set_attribute!(result, :inverse, result)
  return result
end



########################################################################
# Closed embeddings                                                    #
########################################################################

### user facing constructors
function CoveredClosedEmbedding(X::AbsCoveredScheme, I::AbsIdealSheaf;
        covering::Covering=default_covering(X), check::Bool=true)
  space(I) === X || error("ideal sheaf is not defined on the correct scheme")
  mor_dict = IdDict{AbsAffineScheme, ClosedEmbedding}() # Stores the morphism fᵢ : Uᵢ → Vᵢ for some covering Uᵢ ⊂ Z(I) ⊂ X.
  rev_dict = IdDict{AbsAffineScheme, AbsAffineScheme}() # Stores an inverse list to also go back from Vᵢ to Uᵢ for those Vᵢ which are actually hit.
  patch_list = Vector{AbsAffineScheme}()
  for U in patches(covering)
    inc = ClosedEmbedding(U, I(U))
    V = domain(inc)
    if !isempty(V)
      mor_dict[V] = inc
      push!(patch_list, V)
      rev_dict[U] = V
    end
  end
  gluing_dict = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}()
  for Unew in keys(mor_dict)
    U = codomain(mor_dict[Unew])
    for Vnew in keys(mor_dict)
      V = codomain(mor_dict[Vnew])
      gluing_dict[(Unew, Vnew)] = LazyGluing(Unew, Vnew, 
                                               RestrictionDataClosedEmbedding(covering[U, V], Unew, Vnew)
                                              )
    end
  end

  Z = isempty(patch_list) ? CoveredScheme(base_ring(X)) : CoveredScheme(Covering(patch_list, gluing_dict, check=false))
  cov_inc = CoveringMorphism(default_covering(Z), covering, mor_dict, check=false)
  return CoveredClosedEmbedding(Z, X, cov_inc, ideal_sheaf=I, check=false)
end

function base_change(phi::Any, f::CoveredClosedEmbedding;
    domain_map::AbsCoveredSchemeMorphism=base_change(phi, domain(f))[2],
    codomain_map::AbsCoveredSchemeMorphism=base_change(phi, codomain(f))[2]
  )
  g = underlying_morphism(f)
  _, bc_g, _ = base_change(phi, g; domain_map, codomain_map)
  II = image_ideal(f)
  bc_II = pullback(codomain_map, II)
  return domain_map, CoveredClosedEmbedding(domain(bc_g), codomain(bc_g), covering_morphism(bc_g); check=false, ideal_sheaf=bc_II), codomain_map
end


########################################################################
# The standard constructors
########################################################################
@doc raw"""
    composite_map(f::AbsCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)

Realize the composition ``x → g(f(x))`` as a composite map, i.e. an
instance of `CompositeCoveredSchemeMorphism`.

# Examples
```jldoctest
julia> IA2 = affine_space(QQ, [:x, :y])
Affine space of dimension 2
  over rational field
with coordinates [x, y]

julia> (x, y) = gens(OO(IA2));

julia> I = ideal(OO(IA2), [x, y]);

julia> pr = blow_up(IA2, I);

julia> JJ = ideal_sheaf(exceptional_divisor(pr));

julia> inc_E = Oscar.CoveredClosedEmbedding(domain(pr), JJ);

julia> comp = Oscar.composite_map(inc_E, pr)
Composite morphism of
  Hom: scheme over QQ covered with 2 patches -> scheme over QQ covered with 2 patches
  Blow-up: scheme over QQ covered with 2 patches -> scheme over QQ covered with 1 patch

julia> Oscar.maps(comp)[1] === inc_E
true

julia> Oscar.maps(comp)[2] === pr
true

```
"""
function composite_map(f::AbsCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)
  return CompositeCoveredSchemeMorphism([f, g])
end

function composite_map(f::AbsCoveredSchemeMorphism, g::CompositeCoveredSchemeMorphism)
  return CompositeCoveredSchemeMorphism(pushfirst!(Vector{AbsCoveredSchemeMorphism}(copy(maps(g))), f))
end

function composite_map(f::CompositeCoveredSchemeMorphism, g::CompositeCoveredSchemeMorphism)
  return CompositeCoveredSchemeMorphism(vcat(maps(f), maps(g)))
end

function composite_map(f::CompositeCoveredSchemeMorphism, g::AbsCoveredSchemeMorphism)
  return CompositeCoveredSchemeMorphism(push!(Vector{AbsCoveredSchemeMorphism}(copy(maps(f))), g))
end
