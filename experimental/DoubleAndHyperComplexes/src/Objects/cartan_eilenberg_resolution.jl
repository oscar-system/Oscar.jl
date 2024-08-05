#= Cartan Eilenberg resolutions of 1-dimensional complexes
#
#  Suppose 
#
#    0 ← C₀  ← C₁  ← C₂  ← …
#
#  is a bounded below complex. We compute a double complex
#
#         0     0     0
#         ↑     ↑     ↑
#    0 ← P₀₀ ← P₀₁ ← P₀₂ ← …
#         ↑     ↑     ↑
#    0 ← P₁₀ ← P₁₁ ← P₁₂ ← …
#         ↑     ↑     ↑
#    0 ← P₂₀ ← P₂₁ ← P₂₂ ← …
#         ↑     ↑     ↑
#         ⋮     ⋮     ⋮
#
#  which is quasi-isomorphic to C via some augmentation map 
#
#    ε = (εᵢ : P₀ᵢ → Cᵢ)ᵢ
#
#  The challenge is that if we were only computing resolutions of the Cᵢ's 
#  and lifting the maps, then the rows of the resulting diagrams would 
#  not necessarily form complexes. To accomplish that, we split the original 
#  complex into short exact sequences
#
#    0 ← Bᵢ ← Cᵢ ← Zᵢ ← 0
#
#  and apply the Horse shoe lemma to these. Together with the induced maps 
#  from Bᵢ ↪  Zᵢ₋₁ we get the desired double complex.
#
#  If the original complex C is known to be exact, then there is no need 
#  to compute the resolutions of both Bᵢ and Zᵢ and we can shorten the procedure.
=#
### Production of the chains
struct CEChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  c::AbsHyperComplex
  is_exact::Bool
  kernel_resolutions::Dict{Int, <:AbsHyperComplex} # the kernels of  Cᵢ → Cᵢ₋₁
  boundary_resolutions::Dict{Int, <:AbsHyperComplex} # the boundaries of Cᵢ₊₁ → Cᵢ
  induced_maps::Dict{Int, <:AbsHyperComplexMorphism} # the induced maps from the free 
                                                     # resolutions of the boundary and kernel

  function CEChainFactory(c::AbsHyperComplex; is_exact::Bool=false)
    @assert dim(c) == 1 "complex must be 1-dimensional"
    #@assert has_lower_bound(c, 1) "complex must be bounded from below"
    return new{chain_type(c)}(c, is_exact, Dict{Int, AbsHyperComplex}(), Dict{Int, AbsHyperComplex}(), Dict{Int, AbsHyperComplexMorphism}())
  end
end

function kernel_resolution(fac::CEChainFactory, i::Int)
  if !haskey(fac.kernel_resolutions, i)
    Z, _ = kernel(fac.c, i)
    fac.kernel_resolutions[i] = free_resolution(SimpleFreeResolution, Z)[1]
  end
  return fac.kernel_resolutions[i]
end

function boundary_resolution(fac::CEChainFactory, i::Int)
  if !haskey(fac.boundary_resolutions, i)
    Z, _ = boundary(fac.c, i)
    fac.boundary_resolutions[i] = free_resolution(SimpleFreeResolution, Z)[1]
  end
  return fac.boundary_resolutions[i]
end

function induced_map(fac::CEChainFactory, i::Int)
  if !haskey(fac.induced_maps, i)
    Z, inc = kernel(fac.c, i)
    B, pr = boundary(fac.c, i)
    @assert ambient_free_module(Z) === ambient_free_module(B)
    img_gens = elem_type(Z)[Z(g) for g in ambient_representatives_generators(B)]
    res_Z = kernel_resolution(fac, i)
    res_B = boundary_resolution(fac, i)
    aug_Z = augmentation_map(res_Z)
    aug_B = augmentation_map(res_B)
    img_gens = gens(res_B[0])
    img_gens = aug_B[0].(img_gens)
    img_gens = elem_type(res_Z[0])[preimage(aug_Z[0], Z(repres(aug_B[0](g)))) for g in gens(res_B[0])]
    psi = hom(res_B[0], res_Z[0], img_gens; check=true) # TODO: Set to false
    @assert domain(psi) === boundary_resolution(fac, i)[0]
    @assert codomain(psi) === kernel_resolution(fac, i)[0]
    fac.induced_maps[i] = lift_map(boundary_resolution(fac, i), kernel_resolution(fac, i), psi; start_index=0)
  end
  return fac.induced_maps[i]
end

function (fac::CEChainFactory)(self::AbsHyperComplex, I::Tuple)
  (i, j) = I # i the resolution index, j the index in C

  res_Z = kernel_resolution(fac, j)

  if can_compute_map(fac.c, 1, (j,))
    if fac.is_exact # Use the next kernel directly
      res_B = kernel_resolution(fac, j-1)
      return direct_sum(res_B[i], res_Z[i])[1]
    else
      res_B = boundary_resolution(fac, j-1)
      return direct_sum(res_B[i], res_Z[i])[1]
    end
  end
  # We may assume that the next map can not be computed and is, hence, zero.
  return res_Z[i]
end

function can_compute(fac::CEChainFactory, self::AbsHyperComplex, I::Tuple)
  (i, j) = I
  can_compute_index(fac.c, (j,)) || return false
  return i >= 0
end

### Production of the morphisms 
struct CEMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType} end

function (fac::CEMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  (i, j) = I
  cfac = chain_factory(self)
  if p == 1 # vertical upwards maps
    if can_compute_map(cfac.c, 1, (j,))
      # both dom and cod are direct sums in this case
      dom = self[I]
      cod = self[(i-1, j)]
      pr1 = canonical_projection(dom, 1)
      pr2 = canonical_projection(dom, 2)
      @assert domain(pr1) === domain(pr2) === dom
      inc1 = canonical_injection(cod, 1)
      inc2 = canonical_injection(cod, 2)
      @assert codomain(inc1) === codomain(inc2) === cod
      res_Z = kernel_resolution(cfac, j)
      @assert domain(map(res_Z, i)) === codomain(pr2)
      @assert codomain(map(res_Z, i)) === domain(inc2)
      res_B = boundary_resolution(cfac, j-1)
      @assert domain(map(res_B, i)) === codomain(pr1)
      @assert codomain(map(res_B, i)) === domain(inc1)
      return compose(pr1, compose(map(res_B, i), inc1)) + compose(pr2, compose(map(res_Z, i), inc2))
    else
      res_Z = kernel_resolution(cfac, j)
      return map(res_Z, i)
    end
    error("execution should never reach this point")
  elseif p == 2 # the horizontal maps
    dom = self[I]
    cod = self[(i, j-1)]
    if can_compute_map(cfac.c, 1, (j-1,))
      # the codomain is also a direct sum
      if !cfac.is_exact
        psi = induced_map(cfac, j-1)
        phi = psi[i]
        inc = canonical_injection(cod, 2)
        pr = canonical_projection(dom, 1)
        @assert codomain(phi) === domain(inc)
        @assert codomain(pr) === domain(phi)
        return compose(pr, compose(phi, inc))
      else
        inc = canonical_injection(cod, 2)
        pr = canonical_projection(dom, 1)
        return compose(pr, inc)
      end
      error("execution should never reach this point")
    else
      # the codomain is just the kernel
      if !cfac.is_exact
        psi = induced_map(cfac, j-1)
        phi = psi[i]
        pr = canonical_projection(dom, 1)
        return compose(pr, phi)
      else
        pr = canonical_projection(dom, 1)
        return pr
      end
      error("execution should never reach this point")
    end
    error("execution should never reach this point")
  end

    #=
      C = cfac.c[(i,)]
      B, pr_B = boundary(cfac.c, 1, (i-1,))
      if cfac.is_exact
        # use the kernel instead; we have to adjust the map
        Z, inc_Z = kernel(cfac.c, 1, (i-1,))
        B = Z
        phi = map(cfac.c, 1, (i,))
        img_gens = [preimage(inc_Z, phi(g)) for g in gens(C)]
        pr_B = hom(C, B, img_gens; check=true) # TODO: Set to false
      end
      Z, inc_Z = kernel(cfac.c, 1, (i,))
      @assert domain(pr_B) === C
      @assert codomain(inc_Z) === C
      if j == 1
        img_gens = elem_type(C)[preimage(pr_B, g) for g in gens(B)]
        img_gens = vcat(img_gens, elem_type(C)[inc_Z(g) for g in gens(Z)])
        return hom(B_plus_Z, C, img_gens; check=true) # TODO: Set to false
      else
        # may assume j > 1
        cod = self[(i, j-1)]
        res_Z = getindex!(cfac.kernel_resolutions, i) do
          Z, _ = kernel(fac.c, i)
          return free_resolution(SimpleFreeResolution, Z)
        end
        b = map(res_Z, 1, (j,))
        inc2 = canonical_injection(cod, 2)
        @assert domain(inc2) === res_Z[j-1] === codomain(b)
        @assert codomain(inc2) === cod
        
        if cfac.isexact
          res_Z1 = getindex!(cfac.kernel_resolutions, i-1) do
            Z, _ = kernel(cfac.c, i)
            return free_resolution(SimpleFreeResolution, Z)
          end
          a = map(res_Z1, 1, (j,))
          inc1 = canonical_injection(cod, 1)
          @assert domain(inc1) === res_Z1[j-1] === codomain(a)
          @assert codomain(inc1) === cod
          img_gens = elem_type(cod)[inc1(a(g)) for g in gens(res_Z1[j])]
          img_gens = vcat(img_gens, elem_type(cod)[inc2(b(g)) for g in gens(res_Z[j])])
          return hom(B_plus_Z, cod, img_gens; check=true) # TODO: Set to false
        else
          res_B = getindex!(cfac.boundary_resolutions, i-1) do
            B, _ = boundary(cfac.c, i)
            return free_resolution(SimpleFreeResolution, B)
          end
          a = map(res_B, 1, (j,))
          inc1 = canonical_injection(cod, 1)
          @assert domain(inc1) === res_B[j-1] === codomain(a)
          @assert codomain(inc1) === cod
          img_gens = elem_type(cod)[inc1(a(g)) for g in gens(res_Z1[j])]
          img_gens = vcat(img_gens, elem_type(cod)[inc2(b(g)) for g in gens(res_Z[j])])
          return hom(B_plus_Z, cod, img_gens; check=true) # TODO: Set to false
        end
        error("execution should never reach this point")
      end
    else
      res = getindex!(cfac.kernel_resolutions, i) do
        Z, _ = kernel(fac.c, i)
        return free_resolution(SimpleFreeResolution, Z)
      end
      return map(res, 1, (j,))
    end
  end
  =#
  error("direction $p out of bounds")
end

function can_compute(fac::CEMapFactory, self::AbsHyperComplex, p::Int, I::Tuple)
  (i, j) = I
  if p == 1 # vertical maps
    return i > 0 && can_compute(chain_factory(self).c, j)
  elseif p == 2 # horizontal maps
    return i >= 0 && can_compute_map(chain_factory(self).c, j)
  end
  return false
end

### The concrete struct
@attributes mutable struct CartanEilenbergResolution{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function CartanEilenbergResolution(
      c::AbsHyperComplex{ChainType, MorphismType};
      is_exact::Bool=false
    ) where {ChainType, MorphismType}
    @assert dim(c) == 1 "complexes must be 1-dimensional"
    @assert has_lower_bound(c, 1) "complexes must be bounded from below"
    @assert direction(c, 1) == :chain "resolutions are only implemented for chain complexes"
    chain_fac = CEChainFactory(c; is_exact)
    map_fac = CEMapFactory{MorphismType}() # TODO: Do proper type inference here!

    # Assuming d is the dimension of the new complex
    internal_complex = HyperComplex(2, chain_fac, map_fac, [:chain, :chain]; lower_bounds = Union{Int, Nothing}[0, lower_bound(c, 1)])
    # Assuming that ChainType and MorphismType are provided by the input
    return new{ChainType, MorphismType}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::CartanEilenbergResolution) = c.internal_complex

