@doc raw"""
    chern_class(m::AbstractFTheoryModel, k::Int)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
this very cohomology class ``h``, otherwise it raises an error.

!!! note "Completeness check"
    The implemented algorithm is guaranteed to work only for toric ambient spaces
    that are smooth and **complete**. Verifying completeness can be very time 
    consuming. To skip this check, pass the optional keyword argument 
    `completeness_check=false`.

!!!warning
    This method works ONLY for F-theory models which are hypersurfaces in a toric ambient space.
  
!!!warning
    This method represents the Chern classes of said hypersurface by cohomology classes on the toric ambient space.
    These classes counterparts must be restricted to the hypersurface to truly represent the Chern class in question.
    Internally, we integrate those ambient space classes against the class of the hypersurface, which automatically
    executes the restriction to the hypersurface.

# Examples
```jldoctest; setup = :(Oscar.ensure_qsmdb_installed())
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> h = chern_class(qsm_model, 4; completeness_check = false);

julia> is_trivial(h)
false
```
"""
function chern_class(m::AbstractFTheoryModel, k::Int; completeness_check::Bool=true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric ambient space"
  @req k >= 0 "Chern class index must be non-negative"
  @req k <= dim(ambient_space(m)) - 1 "Chern class index must not exceed dimension of the space"
  @req is_smooth(ambient_space(m)) "Chern classes require a smooth and complete ambient space"
  if completeness_check
    @req is_complete(ambient_space(m)) "Chern classes require a smooth and complete ambient space"
  end

  # If thus far, no non-trivial Chern classes have been computed for this toric variety, add an "empty" vector
  coho_R = cohomology_ring(ambient_space(m); completeness_check)
  if !has_attribute(m, :chern_classes)
    cs = Dict{Int64,CohomologyClass}()
    cs[0] = cohomology_class(ambient_space(m), one(coho_R); completeness_check)
    diff =
      degree(leading_term(hypersurface_equation(m))) -
      sum(coordinate_ring(ambient_space(m)).d)
    cs[1] = cohomology_class(
      toric_divisor_class(ambient_space(m), diff); completeness_check
    )
    set_attribute!(m, :chern_classes, cs)
    if k == 0
      return cs[0]
    elseif k == 1
      return cs[1]
    end
  end

  # Check if the Chern class in question is known
  cs = get_attribute(m, :chern_classes)::Dict{Int64,CohomologyClass}
  if haskey(cs, k)
    return cs[k]
  end

  # Chern class is not known, so compute and return it...
  tdc = toric_divisor_class(
    ambient_space(m), degree(leading_term(hypersurface_equation(m)))
  )
  cy = cohomology_class(tdc; completeness_check)
  ck_ambient = chern_class(ambient_space(m), k; completeness_check)
  ckm1 = chern_class(m, k - 1; completeness_check)
  new_poly = lift(polynomial(ck_ambient)) - lift(polynomial(cy)) * lift(polynomial(ckm1))
  cs[k] = cohomology_class(ambient_space(m), coho_R(new_poly); completeness_check)
  set_attribute!(m, :chern_classes, cs)
  return cs[k]
end

@doc raw"""
    chern_classes(m::AbstractFTheoryModel)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute a cohomology class ``h`` on the
toric ambient space ``X_\Sigma``, such that its restriction to ``Y_n`` is the k-th Chern class
``c_k`` of the tangent bundle of ``Y_n``. If those assumptions are satisfied, this method returns
a vector with the cohomology classes corresponding to all non-trivial Chern classes ``c_k`` of
``Y_n``. Otherwise, this methods raises an error.

As of right now, this method is computationally expensive for involved toric ambient spaces,
such as in the example below.

!!! note "Completeness check"
    The implemented algorithm is guaranteed to work only for toric ambient spaces
    that are simplicial and **complete**. Verifying completeness can be very time 
    consuming. To skip this check, pass the optional keyword argument 
    `completeness_check=false`.

# Examples
```jldoctest; setup = :(Oscar.ensure_qsmdb_installed())
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> h = chern_classes(qsm_model; completeness_check = false);

julia> is_one(polynomial(h[0]))
true

julia> is_trivial(h[1])
true

julia> is_trivial(h[2])
false
```
"""
function chern_classes(m::AbstractFTheoryModel; completeness_check::Bool=true)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Chern class of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Chern class of F-theory model currently supported only for toric ambient space"
  for k in 0:(dim(ambient_space(m)) - 1)
    chern_class(m, k; completeness_check)
  end
  return get_attribute(m, :chern_classes)::Dict{Int64,CohomologyClass}
end

@doc raw"""
    euler_characteristic(m::AbstractFTheoryModel)

If the elliptically fibered n-fold ``Y_n`` underlying the F-theory model in question is given
as a hypersurface in a toric ambient space, we can compute the Euler characteristic. If this
assumptions is satisfied, this method returns the Euler characteristic, otherwise it raises an
error.

!!! note "Completeness check"
    The implemented algorithm is guaranteed to work only for toric ambient spaces
    that are smooth and **complete**. Verifying completeness can be very time 
    consuming. To skip this check, pass the optional keyword argument 
    `completeness_check=false`.

# Examples
```jldoctest; setup = :(Oscar.ensure_qsmdb_installed())
julia> using Random;

julia> qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => 4), rng = Random.Xoshiro(1234))
Hypersurface model over a concrete base

julia> h = euler_characteristic(qsm_model; completeness_check = false)
378
```
"""
@attr Int function euler_characteristic(
  m::AbstractFTheoryModel; completeness_check::Bool=true
)
  @req (m isa WeierstrassModel || m isa GlobalTateModel || m isa HypersurfaceModel) "Euler characteristic of F-theory model supported for Weierstrass, global Tate and hypersurface models only"
  @req base_space(m) isa NormalToricVariety "Euler characteristic of F-theory model currently supported only for toric base"
  @req ambient_space(m) isa NormalToricVariety "Euler characteristic of F-theory model currently supported only for toric ambient space"

  # Trigger potential short-cut computation of cohomology ring
  cohomology_ring(ambient_space(m); completeness_check)

  # Compute the cohomology class corresponding to the hypersurface equation
  cy = cohomology_class(
    toric_divisor_class(ambient_space(m), degree(hypersurface_equation(m)))
  )

  # Compute the Euler characteristic
  return Int(integrate(chern_class(m, 4; completeness_check) * cy; completeness_check))
end
