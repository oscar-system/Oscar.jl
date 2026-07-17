# Since coeff_map is at the third position, we wannot use a variable
# with default argument


################################################################################
#
#  Domains of type MPolyRing
#
################################################################################

################################################################################
#
#  Constructors
#
################################################################################

function _check_imgs(S::NCRing, imgs)
  n = length(imgs)
  for i in 2:n, j in 1:(i - 1)
    @req imgs[i] * imgs[j] == imgs[j] * imgs[i] "Images $i and $j do not commute"
  end
  return nothing
end

# no check for commutative codomains
_check_imgs(S::Ring, imgs) = nothing 


@doc raw"""
    hom(R::MPolyRing, S::NCRing, coeff_map, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing)

    hom(R::MPolyRing, S::NCRing, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing)

    hom(::Type{K}, R::MPolyRing, S::NCRing, coeff_map, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing) where {K}

    hom(::Type{K}, R::MPolyRing, S::NCRing, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing) where {K}

Create the ring homomorphism `f : R -> S` sending the `i`-th generator of `R`
to `images[i]`. If `coeff_map` is supplied, coefficients are mapped using
`coeff_map`; otherwise Oscar uses `nothing`.

# Map kinds and grading

Oscar attaches a "kind" tag to every map, one of

- `HomUngraded`
- `HomGraded`
- `HomGradedToUngraded`
- `HomUngradedToGraded`

Untyped constructors `hom(R,S,...)` choose the kind from the types of `R` and `S`.
Typed constructors `hom(HomUngraded, R,S,...)`, `hom(HomGraded, R,S,...)`, etc.,
let the user force the kind.

## Graded maps (HomGraded)

If the chosen kind is `HomGraded`, both `R` and `S` must be graded rings.
Let `GR = grading_group(R)` and `GS = grading_group(S)`. A graded map is
determined by

- a group homomorphism `phi : GR -> GS` (in `grading_group_hom`)
- a shift element `deg_shift` in `GS` (in `degree_shift`)

such that for every generator `x_i = gen(R,i)`:

    degree(images[i]) == phi(degree(x_i)) + deg_shift

and each `images[i]` is homogeneous in `S`.

Defaults:

- If `grading_group_hom === nothing` and `GR == GS`, Oscar uses the canonical
  identity-on-generators map `hom(GR, GS, gens(GS))`.
- If `GR != GS` and `grading_group_hom === nothing`, an error is thrown.
- If `degree_shift === nothing`, `deg_shift` is inferred from the image of
  the first generator.

The keywords `grading_group_hom` and `degree_shift` are only meaningful for
kind `HomGraded`. Passing them for other kinds is an error.

## Checks

If `check=true` (default):

- If `S` is noncommutative, the images must pairwise commute.
- If the kind is `HomGraded`, homogeneity and the degree condition above are
  verified.

If `check=false`, validation is skipped. In the `HomGraded` case the map still
stores the computed `grading_group_hom` and `degree_shift` attributes so that
composition of graded maps remains well-defined.

# Examples
```jldoctest; filter = r"\#\d+" => "#"
julia> K, a = finite_field(2, 2, "a");

julia> R, (x, y) = polynomial_ring(K, [:x, :y]);

julia> F = hom(R, R, z -> z^2, [y, x])
Ring homomorphism
  from multivariate polynomial ring in 2 variables over K
  to multivariate polynomial ring in 2 variables over K
defined by
  x -> y
  y -> x
with map on coefficients
  #1

julia> F(a * y)
(a + 1)*x

julia> Qi, i = quadratic_field(-1)
(Imaginary quadratic field defined by x^2 + 1, sqrt(-1))

julia> S, (x, y) = polynomial_ring(Qi, [:x, :y]);

julia> G = hom(S, S, hom(Qi, Qi, -i), [x^2, y^2])
Ring homomorphism
  from multivariate polynomial ring in 2 variables over Qi
  to multivariate polynomial ring in 2 variables over Qi
defined by
  x -> x^2
  y -> y^2
with map on coefficients
  Map: Qi -> Qi

julia> G(x+i*y)
x^2 - sqrt(-1)*y^2

julia> R, (x, y) = polynomial_ring(ZZ, [:x, :y]);

julia> f = 3*x^2+2*x+1;

julia> S, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> H = hom(R, S, gens(S))
Ring homomorphism
  from multivariate polynomial ring in 2 variables over ZZ
  to multivariate polynomial ring in 2 variables over GF(2)
defined by
  x -> x
  y -> y

julia> H(f)
x^2 + 1
```
"""
function hom(::Type{K}, R::MPolyRing, S::NCRing, coeff_map, images::Vector;
             check::Bool = true,
             grading_group_hom = nothing,
             degree_shift = nothing) where {K <: MPolyHomKind}

  n = ngens(R)
  @req n == length(images) "Number of images must be $n"

  imgs = _coerce(S, images)

  K0 = _kind_for(R, S)
  if K === HomGraded
    @req K0 == HomGraded "HomGraded requested but domain/codomain are not both graded."
  elseif K === HomGradedToUngraded
    @req K0 == HomGradedToUngraded "HomGradedToUngraded requested but rings are not according to graded->ungraded."
  elseif K === HomUngradedToGraded
    @req K0 == HomUngradedToGraded "HomUngradedToGraded requested but rings are not according to ungraded->graded."
  elseif K === HomUngraded
    @req (K0 == HomUngraded || K0 == HomGraded) "HomUngraded requested but rings are mixed graded/ungraded."
  else
    error("Unknown map kind $K")
  end

  if K !== HomGraded
    @req grading_group_hom === nothing "grading_group_hom is only allowed for HomGraded"
    @req degree_shift === nothing "degree_shift is only allowed for HomGraded"
  end

  phi = nothing
  deg_shift = nothing
  if K === HomGraded
    phi, deg_shift = _graded_data(R, S, imgs;
                                  grading_group_hom = grading_group_hom,
                                  degree_shift = degree_shift)
  end

  @check begin
    _check_imgs(S, imgs)
    if K === HomGraded
      _check_homo(S, imgs)
      for i in 1:n
        @req degree(imgs[i]) == phi(degree(gen(R, i))) + deg_shift "Degree mismatch for image of generator $i"
      end
    end
  end

  f = MPolyAnyMap(K, R, S, coeff_map, copy(imgs)) # copy because of #655

  if K === HomGraded
    set_attribute!(f, :grading_group_hom => phi)
    set_attribute!(f, :degree_shift => deg_shift)
  end

  return f
end

function hom(::Type{K}, R::MPolyRing, S::NCRing, images::Vector;
             check::Bool = true,
             grading_group_hom = nothing,
             degree_shift = nothing) where {K <: MPolyHomKind}
  return hom(K, R, S, nothing, images;
             check = check,
             grading_group_hom = grading_group_hom,
             degree_shift = degree_shift)
end

# Untyped constructors for backwards compatibility.
function hom(R::MPolyRing, S::NCRing, coeff_map, images::Vector;
             check::Bool = true,
             grading_group_hom = nothing,
             degree_shift = nothing)
  K = _kind_for(R, S)
  return hom(K, R, S, coeff_map, images;
             check = check,
             grading_group_hom = grading_group_hom,
             degree_shift = degree_shift)
end

function hom(R::MPolyRing, S::NCRing, images::Vector;
             check::Bool = true,
             grading_group_hom = nothing,
             degree_shift = nothing)
  return hom(R, S, nothing, images;
             check = check,
             grading_group_hom = grading_group_hom,
             degree_shift = degree_shift)
end


################################################################################
#
#  Evaluation functions
#
################################################################################

# Some additional methods needed for the test in the constructor for MPolyAnyMap
_is_gen(x::MPolyQuoRingElem) = _is_gen(lift(x))
_is_gen(x::MPolyDecRingElem) = is_gen(forget_grading(x))
_is_gen(x::MPolyRingElem) = is_gen(x)
# default method; overwrite if you want this to work for your rings.
_is_gen(x::NCRingElem) = false

# In case there is a type of ring elements for which hashing is correctly implemented 
# and does not throw an error, this gives the opportunity to overwrite the `allunique`
# to be used within the constructor for maps. 
function _allunique(lst::Vector{T}) where {T<:MPolyRingElem}
  return allunique(lst)
end

# We have a lot of rings which do/can not implement correct hashing. 
# So we make the following the default.
function _allunique(lst::Vector{T}) where {T<:RingElem}
  return all(!(x in lst[i+1:end]) for (i, x) in enumerate(lst))
end

function _allunique(lst::Vector{T}) where {T<:MPolyQuoRingElem}
  rep_list = lift.(lst)
  return _allunique(rep_list)
end

function _build_poly(u::MPolyRingElem, indices::Vector{Int}, S::MPolyRing)
  kk = coefficient_ring(S)
  r = ngens(S)
  ctx = MPolyBuildCtx(S)
  for (c, e) in zip(AbstractAlgebra.coefficients(u), AbstractAlgebra.exponent_vectors(u))
    ee = [0 for _ in 1:r]
    for (i, k) in enumerate(e)
      ee[indices[i]] = k
    end
    push_term!(ctx, kk(c), ee)
  end
  return finish(ctx)
end

function _evaluate_plain(F::MPolyAnyMap{<:MPolyRing, <:MPolyRing}, u)
  fl, var_ind = _maps_variables_to_variables(F)
  if fl
    return _build_poly(u, var_ind, codomain(F))
  end
  return evaluate(u, F.img_gens)
end

function _evaluate_plain(F::MPolyAnyMap{<: MPolyRing}, u)
  return evaluate(u, F.img_gens)
end

# See the comment in MPolyQuo.jl
function _evaluate_plain(F::MPolyAnyMap{<:MPolyRing, <:MPolyQuoRing}, u)
  A = codomain(F)
  R = base_ring(A)
  fl, var_ind = _maps_variables_to_variables(F)
  if fl
    return A(_build_poly(u, var_ind, R))
  end
  v = evaluate(lift(u), lift.(_images(F)))
  return simplify(A(v))
end

# The following assumes `p` to be in `S[x₁,…,xₙ]` where `S` is the 
# actual codomain of the map.
function _evaluate_with_build_ctx(
    p::MPolyRingElem, ind::Vector{Int}, cod_ring::MPolyRing
  )
  @assert cod_ring === coefficient_ring(parent(p))
  r = ngens(cod_ring)
  kk = coefficient_ring(cod_ring)
  ctx = MPolyBuildCtx(cod_ring)
  for (q, e) in zip(AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
    ee = [0 for _ in 1:r]
    for (i, k) in enumerate(e)
      ee[ind[i]] = k
    end
    for (c, d) in zip(AbstractAlgebra.coefficients(q), AbstractAlgebra.exponent_vectors(q))
      push_term!(ctx, kk(c), ee+d)
    end
  end
  return finish(ctx)
end

function _evaluate_general(F::MPolyAnyMap{<:MPolyRing, <:MPolyRing}, u)
  if domain(F) === codomain(F) && coefficient_map(F) === nothing
    return evaluate(map_coefficients(coefficient_map(F), u,
                                     parent = domain(F)), F.img_gens)
  else
    S = temp_ring(F)
    fl, var_ind = _maps_variables_to_variables(F)
    if S !== nothing
      if !fl || coefficient_ring(S) !== codomain(F)
        return evaluate(map_coefficients(coefficient_map(F), u,
                                         parent = S), F.img_gens)
      else
        tmp_poly = map_coefficients(coefficient_map(F), u, parent = S)
        return _evaluate_with_build_ctx(tmp_poly, var_ind, codomain(F))
      end
    else
      if !fl
        return evaluate(map_coefficients(coefficient_map(F), u), F.img_gens)
      else
        # For the case where we can recycle the method above, do so.
        tmp_poly = map_coefficients(coefficient_map(F), u)
        coefficient_ring(parent(tmp_poly)) === codomain(F) && return _evaluate_with_build_ctx(
                   tmp_poly,
                   var_ind,
                   codomain(F)
                 )
        # Otherwise default to the standard evaluation for the time being.
        return evaluate(tmp_poly, F.img_gens)
      end
    end
  end
end

function _evaluate_general(F::MPolyAnyMap{<: MPolyRing}, u)
  if domain(F) === codomain(F) && coefficient_map(F) === nothing
    return evaluate(map_coefficients(coefficient_map(F), u,
                                     parent = domain(F)), F.img_gens)
  else
    S = temp_ring(F)
    if S !== nothing
      return evaluate(map_coefficients(coefficient_map(F), u,
                                       parent = S), F.img_gens)
    else
      return evaluate(map_coefficients(coefficient_map(F), u), F.img_gens)
    end
  end
end

# one more intermediate function

function _evaluate_help(F::MPolyAnyMap{<: MPolyRing, <: Any, Nothing}, g)
  return _evaluate_plain(F, g)
end

function _evaluate_help(F::MPolyAnyMap{<: MPolyRing}, g)
  return _evaluate_general(F, g)
end

function (F::MPolyAnyMap{<: MPolyRing})(g)
  parent(g) === domain(F) || return F(domain(F)(g))
  if g isa elem_type(domain(F))
    if coefficient_map(F) === nothing
      return _evaluate_plain(F, g)
    else 
      return _evaluate_general(F, g)
    end
  else 
    gg = domain(F)(g)
    @assert parent(gg) === domain(F)
    return F(gg)
  end
end
