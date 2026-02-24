################################################################################
#
#  Domains of type MPolyQuoRing
#
################################################################################

function _check_imgs_quo(R, S::NCRing, imgs, coeff_map = nothing)
  n = length(imgs)
  for i in 1:n, j in 1:(i - 1)
    @req imgs[i] * imgs[j] == imgs[j] * imgs[i] "Images $i and $j do not commute"
  end
  gensI = gens(modulus(R))
  for g in gensI
    if coeff_map !== nothing
      @req is_zero(evaluate(map_coefficients(coeff_map, g), imgs)) "Morphism is not well-defined"
    else
      @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
    end
  end
  return nothing
end

# no check for commutative codomains
function _check_imgs_quo(R, S::Ring, imgs, coeff_map = nothing)
  gensI = gens(modulus(R))
  for g in gensI
    if coeff_map !== nothing
      @req is_zero(evaluate(map_coefficients(coeff_map, g), imgs)) "Morphism is not well-defined"
    else
      @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
    end
  end
end


@doc raw"""
    hom(A::MPolyQuoRing, S::NCRing, coeff_map, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing)

    hom(A::MPolyQuoRing, S::NCRing, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing)

    hom(::Type{K}, A::MPolyQuoRing, S::NCRing, coeff_map, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing) where {K}

    hom(::Type{K}, A::MPolyQuoRing, S::NCRing, images::Vector;
        check::Bool = true,
        grading_group_hom = nothing,
        degree_shift = nothing) where {K}

Create the ring homomorphism `f : A -> S` sending the `i`-th generator of `A`
to `images[i]`. If `coeff_map` is supplied, coefficients are mapped using
`coeff_map`; otherwise Oscar uses `nothing`.

When `A = R/I` is a quotient, the defining relations of `I` must map to zero
in `S` (this is checked when `check=true`).

Grading behaviour and map kinds are as for `hom(::MPolyRing, ...)`:
`HomUngraded`, `HomGraded`, `HomGradedToUngraded`, `HomUngradedToGraded`.
Typed constructors `hom(::Type{K}, ...)` allow forcing the kind.

For `HomGraded` maps, there is a group homomorphism `phi` and a shift
`deg_shift` such that

    degree(images[i]) == phi(degree(gen(A,i))) + deg_shift

and each `images[i]` is homogeneous.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A, _ = quo(R, ideal(R, [y-x^2, z-x^3]));

julia> S, (s, t) = polynomial_ring(QQ, [:s, :t]);

julia> F = hom(A, S, [s, s^2, s^3])
Ring homomorphism
  from quotient of multivariate polynomial ring by ideal (-x^2 + y, -x^3 + z)
  to multivariate polynomial ring in 2 variables over QQ
defined by
  x -> s
  y -> s^2
  z -> s^3
```
"""
function hom(::Type{K}, R::MPolyQuoRing, S::NCRing, coeff_map, images::Vector;
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
    _check_imgs_quo(R, S, imgs, coeff_map)
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

function hom(::Type{K}, R::MPolyQuoRing, S::NCRing, images::Vector;
             check::Bool = true,
             grading_group_hom = nothing,
             degree_shift = nothing) where {K <: MPolyHomKind}
  return hom(K, R, S, nothing, images;
             check = check,
             grading_group_hom = grading_group_hom,
             degree_shift = degree_shift)
end

# Untyped constructors
function hom(R::MPolyQuoRing, S::NCRing, coeff_map, images::Vector;
             check::Bool = true,
             grading_group_hom = nothing,
             degree_shift = nothing)
  K = _kind_for(R, S)
  return hom(K, R, S, coeff_map, images;
             check = check,
             grading_group_hom = grading_group_hom,
             degree_shift = degree_shift)
end

function hom(R::MPolyQuoRing, S::NCRing, images::Vector;
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
#  Evaluation
#
################################################################################

function _evaluate_plain(F::MPolyAnyMap{<: MPolyQuoRing}, u)
  return evaluate(lift(u), _images(F))
end

function _evaluate_plain(F::MPolyAnyMap{<:MPolyQuoRing, <:MPolyQuoRing}, u)
  # This workaround deals with the fact that arithmetic in quotient rings is TERRIBLY slow.
  # All the simplify calls make it unusable in this case, probably due to the fact that
  # setting `is_reduced` flags does not pay off in such iterative procedures.
  A = codomain(F)
  # Shortcut needed because the `evaluate` code can not deal with empty lists
  is_empty(_images(F)) && return codomain(F)(_constant_coefficient(u)) 
  v = evaluate(lift(u), lift.(_images(F)))
  return simplify(A(v))
end

_constant_coefficient(f::MPolyRingElem) = constant_coefficient(f)
_constant_coefficient(f::MPolyQuoRingElem) = constant_coefficient(lift(f))

function _evaluate_general(F::MPolyAnyMap{<: MPolyQuoRing}, u)
  if domain(F) === codomain(F) && coefficient_map(F) === nothing
    return evaluate(map_coefficients(coefficient_map(F), lift(u),
                                     parent = domain(F)), F.img_gens)
  else
    S = temp_ring(F)
    if S !== nothing
      return evaluate(map_coefficients(coefficient_map(F), lift(u),
                                                        parent = S), F.img_gens)
    else
      return evaluate(map_coefficients(coefficient_map(F), lift(u)), F.img_gens)
    end
  end
end

# one more intermediate function

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuoRing, <: Any, Nothing}, g)
  return _evaluate_plain(F, g)
end

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuoRing}, g)
  return _evaluate_general(F, g)
end

# The two main evaluation methods
function (F::MPolyAnyMap{<: MPolyQuoRing})(g)
  if g isa elem_type(domain(F))
    @req parent(g) === domain(F) "Element not in domain"
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
