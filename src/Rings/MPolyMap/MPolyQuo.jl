################################################################################
#
#  Domains of type MPolyQuo
#
################################################################################

function _check_imgs_quo(R, S::NCRing, imgs)
  n = length(imgs)
  for i in 1:n, j in 1:(i - 1)
    @req imgs[i] * imgs[j] == imgs[j] * imgs[i] "Images $i and $j do not commute"
  end
  gensI = gens(modulus(R))
  for g in gensI
    @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
  end
  return nothing
end

# no check for commutative codomains
function _check_imgs_quo(R, S::Ring, imgs) 
  gensI = gens(modulus(R))
  for g in gensI
    @req iszero(evaluate(g, imgs)) "Morphism is not well-defined"
  end
end

function hom(R::MPolyQuo, S::NCRing, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  if check
    _check_imgs_quo(R, S, imgs)
    _check_homo(S, imgs)
  end
  return MPolyAnyMap(R, S, nothing, copy(imgs)) # copy because of #655
end

function hom(R::MPolyQuo, S::NCRing, coeff_map, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  if check
    _check_imgs_quo(R, S, imgs)
    _check_homo(S, imgs)
  end

  return MPolyAnyMap(R, S, coeff_map, copy(imgs)) # copy because of #655
end

################################################################################
#
#  Evaluation
#
################################################################################

function _evaluate_plain(F::MPolyAnyMap{<: MPolyQuo}, u)
  return evaluate(lift(u), _images(F))
end

function _evaluate_general(F::MPolyAnyMap{<: MPolyQuo}, u)
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

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuo, <: Any, Nothing}, g)
  return _evaluate_plain(F, g)
end

function _evaluate_help(F::MPolyAnyMap{<: MPolyQuo}, g)
  return _evaluate_general(F, g)
end

# The two main evaluation methods
function (F::MPolyAnyMap{<: MPolyQuo})(g)
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
