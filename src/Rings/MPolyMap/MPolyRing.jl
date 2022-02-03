################################################################################
#
#  Domains of type MPolyRing
#
################################################################################

################################################################################
#
#  Constructos
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

function hom(R::MPolyRing, S::NCRing, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  if check
    _check_imgs(S, imgs)
    _check_homo(S, imgs) # defined in MPolyAnyMap.jl
  end
  return MPolyAnyMap(R, S, nothing, copy(imgs)) # copy because of #655
end

# Since coeff_map is at the third position, we wannot use a variable
# with default argument
function hom(R::MPolyRing, S::Ring, coeff_map, images::Vector; check::Bool = true)
  n = ngens(R)
  @req n == length(images) "Number of images must be $n"
  # Now coerce into S or throw an error if not possible
  imgs = _coerce(S, images)
  if check
    _check_imgs(S, imgs)
    _check_homo(S, imgs) # defined in MPolyAnyMap.jl
  end

  return MPolyAnyMap(R, S, coeff_map, copy(imgs)) # copy because of #655
end

################################################################################
#
#  Evaluation functions
#
################################################################################

function _evaluate_plain(F::MPolyAnyMap{<: MPolyRing}, u)
  return evaluate(u, F.img_gens)
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
