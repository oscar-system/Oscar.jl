function (h::QuantumGroupHom)(x::QuantumGroupElem)
  @req parent(x) == h.domain "parent mismatch"

  val = zero(h.codomain.algebra)
  t = one(U.alg)
  barred = zero(coefficient_ring(U))

  exp = Memory{Int}(undef, ngens(h.domain.algebra))
  for i in 1:length(x)
    exponent_vector!(exp, x.elem, i)
    for j in 1:length(exp)
      for _ in 1:exp[j]
        t = mul!(t, h.img[j])
      end
    end

    bar!(barred, coeff(x.elem, i))
    val = addmul!(val, t, barred)
    t = one!(t)
  end
end

@doc raw"""
    bar_automorphism(U::QuantumGroup) -> QuantumGroupHom
    
Return the bar automorphism of the quantum group `U`.
"""
function bar_automorphism(U::QuantumGroup)
  return _bar_involution(U)
end

function _bar_involution(U::QuantumGroup)
  nsim = number_of_simple_roots(root_system(U))
  npos = number_of_positive_roots(root_system(U))
  cvx = U.cvx

  img = Vector{PBWAlgebraElem{QuantumFieldElem}}(undef, npos)
  for i in 1:nsim
    img[cvx[i]] = gen(U.algebra, cvx[i])
  end

  # we construct the images of the PBW generators inductively by height
  R = root_system(U)
  refl = weyl_group(root_system(U)).refl
  for m in (nsim + 1):npos
    n = 0
    s = 0
    # find positive root with smaller height
    for i in 1:nsim
      n = Int(refl[i, m])
      if n < m
        s = i
        break
      end
    end

    pow = Int(height(positive_root(R, m)) - height(positive_root(R, n)))
    if cvx[s] < cvx[n]
      rel = gen(U.algebra, cvx[n]) * gen(U.algebra, cvx[s])^pow
      img[cvx[m]] = img[cvx[n]] * img[cvx[s]]^pow
    else
      rel = gen(U.algebra, cvx[s])^pow * gen(U.algebra, cvx[n])
      img[cvx[m]] = img[cvx[s]]^pow * img[cvx[n]]
    end

    b = one(U.algebra)
    barred = zero(coefficient_ring(U))
    c = zero(coefficient_ring(U))

    exp = Memory{Int}(undef, ngens(U.algebra))
    for i in 1:length(rel)
      exponent_vector!(exp, rel, i)
      if exp[cvx[m]] != 0
        c = inv!(bar!(c, coeff(rel, i)))
        continue
      end

      for n in eachindex(exp)
        for _ in 1:exp[n]
          b = mul!(b, img[n])
        end
      end

      bar!(barred, coeff(rel, i))
      img[cvx[m]] = submul!(img[cvx[m]], b, barred)
      b = one!(b)
    end
    img[cvx[m]] = mul!(img[cvx[m]], c)
  end

  return function (x::QuantumGroupElem)
    @req parent(x) == U "parent mismatch"

    val = zero(U.alg)
    t = one(U.alg)
    barred = zero(coefficient_ring(U))
    for term in terms(x.elem)
      exp = leading_exponent_vector(term)
      for i in eachindex(exp)
        for _ in 1:exp[i]
          t = mul!(t, img[i])
        end
      end

      bar!(barred, leading_coefficient(term))
      val = addmul!(val, t, barred)
      t = one!(t)
    end

    return QuantumGroupElem(U, val)
  end
end

function braid_automorphism(U::QuantumGroup, i::Int)
  error("not implemented")
end

###############################################################################
#
#   Internal helpers
#
###############################################################################

function _image(imgs::Vector{})
  # we construct the images of the PBW generators inductively by height
  for m in (nsim + 1):npos
    n = 0
    s = 0
    # find positive root with smaller height
    for i in 1:nsim
      n = Int(refl[i, m])
      if n < m
        s = i
        break
      end
    end

    pow = Int(height(positive_root(R, m)) - height(positive_root(R, n)))
    if cvx[s] < cvx[n]
      rel = gen(U.algebra, cvx[n]) * gen(U.algebra, cvx[s])^pow
      img[cvx[m]] = img[cvx[n]] * img[cvx[s]]^pow
    else
      rel = gen(U.algebra, cvx[s])^pow * gen(U.algebra, cvx[n])
      img[cvx[m]] = img[cvx[s]]^pow * img[cvx[n]]
    end

    b = one(U.algebra)
    barred = zero(coefficient_ring(U))
    coeff = zero(coefficient_ring(U))

    exp = Memory{Int}(undef, ngens(U.algebra))
    for i in 1:length(rel)
      exponent_vector!(exp, rel, i)
      if exp[cvx[m]] != 0
        coeff = inv!(bar!(coeff, coeff(rel, i)))
        continue
      end

      for n in eachindex(exp)
        for _ in 1:exp[n]
          b = mul!(b, img[n])
        end
      end

      bar!(barred, coeff(rel, i))
      img[cvx[m]] = submul!(img[cvx[m]], b, barred)
      b = one!(b)
    end
    img[cvx[m]] = mul!(img[cvx[m]], coeff)
  end
end
