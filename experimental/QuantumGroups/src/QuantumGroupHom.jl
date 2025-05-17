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

function bar_involution(U::QuantumGroup)
  if isdefined(U, :bar_involution)
    return U.bar_involution
  end

  cvx = U.cvx

  img = zeros(U.algebra, ngens(U.algebra))
  for i in 1:nsim
    img[cvx[i]] = add!(img[cvx[i]], gen(U.algebra, cvx[i]))
  end

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

  return nothing
end

function braid_automorphism(U::QuantumGroup, i::Int)
  error("not implemented")
end
