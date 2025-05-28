###############################################################################
#
#   QuantumGroupHom
#
###############################################################################

function domain(hom::QuantumGroupHom)
  return hom.domain
end

function codomain(hom::QuantumGroupHom)
  return hom.codomain
end

function image!(z::QuantumGroupElem, hom::QuantumGroupHom, x::QuantumGroupElem)
  z.elem = image!(hom.hom, z.elem, x.elem)
  return z
end

function image(hom::QuantumGroupHom, x::QuantumGroupElem)
  @req parent(x) === hom.domain "parent mismatch"
  return image!(hom, zero(hom.codomain), x)
end

function (hom::QuantumGroupHom)(x::QuantumGroupElem)
  return image(hom, x)
end

###############################################################################
#
#   Printing
#
###############################################################################

function Base.show(io::IO, ::QuantumGroupHom)
  print(io, "Homomorphism of quantum groups")
end

###############################################################################
#
#   
#
###############################################################################

@doc raw"""
    bar_automorphism(U::QuantumGroup) -> QuantumGroupHom
    
Return the bar automorphism of the quantum group `U`.
"""
function bar_automorphism(U::QuantumGroup)
  return QuantumGroupHom(U, U, U.bar_automorphism)
end

function _bar_automorphism(A::PBWAlgebra{QuantumFieldElem}, R::RootSystem, cvx::Vector{Int})
  nsim = number_of_simple_roots(R)
  npos = number_of_positive_roots(R)

  img = Vector{PBWAlgebraElem{QuantumFieldElem}}(undef, npos)
  for i in 1:nsim
    img[cvx[i]] = gen(A, cvx[i])
  end

  refl = weyl_group(R).refl
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
      rel = gen(A, cvx[n]) * gen(A, cvx[s])^pow
      img[cvx[m]] = img[cvx[n]] * img[cvx[s]]^pow
    else
      rel = gen(A, cvx[s])^pow * gen(A, cvx[n])
      img[cvx[m]] = img[cvx[s]]^pow * img[cvx[n]]
    end

    b = one(A)
    barred = zero(coefficient_ring(A))
    c = zero(coefficient_ring(A))

    exp = Memory{Int}(undef, ngens(A))
    for i in 1:length(rel)
      exponent_vector!(exp, rel, i)
      if exp[cvx[m]] != 0
        c = inv!(image!(c, _BarAutomorphism, coeff(rel, i)))
        continue
      end

      for n in eachindex(exp)
        for _ in 1:exp[n]
          b = mul!(b, img[n])
        end
      end

      image!(barred, _BarAutomorphism, coeff(rel, i))
      img[cvx[m]] = submul!(img[cvx[m]], b, barred)
      b = one!(b)
    end

    img[cvx[m]] = mul!(img[cvx[m]], c)
  end

  return PBWAlgebraHom(A, A, _BarAutomorphism(), img)
end

function braid_automorphism(U::QuantumGroup, i::Int)
  error("not implemented")
end

###############################################################################
#
#   Internal helpers
#
###############################################################################
