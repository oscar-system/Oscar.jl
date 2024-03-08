# Function to check if an element of an algebraic extension of Q is a root of unity
# and if so determine its order.
#
# Ulrich Thiel, 2024


function is_root_of_unity_with_data(x::QQFieldElem)

  if x == 0
    return false, 0
  elseif x == 1
    return true, 1
  elseif x == -1
    return true, 2
  else
    return false, 0
  end

end

function is_root_of_unity_with_data(x::QQAlgFieldElem)

  if x == 0
    return false, 0
  elseif x == 1
    return true, 1
  elseif x == -1
    return true, 2
  else
    p = absolute_minpoly(x)
    if !is_cyclotomic_polynomial(p)
      return false, 0
    end
    candidates = euler_phi_inv(degree(p))
    for n in sort(candidates)
      if x^n == 1
        return true, n
      end
    end
  end
end

function is_root_of_unity(x::QQAlgFieldElem)
  b,n = is_root_of_unity_with_data(x)
  return b
end

function order(x::QQAlgFieldElem)
  b,n = is_root_of_unity_with_data(x)
  if !b
    throw(ArgumentError("Element is not of finite order."))
  else
    return n
  end
end
