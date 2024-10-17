module KashiwaraCrystal

function e end

function e! end

function f end

function f! end

function eps end

function phi end

function weight end

end

function deepcopy(::KashiwaraCrystalElem)
  error("not implemented")
end

function KashiwaraCrystal.e(::KashiwaraCrystalElem, ::Int)
  error("not implemented")
end

function KashiwaraCrystal.f(::KashiwaraCrystalElem, ::Int)
  error("not implemented")  
end

@doc raw"""
    KashiwaraCrystal.eps(b::KashiwaraCrystalElem, i::Int) -> Int

Return $\epsilon_i(b)$.
"""
function KashiwaraCrystal.eps(::KashiwaraCrystalElem, ::Int)
  error("not implemented")  
end

@doc raw"""
    KashiwaraCrystal.phi(b::KashiwaraCrystalElem, i::Int) -> Int

Return $\varphi_i(b)$.
"""
function KashiwaraCrystal.phi(::KashiwaraCrystalElem, ::Int)
  error("not implemented")  
end

@doc raw"""
    weight(b::KashiwaraCrystalElem) -> WeightLatticeElem

Return the weight of `b`.
"""
function weight(::KashiwaraCrystalElem)
  error("not implemented")
end

# KashiwaraCrystal optional methods

@doc raw"""
    is_upper_normal(B::KashiwaraCrystal) -> Bool

Return `true` if `B` is upper normal, and `false` if `B` is not upper normal or it is unknown.
In the case of `false` consult the specific implemenation for details.
"""
function is_upper_normal(::KashiwaraCrystal)
  return true
end

@doc raw"""
    is_lower_normal(B::KashiwaraCrystal) -> Bool

Return `true` if `B` is lower normal, and `false` if `B` is not lower normal or it is unknown.
In the case of `false` consult the specific implemenation for details.
"""
function is_lower_normal(::KashiwaraCrystal)
  return true
end

@doc raw"""
    is_normal(B::KashiwaraCrystal) -> Bool

Return `true` if `B` is normal, and `false` if `B` is not normal or it is unknown,
whether `B` is upper normal or lower normal.
In the case of `false` consult the specific implemenation for details.
"""
function is_normal(B::KashiwaraCrystal)
  return is_upper_normal(B) && is_lower_normal(B)
end

# KashiwaraCrystalElem provided methods

# e_i variants

@doc raw"""
    KashiwaraCrystal.e(b::KashiwaraCrystalElem, i::Int) -> KashiwaraCrystalElem

Return the result of applying $\tilde e_i$ to `b`.
"""
function KashiwaraCrystal.e(b::KashiwaraCrystalElem, i::Int)
  return KashiwaraCrystal.e!(deepcopy(b), i)
end

@doc raw"""
    KashiwaraCrystal.e(b::KashiwaraCrystalElem, i::Int, n::Int) -> KashiwaraCrystalElem

Shorthand for `n`fold appliction of `KashiwaraCrystal.e(b, i)`.
"""
function KashiwaraCrystal.e(b::KashiwaraCrystalElem, i::Int, n::Int)
  return KashiwaraCrystal.e!(deepcopy(b), i, n)
end

@doc raw"""
    KashiwaraCrystal.e!(b::KashiwaraCrystalElem, i::Int, n::Int) -> KashiwaraCrystalElem
    
Shorthand for `n`fold appliction of `KashiwaraCrystal.e!(b, i)`.
"""
function KashiwaraCrystal.e!(b::KashiwaraCrystalElem, i::Int, n::Int)
  for _ in 1:n
    KashiwaraCrystal.e!(b, i)
  end
  return b
end

@doc raw"""
    KashiwaraCrystal.e!(b::KashiwaraCrystalElem, i::Int, n::Int) -> KashiwaraCrystalElem
    
Shorthand for `n`fold appliction of `KashiwaraCrystal.e!(b, i)`.
"""
function KashiwaraCrystal.e(b::KashiwaraCrystalElem, i::Vector{Int}, n::Vector{Int})
  KashiwaraCrystal.e!(deepcopy(b), i, n)
end

@doc raw"""
    KashiwaraCrystal.e!(b::KashiwaraCrystalElem, i::Vector{Int}, n::Vector{Int}) -> KashiwaraCrystalElem

Return the result of $\tilde e_{i_1} \dots \tilde e_{i_r}b$.
"""
function KashiwaraCrystal.e!(b::KashiwaraCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "i and n must have same length"

  for l in length(i):-1:1
    KashiwaraCrystal.e!(b, i[l], n[l])
  end
  return b
end

# f_i variants

@doc raw"""
    KashiwaraCrystal.f(b::KashiwaraCrystalElem, i::Int) -> KashiwaraCrystalElem

Return the result of applying $\tilde f_i$ to `b`.
"""
function KashiwaraCrystal.f(b::KashiwaraCrystalElem, i::Int)
  return KashiwaraCrystal.f!(deepcopy(b), i)
end

@doc raw"""
    KashiwaraCrystal.e(b::KashiwaraCrystalElem, i::Int, n::Int) -> KashiwaraCrystalElem

Shorthand for `n`fold appliction of `KashiwaraCrystal.f(b, i)`.
"""
function KashiwaraCrystal.f(b::KashiwaraCrystalElem, i::Int, n::Int)
  return KashiwaraCrystal.f!(deepcopy(b), i, n)
end

@doc raw"""
    KashiwaraCrystal.f!(b::KashiwaraCrystalElem, i::Int, n::Int) -> KashiwaraCrystalElem
    
Shorthand for `n`fold appliction of `KashiwaraCrystal.f!(b, i)`.
"""
function KashiwaraCrystal.f!(b::KashiwaraCrystalElem, i::Int, n::Int)
  for _ in 1:n
    KashiwaraCrystal.f!(b, i)
  end
  return b
end

@doc raw"""
    KashiwaraCrystal.f(b::KashiwaraCrystalElem, i::Vector{Int}, n::Vector{Int}) -> KashiwaraCrystalElem

"""
function KashiwaraCrystal.f(b::KashiwaraCrystalElem, i::Vector{Int}, n::Vector{Int})
  KashiwaraCrystal.f!(deepcopy(b), i, n)
end

@doc raw"""
    KashiwaraCrystal.f!(b::KashiwaraCrystalElem, i::Vector{Int}, n::Vector{Int}) -> KashiwaraCrystalElem

Return the result of $\tilde e_{i_1} \dots \tilde e_{i_r}b$.
"""
function KashiwaraCrystal.f!(b::KashiwaraCrystalElem, i::Vector{Int}, n::Vector{Int})
  @req length(i) == length(n) "i and n must have same length"

  for l in length(i):-1:1
    KashiwaraCrystal.f!(b, i[l], n[l])
  end
  return b
end
