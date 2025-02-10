###############################################################################
# more functions

# computes the maximum `p` such that `beta - p*alpha` is still a root
# beta is assumed to be a root
function _root_string_length_down(alpha::RootSpaceElem, beta::RootSpaceElem)
  p = 0
  beta_sub_p_alpha = beta - alpha
  while is_root(beta_sub_p_alpha)
    p += 1
    beta_sub_p_alpha = sub!(beta_sub_p_alpha, alpha)
  end
  return p
end

@doc raw"""
    dim_of_simple_module([T = Int], R::RootSystem, hw::WeightLatticeElem) -> T
    dim_of_simple_module([T = Int], R::RootSystem, hw::Vector{<:IntegerUnion}) -> T

Compute the dimension of the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw` using Weyl's dimension formula.
The return value is of type `T`.

# Examples
```jldoctest
julia> R = root_system(:B, 2);

julia> dim_of_simple_module(R, [1, 0])
5
```
"""
function dim_of_simple_module(R::RootSystem, hw::WeightLatticeElem)
  return dim_of_simple_module(Int, R, hw)
end

function dim_of_simple_module(T::Type, R::RootSystem, hw::WeightLatticeElem)
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"
  rho = weyl_vector(R)
  hw_plus_rho = hw + rho
  num = one(ZZ)
  den = one(ZZ)
  for alpha in positive_roots(R)
    num *= ZZ(dot(hw_plus_rho, alpha))
    den *= ZZ(dot(rho, alpha))
  end
  return T(div(num, den))
end

function dim_of_simple_module(T::Type, R::RootSystem, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(T, R, WeightLatticeElem(R, hw))
end

function dim_of_simple_module(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(R, WeightLatticeElem(R, hw))
end

@doc raw"""
    dominant_weights(R::RootSystem, hw::WeightLatticeElem) -> Vector{WeightLatticeElem}
    dominant_weights(R::RootSystem, hw::Vector{<:IntegerUnion}) -> Vector{WeightLatticeElem}

Computes the dominant weights occurring in the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw`,
sorted ascendingly by the total height of roots needed to reach them from `hw`.

See [MP82](@cite) for details and the implemented algorithm.

# Examples
```jldoctest
julia> R = root_system(:B, 3);

julia> dominant_weights(R, [3, 0, 1])
7-element Vector{WeightLatticeElem}:
 3*w_1 + w_3
 w_1 + w_2 + w_3
 2*w_1 + w_3
 3*w_3
 w_2 + w_3
 w_1 + w_3
 w_3
```
"""
function dominant_weights(R::RootSystem, hw::WeightLatticeElem)
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"

  pos_roots = positive_roots(R)
  pos_roots_w = WeightLatticeElem.(positive_roots(R))

  ws_with_level = Dict(hw => 0)
  todo = [hw]
  while !isempty(todo)
    new_todo = empty(todo)
    for w in todo
      for (alpha, alpha_w) in zip(pos_roots, pos_roots_w)
        w_sub_alpha = w - alpha_w
        if is_dominant(w_sub_alpha) && !haskey(ws_with_level, w_sub_alpha)
          push!(new_todo, w_sub_alpha)
          push!(ws_with_level, w_sub_alpha => ws_with_level[w] + Int(height(alpha)))
        end
      end
    end
    todo = new_todo
  end
  return first.(sort!(collect(ws_with_level); by=last)) # order by level needed for dominant_character
end

function dominant_weights(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_weights(R, WeightLatticeElem(R, hw))
end

function _action_matrices_on_weights(W::WeylGroup)
  R = root_system(W)
  return map(1:rank(R)) do i
    x = gen(W, i)
    matrix(
      ZZ, reduce(vcat, coefficients(fundamental_weight(R, j) * x) for j in 1:rank(R))
    )
  end
end

@doc raw"""
    dominant_character([T = Int], R::RootSystem, hw::WeightLatticeElem) -> Dict{WeightLatticeElem, T}
    dominant_character([T = Int], R::RootSystem, hw::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}

Computes the dominant weights occurring in the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw`, together with their multiplicities.

This function uses an optimized version of the Freudenthal formula, see [MP82](@cite) for details.

# Examples
```jldoctest
julia> R = root_system(:B, 3);

julia> dominant_character(R, [2, 0, 1])
Dict{WeightLatticeElem, Int64} with 4 entries:
  w_2 + w_3   => 1
  2*w_1 + w_3 => 1
  w_1 + w_3   => 3
  w_3         => 6
```
"""
function dominant_character(R::RootSystem, hw::WeightLatticeElem)
  return dominant_character(Int, R, hw)
end

function dominant_character(T::DataType, R::RootSystem, hw::WeightLatticeElem)
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"
  W = weyl_group(R)
  rho = weyl_vector(R)
  hw_plus_rho = hw + rho
  dot_how_plus_rho = dot(hw_plus_rho, hw_plus_rho)

  pos_roots = positive_roots(R)
  pos_roots_w = WeightLatticeElem.(positive_roots(R))
  pos_roots_w_coeffs = coefficients.(pos_roots_w)

  char = Dict(hw => T(1))

  todo = dominant_weights(R, hw)
  all_orbs = Dict{Vector{Int},Vector{Tuple{WeightLatticeElem,Int}}}()
  action_matrices_on_weights = _action_matrices_on_weights(W)

  for w in Iterators.drop(todo, 1) # drop hw as its multiplicity is already known
    stab_inds = [i for (i, ci) in enumerate(coefficients(w)) if iszero(ci)]
    orbs = get!(all_orbs, stab_inds) do
      gens = action_matrices_on_weights[stab_inds]
      push!(gens, -identity_matrix(ZZ, rank(R)))
      G = matrix_group(gens)
      O = orbits(gset(G, *, [pos_roots_w_coeffs; -pos_roots_w_coeffs]))
      [
        (
          WeightLatticeElem(
            R,
            first(intersect(elements(o), pos_roots_w_coeffs)),
          ),
          length(o),
        ) for o in O
      ]
    end

    accum = sum(
      data -> begin
        rep, len = data
        accum2 = 0
        w_plus_i_rep = w + rep
        while true
          w_plus_i_rep_conj = conjugate_dominant_weight(w_plus_i_rep)
          haskey(char, w_plus_i_rep_conj) || break
          accum2 += char[w_plus_i_rep_conj] * dot(w_plus_i_rep, rep)
          add!(w_plus_i_rep, rep)
        end
        len * accum2
      end, orbs; init=zero(QQ))
    if !iszero(accum)
      w_plus_rho = w + rho
      denom = dot_how_plus_rho - dot(w_plus_rho, w_plus_rho)
      if !iszero(denom)
        char[w] = T(ZZ(div(accum, denom)))
      end
    end
  end
  return char
end

function dominant_character(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_character(R, WeightLatticeElem(R, hw))
end

function dominant_character(T::DataType, R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_character(T, R, WeightLatticeElem(R, hw))
end

@doc raw"""
    character([T = Int], R::RootSystem, hw::WeightLatticeElem) -> Dict{WeightLatticeElem, T}
    character([T = Int], R::RootSystem, hw::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}

Computes all weights occurring in the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw`, together with their multiplicities.
This is achieved by acting with the Weyl group on the [`dominant_character`](@ref dominant_character(::RootSystem, ::WeightLatticeElem)).

# Examples
```jldoctest
julia> R = root_system(:B, 3);

julia> character(R, [0, 0, 1])
Dict{WeightLatticeElem, Int64} with 8 entries:
  -w_1 + w_2 - w_3 => 1
  w_1 - w_3        => 1
  -w_2 + w_3       => 1
  w_3              => 1
  w_2 - w_3        => 1
  -w_3             => 1
  w_1 - w_2 + w_3  => 1
  -w_1 + w_3       => 1
```
"""
function character(R::RootSystem, hw::WeightLatticeElem)
  return character(Int, R, hw)
end

function character(T::DataType, R::RootSystem, hw::WeightLatticeElem)
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"
  dom_char = dominant_character(T, R, hw)
  char = Dict{WeightLatticeElem,T}()

  for (w, m) in dom_char
    for w_conj in weyl_orbit(w)
      push!(char, w_conj => m)
    end
  end

  return char
end

function character(R::RootSystem, hw::Vector{<:IntegerUnion})
  return character(R, WeightLatticeElem(R, hw))
end

function character(T::DataType, R::RootSystem, hw::Vector{<:IntegerUnion})
  return character(T, R, WeightLatticeElem(R, hw))
end

@doc raw"""
    tensor_product_decomposition(R::RootSystem, hw1::WeightLatticeElem, hw2::WeightLatticeElem) -> MSet{Vector{Int}}
    tensor_product_decomposition(R::RootSystem, hw1::Vector{<:IntegerUnion}, hw2::Vector{<:IntegerUnion}) -> MSet{Vector{Int}}

Computes the decomposition of the tensor product of the simple modules of the Lie algebra defined by the root system `R`
with highest weights `hw1` and `hw2` into simple modules with their multiplicities.
This function uses Klymik's formula.

The return type may change in the future.

# Examples
```jldoctest
julia> R = root_system(:B, 2);

julia> tensor_product_decomposition(R, [1, 0], [0, 1])
MSet{Vector{Int64}} with 2 elements:
  [1, 1]
  [0, 1]

julia> tensor_product_decomposition(R, [1, 1], [1, 1])
MSet{Vector{Int64}} with 10 elements:
  [0, 0]
  [0, 4]
  [0, 2] : 2
  [2, 0]
  [1, 0]
  [2, 2]
  [3, 0]
  [1, 2] : 2
```
"""
function tensor_product_decomposition(
  R::RootSystem, hw1::WeightLatticeElem, hw2::WeightLatticeElem
)
  @req root_system(hw1) === R "parent root system mismatch"
  @req root_system(hw2) === R "parent root system mismatch"
  @req is_dominant(hw1) "not a dominant weight"
  @req is_dominant(hw2) "not a dominant weight"

  rho = weyl_vector(R)
  hw2_plus_rho = hw2 + rho

  mults = multiset(WeightLatticeElem)
  for (w_, m) in dominant_character(R, hw1)
    for w in weyl_orbit(w_)
      add!(w, hw2_plus_rho)
      w_dom, x = conjugate_dominant_weight_with_elem!(w)
      if all(!iszero, coefficients(w_dom))
        sub!(w_dom, rho)
        coeff = m * (-1)^length(x)
        push!(mults, w_dom, coeff)
      end
    end
  end

  # return mults
  return multiset(
    Dict(Int.(_vec(coefficients(w))) => multiplicity(mults, w) for w in unique(mults))
  )
end

function tensor_product_decomposition(
  R::RootSystem, hw1::Vector{<:IntegerUnion}, hw2::Vector{<:IntegerUnion}
)
  return tensor_product_decomposition(
    R, WeightLatticeElem(R, hw1), WeightLatticeElem(R, hw2)
  )
end

###############################################################################
# demazures character formula
function _demazure_operator(r::RootSpaceElem, w::WeightLatticeElem)
  @req is_positive_root(r) "r is not a positive root"

  d = 2 * dot(w, r)//dot(r, r)
  list_of_occuring_weights = WeightLatticeElem[]

  refl = reflect(w, r)

  wlelem_r = WeightLatticeElem(r)
  if d > -1
    while w != refl
      push!(list_of_occuring_weights, w)
      w -= wlelem_r
    end
    push!(list_of_occuring_weights, w)
    return 1, list_of_occuring_weights
  elseif d < -1
    w += wlelem_r
    push!(list_of_occuring_weights, w)
    while w != refl - wlelem_r
      w += wlelem_r
      push!(list_of_occuring_weights, w)
    end
    return -1, list_of_occuring_weights
  else
    return 0, list_of_occuring_weights
  end
end

@doc raw"""
    demazure_operator(r::RootSpaceElem, w::WeightLatticeElem) -> Dict{WeightLatticeElem,<:IntegerUnion}
    demazure_operator(r::RootSpaceElem, groupringelem::Dict{WeightLatticeElem,<:IntegerUnion}) -> Dict{WeightLatticeElem,<:IntegerUnion}

Computes the action of the Demazure operator associated to the positive root `r` on the given element of the group ring $\mathbb{Z}[P]$.

If a single Weight lattice element `w` is supplied, this is interpreted as `Dict(w => 1)`.

# Examples
```jldoctest
julia> R = root_system(:A, 3);

julia> pos_r = positive_root(R, 4)
a_1 + a_2

julia> w = fundamental_weight(R, 1)
w_1

julia> demazure_operator(pos_r, w)
Dict{WeightLatticeElem, Int64} with 2 entries:
  -w_2 + w_3 => 1
  w_1        => 1
```
"""
function demazure_operator(
  r::RootSpaceElem, groupringelem::Dict{WeightLatticeElem,<:IntegerUnion}
)
  dict = empty(groupringelem)
  for (w, dim) in groupringelem
    sign, weights = _demazure_operator(r, w)
    for w_ in weights
      val = get(dict, w_, 0) + sign * dim
      if is_zero(val)
        delete!(dict, w_)
      else
        dict[w_] = val
      end
    end
  end
  return dict
end

function demazure_operator(r::RootSpaceElem, w::WeightLatticeElem)
  return demazure_operator(r, Dict(w => 1))
end

@doc raw"""
    demazure_character([T = Int], R::RootSystem, w::WeightLatticeElem, x::WeylGroupElem) -> Dict{WeightLatticeElem, T}
    demazure_character([T = Int], R::RootSystem, w::Vector{<:IntegerUnion}, x::WeylGroupElem) -> Dict{WeightLatticeElem, T}
    demazure_character([T = Int], R::RootSystem, w::WeightLatticeElem, reduced_expr::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}
    demazure_character([T = Int], R::RootSystem, w::Vector{<:IntegerUnion}, reduced_expr::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}

Computes all weights occurring in the Demazure module of the Lie algebra defined by the root system `R`
with extremal weight `w * x`, together with their multiplicities.

Instead of a Weyl group element `x`, a reduced expression for `x` can be supplied.
This function may return arbitrary results if the provided expression is not reduced.

# Examples
```jldoctest
julia> R = root_system(:B, 3);

julia> demazure_character(R, [0, 1, 0], [1, 2, 3])
Dict{WeightLatticeElem, Int64} with 4 entries:
  w_1 + w_2 - 2*w_3 => 1
  w_1               => 1
  w_1 - w_2 + 2*w_3 => 1
  w_2               => 1
```
"""
function demazure_character(R::RootSystem, w::WeightLatticeElem, x::WeylGroupElem)
  @req root_system(parent(x)) === R "parent root system mismatch"
  return demazure_character(R, w, word(x))
end

function demazure_character(
  T::DataType, R::RootSystem, w::WeightLatticeElem, x::WeylGroupElem
)
  @req root_system(parent(x)) === R "parent root system mismatch"
  return demazure_character(T, R, w, word(x))
end

function demazure_character(
  R::RootSystem, w::WeightLatticeElem, reduced_expression::Vector{<:IntegerUnion}
)
  return demazure_character(Int, R, w, reduced_expression)
end

function demazure_character(
  T::DataType,
  R::RootSystem,
  w::WeightLatticeElem,
  reduced_expression::Vector{<:IntegerUnion},
)
  @req root_system(w) === R "parent root system mismatch"
  @req is_dominant(w) "not a dominant weight"
  char = Dict{WeightLatticeElem,T}(w => T(1))
  for i in reduced_expression
    char = demazure_operator(simple_root(root_system(w), Int(i)), char)
  end
  return char
end

function demazure_character(R::RootSystem, w::Vector{<:IntegerUnion}, x::WeylGroupElem)
  return demazure_character(R, WeightLatticeElem(R, w), x)
end

function demazure_character(
  T::DataType, R::RootSystem, w::Vector{<:IntegerUnion}, x::WeylGroupElem
)
  return demazure_character(T, R, WeightLatticeElem(R, w), x)
end

function demazure_character(
  R::RootSystem, w::Vector{<:IntegerUnion}, reduced_expression::Vector{<:IntegerUnion}
)
  return demazure_character(R, WeightLatticeElem(R, w), reduced_expression)
end

function demazure_character(
  T::DataType,
  R::RootSystem,
  w::Vector{<:IntegerUnion},
  reduced_expression::Vector{<:IntegerUnion},
)
  return demazure_character(T, R, WeightLatticeElem(R, w), reduced_expression)
end
