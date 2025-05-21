###############################################################################
#
#   Representation theory for semisimple Lie algebras
#
###############################################################################

###############################################################################
#
#   Simple modules (via highest weight) of semisimple Lie algebras
#
###############################################################################

@doc raw"""
    simple_module(L::LieAlgebra{C}, hw::WeightLatticeElem) -> LieAlgebraModule{C}
    simple_module(L::LieAlgebra{C}, hw::Vector{<:IntegerUnion}) -> LieAlgebraModule{C}

Construct the simple module of the Lie algebra `L` with highest weight `hw`.

`L` needs to be a semisimple Lie algebra of characteristic $0$.
"""
function simple_module(L::LieAlgebra, hw::WeightLatticeElem)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  @req is_dominant(hw) "Not a dominant weight"
  struct_consts = lie_algebra_simple_module_struct_consts_gap(L, hw)
  dimV = size(struct_consts, 2)
  V = abstract_module(L, dimV, struct_consts; check=false)
  # TODO: set appropriate attributes
  return V
end

function simple_module(L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return simple_module(L, WeightLatticeElem(root_system(L), hw))
end

@doc raw"""
    dim_of_simple_module(
      [T = Int],
      L::LieAlgebra *or* R::RootSystem,
      hw::WeightLatticeElem *or* hw::Vector{<:IntegerUnion},
    ) -> T

Compute the dimension of the simple module with highest weight `hw`
using Weyl's dimension formula.

One can either provide a semisimple Lie algebra of characteristic $0$,
or its root system.

The return value is of type `T`.

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :A, 3);

julia> dim_of_simple_module(L, [1, 1, 1])
64
```

```jldoctest
julia> R = root_system(:B, 2);

julia> dim_of_simple_module(R, fundamental_weight(R, 1))
5
```
"""
function dim_of_simple_module(
  L::LieAlgebra, hw::Union{WeightLatticeElem,Vector{<:IntegerUnion}}
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return dim_of_simple_module(root_system(L), hw)
end

function dim_of_simple_module(
  T::Type, L::LieAlgebra, hw::Union{WeightLatticeElem,Vector{<:IntegerUnion}}
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return dim_of_simple_module(T, root_system(L), hw)
end

function dim_of_simple_module(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(R, WeightLatticeElem(R, hw))
end

function dim_of_simple_module(T::Type, R::RootSystem, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(T, R, WeightLatticeElem(R, hw))
end

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

@doc raw"""
    dominant_weights(
      L::LieAlgebra *or* R::RootSystem, 
      hw::WeightLatticeElem *or* hw::Vector{<:IntegerUnion},
    ) -> Vector{WeightLatticeElem}

Compute the dominant weights occurring in the simple module with highest weight `hw`,
sorted ascendingly by the total height of roots needed to reach them from `hw`.

One can either provide a semisimple Lie algebra of characteristic $0$,
or its root system.

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :B, 3);

julia> dominant_weights(L, [1, 0, 3])
7-element Vector{WeightLatticeElem}:
 w_1 + 3*w_3
 w_1 + w_2 + w_3
 2*w_1 + w_3
 3*w_3
 w_2 + w_3
 w_1 + w_3
 w_3
```

```jldoctest
julia> R = root_system(:B, 3);

julia> dominant_weights(R, 3 * fundamental_weight(R, 1) + fundamental_weight(R, 3))
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
function dominant_weights(
  L::LieAlgebra, hw::Union{WeightLatticeElem,Vector{<:IntegerUnion}}
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return dominant_weights(root_system(L), hw)
end

function dominant_weights(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_weights(R, WeightLatticeElem(R, hw))
end

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

@doc raw"""
    dominant_character(
      [T = Int],
      L::LieAlgebra *or* R::RootSystem,
      hw::WeightLatticeElem *or* hw::Vector{<:IntegerUnion}
    ) -> Dict{WeightLatticeElem, T}
    

Compute the dominant weights occurring in the simple module with highest weight `hw`
together with their multiplicities.

One can either provide a semisimple Lie algebra of characteristic $0$,
or its root system.

This function uses an optimized version of the Freudenthal formula, see [MP82](@cite) for details.

# Examples
```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> L = lie_algebra(QQ, :A, 3);

julia> dominant_character(L, [2, 1, 0])
Dict{WeightLatticeElem, Int64} with 4 entries:
  0           => 3
  2*w_2       => 1
  w_1 + w_3   => 2
  2*w_1 + w_2 => 1
```

```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> R = root_system(:B, 3);

julia> dominant_character(R, 2 * fundamental_weight(R, 1) + fundamental_weight(R, 3))
Dict{WeightLatticeElem, Int64} with 4 entries:
  w_2 + w_3   => 1
  2*w_1 + w_3 => 1
  w_1 + w_3   => 3
  w_3         => 6
```
"""
function dominant_character(
  L::LieAlgebra, hw::Union{WeightLatticeElem,Vector{<:IntegerUnion}}
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return dominant_character(root_system(L), hw)
end

function dominant_character(
  T::DataType, L::LieAlgebra, hw::Union{WeightLatticeElem,Vector{<:IntegerUnion}}
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return dominant_character(T, root_system(L), hw)
end

function dominant_character(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_character(R, WeightLatticeElem(R, hw))
end

function dominant_character(T::DataType, R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_character(T, R, WeightLatticeElem(R, hw))
end

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
    character(
      [T = Int],
      L::LieAlgebra *or* R::RootSystem,
      hw::WeightLatticeElem *or* hw::Vector{<:IntegerUnion},
    ) -> Dict{WeightLatticeElem, T}

Compute all dominant weights occurring in the simple module with highest weight `hw`
together with their multiplicities.

One can either provide a semisimple Lie algebra of characteristic $0$,
or its root system.

This is achieved by acting with the Weyl group on the [`dominant_character`](@ref dominant_character(::LieAlgebra, ::WeightLatticeElem)).

# Examples
```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> L = lie_algebra(QQ, :A, 3);

julia> character(L, [2, 0, 0])
Dict{WeightLatticeElem, Int64} with 10 entries:
  -2*w_3           => 1
  -2*w_2 + 2*w_3   => 1
  2*w_1            => 1
  -2*w_1 + 2*w_2   => 1
  -w_1 + w_3       => 1
  w_2              => 1
  w_1 - w_2 + w_3  => 1
  w_1 - w_3        => 1
  -w_1 + w_2 - w_3 => 1
  -w_2             => 1
```

```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> R = root_system(:B, 3);

julia> character(R, fundamental_weight(R, 3))
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
function character(L::LieAlgebra, hw::Union{WeightLatticeElem,Vector{<:IntegerUnion}})
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return character(root_system(L), hw)
end

function character(
  T::DataType, L::LieAlgebra, hw::Union{WeightLatticeElem,Vector{<:IntegerUnion}}
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return character(T, root_system(L), hw)
end

function character(R::RootSystem, hw::Vector{<:IntegerUnion})
  return character(R, WeightLatticeElem(R, hw))
end

function character(T::DataType, R::RootSystem, hw::Vector{<:IntegerUnion})
  return character(T, R, WeightLatticeElem(R, hw))
end

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

@doc raw"""
    tensor_product_decomposition(
      L::LieAlgebra *or* R::RootSystem,
      hw1::WeightLatticeElem *or* hw1::Vector{<:IntegerUnion},
      hw2::WeightLatticeElem *or* hw1::Vector{<:IntegerUnion},
    ) -> MSet{Vector{Int}}

Compute the decomposition of the tensor product of the simple modules with highest weights `hw1` and `hw2`
into simple modules with their multiplicities.

One can either provide a semisimple Lie algebra of characteristic $0$,
or its root system.

This function uses Klimyk's formula (see [Hum72; Exercise 24.9](@cite)).

The return type may change in the future.

# Examples
```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> L = lie_algebra(QQ, :A, 2);

julia> tensor_product_decomposition(L, [1, 0], [0, 1])
MSet{Vector{Int64}} with 2 elements:
  [0, 0]
  [1, 1]

julia> tensor_product_decomposition(L, [1, 1], [1, 1])
MSet{Vector{Int64}} with 6 elements:
  [0, 0]
  [1, 1] : 2
  [2, 2]
  [3, 0]
  [0, 3]
```

```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> R = root_system(:B, 2);

julia> tensor_product_decomposition(R, fundamental_weight(R, 1), fundamental_weight(R, 2))
MSet{Vector{Int64}} with 2 elements:
  [1, 1]
  [0, 1]

julia> tensor_product_decomposition(R, weyl_vector(R), weyl_vector(R))
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
  L::LieAlgebra,
  hw1::Union{WeightLatticeElem,Vector{<:IntegerUnion}},
  hw2::Union{WeightLatticeElem,Vector{<:IntegerUnion}},
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return tensor_product_decomposition(root_system(L), hw1, hw2)
end

function tensor_product_decomposition(
  R::RootSystem,
  hw1::Union{WeightLatticeElem,Vector{<:IntegerUnion}},
  hw2::Vector{<:IntegerUnion},
)
  return tensor_product_decomposition(R, hw1, WeightLatticeElem(R, hw2))
end

function tensor_product_decomposition(
  R::RootSystem, hw1::Vector{<:IntegerUnion}, hw2::WeightLatticeElem
)
  return tensor_product_decomposition(R, WeightLatticeElem(R, hw1), hw2)
end

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

###############################################################################
#
#   Demazure modules (via highest weight) of semisimple Lie algebras
#
###############################################################################

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
    demazure_operator(r::RootSpaceElem, w::WeightLatticeElem) -> Dict{WeightLatticeElem,Int}
    demazure_operator(r::RootSpaceElem, groupringelem::Dict{WeightLatticeElem,T}) where {T<:IntegerUnion} -> Dict{WeightLatticeElem,T}

Compute the action of the Demazure operator (see [Dem74](@cite)) associated to the positive root `r` on the given element of the group ring $\mathbb{Z}[P]$,
where $P$ denotes the additive group of the weight lattice.

If a single weight lattice element `w` is supplied, this is interpreted as `Dict(w => 1)`.

# Examples
```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
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
    demazure_character(
      [T = Int],
      L::LieAlgebra *or* R::RootSystem,
      w::WeightLatticeElem *or* w::Vector{<:IntegerUnion},
      x::WeylGroupElem *or* reduced_expr::Vector{<:IntegerUnion},
    ) -> Dict{WeightLatticeElem, T}

Compute all weights occurring in the Demazure module with extremal weight `w * x`
together with their multiplicities,
using the Demazure dimension formula (see [Dem74](@cite)).

One can either provide a semisimple Lie algebra of characteristic $0$,
or its root system.

Instead of a Weyl group element `x`, a reduced expression for `x` can be supplied.
This function may return arbitrary results if the provided expression is not reduced.

For Demazure characters of generalized flag manifolds, as in [PS09](@cite),
see [`demazure_character(::AbstractVector, ::PermGroupElem)`](@ref).

# Examples
```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> L = lie_algebra(QQ, :A, 2);

julia> demazure_character(L, [1, 1], [2, 1])
Dict{WeightLatticeElem, Int64} with 5 entries:
  2*w_1 - w_2  => 1
  w_1 + w_2    => 1
  0            => 1
  -w_1 + 2*w_2 => 1
  -2*w_1 + w_2 => 1
```

```jldoctest; filter = :(Oscar.doctestfilter_hash_changes_in_1_13())
julia> R = root_system(:B, 3);

julia> demazure_character(R, fundamental_weight(R, 2), weyl_group(R)([1, 2, 3]))
Dict{WeightLatticeElem, Int64} with 4 entries:
  w_1 + w_2 - 2*w_3 => 1
  w_1               => 1
  w_1 - w_2 + 2*w_3 => 1
  w_2               => 1
```
"""
function demazure_character(
  L::LieAlgebra,
  w::Union{WeightLatticeElem,Vector{<:IntegerUnion}},
  x::Union{WeylGroupElem,Vector{<:IntegerUnion}},
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return demazure_character(root_system(L), w, x)
end

function demazure_character(
  T::DataType,
  L::LieAlgebra,
  w::Union{WeightLatticeElem,Vector{<:IntegerUnion}},
  x::Union{WeylGroupElem,Vector{<:IntegerUnion}},
)
  @req is_zero(characteristic(L)) "Characteristic must be zero"
  @req is_semisimple(L) "Lie algebra not semisimple"
  return demazure_character(T, root_system(L), w, x)
end

function demazure_character(
  R::RootSystem, w::Vector{<:IntegerUnion}, x::Union{WeylGroupElem,Vector{<:IntegerUnion}}
)
  return demazure_character(R, WeightLatticeElem(R, w), x)
end

function demazure_character(
  T::DataType,
  R::RootSystem,
  w::Vector{<:IntegerUnion},
  x::Union{WeylGroupElem,Vector{<:IntegerUnion}},
)
  return demazure_character(T, R, WeightLatticeElem(R, w), x)
end

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
