
@doc raw"""
    newell_patch(k::Union{QQField, QQBarFieldElem}, n::Int=1)

Let $n$ be an integer between 1 and 32. Returns the ideal corresponding to 
the implicitization of the $n$-th bi-cubic patch describing
the Newell's teapot as a parametric surface.

The specific generators for each patch have been taken from [Tra04](@cite).
"""
function newell_patch(k::Union{QQField, QQBarFieldElem}, n::Int=1)
  return get_newell_patch_generators(n) |> ideal
end

@doc raw"""
    newell_patch(k::Field, n::Int=1)

Let $n$ be an integer between 1 and 32. Returns the ideal corresponding to 
the implicitization of the $n$-th bi-cubic patch describing
the Newell's teapot as a parametric surface.

The specific generators for each patch have been taken from [Tra04](@cite).

For fields $k\neq\mathbb{Q},\bar{\mathbb{Q}}$, this gives a variant of the ideal with 
integer coefficients.
"""
function newell_patch(k::Field, n::Int=1)
  F = get_newell_patch_generators(n)

  lcm_denom = [
      lcm(numerator.(coefficients(f))) for f in F
  ]
  integral_F = [
      change_coefficient_ring(
          k, l*f
      ) for (l, f) in zip(lcm_denom, F)
  ]

  return ideal(integral_F)
end

@doc raw"""
    newell_patch_with_orderings(k::Field, n::Int=1)

Let $n$ be an integer between 1 and 32. Returns the ideal corresponding to 
the implicitization of the $n$-th bi-cubic patch describing
the Newell's teapot as a parametric surface.
Additionally returns suitable start and target orderings, e.g. for use with the Gröbner walk.

The specific generators for each patch have been taken from [Tra04](@cite).

For fields $k\neq\mathbb{Q},\bar{\mathbb{Q}}$, this gives a variant of the ideal with 
integer coefficients.
"""
function newell_patch_with_orderings(k::Field, n::Int=1)
  I = newell_patch(k, n)
  R = base_ring(I)
  o1 = matrix_ordering(R, [1 1 1 0 0; 0 0 0 1 1; 0 0 0 1 0; 1 1 0 0 0; 1 0 0 0 0])
  o2 = matrix_ordering(R, [0 0 0 1 1; 1 1 1 0 0; 1 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0])
  return I, o2, o1
end

function get_newell_patch_generators(n::Int)
  R, (x,y,z,u,v) = polynomial_ring(QQ, [:x,:y,:z,:u,:v])

  if n == 1
    return [
        -x + 7//5 - 231//125 * v^2 + 39//80 * u^2 - 1//5 * u^3 + 99//400 * u * v^2 - 1287//2000 * u^2 * v^2 + 33//125 * u^3 * v^2 - 3//16 * u + 56//125 * v^3 - 3//50 * u * v^3 + 39//250 * u^2 * v^3 - 8//125 * u^3 * v^3,
        -y + 63//125 * v^2 - 294//125 * v + 56//125 * v^3 - 819//1000 * u^2 * v + 42//125 * u^3 * v - 3//50 * u * v^3 + 351//2000 * u^2 * v^2 + 39//250 * u^2 * v^3 - 9//125 * u^3 * v^2 - 8//125 * u^3 * v^3,
        -z + 12//5 - 63//160 * u^2 + 63//160 * u
    ]
  else # TODO: Add the other patches
    # TODO: Throw error
  end
end

