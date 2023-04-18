@doc raw"""
    su5_tate_model_over_arbitrary_3d_base()

Return the SU(5) Tate model over an arbitrary
3-dimensional base space. For more details see
e.g. [Wei18](@cite) and references therein.

```jldoctest
julia> tm = su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base

julia> v = toric_ambient_space(tm)
Normal, simplicial toric variety

julia> a10,a21,a32,a43,a65,w,x,y,z = gens(cox_ring(v));

julia> I = ideal([x,y,w]);

julia> v2 = blow_up(v,I)
Normal toric variety

julia> cox_ring(v2)
Multivariate Polynomial Ring in 10 variables a10, a21, a32, x, ..., e over Rational Field graded by 
  a10 -> [0 0]
  a21 -> [0 0]
  a32 -> [0 0]
  x -> [1 0]
  y -> [0 1]
  z -> [-1 1]
  a43 -> [0 0]
  a65 -> [0 0]
  w -> [3 -2]
  e -> [-3 2]
```
"""
function su5_tate_model_over_arbitrary_3d_base()
    auxiliary_base_ring, (a10, a21, a32, a43, a65, w) = QQ["a10", "a21", "a32", "a43", "a65", "w"];
    a1 = a10;
    a2 = a21 * w;
    a3 = a32 * w^2;
    a4 = a43 * w^3;
    a6 = a65 * w^5;
    ais = [a1, a2, a3, a4, a6];
    return global_tate_model(ais, auxiliary_base_ring, 3)
end
