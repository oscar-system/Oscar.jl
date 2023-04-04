@doc Markdown.doc"""
    su5_tate_model_over_arbitrary_3d_base()

Return the SU(5) Tate model over an arbitrary
3-dimensional base space. For more details see
e.g. [Wei18](@cite) and references therein.

```jldoctest
julia> su5_tate_model_over_arbitrary_3d_base()
Global Tate model over a not fully specified base
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
