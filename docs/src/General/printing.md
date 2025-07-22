```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Printing Options

## The 2 + 1 print modes of Oscar
Oscar has two user print modes `detailed` and `one line` and one internal
print mode `terse`. The latter is for use during recursion,
e.g. to print the `base_ring(X)` when in `one line` mode.
It exists to make sure that `one line` stays compact and human readable.

Top-level REPL printing of an object will use `detailed` mode by default
```julia-repl
julia> X
detailed
```
Inside nested structures, e.g. inside a `Vector`, the `one line` mode is used.
```julia-repl
julia> [X,X]
3-element Vector{TypeofX{T}}
 one line
 one line
 one line
```

#### An Example for the 2 + 1 print modes

detailed mode:
```jldoctest printing_example
julia> E = elliptic_curve(QQ, [-82, 0])
Elliptic curve
  over rational field
with equation
  y^2 = x^3 - 82*x
```

one line mode:
```jldoctest printing_example
julia> println(E)
Elliptic curve over QQ with equation y^2 = x^3 - 82*x
```

terse mode:
```jldoctest printing_example
julia> println(Oscar.terse(stdout),E)
Elliptic curve
```

## Unicode printing
Per default the output will be ASCII only (no Unicode). Implementors of
`Base.show` and related functions can branch on the output of
`is_unicode_allowed()` to display objects using non-ASCII characters.
This will then be used, if users have enabled Unicode using
`allow_unicode(true)`.

```@docs
allow_unicode(allowed::Bool; temporary::Bool=false)
with_unicode(f::Function, allowed::Bool=true)
is_unicode_allowed()
```

Objects may follow the value of `is_unicode_allowed()` at the time of their
creation for their printing, i.e. ignore later changes of the setting.
This is useful for objects storing a string representation of themselves, e.g.
generators of a module.

Here is an example with and without output using Unicode:
```julia-repl
julia> allow_unicode(false);

julia> R, (x,y) = QQ[:x,:y];

julia> Q,_ = quo(R,ideal([x,y]));

julia> variety(Q)
Affine variety
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (x, y)

julia> allow_unicode(true);

julia> variety(Q)
Affine variety
  in 𝔸² over QQ with coordinates [x, y]
defined by ideal (x, y)
```

## LaTeX and HTML printing

Some types support LaTeX or HTML output:
```julia-repl
julia> Qx, x = QQ[:x];

julia> show(stdout, MIME"text/latex"(), x^2 + 2x + x^10)
x^{10} + x^{2} + 2 x

julia> show(stdout, MIME"text/html"(), x^2 + 2x + x^10)
x^10 + x^2 + 2*x

julia> show(stdout, MIME"text/latex"(), Qx[x x^2; 1 1])
\begin{array}{cc}
x & x^{2} \\
1 & 1
\end{array}
```
Note: While using IJulia, you can enable LaTeX rendering in Jupyter notebooks using the following function:
```julia-repl
julia> Qx, x = QQ[:x];

julia> AbstractAlgebra.set_html_as_latex(true);

julia> show(stdout, "text/html", x^2 + 2x + x^10)
$x^{10} + x^{2} + 2 x$
```
