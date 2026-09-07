# Closure boxes

Julia allocates a `Core.Box` for a local variable that is captured by a closure
and is either assigned more than once, or assigned at a point where the
compiler cannot prove it is already defined. A box is a mutable cell with an
untyped field, so *every* read of that variable -- inside the closure and
outside it -- turns into a dynamic lookup, and type inference gives up on the
surrounding code as well. The cost is therefore not confined to the closure.

This is easy to trigger by accident:

```julia
acc = T[]
for i in 1:n
  acc = vcat(acc, [f(y, acc) for y in g(i)])  # the comprehension captures `acc`
end
```

and usually cheap to avoid, by giving the value the closure should see a name
of its own that is assigned exactly once:

```julia
acc = T[]
for i in 1:n
  seen = acc                                  # single-assigned, not boxed
  acc = vcat(seen, [f(y, seen) for y in g(i)])
end
```

Such an alias looks redundant, and OSCAR carries a fair number of them. Please
do not remove one without checking: doing so reintroduces the box.

## Shapes that box a captured variable

| Situation | Remedy |
|:--------- |:------ |
| assigned in both branches of an `if`/`else` | use an `if` *expression* |
| a local function calling itself, or two calling each other | move them to the top level |
| a local function with more than one method, including via default arguments | give it a single method, or move it to the top level |
| two sibling blocks reusing the same names | wrap *both* in `let` |
| a value that genuinely has to change | capture a `Ref` |
| a variable assigned after the closure that captures it | assign it before |

Note that a comprehension or generator captures only what its *body* uses; the
iterable is evaluated eagerly and is not captured. So in
`[h(x, v) for x in v]` only the `v` inside the body matters.

Watch out for macros that expand their argument more than once. In particular
`@vtime :Scope n x = expr` expands `expr` twice and hence assigns `x` twice;
write `x = @vtime :Scope n expr` instead.

## Boxes as a symptom

An assignment inside a nested function writes the enclosing local of the same
name rather than introducing a new one. A reported box therefore often marks a
place where the two uses were never meant to be the same variable, and fixing
it fixes a bug. Look at each report before mechanically renaming.

## Detection

[`Test.detect_closure_boxes`](https://github.com/JuliaLang/julia/pull/60478),
available since Julia 1.14, reports every method in a module whose lowered code
allocates a box, together with the variables responsible:

```julia
using Test, Oscar
for (m, vars) in Test.detect_closure_boxes(Oscar)
  println("Boxed variable(s) ", join(vars, ", "), " in ", m)
end
```

`test/ClosureBoxes.jl` runs this over all of OSCAR and requires the result to
be empty, so new boxes are caught by CI. The test body is skipped on older
Julia versions.
