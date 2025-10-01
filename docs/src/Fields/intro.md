```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id fields)

OSCAR provides functionality for working with a wide variety of different fields.
The methods applicable to any field and its element can be found in the [interface](@ref fields_interface) section.
For functionality available only for specific fields, consult the corresponding section of the manual.

## Available fields

Here is a list of the fields available in OSCAR:

| Field        | How to create | Remark | Reference |
| ------------ | ----------- | --------- |------|
| $\mathbb{Q}$ | `rational_field()` | Also available as `QQ` | [Rationals](@ref rationals_section) |
| $\mathbb{F}_q$ | `GF(q)` | See also `finite_field` | [Finite fields](@ref) |
| $\mathbb{F}_q[X]/(f)$ | `finite_field(f)` |  | [Finite fields](@ref)|
| $\overline{\mathbb{Q}}$ | `algebraic_closure(QQ)` | | [Algebraic closure of the rational numbers](@ref) |
| $\overline{\mathbb{F}}_q$ | `algebraic_closure(F)` | | [Algebraic closure of finite prime fields](@ref)
| $\mathbb{Q}^{\mathrm{ab}}$ | `abelian_closure(QQ)` | | [Abelian closure of the rationals](@ref)
| $\mathbb{Q}[X]/(f)$ | `number_field(f)` |
| $\mathbb{Q}(\alpha) \subseteq \R$ | `embedded_number_field` | Ordered field
| $\mathbb{R}$ | `real_field()` | Ball arithmetic | [Arbitrary precision real balls](@ref)
| $\mathbb{C}$ | `complex_field()` | Ball arithmetic | [Arbitrary precision complex balls](@ref)
| $\mathbb{Q}_p$ | `padic_field(p)` |  | [Padics](@ref)
| $\mathbb{Q}_{p^n}$ | `qadic_field(p, n)` | Unramified extensions of $\mathbb{Q}_p$ | [Qadics](@ref)
| $R/(f)$ | `residue_field(R, f)` | $R$ must be a principal ideal domain
| $\mathrm{Quot}(R)$ | `fraction_field(R)` | $R$ must be an integral domain | [Generic fraction fields](@ref)

## Converting between fields

For fields $K$ and $L$ that admit a "canonical" embedding $K \to L$, elements from $K$ can be converted to elements from $L$ using "coercion" as in the following example:

```jldoctest
julia> a = QQ(2)
2

julia> Qbar = algebraic_closure(QQ);

julia> b = Qbar(a)
{a1: 2.00000}

julia> C = complex_field();

julia> sqrt(C(b))
[1.414213562373095049 +/- 3.45e-19]

julia> A = QQ[1 2; 3 4]
[1   2]
[3   4]

julia> B = change_base_ring(C, A)
[1.0000000000000000000   2.0000000000000000000]
[3.0000000000000000000   4.0000000000000000000]
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Claus Fieker](https://math.rptu.de/en/wgs/agag/people/head/fieker),
* [Tommy Hofmann](https://www.thofma.com/),
* [Max Horn](https://math.rptu.de/en/wgs/agag/people/head/prof-dr-max-horn).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
