```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Field interface](@id fields_interface)

We list the methods that work for any field and its elements.

## Fields

| Method | Remark |
|:------ |:------ |
| `is_finite(K)` | Return whether $K$ is finite |
| `order(K)` | $\lvert K \rvert$ or error if the field is infinite |
| `characteristic(K)` | $\operatorname{char}(K)$ |

## Field elements

| Method | Remark |
|:------ |:------ |
| `is_one(a)` | |
| `is_zero(a)` | |
| `is_invertible(a)` | |
| `is_unit(a)` | |
| `one(F)` | |
| `zero(F)` | |
| `a + b` | |
| `a - b` | |
| `a * b` | |
| `a / b` | Errors if `b` is zero |
| `a^n` | |
| `is_power(a, n)` | |
| `is_square(a)` | | 
| `is_square_with_sqrt(a)` | |
