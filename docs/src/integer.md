```@meta
CurrentModule = Oscar
```

# Integers

An important design decision in Oscar.jl is to use Julia as the user language
by default. This means that integers typed at the
[REPL](https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop)
are Julia integers.

For performance reasons, Oscar has its own integer format. These are
entered using the `ZZ` constructor.

```@repl
a = ZZ(2)^100
```

For convenience, many Oscar functions also accept Julia integers as inputs by
converting them to Oscar integers, especially if they do not fit in a machine
word. For example:

```@repl
R, x = ZZ["x"] # create a polynomial ring over the integers
f = 2x
```

In this example, `2` is a Julia integer but is still valid in the
call to the Oscar polynomial multiplication function that is implicit in the
expression `2x`.

!!! note

	In Julia, `2^64` will return `0`, as the Julia integer `2` is a machine
	word. In Oscar, the expression `ZZ(2)^64` will return the correct
	result (as an Oscar integer). In general, Oscar can only do an
	automatic conversion to an Oscar integer if the Julia integer is
	combined with another Oscar expression, or passed to an Oscar function.




