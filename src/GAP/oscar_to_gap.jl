## conversions of Oscar objects to GAP objects
## (extends the conversions from GAP.jl's `src/julia_to_gap.jl`,
## where low level Julia objects are treated)

import GAP: julia_to_gap

## WORKAROUND: the following method is also defined in GAP.jl, but apparently
## we need to redefine it here, otherwise it doesn't work
GAP.julia_to_gap(obj::Any; recursive::Bool) = GAP.julia_to_gap(obj)

## `fmpz` to GAP integer
GAP.julia_to_gap(obj::fmpz) = GAP.julia_to_gap(BigInt(obj))

## `fmpq` to GAP rational
GAP.julia_to_gap(obj::fmpq) = GAP.Globals.QUO(GAP.julia_to_gap(numerator(obj)), GAP.julia_to_gap(denominator(obj)))

## `fmpz_mat` to matrix of GAP integers
GAP.julia_to_gap(obj::fmpz_mat) = GAP.julia_to_gap(Matrix(obj), recursive = true)

## `fmpq_mat` to matrix of GAP rationals or integers
GAP.julia_to_gap(obj::fmpq_mat) = GAP.julia_to_gap(Matrix(obj), recursive = true)
