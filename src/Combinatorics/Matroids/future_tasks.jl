#This file is not meant to be included but rather serves as a list of our future tasks

is_graphic(M::Matroid) = true #TODO

is_cographic(M::Matroid) = is_graphic(dual_matroid(M))

is_quotient() = 0