export zero_object

zero_object(M::FreeMod) = free_module(base_ring(M), 0)
