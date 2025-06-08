# Toric Ideal Sheaves (Experimental)

Ideal sheaves on toric varieties are currently in experimental state.
Currently, we support the following functionality. Note that, as of
October 2023, this is limited to smooth toric varieties.
```@docs
ideal_sheaf(td::ToricDivisor)
ideal_sheaf(X::NormalToricVariety, I::MPolyIdeal)
```
