exclude = [:Nemo, :AbstractAlgebra, :Rational, :change_uniformizer, :genus_symbol,
    :isdefintie, :narrow_class_group]

for i in names(Hecke)
  i in exclude && continue
  eval(Meta.parse("import Hecke." * string(i)))
  eval(Expr(:export, i))
end

