exclude = [:Nemo, :AbstractAlgebra, :Rational]

for i in names(Hecke)
  i in exclude && continue
  eval(Meta.parse("import Hecke." * string(i)))
  eval(Expr(:export, i))
end

