exclude = [:Nemo, :AbstractAlgebra, :Rational, :change_uniformizer, :genus_symbol, :data,
    :isdefintie, :narrow_class_group]

for i in names(Hecke)
  i in exclude && continue
  eval(Meta.parse("import Hecke." * string(i)))
  eval(Expr(:export, i))
end

# TODO: remove the following once Hecke or Nemo have IntegerUnion
# (and the version adding it is required in Project.toml)
if isdefined(Hecke, :IntegerUnion)
  import Hecke.IntegerUnion
elseif isdefined(Nemo, :IntegerUnion)
  import Nemo.IntegerUnion
else
  const IntegerUnion = Union{Integer, fmpz}
end
