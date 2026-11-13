################################################################################
# Upgrade Summary
# Types that are part of the mrdi namespace are renamed to names that do not
# depend on Julia, so that they can be shared with other systems.
# Univariate and multivariate polynomial rings share the name `PolynomialRing`
# and are told apart by their `_instance`.

const renamings_1_9_0 = Dict{String, Any}([
  ("Bool", "Boolean"),
  ("ZZRing", "Integers"),
  ("ZZRingElem", "Integer"),
  ("QQField", "Rationals"),
  ("QQFieldElem", "Rational"),
  ("Vector", "HomogeneousList"),
  ("Tuple", "InhomogeneousList"),
  ("FqFieldElem", "FiniteFieldElement"),
  ("PolyRing", Dict("name" => "PolynomialRing", "_instance" => "PolyRing")),
  ("MPolyRing", Dict("name" => "PolynomialRing", "_instance" => "MPolyRing")),
  ("MPolyRingElem", "Polynomial"),
])

push!(upgrade_scripts_set, UpgradeScript(
  v"1.9.0",
  function upgrade_1_9_0(s::UpgradeState, dict::AbstractDict{String, Any})
    upgrade_recursive(upgrade_1_9_0, s, dict)
    return rename_types(dict, renamings_1_9_0)
  end
))
