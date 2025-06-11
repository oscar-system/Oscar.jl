#
# Patch the Alnuth package to call OSCAR instead of Pari
#
LoadPackage("alnuth"); # HACK

# store isomorphism between GAP and OSCAR univariate polynomial ring over the
# rationals as we need it a lot
BindGlobal("_PolyRingIso", IsoGapOscar(PolynomialRing(Rationals)));

# given a GAP number field, compute an isomorphic Oscar number field f
# TODO: also store the isomorphism
# TODO: turn this into an attribute resp. a method for JuliaWrapper / JuliaData ???
BindGlobal("_OscarField", function(F)
  local f;
  if not IsBound(F!.oscarField) then
    f := _PolyRingIso(IntegerDefiningPolynomial(F));
    F!.oscarField := Oscar_jl.number_field(f)[1];
  fi;
  return F!.oscarField;
end);

BindGlobal("MaximalOrderDescriptionOscar", function(F)
  local K, O, basis;

  K := _OscarField(F);
  O := Oscar_jl.maximal_order(K);
  basis := Julia.map(Oscar_jl.coordinates,Oscar_jl.basis(O,K));

  return JuliaToGAP(IsList, basis, true);
end);

BindGlobal("UnitGroupDescriptionOscar", function(F)
  local K, O, U_m, U, m, basis;

  K := _OscarField(F);
  O := Oscar_jl.maximal_order(K);
  U_m := Oscar_jl.unit_group(O);
  U := U_m[1];
  m := U_m[2];

  basis := Julia.map(Oscar_jl.coordinates, Julia.map(K, Julia.map(m, Oscar_jl.gens(U))));

  return JuliaToGAP(IsList, basis, true);
end);

BindGlobal("ExponentsOfUnitsDescriptionWithRankOscar", function(F, elms)
  local K, O, U_m, U, m, basis, units, rank, expns, x;

  K := _OscarField(F);
  O := Oscar_jl.maximal_order(K);
  U_m := Oscar_jl.unit_group(O);
  U := U_m[1];
  m := U_m[2];

  basis := Julia.map(Oscar_jl.coordinates, Julia.map(K, Julia.map(m, Oscar_jl.gens(U))));
  units := JuliaToGAP(IsList, basis, true);

  # the order of the torsion part of the full unit group
  rank := GAP_jl.GapObj(Oscar_jl.order(U_m[1][1]));

  expns := [];
  for x in elms do
     # map GAP int vector to Oscar element
     Add(expns, GAP_jl.GapObj(Oscar_jl.preimage(m, K(GAPToJulia(x))).coeff)[1]);
  od;

  # return result
  return rec(units := units, expns := expns, rank := rank);

end);

BindGlobal("ExponentsOfFractionalIdealDescriptionOscar", function(F, elms)
  local K, O, tmp, ideals, result, blah, vec, I;

  K := _OscarField(F);
  O := Oscar_jl.maximal_order(K);

  tmp := Julia.map(x -> Oscar_jl.factor(K(GAPToJulia(x)) * O), elms);
  # take the union of the keys
  ideals := Julia.collect(Julia.reduce(Oscar_jl.union, Julia.map(Oscar_jl.Set, Julia.map(Oscar_jl.keys, tmp))));
  if Julia.isempty(ideals) then
    return [];
  fi;

  result := [];
  for blah in JuliaToGAP(IsList, tmp) do
    vec := [];
    for I in JuliaToGAP(IsList, ideals) do
      Add(vec, Julia.get(blah, I, 0));
    od;
    Add(result, vec);
  od;

  return result;
end);

BindGlobal("NormCosetsDescriptionOscar", function(F, norm)
  local K, O, U_m, U, m, basis, units, eqn, creps;

  K := _OscarField(F);
  O := Oscar_jl.maximal_order(K);
  U_m := Oscar_jl.unit_group(O);
  U := U_m[1];
  m := U_m[2];

  basis := Julia.map(Oscar_jl.coordinates, Julia.map(K, Julia.map(m, Oscar_jl.gens(U))));
  units := JuliaToGAP(IsList, basis, true);

  eqn := Oscar_jl.norm_equation(O, norm);

  creps := JuliaToGAP(IsList, Julia.map(Oscar_jl.coordinates, Julia.map(K, eqn)), true);

  # return result
  return rec(units := units, creps := creps);
end);

BindGlobal("PolynomialFactorsDescriptionOscar", function(F, coeffs)
  local K, cf, poly, facs, result, f, g, i;

  K := _OscarField(F);

  cf := Oscar_jl.OscarInterfaceConstants._Vector_nf_elem(Reversed(List(coeffs, x -> K(GAPToJulia(x)))));
  poly := Oscar_jl.polynomial(K, cf);
  facs := Oscar_jl.factor(poly);
  Assert(0, Oscar_jl.is_one(facs.unit));

  result := [];
  for f in JuliaToGAP(IsList, Oscar_jl.collect(facs.fac)) do
    # Convert factor to GAP
    g := Julia.map(Oscar_jl.coordinates, f[1].coeffs);
    g := JuliaToGAP(IsList, g, true);
    g := Reversed(g);

    # add as many copies as necessary
    for i in [1..f[2]] do
      Add(result, g);
    od;
  od;

  # Sort result by ascending degrees
  SortBy(result, Length);

  return result;
end);

# Now make the functions also available under the *Pari name, so that the calling
# code does not have to be changed.
ReplaceGapFunc("MaximalOrderDescriptionPari", MaximalOrderDescriptionOscar);
ReplaceGapFunc("UnitGroupDescriptionPari", UnitGroupDescriptionOscar);
ReplaceGapFunc("ExponentsOfUnitsDescriptionWithRankPari", ExponentsOfUnitsDescriptionWithRankOscar);
ReplaceGapFunc("ExponentsOfFractionalIdealDescriptionPari", ExponentsOfFractionalIdealDescriptionOscar);
ReplaceGapFunc("NormCosetsDescriptionPari", NormCosetsDescriptionOscar);
ReplaceGapFunc("PolynomialFactorsDescriptionPari", PolynomialFactorsDescriptionOscar);

# To ensure we really don't call PARI anymore:
ReplaceGapFunc("ProcessPariGP", function(arg...) Error("ProcessPariGP not supported"); end);
