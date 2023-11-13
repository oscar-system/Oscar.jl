BindGlobal("ReplaceGapFunc", function(name, func)
  local orig_name;
  # Store a copy of the original value, but only if that copy does not yet
  # exist. This ensures we don't overwrite it during a call to `Reread`.
  orig_name := Concatenation("_ORIG_", name);
  if not IsBoundGlobal(orig_name) then
    BindGlobal(orig_name, ValueGlobal(name));
  fi;
  MakeReadWriteGlobal(name);
  UnbindGlobal(name);
  BindGlobal(name, func);
end);

############################################################################

InstallOtherMethod( ImagesRepresentative,
    [ IsActionHomomorphism and HasJuliaData, IsMultiplicativeElementWithInverse ],
function( hom, elm )
    local data;
    data:= JuliaData( hom );
    return Julia.Oscar.permutation(data[1], Julia.Oscar.group_element(data[2], elm)).X;
end );

InstallMethod( RestrictedMapping,
    CollFamSourceEqFamElms,
    [ IsActionHomomorphism and HasJuliaData, IsGroup ],
function( hom, H )
    local data, OscarG, xset, Omega, Hgens, Hacts, OscarH, res;
    data:= JuliaData( hom ); # the Oscar G-set and the acting Oscar group G
    OscarG:= data[2]; # the acting Oscar group G
    xset:= UnderlyingExternalSet( hom );
    Omega:= HomeEnumerator( xset ); # the set of Oscar objects
    Hgens:= GeneratorsOfGroup( H ); # GAP generators of H
    Hacts:= List( Hgens, x -> Julia.Oscar.group_element( OscarG, x ) ); # corresponding Oscar generators of H
    OscarH:= Julia.Oscar._as_subgroup_bare( OscarG, H );
    res:= ActionHomomorphism( H, Omega, Hgens, Hacts, FunctionAction( xset ) );
    SetJuliaData( res, [ data[1], OscarH ] );
    return res;
end );

############################################################################

# Install the methods for `IsoGapOscar`,
# in order to make Oscar's `iso_gap_oscar` work.
Perform( Oscar._iso_gap_oscar_methods,
         function( pair )
           InstallMethod( IsoGapOscar, [ JuliaToGAP( IsString, pair[1] ) ],
               Oscar.GAP.WrapJuliaFunc( pair[2] ) );
         end );

############################################################################

BindGlobal("_OSCAR_GroupElem", Oscar.AbstractAlgebra.GroupElem);

# the following code ensures that `GAP.Globals.Group` accepts a GAP list
# of Oscar group elements, as there is a GAP `OrbitStabilizerAlgorithm` method
# which sometimes does that when called by our `stabilize` method.
BindGlobal("JuliaObjectCollectionsFamily", CollectionsFamily(JuliaObjectFamily));
InstallMethod( IsGeneratorsOfMagmaWithInverses,
    "for a list or collection of Julia objects",
    f -> f= JuliaObjectCollectionsFamily,
    [ IsListOrCollection ],
    gens -> IsCollection( gens ) and
       ForAll( gens, x -> Julia.isa(x, _OSCAR_GroupElem) ) );

############################################################################

Perform( Oscar._GAP_serializations,
         function( entry )
           InstallMethod( SerializeInOscar,
             [ JuliaToGAP( IsString, entry[1] ), "IsObject" ],
             Julia.GAP.WrapJuliaFunc( entry[2] ) );
         end );

Perform( Oscar._GAP_deserializations,
         function( entry )
           InstallMethod( DeserializeInOscar,
             [ JuliaToGAP( IsString, entry[1] ), "IsObject", "IsObject", "IsObject" ],
             Julia.GAP.WrapJuliaFunc( entry[2] ) );
         end );
