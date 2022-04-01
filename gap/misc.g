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
