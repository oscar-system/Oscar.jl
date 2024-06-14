DeclareAttribute( "JuliaData", IsObject );

############################################################################

# Use a GAP attribute for caching the mapping.
DeclareAttribute( "IsoGapOscar", IsDomain );

############################################################################

# Create GAP filters that describe
# - the union of `IsMultiplicativeElementWithInverseByPolycyclicCollector`
#   and `IsPcpElement` and
# - the union of `IsPcGroup` and `IsPcpGroup`.
DeclareFilter("IsPcElementOrPcpElement");
InstallTrueMethod(IsPcElementOrPcpElement, IsMultiplicativeElementWithInverseByPolycyclicCollector);
InstallTrueMethod(IsPcElementOrPcpElement, IsPcpElement);
BindGlobal("IsPcGroupOrPcpGroup", IsGroup and CategoryCollections(IsPcElementOrPcpElement));

############################################################################

# Use a GAP property for caching whether a fp/pc/pcp group is a full group
# and its stored generators are the generators for the defining presentation.
DeclareProperty( "GroupGeneratorsDefinePresentation", IsGroup );

############################################################################

# Use GAP operations for the serialization of GAP objects.
# (The methods will be Julia functions.)
DeclareOperation( "SerializeInOscar", [ IsObject, IsObject ] );
DeclareConstructor( "DeserializeInOscar", [ IsObject, IsObject, IsObject ] );
