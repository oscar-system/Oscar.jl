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
DeclareOperation( "SerializationInOscarDependentObjects", [ IsObject ] );
DeclareOperation( "SerializeInOscar", [ IsObject, IsObject ] );
DeclareConstructor( "DeserializeInOscar", [ IsObject, IsObject, IsObject ] );
DeclareConstructor( "DeserializeInOscar", [ IsObject, IsObject, IsObject, IsObject ] );

############################################################################

# In Oscar 1.3.0, `Oscar_jl` was introduced to replace `Oscar` as the
# global GAP variable for the Oscar module. We keep the old name for
# compatibility.
DeclareSynonym( "Oscar", Oscar_jl );
