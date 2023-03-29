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
