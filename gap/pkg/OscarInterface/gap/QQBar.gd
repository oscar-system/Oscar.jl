DeclareCategory( "IsQQBarFieldElement",
                 IsMultiplicativeElementWithInverse
                 and IsAdditiveElementWithInverse and IsZDFRE );

DeclareCategoryFamily( "IsQQBarFieldElement" );
DeclareCategoryCollections( "IsQQBarFieldElement" );
DeclareCategoryFamily( "IsQQBarFieldElementCollection" );
DeclareCategoryCollections( "IsQQBarFieldElementCollection" );

DeclareCategory( "IsQQBarField", IsField );

DeclareOperation( "_QQBarFieldElement", [ IsObject ] );
