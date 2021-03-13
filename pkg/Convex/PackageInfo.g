#############################################################################
##
##  PackageInfo.g       Convex package
##                      Martin Bies
##
##  Copyright 2021      University of Pennsylvania
##
##  A Gap package to do convex geometry by Polymake, Cdd and Normaliz
##
#############################################################################

SetPackageInfo( rec(

PackageName := "Convex",

Subtitle := "A Gap package to do convex geometry by Polymake, Cdd and Normaliz",

Version :=  Maximum( [
  "2021.03.13", # Martins version
] ),

Date := ~.Version{[ 1 .. 10 ]},
Date := Concatenation( ~.Date{[ 9, 10 ]}, "/", ~.Date{[ 6, 7 ]}, "/", ~.Date{[ 1 .. 4 ]} ),

License := "GPL-2.0-or-later",

Persons := [
rec(
    LastName      := "Bies",
    FirstNames    := "Martin",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email := "martin.bies@alumni.uni-heidelberg.de",
    WWWHome := "https://martinbies.github.io/",
    PostalAddress := Concatenation(
                 "Department of Mathematics\n",
                 "University of Pennsylvania\n",
                 "David Rittenhouse Laboratory\n",
                 "209 S 33rd St\n",
                 "Philadelphia\n",
                 "PA 19104\n",
                 "United States of America" ),
    Place         := "Pennsylvania",
    Institution   := "University of Pennsylvania"
  ),
],

Status := "dev",
PackageWWWHome := "https://github.com/homalg-project/ToricVarieties_project/tree/Convex/Convex/",
ArchiveFormats := ".zip",
ArchiveURL     := "https://github.com/homalg-project/ToricVarieties_project/releases/download/TO-COME/Convex.zip",
README_URL     := Concatenation( ~.PackageWWWHome, "README" ),
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),

AbstractHTML := "Convex provides functionality to perform convex geometry with Polymake, Cdd and Normaliz.",

PackageDoc := rec(
  BookName  := "Convex",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "A Gap package to do convex geometry by Polymake, Cdd and Normaliz",
  Autoload  := false
),


Dependencies := rec(
  GAP := ">=4.9",
  NeededOtherPackages := [ [ "AutoDoc", ">= 2019.05.20" ],
                           [ "NConvex", ">= 2020.11-04" ]
                         ],
  SuggestedOtherPackages := [ [ "JuliaInterface", ">= 0.5.2" ] ],
  ExternalConditions := [ ],
),


AvailabilityTest := function()

    return true;
  end,

Autoload := false,


Keywords := [ "Convex geometry", "Polymake", "Cdd", "Normaliz" ],

AutoDoc := rec(
    TitlePage := rec(
        Copyright := """
This package may be distributed under the terms and conditions
of the GNU Public License Version 2 or (at your option) any later version.
"""
    ),
),

));
