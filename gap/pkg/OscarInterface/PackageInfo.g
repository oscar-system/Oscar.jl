#
# OscarInterface: GAP interface to OSCAR
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "OscarInterface",
Subtitle := "GAP interface to OSCAR",
Version := "1.5.0-DEV",
Date := "28/05/2025", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    LastName := "Breuer",
    FirstNames := "Thomas",
    IsAuthor := true,
    IsMaintainer := true,
    Email := "sam@math.rwth-aachen.de",
    WWWHome := "https://www.math.rwth-aachen.de/~Thomas.Breuer",
    PostalAddress := Concatenation(
               "Thomas Breuer\n",
               "Lehrstuhl für Algebra und Zahlentheorie\n",
                       "RWTH Aachen\n",
               "Pontdriesch 14/16\n",
               "52062 Aachen\n",
               "Germany" ),
    Place := "Aachen",
    Institution := "RWTH Aachen",
  ),
  rec(
    LastName := "Horn",
    FirstNames := "Max",
    IsAuthor := true,
    IsMaintainer := true,
    Email := "mhorn@rptu.de",
    WWWHome := "https://www.quendi.de/math",
    PostalAddress := Concatenation(
               "Fachbereich Mathematik\n",
               "RPTU Kaiserslautern-Landau\n",
               "Gottlieb-Daimler-Straße 48\n",
               "67663 Kaiserslautern\n",
               "Germany" ),
    Place := "Kaiserslautern, Germany",
    Institution := "RPTU Kaiserslautern-Landau",
  ),
],

SourceRepository := rec( Type := "git", URL := "https://github.com/oscar-system/Oscar.jl" ),
IssueTrackerURL := "https://github.com/oscar-system/Oscar.jl/issues",
PackageWWWHome := "https://www.oscar-system.org/",

ArchiveURL     := Concatenation( ~.PackageWWWHome, "OscarInterface-", ~.Version ),
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "README.md" ),

ArchiveFormats := ".tar.gz",

Status := "other",

AbstractHTML   :=  "",

PackageDoc := [],

Dependencies := rec(
  NeededOtherPackages := [
    ["JuliaInterface", ">=0.9"],
    ["Polycyclic", ">=2.16"],
  ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := [ "Oscar", "Julia", "Interface" ],

));
