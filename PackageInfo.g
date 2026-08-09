#############################################################################
##
##  PackageInfo.g for the SpaceGroupCohomology GAP package
##
##  Reference: Chunxiao Liu and Weicheng Ye,
##    "Crystallography, group cohomology, and Lieb-Schultz-Mattis constraints",
##    SciPost Phys. 18, 161 (2025).  arXiv:2410.03607v3
##
SetPackageInfo( rec(

  PackageName    := "SpaceGroupCohomology",
  Subtitle       := "Mod-2 cohomology rings and LSM anomaly classes for the 230 space groups",
  Version        := "2.1.0",
  Date           := "2026-08-09",
  License        := "GPL-2.0-or-later",

  ##  The four URL fields below are required by GAP's package validator and
  ##  must each start with http:// , https:// or ftp:// .  The defaults
  ##  point at the published article and the arXiv preprint; replace them
  ##  with your own repository URLs once the package is hosted.
  PackageWWWHome  := "https://scipost.org/10.21468/SciPostPhys.18.5.161",
  README_URL      := "https://arxiv.org/abs/2410.03607",
  PackageInfoURL  := "https://arxiv.org/abs/2410.03607",
  ArchiveURL      := "https://arxiv.org/abs/2410.03607",
  ArchiveFormats  := ".tar.gz",

  Persons := [
    rec(
      LastName       := "Liu",
      FirstNames     := "Chunxiao",
      IsAuthor       := true,
      IsMaintainer   := true,
      Email          := "",
      WWW            := "",
      PostalAddress  := "",
      Place          := "",
      Institution    := "",
    ),
    rec(
      LastName       := "Ye",
      FirstNames     := "Weicheng",
      IsAuthor       := true,
      IsMaintainer   := true,
      Email          := "",
      WWW            := "",
      PostalAddress  := "",
      Place          := "",
      Institution    := "",
    ),
  ],

  Status         := "other",

  AbstractHTML   := Concatenation(
    "<span class=\"pkgname\">SpaceGroupCohomology</span> computes the ",
    "mod-2 cohomology ring and Lieb-Schultz-Mattis anomaly classes for ",
    "each of the 230 crystallographic space groups, using HAP as the ",
    "underlying resolution engine.  The package includes an optional C++17 ",
    "helper (<code>cpp/</code>, built as <code>bin/sgclinalg</code>) for ",
    "some GF(2) linear algebra, but it is not required; without it the ",
    "package uses GAP-native routines." ),

  PackageDoc := [
    rec(
      BookName  := "SpaceGroupCohomology",
      ArchiveURLSubset := [ "doc" ],
      HTMLStart := "doc/chap0.html",
      PDFFile   := "doc/manual.pdf",
      SixFile   := "doc/manual.six",
      LongTitle := "Mod-2 cohomology and LSM constraints for space groups",
    ),
  ],

  Dependencies := rec(
    GAP                    := ">= 4.15.0",
    NeededOtherPackages    := [ [ "HAP", ">= 1.30" ] ],
    SuggestedOtherPackages := [ ],
    ExternalConditions     := [ "polymake (https://polymake.org) must be on the PATH: HAP's CrystallographicComplex shells out to it to compute fundamental-domain vertices and face lattices" ],
  ),

  AvailabilityTest := ReturnTrue,

  TestFile := "tst/smoke.tst",

  Keywords := [ "space group", "cohomology", "Lieb-Schultz-Mattis",
                "anomaly", "crystallography", "HAP" ],
));
