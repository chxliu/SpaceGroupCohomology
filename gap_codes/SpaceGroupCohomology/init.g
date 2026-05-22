# init.g — read at package load time, *before* read.g.
# Purpose: declare the names of the package's data globals so that the
# parser does not warn "Unbound global variable" when it later sees
# them referenced in gap/functions.gi.
#
# Each name is filled in at read time by InstallValue(...) in gap/data.gi.

DeclareGlobalVariable( "PGGens230",
    "Generators of the 230 space groups in 4x4 matrix form." );
DeclareGlobalVariable( "IWP",
    "Irreducible Wyckoff positions for each of the 230 space groups." );
DeclareGlobalVariable( "GENNAMES",
    "Display names of mod-2 cohomology ring generators per space group." );
DeclareGlobalVariable( "funcs230",
    "Explicit cocycle function names per space group (up to degree 3)." );
