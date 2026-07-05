# read.g — read at package load time, *after* init.g.
# Load the data tables first (they populate the globals declared in init.g),
# then load the function library that uses them.


ReadPackage( "SpaceGroupCohomology", "gap/data.gi" );
ReadPackage( "SpaceGroupCohomology", "gap/linalg.gi" );
ReadPackage( "SpaceGroupCohomology", "gap/homotopy.gi" );
ReadPackage( "SpaceGroupCohomology", "gap/functions.gi" );
