# DEPRECATED entry point.
#
# This file used to contain the cohomology-ring routines for the
# SpaceGroupCohomology project.  It has been moved into a proper GAP
# package; see ./SpaceGroupCohomology/ in this directory.
#
# Replace any:
#     Read("SpaceGroupCohomologyData.gi");
#     Read("SpaceGroupCohomologyFunctions.gi");
# with the single call:
#     LoadPackage("SpaceGroupCohomology");
#
# Reference: Chunxiao Liu and Weicheng Ye, SciPost Phys. 18, 161 (2025).
# arXiv:2410.03607

Print( "SpaceGroupCohomologyFunctions.gi is deprecated.\n",
       "Use  LoadPackage(\"SpaceGroupCohomology\");  instead.\n",
       "See gap_codes/SpaceGroupCohomology/README.md for installation.\n" );
