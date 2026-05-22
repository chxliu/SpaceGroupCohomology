LoadPackage("SpaceGroupCohomology");

Print("Starting tests...\n\n");

# Test 1: Simple triclinic group
Print("Testing Space Group 1 (P1)...\n");
t1 := Runtime();
res1 := SpaceGroupCohomologyRingGapInterface(1);
t2 := Runtime();
if res1 = true then
    Print("Test passed for IT=1 (", t2-t1, " ms)\n\n");
else
    Print("Test failed for IT=1\n\n");
fi;

# Test 2: Orthorhombic group
Print("Testing Space Group 16 (P222)...\n");
t1 := Runtime();
res2 := SpaceGroupCohomologyRingGapInterface(16);
t2 := Runtime();
if res2 = true then
    Print("Test passed for IT=16 (", t2-t1, " ms)\n\n");
else
    Print("Test failed for IT=16\n\n");
fi;

# Test 3: Tetragonal group
Print("Testing Space Group 89 (P422)...\n");
t1 := Runtime();
res3 := SpaceGroupCohomologyRingGapInterface(89);
t2 := Runtime();
if res3 = true then
    Print("Test passed for IT=89 (", t2-t1, " ms)\n\n");
else
    Print("Test failed for IT=89\n\n");
fi;


Print("All tests finished.\n");
QUIT;
