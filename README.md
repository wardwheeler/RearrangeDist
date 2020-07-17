# rearrangeDist
Haskell program to create genome rearrangement distances from fastc input file.

Creates new POY-type "custom" alphabet file with associated "tcm" file containing 
distances.

Source files in src directory
Test data files in testData 
Documentation in doc (to come)

Requires cabal and ghc

To build:  cabal build 

Usage :  rearrangeDist inputFile method (breakpoint/inversion) topology (linear/circular) rearrangecost (integer) locusIndelcost (integer) fileStub 

