# uegCCD
Code for calculations on the uniform electron gas

# licence
MIT License
Copyright (c) 2018 Tom Henderson, James J. Shepherd, Gustavo Scuseria
For more details, see the LICENCE file

# references
Please cite:
(1) J. Chem. Phys. 140, 124102 (2014) 10.1063/1.4867783
(2) Phys. Rev. Lett. 112, 133002 (2014) 10.1103/PhysRevLett.112.133002

# how to use: Mac
1. On a mac with xCode (& command line tools) installed and the GNU Fortran (GCC) 6.3.0 compiler, run the make script (./make) to compile. 
2. ZCode inside directory ZRun will execute with a file labelled Input in the same directory.
3. Output will be produced including MP2 and CCD calculations

On a different system, you need to modify the make script to include
a reference to Lapack. This code will compile with Lapack, but by
default it's just set up to compile on a mac.

# input syntax

## self-explanatory UEG parameters:
\# Electrons:      14

\# rS Points:      1

Momentum Cutoff:  5

rSMin:            1.0

rSMax:            1.0

## these refer to the channels in JCP 140, 124102 (2014) and the pattern of inputs below corresponds to CCD:
Do Rings:         T

Do XRings:        T

Do Ladders:       T

Do Mosaics:       T

Rings Range:      0

XRings Range:     0

Ladders Range:    0

Mosaics Range:    0

DriverDir Range:  0

DriverEx Range:   0

Energy Range:     0

LinRings Range:   0

QuadRings Range:  0

DRings Range:     0

ExRings Range:    0

LinLadd Range:    0

QuadLadd Range:   0

DLadders Range:   0

ExLadders Range:  0
