#!/bin/bash

ifort MCSCO.f	\
readinput.f	\
corepart.f	\
MonteCarlo.f	\
Energy.f	\
dEnergy.f	\
latinit.f	\
nexttype.f	\
lattice.f	\
random.f
