#!/bin/bash

HEPMCDIR="/home/hep/kjd110/hepmc_visualisation/HepMC-2.06.09/test/"

echo copy source code to \"HepMC-2.06.09/test/\" folder and complie
cp testHepMC.cc $HEPMCDIR
cd $HEPMCDIR
make testHepMC
cd -
