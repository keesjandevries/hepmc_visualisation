#!/bin/bash

HEPMCDIR="/home/hep/kjd110/hepmc_visualisation/HepMC-2.06.09/test/"

cd $HEPMCDIR
./testHepMC -i ./pythia-output.hepmc  -o /home/hep/kjd110/hepmc_visualisation/HepMC_Visualisation/test.tex -n 60
cd -
