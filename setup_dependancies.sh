#!/usr/bin/env bash

#This script is to setup all dependancies for module_c
#You will need to install PyRosetta independantly

echo "PyRosetta must be installed separately for certain modules and applications.  Install Pip via Curl."

pip install biopython
pip install numpy
pip install scipy
pip install scikit-learn
pip install pandas
pip install seaborn
