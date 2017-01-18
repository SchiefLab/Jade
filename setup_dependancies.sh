#!/usr/bin/env bash

#This script is to setup all dependancies for module_c
#You will need to install PyRosetta independantly

echo "PyRosetta must be installed separately for certain modules and applications.  Install Pip via Curl."
echo "PyMol must be installed separately for certain modules and applications"
echo "Clustal Omega must be installed separately for certain modules and application and in your $PATH"

pip install biopython
pip install numpy
pip install scipy
pip install scikit-learn
pip install pandas
pip install seaborn
pip install overrides
