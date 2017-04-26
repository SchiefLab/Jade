#!/bin/bash



#Will make PDBLIST and run features on a particular benchmark.


#/make_pdblists_run_features.sh /home/jadolfbr/Documents/rosetta/projects/antibody_design/benchmarks/baseline_decoys $exp.$pdbs.$type.$benchmark cluster_features talaris2013 30ab.final False


##### Variable Setup #############

decoy_list=$1
name=$2
feature_type=$3
xml=$3.xml

score=$4
batch=$5
ens=$6

####### FEATURES ##################


#1) Top Designs
l=$decoys/$benchmark.top.txt

name=$benchmark.top.$feature_type

rosetta_scripts.linuxclangrelease -database $ROSETTA3/database -parser:protocol $xml @ features.flag -l $decoy_list -in:path $decoys -out:path:all $out_dir -parser:script_vars name=$name score=$score batch=$batch
mv $name.$score.db3 features/databases





