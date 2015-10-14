#!/bin/bash



#Will make PDBLIST and run features on a particular benchmark.


#/make_pdblists_run_features.sh /home/jadolfbr/Documents/rosetta/projects/antibody_design/benchmarks/baseline_decoys $exp.$pdbs.$type.$benchmark cluster_features talaris2013 30ab.final False


##### Variable Setup #############

decoys=$1
benchmark=$2
feature_type=$3
xml=features/$3.xml

score=$4
batch=$5
ens=$6
ext=$7

out_dir=features/databases

cd $decoys


####### PDBLISTS ##################

#1) Top Designs
ls $benchmark*.pdb.gz > $benchmark.top.txt


#2) All Ensembles
ls ensemble*$benchmark*.pdb.gz > $benchmark.ens.txt

cd -



####### FEATURES ##################


#1) Top Designs
l=$decoys/$benchmark.top.txt

name=$benchmark.top.$feature_type

rosetta_scripts.$ext -parser:protocol $xml @ features/features.flag -l $l -in:path $decoys -out:path:all $out_dir -parser:script_vars name=$name score=$score batch=$batch
mv $name.$score.db3 databases




#2) All Ensembles if passed

if  [ "$ens" = "True" ]; then
l=$decoys/$benchmark.ens.txt
name=$benchmark.ens.$feature_type

rosetta_scripts.linuxclangrelease -database $ROSETTA3/database -parser:protocol $xml @ features/features.flag -l $l -in:path $decoys -out:path:all $out_dir -parser:script_vars name=$name score=$score batch=$batch
mv $name.$score.db3 features/databases

fi


