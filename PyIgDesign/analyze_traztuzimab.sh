#!/bin/sh


out_dir=traztuzimab/strategy_analysis


#strategyA_first_run strategyA_no_relax strategyB strategyC
#include_clusters

########################################################################################################################
#####  First Get PDBLISTs and score everything using 8 procs
########################################################################################################################
for strategy in strategyA strategyB strategyC
do
    for include in include_clusters exclude_clusters
    do

    in_dir=traztuzimab/decoys/$include/$strategy
    outname=$include.$strategy

    #./analyze_antibody_design_strategy.py --do_rescore --outdir $out_dir --indir $in_dir --out_name $outname --out_db_batch traz_6_strat --np 8 --native_path traztuzimab/1N8Z_coord_constrain_rel_1.pdb


    done
done

### Accidental runs:
for strategy in strategyA_no_relax
do
    for include in exclude_clusters include_clusters
    do

    in_dir=traztuzimab/decoys/$include/$strategy
    outname=$include.$strategy

    #./analyze_antibody_design_strategy.py --do_rescore --outdir $out_dir --indir $in_dir --out_name $outname --out_db_batch traz_6_strat --np 8 --native_path traztuzimab/1N8Z.pdb


    done
done

in_dir=traztuzimab/decoys/exclude_clusters/strategyA_combined_rel
outname=exclude_clusters.strategyA_combined_rel
#./analyze_antibody_design_strategy.py --do_rescore --outdir $out_dir --indir $in_dir --out_name $outname --out_db_batch traz_6_strat --np 8 --native_path traztuzimab/1N8Z_combined_rel_1.pdb




########################################################################################################################
#####   Then, get features databases using all procs.
########################################################################################################################

for strategy in strategyA strategyB strategyC strategyA_no_relax
do
    for include in include_clusters exclude_clusters
    do

    in_dir=traztuzimab/decoys/$include/$strategy
    outname=$include.$strategy

    #./analyze_antibody_design_strategy.py --do_run_antibody_features_all --do_run_cluster_features_all --do_run_cluster_features_top --do_run_antibody_features_top --outdir $out_dir --indir $in_dir --out_name $outname --out_db_batch traz_6_strat --np 8 &


    done
done



in_dir=traztuzimab/decoys/exclude_clusters/strategyA_combined_rel
outname=exclude_clusters.strategyA_combined_rel
#./analyze_antibody_design_strategy.py --do_run_antibody_features_all --do_run_cluster_features_all --do_run_cluster_features_top --do_run_antibody_features_top --outdir $out_dir --indir $in_dir --out_name $outname --out_db_batch traz_6_strat --np 8


########################################################################################################################
#####   Antibody Features on Relaxed Natives
########################################################################################################################

in_dir=traztuzimab/decoys/relaxed_natives
outname=native.relaxed

./analyze_antibody_design_strategy.py --do_run_antibody_features_all --outdir $out_dir --indir $in_dir --out_name $outname --out_db_batch traz_6_strat
