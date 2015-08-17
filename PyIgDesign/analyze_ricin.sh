#!/bin/sh


outdir=ricin_camelid/strategy_analysis
JAB_dir=ricin_camelid/decoys/JAB


weights=talaris2013_occ_sol

#Natives pareto optimal:
#./analyze_antibody_design_strategy.py --do_run_antibody_features_all --score_weights $weights --outdir $outdir --indir ricin_camelid/decoys/JAB/pareto_pre_min --out_name natives.pareto_optimal --out_db_batch JAB_SET_1

#Run1 designs
for strat in strat1.occ_sol strat1A.occ_sol strat1B.occ_sol strat2.occ_sol
do
    indir=$JAB_dir/$strat
    outname=$strat
    #./analyze_antibody_design_strategy.py --do_run_antibody_features_all --score_weights $weights --outdir $outdir --indir $indir --out_name $outname --out_db_batch JAB_SET_1 &
done


#Run Andrea's runs:
pre_backrub_out=ricin_camelid/backrub/output/pre-backrub/min
#./analyze_antibody_design_strategy.py --do_run_antibody_features_all --score_weights $weights --outdir $outdir --indir $pre_backrub_out --out_name pre_backrub_min --out_db_batch ANDREA_BR_MIN


#Dualspace Relaxed runs for Decoy Descrimination:

for strat in strat1.occ_sol strat1A.occ_sol strat1B.occ_sol strat2.occ_sol
do
    indir=$JAB_dir/$strat.rel
    outname=rel.$strat
    ./analyze_antibody_design_strategy.py --do_run_antibody_features_all --score_weights $weights --outdir $outdir --indir $indir --out_name $outname --out_db_batch JAB_SET_1_REL &
done
