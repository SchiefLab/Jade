# Antibody Design Strategy Analysis

 - 1) Run generate_rabd_features_dbs.py on each strategy.  Relative decoy paths are output.  Keep the decoy paths the same for further analysis.  This script uses the RunRosettaMPI code base to run the features reporters on a cluster.  Rosetta paths should be setup in your bashrc.  This has only been tested using Slurm, but should work for qsub as well.  Use ```--help``` for a full list of options.  A typical run looks like this:
    ```
    generate_rabd_features_dbs.py --indir L1L2 --np 101 --one_file_mpi --nstruct 1 --db_prefix L1L2_
    ```
 The ```--print_only``` option can be used to double check the run before actually running it. 5k decoys make about 15gb dbs at the moment due to extra h-bonding information.
 
 - 2) Move databases into a single directory.
 - 3) Load RAbD_Jade GUI.  See --help for more.  A typical command is this:
 
    ```
    RAbD_Jade.py --db_dir databases --analysis_name all_strat --native input_pdbs/19.99_incl_stratA.pdb \
    --pyigclassify_dir /Users/jadolfbr/Documents/projects/PyIgClassify/DBOUT
    ```
    
 - 4) Choose which CDRs you will analyze.  Give a name to your analysis.  Double click on strategies from the leftmost Listbox, which will select them for analysis.  
 
    Checkboxes are at the bottom for individual analysis and combined analysis.  Individual analysis will get the best decoys for each  strategy.  Combined will combine all of the strategies, and select the best decoys. 
 
    In the file menu, one can pick the type of scores to use.  By default, we output the best dG, dSASA, total_score, and dG by top 10% of the total score.  We can also filter models in the file menu, as well as set the top N, among other things.
 
    Most of your analysis can be complete in the Strategy Data and Model Data menus.  
    
   
 - 5) Select output CSV(Top) for both Strategy and Model data to get all of the data concatonated.  The Model data will list all score terms, aligned CDRs, lengths, and clusters with the relative path to the decoy.  The Strategy data will give you a nice summary of the strategy for each score term. 
 
 - 6) Select the PyMol menu item.  Select create_pymol_sessions (Top).  This will create pymol sessions for all top decoys for each strategy (individual and/or combined) for each score term we are using for descrimination.  It will run in parellel and should be quick.  If this does not work for you (lets say you are on a mac), make sure to alias the pymol executable so that it can be called from the command line.
 
 - 7) Begin working on an excel file.  Copy the CSV files into the excel file.  Have a tab for final expressed designs.  If you are going to choose designs based on a strategic rule, instead of by eye, you can use the ```copy_top_each_strategy.py``` script.  Use ```--help``` for more info.  CD into the individual or combined analysis directory.  Use this script to copy X top models for each strategy and for each scoring strategy into a final directory.  This requires you to already run the output PyMol sessions in our GUI. Here is an example of copying the best two models into our final directory:
 
    ```
    copy_top_each_strategy.py -i pdbs_sessions -o ../../expressed -n 2
    ```
 
    