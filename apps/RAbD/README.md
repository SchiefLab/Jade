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
 
 - 8) Have a page in Excel made of up three columns.  new_name expressed rel_path.  The expressed column should be * with ones you will actually express.  New name is what you will now call the designs.  Copy and paste these three columns in a text file.  Save the file.  Run this script within the directory with all the PDBs (I recommend copying them into a new directory and then CDing into that first):
 
    ```
    rename_designs.py -i old_to_new.txt
    ```
 
 - 9) Now we get the sequences ready for ordering.  In this case, we have antibodies.  We will have to put them into a construct.  Currently, we have a script to create the orders for the Schief site, but we can generally make one using the code.   Open up the file and see what we are currently using as the vector, with an example sequence.  Use clustal omega to align the example sequence with one of your sequences. See what will need to be added or stripped.  Then use the following commands:
    
    ```
    get_seq.py --help
    ```
 Take a look at our options.  Now, for an example of an antibody order for IgG:
  
    ```
    get_seq.py --pdblist PDBLIST.txt --chain H --format IgG_order_heavy --outpath 1999_heavy_order.txt --pad_c_term S
    get_seq.py --pdblist PDBLIST.txt --chain L --format IgG_order_kappa --outpath 1999_light_order.txt
    ```
    
 Add any control sequences you want to order and you are now ready to send your sequences to genscript!
 
    