#!/usr/bin/env Rscript
#Jared Adolf-Bryfogle

#This script accepts a single command-line argument with a JSON configuration file that specifies how to run the analysis.
#You must install the RosettaFeatures library first!

library(RosettaFeatures)


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("At single argument, the JSON configuration file, must be supplied (analysis_configuration).json", call.=FALSE)
} 

#Run the analysis
compare_sample_sources( config_filename=args[1] )

