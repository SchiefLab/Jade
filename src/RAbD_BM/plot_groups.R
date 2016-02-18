#!/usr/bin/env Rscript

#Used to plot data from pooled database.

#Args:
  #1) Outpath - Use Python to create temporary databases of given order.
  #2 - Databases to plot

require(RSQLite)
library(ggplot2)
require(plyr)

#library(Cairo)

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied!")
}


drv <- dbDriver("SQLite")


outdir = toString(args[1])
dbs = args[2:length(args)]

if (! file.exists(outdir)){
  dir.create(outdir)
}

output_pdf = function(p, outpath){

    outpath = paste(outpath, ".pdf", sep = "")
    #pdf(file = outpath, height=18, width=10.5)
    pdf(file = outpath)
    print(p)
    dev.off()

}

output_png = function(p, outpath){
    outpath = paste(outpath, ".png", sep = "")
    png(file = outpath, height=8, width=10.5, res=300, units="in")
    print(p)
    dev.off()
}

output_tiff = function(p, outpath){
    outpath = paste(outpath, ".tiff", sep = "")
    tiff(outpath, res = 600, compression = "lzw", height = 8, width = 10.5, units = "in")
    print(p)
    dev.off()
}

for (rec_type in c("length", "cluster")){

    plot_field_grid = function(p, plot_name, grid = NULL){
  
        if (! is.null(grid)){
            p <- p+ facet_grid(facets=grid)
        }
        outname = paste(plot_name, rec_type, sep="_")
        outpath = paste(outdir, outname, sep="/")

        output_pdf(p, outpath)
        output_png(p, outpath)
        output_tiff(p, outpath)


    }

    plot_field_wrap = function(p, plot_name, grid = NULL){
    
        if (! is.null(grid)){
          p <- p+ facet_wrap(facets=grid, ncol=3)
        }
  
        outname = paste(plot_name, rec_type, sep="_")
        outpath = paste(outdir, outname, sep="/")

        output_pdf(p, outpath)
        output_png(p, outpath)
        output_tiff(p, outpath)
  
    #Keeps saying plot should be a ggplot2 plot - which makes no damn sense!
    #for (plot_type in plot_types){
    #  out = paste(outname, plot_type, sep=".")
    #  outpath = paste(outdir, out, sep="/")
    #  print(outpath)
    #  ggsave(self, "test_path.png", p)
    #}
  
  
    }
    
    #Get selection according to name = assumes FROM xxx is last thing in string
    get_sele <- function(sele){
      
        exp_data = adply(dbs, 1, function(db_path){
            print(db_path)
            sele = paste(sele, rec_type, sep="_")
            con <- dbConnect(drv, toString(db_path))
            exp_data = dbGetQuery(con, sele)
            print(exp_data)
  
            #Rename docked here for clarity:
            exp_data[exp_data$exp_group=="docked", "exp_group"] = "baseline_w_dock"
            exp_data[exp_data$exp_group=="normal", "exp_group"] = "baseline"
  
            #Factor to retain order of creation - should actually be from the input list.  But this does not happen.
            exp_data$exp = factor(exp_data$exp, as.character(exp_data$exp))
            exp_data$exp_group = factor(exp_data$exp_group, as.character(exp_data$exp_group))
            exp_data$exp_type = factor(exp_data$exp_type, as.character(exp_data$exp_type))

            exp_data
        })
    
        exp_data
    }

    plot_types = c("png", "pdf", "eps")

    #Plot Experimental data

    exp_data = get_sele("SELECT exp, exp_group, exp_type, h3_present, top_rec, top_rr FROM exp")
    cdr_data = get_sele("SELECT exp, exp_group, exp_type, CDR, top_rec, top_rr FROM cdr")

    h3_present = rbind('Exclude H3', 'All CDRs')
    #h3_present_fac = factor(h3_present, as.character(h3_present))

    exp_data[exp_data$h3_present==1, "h3_present"] = h3_present[2]
    exp_data[exp_data$h3_present==0, "h3_present"] = h3_present[1]

    exp_data$h3_present = factor(exp_data$h3_present, as.character(exp_data$h3_present))

    #Plot Recovery
    p <- ggplot(data=exp_data, na.rm=T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_group, y= top_rec, fill=exp_type)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Recovery")) +
        theme_bw() +
        xlab("")+
        ylab("% Recovery")

    plot_field_grid(p, "exp_group_rec_by_exp", ~h3_present)


    #Plot Experimental RR
    p <- ggplot(data=exp_data, na.rm=T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_group, y= top_rr, fill=exp_type)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Risk Ratios")) +
        theme_bw() +
        xlab("")+
        geom_abline(intercept=1, slope=0) +
        ylab("Avg RR")

    plot_field_grid(p, "exp_group_rr_by_exp", ~h3_present)




    #Plot CDR data
    p <- ggplot(data=cdr_data, na.rm = T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_group, y = top_rec, fill=exp_type)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Recovery")) +
        theme_bw() +
        xlab("") +
        ylab("% Recovery")

    plot_field_wrap(p, "exp_group_rec_by_cdr", ~CDR)

    #Plot CDR RR
    p <- ggplot(data=cdr_data, na.rm=T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_group, y= top_rr, fill=exp_type)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Risk Ratios")) +
        theme_bw() +
        xlab("")+
        geom_abline(intercept=1, slope=0) +
        ylab("Avg RR")

    plot_field_wrap(p, "exp_group_rr_by_cdr", ~CDR)





    #############################Group By MinType ###############################################3


    #Plot Recovery
    p <- ggplot(data=exp_data, na.rm=T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_type, y= top_rec, fill=exp_group)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Recovery")) +
        theme_bw() +
        xlab("")+
        ylab("% Recovery")

    plot_field_grid(p, "min_group_rec_by_exp", ~h3_present)


    #Plot Experimental RR
    p <- ggplot(data=exp_data, na.rm=T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_type, y= top_rr, fill=exp_group)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Risk Ratios")) +
        theme_bw() +
        xlab("")+
        geom_abline(intercept=1, slope=0) +
        ylab("Avg RR")

    plot_field_grid(p, "min_group_rr_by_exp", ~h3_present)

    #Plot CDR data
    p <- ggplot(data=cdr_data, na.rm = T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_type, y = top_rec, fill=exp_group)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Recovery")) +
        theme_bw() +
        xlab("") +
        ylab("% Recovery")

    plot_field_wrap(p, "min_group_rec_by_cdr", ~CDR)

    #Plot CDR RR
    p <- ggplot(data=cdr_data, na.rm=T) +
        geom_bar(position="dodge", stat="identity", aes(x=exp_type, y= top_rr, fill=exp_group)) +
        ggtitle(paste("GraftDesign Native", rec_type, "Risk Ratios")) +
        theme_bw() +
        xlab("")+
        geom_abline(intercept=1, slope=0) +
        ylab("Avg RR")

    plot_field_wrap(p, "min_group_rr_by_cdr", ~CDR)

  
} #RecTypes



#Plot BoxPlot of all data - looks very strange
#all_data = get_sele("SELECT exp, exp_group, exp_type, CDR, native, top_rec, top_rr FROM all")

#p <- ggplot(data=all_data, na.rm = T) +
#  geom_point(data=all_data, size=1.2, position="jitter", aes(x=exp_type, y=top_rec, colour=exp_group)) + 
#  scale_y_continuous(limit=c(0, 100)) +
#  xlab("") +
#  ylab("% Recovery")
  
#plot_field_grid(p, "rec_by_exp_jitterplot")
#plot_field_wrap(p, "rec_by_cdr_jitterplot", ~CDR)

