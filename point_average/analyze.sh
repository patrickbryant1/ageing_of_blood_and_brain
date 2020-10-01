#!/usr/bin/env bash
ANNO=/hdd/pbryant/data/Methylation/anno450k.csv
BETAS=/hdd/pbryant/data/Methylation/Hannum_blood/joined_betas_full450k.csv
SAMPLES=/hdd/pbryant/data/Methylation/Hannum_blood/sample_sheet.tsv
AL='blood'
OUTDIR=../results/point_average/
#./analyze_blood.py --gene_annotations $ANNO --joined_betas $BETAS --sample_sheet $SAMPLES --agelabel $AL --outdir $OUTDIR

#Adjust the pvals, select markers and and visualize
RA=/home/pbryant/results/methylation/aging_of_blood/running_averages.npy
FCDF=../results/point_average/marker_max_FC_pval.csv
./adjust_and_vis.py --gene_annotations $ANNO --running_averages $RA --max_fold_change_df $FCDF --outdir $OUTDIR
