#!/usr/bin/env bash
ANNO=/hdd/pbryant/data/Methylation/anno450k.csv
BETAS=/hdd/pbryant/data/Methylation/Hannum_blood/joined_betas_full450k.csv
SAMPLES=/hdd/pbryant/data/Methylation/Hannum_blood/sample_sheet.tsv
AL='blood'
OUTDIR=../../results/point_average/blood/
#./analyze_blood.py --gene_annotations $ANNO --joined_betas $BETAS --sample_sheet $SAMPLES --agelabel $AL --outdir $OUTDIR

#Adjust the pvals, select markers and and visualize
RA=/home/pbryant/results/methylation/aging_of_blood/10/median/running_averages.npy
FCDF=../../results/point_average/blood/marker_max_FC_pval.csv
MV=/home/pbryant/results/methylation/aging_of_blood/10/median/marker_values.npy
AGES=../../results/point_average/blood/ages.csv
AGEPOINTS=../../results/point_average/blood/age_points.npy
HANNUM_MARKERS=./hannum_markers.csv
CORRELATIONRESULTS=../../results/point_average/blood/correlation_results.csv
NCLUSTERS=10
./adjust_and_vis.py --gene_annotations $ANNO --running_averages $RA --max_fold_change_df $FCDF --marker_values $MV --ages $AGES --point_indices $AGEPOINTS --sample_sheet $SAMPLES --hannum_markers $HANNUM_MARKERS --correlation_results $CORRELATIONRESULTS --n_clusters $NCLUSTERS --outdir $OUTDIR

#Plot the GO
RAGO=../../results/point_average/blood/genes/pantherChart.txt
HANNUMGO=../../results/point_average/blood/genes/Hannum/pantherChart.txt
#./gene_ontology.py --ra_GO $RAGO --hannum_GO $HANNUMGO --outdir $OUTDIR
