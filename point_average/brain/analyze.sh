#!/usr/bin/env bash
ANNO=/hdd/pbryant/data/Methylation/anno27k.csv
BETAS=/hdd/pbryant/data/Methylation/Brain/joined_betas.csv #27k samples
SAMPLES36194=/hdd/pbryant/data/Methylation/Brain/sample_sheet_36194.tsv
SAMPLES1575=/hdd/pbryant/data/Methylation/Brain/sample_sheet_15745.tsv

OUTDIR=../../results/point_average/brain/
./analyze_brain.py --gene_annotations $ANNO --joined_betas $BETAS --sample_sheet36194 $SAMPLES36194 --sample_sheet1575 $SAMPLES1575 --outdir $OUTDIR

#Adjust the pvals, select markers and and visualize
RA=/home/pbryant/results/methylation/aging_of_blood/running_averages.npy
FCDF=../../results/point_average/brain/marker_max_FC_pval.csv
MV=/home/pbryant/results/methylation/aging_of_blood/marker_values.npy
AGES=.../../results/point_average/brain/ages.csv
AGEPOINTS=../../results/point_average/brain/age_points.npy
HANNUM_MARKERS=./hannum_markers.csv
CORRELATIONRESULTS=../../results/point_average/brain/correlation_results.csv
#./adjust_and_vis.py --gene_annotations $ANNO --running_averages $RA --max_fold_change_df $FCDF --marker_values $MV --ages $AGES --age_points $AGEPOINTS --sample_sheet $SAMPLES --hannum_markers $HANNUM_MARKERS --correlation_results $CORRELATIONRESULTS --outdir $OUTDIR

#Plot the GO
RAGO=../../results/point_average/brain/genes/gene_chart_ra.txt
HANNUMGO=../../results/point_average/brain/genes/Hannum/gene_chart_hannum.txt
#./gene_ontology.py --ra_GO $RAGO --hannum_GO $HANNUMGO --outdir $OUTDIR
