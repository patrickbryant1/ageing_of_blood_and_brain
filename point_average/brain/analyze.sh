#!/usr/bin/env bash
ANNO=/hdd/pbryant/data/Methylation/anno27k.csv
BETAS=/hdd/pbryant/data/Methylation/Brain/joined_betas.csv #27k samples
SAMPLES36194=/hdd/pbryant/data/Methylation/Brain/sample_sheet_36194.tsv
SAMPLES1575=/hdd/pbryant/data/Methylation/Brain/sample_sheet_15745.tsv
OUTDIR=../../results/point_average/brain/
#./analyze_brain.py --gene_annotations $ANNO --joined_betas $BETAS --sample_sheet36194 $SAMPLES36194 --sample_sheet1575 $SAMPLES1575 --outdir $OUTDIR

#Adjust the pvals, select markers and and visualize
###Cerebellum
RA=../../results/point_average/brain/cerebellum_running_averages.npy
FCDF=../../results/point_average/brain/cerebellum_marker_max_FC_pval.csv
MV=../../results/point_average/brain/cerebellum_marker_values.npy
AGES=../../results/point_average/brain/cerebellum_ages.csv
AGEPOINTS=../../results/point_average/brain/cerebellum_age_points.npy
HORVATH_MARKERS=./horvath_markers.csv
CORRELATIONRESULTS=../../results/point_average/brain/cerebellum_correlation_results.csv
NCLUSTERS=4
MEDIANRANGE='8,95'
OUTDIR=../../results/point_average/brain/cerebellum/
./adjust_and_vis.py --gene_annotations $ANNO --running_averages $RA --max_fold_change_df $FCDF --marker_values $MV --ages $AGES --point_indices $AGEPOINTS --sample_sheet36194 $SAMPLES36194 --sample_sheet1575 $SAMPLES1575 --horvath_markers $HORVATH_MARKERS --correlation_results $CORRELATIONRESULTS --n_clusters $NCLUSTERS --median_range $MEDIANRANGE --outdir $OUTDIR

###frontal_cortex
RA='../../results/point_average/brain/frontal_cortex_running_averages.npy'
FCDF='../../results/point_average/brain/frontal_cortex_marker_max_FC_pval.csv'
MV='../../results/point_average/brain/frontal_cortex_marker_values.npy'
AGES='../../results/point_average/brain/frontal_cortex_ages.csv'
AGEPOINTS='../../results/point_average/brain/frontal_cortex_age_points.npy'
HORVATH_MARKERS='./horvath_markers.csv'
CORRELATIONRESULTS='../../results/point_average/brain/frontal_cortex_correlation_results.csv'
OUTDIR=../../results/point_average/brain/frontal_cortex/
NCLUSTERS=1
MEDIANRANGE='8,96'
#./adjust_and_vis.py --gene_annotations $ANNO --running_averages $RA --max_fold_change_df $FCDF --marker_values $MV --ages $AGES --point_indices $AGEPOINTS --sample_sheet36194 $SAMPLES36194 --sample_sheet1575 $SAMPLES1575 --horvath_markers $HORVATH_MARKERS --correlation_results $CORRELATIONRESULTS --n_clusters $NCLUSTERS --median_range $MEDIANRANGE --outdir $OUTDIR
