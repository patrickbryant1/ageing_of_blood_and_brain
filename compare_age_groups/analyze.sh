#!/usr/bin/env bash
ANNO=/hdd/pbryant/data/Methylation/anno450k.csv
BETAS=/hdd/pbryant/data/Methylation/Hannum_blood/joined_betas_full450k.csv
SAMPLES=/hdd/pbryant/data/Methylation/Hannum_blood/sample_sheet.tsv
AL='blood'
OUTDIR=../results/
#./analyze_blood.py --gene_annotations $ANNO --joined_betas $BETAS --sample_sheet $SAMPLES --agelabel $AL --outdir $OUTDIR

#Compare significant markers
INDIR=../results/
./overlap_sig.py --gene_annotations $ANNO --joined_betas $BETAS --sample_sheet $SAMPLES --agelabel $AL --indir $INDIR --outdir $OUTDIR

#Visualize selected markers
#./vis_selected_markers.py --gene_annotations $ANNO --joined_betas $BETAS --sample_sheet $SAMPLES --agelabel $AL --overlapping_probes ../results/overlapping_probes_df.csv --overlapping_genes ../results/total_gene_df.csv --diff_probes ../results/diff_btw_corr_age_comp.npy --top10_corr ../results/top10_corr.csv --outdir $OUTDIR
